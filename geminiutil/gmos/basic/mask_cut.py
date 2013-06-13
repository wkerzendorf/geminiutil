import numpy as np
import itertools
from astropy import stats, units, table
from astropy.io import fits

from scipy import ndimage
from geminiutil.gmos.gmos_alchemy import GMOSMOSSlice
import logging

logger = logging.getLogger(__name__)

def find_mask_edges(flat_image, use_image_columns=slice(None), gauss_filter_sigma=3, sigma_clip_sigma=2, sigma_clip_iter=5):
    """
        Finding the mask edges with a gradient image

        Parameters
        ----------

        flat_image : ~numpy.ndarray
            Image (preferably a flat) to be used for finding the edges

    """

    gradient_image = np.diff(flat_image, axis=0)

    gradient_profile = ndimage.gaussian_filter(np.median(gradient_image[:, use_image_columns], axis=1),
                                               sigma=gauss_filter_sigma)

    clipped_gradient_profile = stats.sigma_clip(gradient_profile, sig=sigma_clip_sigma, iters=sigma_clip_iter, maout=True)

    peak_mask = clipped_gradient_profile.mask & (clipped_gradient_profile.data >= 0)
    trough_mask = clipped_gradient_profile.mask & (clipped_gradient_profile.data < 0)

    peaks = np.ma.MaskedArray(clipped_gradient_profile.data, mask=~peak_mask)
    troughs = np.ma.MaskedArray(clipped_gradient_profile.data, mask=~trough_mask)

    pixel_index = np.arange(len(peaks))
    lower_edge_groups = [list(group) for key, group in itertools.groupby(pixel_index, lambda x: peaks.mask[x])
                         if key==False]

    lower_edges = []
    for edge_group in lower_edge_groups:
        if len(edge_group) == 1:
            continue

        lower_edges.append(np.average(edge_group, weights=gradient_profile[edge_group]))

    lower_edges = np.array(lower_edges)

    upper_edge_groups = [list(group) for key, group in itertools.groupby(pixel_index, lambda x: troughs.mask[x])
                         if key==False]

    upper_edges = []
    for edge_group in upper_edge_groups:
        if len(edge_group) == 1:
            continue

        upper_edges.append(np.average(edge_group, weights=gradient_profile[edge_group]))

    upper_edges = np.array(upper_edges)

    return peaks, troughs, lower_edges, upper_edges


def cut_slits(chip_data, mdf_table):

    #check if the table has been prepared
    if 'SECX1' not in mdf_table.colnames:
        raise ValueError('The supplied table has not been prepared')

    final_hdu_list = fits.HDUList()

    for chip_id, data in enumerate(chip_data):
        for i, (sec_y1, sec_y2) in enumerate(mdf_table['SECY1', 'SECY2']):
            current_slice = slice(sec_y1, sec_y2)
            data_slice = data[current_slice].copy()


            final_hdu_list.append(fits.ImageHDU(data_slice, name='slice%d.chip%d.data' % (i, chip_id+1)))

    return final_hdu_list


def calculate_slice_geometries(science_frame, shift_bounds=[-20, 20], shift_samples=100, fit_sample=5):


    raw_fits = science_frame.raw_fits
    mdf_table = table.Table(raw_fits.mask.fits.fits_data[1].data)


    #calculating the slitsize and slit position from the MDF and instrument information
    naxis1, naxis2 = raw_fits.prepared_fits.fits.fits_data[1].header['naxis1'] * units.pix, \
                     raw_fits.prepared_fits.fits.fits_data[1].header['naxis2'] * units.pix

    slit_pos_x = mdf_table['slitpos_mx'] * units.mm
    slit_pos_x *= raw_fits.instrument_setup.x_pix_per_mm
    slit_pos_x += (naxis1 / 2) - 1 * units.pix

    slit_pos_y = np.polyval(raw_fits.instrument_setup.y_distortion_coefficients, mdf_table['slitpos_my']) * units.mm
    slit_pos_y *= raw_fits.instrument_setup.y_pix_per_mm
    slit_pos_y += (naxis2 / 2) - 1 * units.pix


    slit_size_x = mdf_table['slitsize_mx'] * units.mm * raw_fits.instrument_setup.x_pix_per_mm
    slit_size_y = mdf_table['slitsize_my'] * units.mm * raw_fits.instrument_setup.y_pix_per_mm


    slice_lower_edge = (slit_pos_y  - slit_size_y/2 + raw_fits.instrument_setup.y_offset).value
    slice_upper_edge = slice_lower_edge + slit_size_y.value


    if science_frame.flat is None:
        raise ValueError('science_frame does not have flat associated with it')

    flat_slices = np.median(science_frame.flat.fits.fits_data[2].data, axis=1)
    slice_model = np.zeros_like(flat_slices)

    for slice_lower, slice_upper in zip(slice_lower_edge, slice_upper_edge):
        lower_idx = np.int(np.round(slice_lower))
        upper_idx = np.int(np.round(slice_upper))
        slice_model[lower_idx:upper_idx] = 1.0

    slice_model *= np.median(flat_slices)

    rms_space = []
    pixel_shifts = np.linspace(shift_bounds[0], shift_bounds[1], shift_samples)
    for shift in pixel_shifts:
        rms_space.append(((ndimage.shift(slice_model, shift) - flat_slices)**2).sum())

    rms_space = np.array(rms_space)

    rms_slice = slice(rms_space.argmin()-5, rms_space.argmin()+5)
    a, b, c = np.polyfit(pixel_shifts[rms_slice], rms_space[rms_slice], 2)
    fitted_shift = -b/(2*a)


    slices = []
    for i in xrange(len(slice_lower_edge)):
        mdf_line = mdf_table[i]
        slices.append(GMOSMOSSlice(list_id=i, object_id=int(mdf_line['ID']), priority=int(mdf_line['priority']),
                            slice_set_id=science_frame.id, lower_edge=slice_lower_edge[i]+fitted_shift,
                            upper_edge=slice_upper_edge[i]+fitted_shift))
    return slices



