import numpy as np
import itertools
from astropy import stats
from astropy import units

from scipy import ndimage

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


def cut_slits(image_array, table):
    pass

def prepare_mdf(mdf_table, naxis1, naxis2, x_scale, y_scale, anamorphic_factor, wavelength_offset, spectral_pixel_scale,
                wavelength_start, wavelength_central, wavelength_end, y_distortion_coefficients=[1, 0, 0], arcsecpermm = 1.611444 * units.Unit('arcsec/mm')):
    """
    Prepares the MDF Table to reflect the pixel sections for the different slits

    Parameters
    ----------

    mdf_table : table like object (MDF)

    x_scale : ~astropy.units.Quantity
        in 'arcsec/pixel'

    y_scale : ~astropy.units.Quantity
        in 'arcsec/pixel'

    y_distortion_coefficients : ~np.ndarray
        y distortion correction coeffiencents in the form of coef[0] * y_pos + coef[1] * y_pos**2 + coef[2] * y_pos**3
        the coefficients are applied to the slits still in 'mm'


    """

    slit_pos_mx, slit_pos_my = mdf_table['slitpos_mx'], mdf_table['slitpos_my']
    slit_size_mx, slit_size_my = mdf_table['slitsize_mx'], mdf_table['slitsize_my']
    corrected_ypos = (y_distortion_coefficients[0] * slit_pos_my + y_distortion_coefficients[1] * slit_pos_my**2 +
                                    y_distortion_coefficients[2] * slit_pos_my**3) * units.Unit('mm')

    x_center =  naxis1 / 2.
    y_center = naxis2 / 2.

    slit_pos_y = corrected_ypos * arcsecpermm / y_scale
    slit_length = slit_size_my * arcsecpermm
    assert slit_pos_y.unit == units.Unit('pix')

    slit_pos_x = slit_pos_mx * arcsecpermm / x_scale
    slit_width = slit_size_mx * arcsecpermm
    assert slit_pos_x.unit == units.Unit('pix')

    slit_pos_x = slit_pos_x.value + x_center
    slit_pos_y = slit_pos_y.value + y_center

    spectrum_width = np.round(1.05 * (slit_length / y_scale).to('pix').value) * units.Unit('pix')
    spectrum_length = np.round(((wavelength_end - wavelength_start) / spectral_pixel_scale).to('pix').value) * units.Unit('pix')


    wavelength_central_pixel = (spectrum_length - (wavelength_central - wavelength_start) / spectral_pixel_scale).\
        to('pix').value
    #Simple correction for distortion in x
    y = (slit_pos_y / naxis1) - 0.05
    distortion_x = naxis1 *  (0.0014 * y - 0.0167 * y**2)



    x1 = np.int(np.round(x_center - ((x_center - slit_pos_x) / anamorphic_factor) - wavelength_central_pixel)) + \
         (wavelength_offset / spectral_pixel_scale).to('pix').value + distortion_x
    x2 = x1 + spectrum_length #Check if -1 is needed

    y1 = nint(yccd-center+l_yoff)
    y2 = y1 + specwid-1





