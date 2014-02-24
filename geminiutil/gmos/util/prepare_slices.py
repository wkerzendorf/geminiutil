import numpy as np
from astropy import units, table
from scipy import ndimage
from geminiutil.gmos.gmos_alchemy import GMOSMOSSlice
import logging

logger = logging.getLogger(__name__)


def calculate_slice_geometries(science_set, shift_bounds=[-20, 20], shift_samples=100, fit_sample=5):


    science_frame = science_set.science
    mdf_table = table.Table(science_frame.mask.fits.fits_data[1].data)


    #calculating the slitsize and slit position from the MDF and instrument information
    naxis1, naxis2 = science_frame.prepared.fits.fits_data[1].header['naxis1'] * units.pix, \
                     science_frame.prepared.fits.fits_data[1].header['naxis2'] * units.pix

    slit_pos_x = mdf_table['slitpos_mx'].astype(np.float64) * units.mm
    slit_pos_x *= science_frame.instrument_setup.x_pix_per_mm
    slit_pos_x *= slit_pos_x + (naxis1 / 2) - 1 * units.pix

    slit_pos_y = np.polyval(science_frame.instrument_setup.y_distortion_coefficients,
                            mdf_table['slitpos_my'].astype(np.float64)) * units.mm
    slit_pos_y *= science_frame.instrument_setup.y_pix_per_mm
    slit_pos_y += (naxis2 / 2) - 1 * units.pix


    slit_size_x = mdf_table['slitsize_mx'].astype(np.float64) * units.mm * science_frame.instrument_setup.x_pix_per_mm
    slit_size_y = mdf_table['slitsize_my'].astype(np.float64) * units.mm * science_frame.instrument_setup.y_pix_per_mm


    slice_lower_edge = (slit_pos_y  - slit_size_y/2 + science_frame.instrument_setup.y_offset).value
    slice_upper_edge = slice_lower_edge + slit_size_y.value


    if science_set.flat is None:
        raise ValueError('science_frame does not have flat associated with it')

    flat_slices = np.median(science_set.flat.fits.fits_data[2].data, axis=1)
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
