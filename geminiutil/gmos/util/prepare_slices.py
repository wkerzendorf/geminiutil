import numpy as np
from astropy import units, table
from scipy import ndimage, optimize

import logging

logger = logging.getLogger(__name__)


def calculate_slice_geometries_from_mdf(gmos_mos_fits):
    """
    Calculate the slice geometries from data given in the MDF file associated
    with the science set

    Parameters
    ----------

    gmos_mos_fits: ~geminiutil.gmos.gmos_alchemy.GMOSMOSRawFITS
        Raw FITS file associated with masks and instrument setups
    """
    mdf_table = gmos_mos_fits.mask.table

    if gmos_mos_fits.prepared is None:
        raise ValueError('given fits {0} is not prepared'.format(gmos_mos_fits))

    #calculating the slitsize and slit position from the MDF and instrument information
    naxis1, naxis2 = gmos_mos_fits.prepared.fits.fits_data[1].header['naxis1'] * units.pix, \
                     gmos_mos_fits.prepared.fits.fits_data[1].header['naxis2'] * units.pix

    slit_pos_x = mdf_table['slitpos_mx'].astype(np.float64) * units.mm
    slit_pos_x *= gmos_mos_fits.instrument_setup.x_pix_per_mm
    slit_pos_x *= slit_pos_x + (naxis1 / 2) - 1 * units.pix

    slit_pos_y = np.polyval(gmos_mos_fits.instrument_setup.y_distortion_coefficients,
                            mdf_table['slitpos_my'].astype(np.float64)) * units.mm
    slit_pos_y *= gmos_mos_fits.instrument_setup.y_pix_per_mm
    slit_pos_y += (naxis2 / 2) - 1 * units.pix


    slit_size_x = mdf_table['slitsize_mx'].astype(np.float64) * units.mm \
                  * gmos_mos_fits.instrument_setup.x_pix_per_mm
    slit_size_y = mdf_table['slitsize_my'].astype(np.float64) * units.mm \
                  * gmos_mos_fits.instrument_setup.y_pix_per_mm

    #Added -1 offset to go to numpy counting and not iraf counting.

    slice_lower_edge = (slit_pos_y  - slit_size_y/2 + gmos_mos_fits.instrument_setup.y_offset).value - 1
    slice_upper_edge = slice_lower_edge + slit_size_y.value - 1

    return slice_lower_edge, slice_upper_edge


def fit_slice_y_distortion(gmos_mos_fits, initial_slice_edges_guess,
                           method='Nelder-Mead', slice_model_chip=2,
                           slice_model_slice=slice(None), degree=5):
    """
    Fit slice

    Parameters
    ----------

    gmos_mos_fits: ~geminiutil.gmos.gmos_alchemy.GMOSMOSRawFITS
    slice_edges_initial_guess: tuple of ~np.ndarray
    method: str

    """
    if gmos_mos_fits.prepared is None:
        raise ValueError('given fits {0} is not prepared'.format(gmos_mos_fits))

    (initial_slice_lower_edges,
     initial_slice_upper_edges) = initial_slice_edges_guess


    raw_data = gmos_mos_fits.prepared.fits.fits_data[slice_model_chip].data
    data_1D_slice = np.median(raw_data[slice_model_slice], axis=1)

    distorted_model = DistortedSliceModel(data_1D_slice,
                                          initial_slice_lower_edges,
                                          initial_slice_upper_edges,
                                          degree=degree)
    initial_p_coeff = distorted_model.p_coef.copy()
    optimize.minimize(distorted_model.fit, initial_p_coeff, method=method)

    return distorted_model.polynomial(initial_slice_lower_edges), \
           distorted_model.polynomial(initial_slice_upper_edges)



#TODO: make this a model
class DistortedSliceModel(object):
    """
    Class to make a distorted slice model for fitting

    Parameters
    ----------

    flat_slices: ~np.ndarray
    slice_lower_edges: ~np.ndarray
        lower edges of the slices
    slice_upper_edges: ~np.ndarray
        upper edges of the slices
    normalize_factor: ~float

    """

    def __init__(self, flat_slices, slice_lower_edges, slice_upper_edges,
                 normalize_factor=None, auto_normalize_threshold=5000,
                 degree=5):
        self.slice_model = np.zeros_like(flat_slices)

        for slice_lower, slice_upper in zip(slice_lower_edges, slice_upper_edges):
            lower_idx = np.int(np.round(slice_lower))
            upper_idx = np.int(np.round(slice_upper))
            self.slice_model[lower_idx:upper_idx] = 1.0


        if normalize_factor is None:
            normalize_factor = np.median(
                flat_slices[flat_slices > auto_normalize_threshold])

        self.slice_model *= normalize_factor

        self.y = np.indices(self.slice_model.shape)[0]

        self.flat_slices = flat_slices
        self.p_coef = np.ones(degree) * 1e-14
        self.p_coef[1] = 1.0

        self.polynomial = np.polynomial.Polynomial(self.p_coef,
                                  domain=[0, len(self.slice_model)],
                                  window=[0, len(self.slice_model)])


    def __call__(self, p_coef):

        self.polynomial.coef = p_coef
        y_distorted = self.polynomial(self.y)

        return np.interp(self.y, y_distorted, self.slice_model)


    def fit(self, p_coef):
        self.p_coef = p_coef
        return np.sum((self(p_coef) - self.flat_slices)**2)