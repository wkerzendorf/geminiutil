from .. base import Base, FITSFile
import astropy.io.fits as fits
from astropy.nddata import NDData, StdDevUncertainty
import numpy as np
from numpy.polynomial.polynomial import polyfit as polyfit
import logging
import os
from collections import OrderedDict

from .basic import prepare, mask_cut

logger = logging.getLogger(__name__)

"""Basic reduction classes Prepare and GMOSPrepare

Typical usage:
Prep = Prepare()   # can give bias_subslice, data_subslice, clip, gains, ...
infits = fits.open(fitsfile)
outfits = Prep(infits)

***TODO*** 
Ensure the out fits headers contain everything required to reproduce the result
"""

class GMOSPrepare(object): # will be base when we know what
    __tablename__ = 'gmos_prepare'

    file_prefix = 'prep'

    def __init__(self, bias_subslice=None,
                 data_subslice=None,
                 bias_clip_sigma=3.):
        self.bias_subslice = bias_subslice
        self.data_subslice = data_subslice
        self.bias_clip_sigma = bias_clip_sigma

    def __call__(self, gmos_raw_object, fname=None, destination_dir='.', write_steps=False, write_cut_image=None):
        """
        Preparing the Image

        Parameters
        ----------

        gmos_raw_object :

        fname :
            output filename

        write_steps :
            write out the individual steps to the individual fits files

        write_cut_image :
            write out the cut_image to a fits file with name specified in variable



        """

        if fname is None:
            fname = '%s-%s' % (self.file_prefix, gmos_raw_object.fits.fname)

        full_path = os.path.join(destination_dir, fname)

        fits_data = gmos_raw_object.fits.fits_data

        final_hdu_list = [fits_data[0].copy()]

        for i in xrange(1, len(fits_data)):
            current_amplifier = fits_data[i]
            detector = gmos_raw_object.instrument_setup.detectors[i - 1]



            #####
            # Subtracting Overscan
            #####
            amplifier_data = prepare.correct_overscan(current_amplifier, 
                                                      self.bias_subslice, 
                                                      self.data_subslice,
                                                      self.bias_clip_sigma)
            #####
            #Correcting the amplifier data gain
            #####

            amplifier_data = prepare.correct_gain(amplifier_data, gain=detector.gain)
            amplifier_data.name = 'DATA_%d' % i
            bias_uncertainty = amplifier_data.header['BIASSTD']


            #####
            #Create uncertainty Frame
            ####
            # ***TODO*** bias uncertainty should not be here, since this is correlated
            # for all pixels.  What is important later is a flat field/sensitivity error.

            amplifier_uncertainty = create_uncertainty_frame(amplifier_data, readout_noise=detector.readout_noise,
                                                             bias_uncertainty=bias_uncertainty)

            amplifier_uncertainty.name = 'UNCERTAINTY_%d' % i

            ####
            #CREATE MASK FRAME
            ####

            # TODO add BPM from gmos_data

            amplifier_mask = create_mask(amplifier_data)
            amplifier_mask.name = 'MASK_%d' %i

            final_hdu_list += [amplifier_data, amplifier_uncertainty, amplifier_mask]
        ######
        #MOSAICING DATA
        #####
        final_hdu_list = fits.HDUList(final_hdu_list)
        final_hdu_list = prepare.mosaic(final_hdu_list, chip_gap=gmos_raw_object.instrument_setup.chip_gap)
        if write_steps:
            mosaic_fname = 'mosaic-%s' % gmos_raw_object.fits.fname
            final_hdu_list.writeto(os.path.join(destination_dir, mosaic_fname), clobber=True)
        #####
        #Cutting the Mask
        #####
        mdf_table = gmos_raw_object.mask.fits.fits_data['MDF'].data
        naxis1 = final_hdu_list['DATA'].header['naxis1']
        naxis2 = final_hdu_list['DATA'].header['naxis2']

        current_instrument_setup = gmos_raw_object.instrument_setup
        x_scale = current_instrument_setup.x_scale
        y_scale = current_instrument_setup.y_scale

        wavelength_offset = current_instrument_setup.grating.wavelength_offset
        spectral_pixel_scale = current_instrument_setup.spectral_pixel_scale
        wavelength_start = current_instrument_setup.wavelength_start
        wavelength_end = current_instrument_setup.wavelength_end
        wavelength_central = current_instrument_setup.grating_central_wavelength
        y_distortion_coefficients = current_instrument_setup.y_distortion_coefficients
        y_offset = current_instrument_setup.y_offset
        anamorphic_factor = current_instrument_setup.anamorphic_factor
        arcsecpermm = current_instrument_setup.arcsecpermm
        prepared_mdf_table = mask_cut.prepare_mdf_table(mdf_table, naxis1, naxis2, x_scale, y_scale, anamorphic_factor,
                                                        wavelength_offset, spectral_pixel_scale, wavelength_start,
                                                        wavelength_central, wavelength_end,
                                                        y_distortion_coefficients=y_distortion_coefficients,
                                                        arcsecpermm = arcsecpermm, y_offset=y_offset)

        if write_cut_image:
            cut_image, cut_hdu_list = mask_cut.cut_slits(final_hdu_list['DATA'].data, prepared_mdf_table,
                                            uncertainty=final_hdu_list['UNCERTAINTY'].data,
                                            mask=final_hdu_list['MASK'].data, return_cut_image=True)

            fits.ImageHDU(data=cut_image).writeto(os.path.join(destination_dir, write_cut_image), clobber=True)
        else:
            cut_hdu_list = mask_cut.cut_slits(final_hdu_list['DATA'].data, prepared_mdf_table,
                                            uncertainty=final_hdu_list['UNCERTAINTY'].data,
                                            mask=final_hdu_list['MASK'].data)




        cut_hdu_list.insert(0, final_hdu_list[0])
        fits.HDUList(cut_hdu_list).writeto(full_path, clobber=True)
        return FITSFile.from_fits_file(full_path)

class Prepare(object):
    def __init__(self, bias_subslice=[slice(None), slice(1,11)], 
                 data_subslice=[slice(None), slice(-1)], clip=3.,
                 gain=None, read_noise=None, combine=True, 
                 overscan_std_threshold=3.):
        """Class that extracts and combines bias- and gain-corrected extensions.

        Wrapper around function ~basic.prepare.prepare

        Returns
        -------
        Class instance which can be applied to input fits files, using call
        method (which calls functions correct_overscan, correct_gain, and
        combine_halves with parameters set here)

        Examples
        --------
        PrepDA = gmos_basics.Prepare(data_subslice=[slice(72,172),slice(-1)])
        dayarc = PrepDA(fits.open('N20120423S0019.fits'))
        """
        self.bias_subslice = bias_subslice
        self.data_subslice = data_subslice
        self.clip = clip
        self.gain = gain
        self.read_noise = read_noise
        self.combine = combine
        self.overscan_std_threshold = overscan_std_threshold

    def __call__(self, image):
        """Measure bias from overscan, correct for gain, combine halves.

        Notes
        -----
        calls ~basic.prepare.prepare
        """
        return prepare.prepare(image, self.bias_subslice, 
                               self.data_subslice, self.clip, self.gain,
                               self.read_noise, self.combine,
                               self.overscan_std_threshold)

def create_uncertainty_frame(amplifier, readout_noise=None, 
                             bias_uncertainty=None):
    """
    Create an standard deviation uncertainty frame for GMOS

    Parameters
    ----------

    amplifier: ~astropy.io.fits.ImageHDU

    readout_noise: ~float
        e-; default (~None), from header['RDNOISE']
    bias_uncertainty: ~float
        e-; default (~None), use header['BIASUNC']*header['GAINUSED']

    ***TODO*** bias uncertainty should not be here, since this is correlated 
    for all pixels.  What is important later is a flat field/sensitivity error.
    """

    if readout_noise is None:
        readout_noise = amplifier.header['RDNOISE']

    if bias_uncertainty is None:
        bias_uncertainty = (amplifier.header['BIASUNC'] * 
                            amplifier.header['GAINUSED'])

    uncertainty_data = np.sqrt(np.abs(amplifier.data) + 
                               readout_noise**2 + bias_uncertainty**2)
    uncertainty = fits.ImageHDU(uncertainty_data, amplifier.header)
    uncertainty.header['uncertainty_type'] = 'stddev'

    return uncertainty

def create_mask(amplifier, min_data=0, max_data=None, template_mask=None):
    if template_mask is None:
        mask_data = np.zeros_like(amplifier.data, dtype=bool)
    else:
        mask_data = template_mask

    #Make this more complex - like loading initial bpm from somewhere else

    if min_data is not None:
        mask_data |= amplifier.data < min_data

    if max_data is not None:
        mask_data |= amplifier.data > max_data

    mask_data = mask_data.astype(dtype=np.uint8)

    return fits.ImageHDU(mask_data, header=amplifier.header)

def gmos_ccd_image_arithmetic(func, *args):
    pass

