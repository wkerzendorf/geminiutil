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

        #make sure we only have data for 3 fits files. This does not yet work for GMOS-N
        assert len(fits_data) - 1 == 3

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
            amplifier_data.name = 'chip%d.data' % i

            final_hdu_list += [amplifier_data]

            fits.HDUList(final_hdu_list).writeto(full_path, clobber=True)
        return FITSFile.from_fits_file(full_path)


class GMOSMOSCutSlits(object):
        __tablename__ = 'gmosmos_cut_slits'

        file_prefix = 'cut'

        def __init__(self, mode='header_cut'):
            pass

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


            return prepared_mdf_table

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

