import numpy as np
from astropy.io import fits
from .. base import FITSFile
import os

class GMOSMOSChipCombine(object):

    operation_prefix = 'c'

    def __init__(self, fits_object, chip1_offset=3, chip2_offset=4, chip3_offset=-5, reduction_dir='.'):
        self.fits_object = fits_object
        self.fits_in_id = fits_object.id
        self.chip1_offset = chip1_offset
        self.chip2_offset = chip2_offset
        self.chip3_offset = chip3_offset
        new_fname = os.path.join(reduction_dir, self.operation_prefix + fits_object.fname)
        self._combine_chips(new_fname)

    def _combine_chips(self, fname):
        chip1_data = self.fits_object.fits_data[1].data * self.chip1_offset
        chip2_data = self.fits_object.fits_data[2].data * self.chip2_offset
        chip3_data = self.fits_object.fits_data[2].data * self.chip3_offset
        new_data = np.zeros((chip1_data.shape[0], chip1_data.shape[1] + chip2_data.shape[1] + chip3_data.shape[1]))
        new_data[:, 0:chip1_data.shape[1]] = chip1_data
        new_data[:, chip1_data.shape[1]:chip1_data.shape[1] + chip2_data.shape[1]] = chip2_data
        new_data[:, chip1_data.shape[1] + chip2_data.shape[1]:] = chip2_data
        new_hdu_list = fits.HDUList(hdus=[self.fits_object.fits_data[0]])
        new_hdu = fits.ImageHDU(new_data)
        new_hdu_list.append(new_hdu)
        new_hdu_list.writeto(fname)
        return FITSFile.from_fits_file(fname)


