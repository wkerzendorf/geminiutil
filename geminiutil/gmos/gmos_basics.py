from .. base import Base, FITSFile
import astropy.io.fits as fits
from astropy.nddata import NDData, StdDevUncertainty
import numpy as np
from numpy.polynomial.polynomial import polyfit as polyfit
import logging
import os
from collections import OrderedDict

from .basic import prepare

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

    def __init__(self, bias_subslice=[slice(None), slice(1, 11)], 
                 data_subslice=[slice(None), slice(-1)],
                 bias_clip_sigma=3.):
        self.bias_subslice = bias_subslice
        self.data_subslice = data_subslice
        self.bias_clip_sigma = bias_clip_sigma

    def __call__(self, gmos_raw_object, fname=None, destination_dir='.'):

        if fname is None:
            fname = '%s-%s' % (self.file_prefix, gmos_raw_object.fits.fname)

        full_path = os.path.join(destination_dir, fname)

        fits_data = gmos_raw_object.fits.fits_data

        final_hdu_list = [fits_data[0].copy()]

        for i in xrange(1, len(fits_data)):
            current_amplifier = fits_data[i]
            detector = gmos_raw_object.instrument_setup.detectors[i - 1]
            if detector.readout_direction.lower() == 'left':
                isleftamp = True
            elif detector.readout_direction.lower() == 'right':
                isleftamp = False
            else:
                raise ValueError('readout direction is unknown string: %s' % detector.readout_direction.lower())

            amplifier_data = prepare.correct_overscan(current_amplifier, 
                                                      self.bias_subslice, 
                                                      self.data_subslice,
                                                      self.bias_clip_sigma)
            amplifier_data = prepare.correct_gain(amplifier_data, gain=detector.gain)
            amplifier_data.name = 'DATA_%d' % i
            bias_uncertainty = amplifier_data.header['BIASSTD']
            # ***TODO*** bias uncertainty should not be here, since this is correlated 
            # for all pixels.  What is important later is a flat field/sensitivity error.
            amplifier_uncertainty = create_uncertainty_frame(amplifier_data, readout_noise=detector.readout_noise,
                                                             bias_uncertainty=bias_uncertainty)

            amplifier_uncertainty.name = 'UNCERTAINTY_%d' % i

            amplifier_mask = create_mask(amplifier_data)
            amplifier_mask.name = 'MASK_%d' %i

            final_hdu_list += [amplifier_data, amplifier_uncertainty, amplifier_mask]

        fits.HDUList(final_hdu_list).writeto(full_path, clobber=True)
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

    mask_data = mask_data.astype(dtype=np.unit8)

    return fits.ImageHDU(mask_data, header=amplifier.header)

def gmos_ccd_image_arithmetic(func, *args):
    pass

class GMOSCCDImage(object):

    @classmethod
    def from_fits_object(cls, fits_object):
        #searching for independent datasets:
        fits_data = fits_object.fits_data
        chips = []
        meta = OrderedDict(fits_data[0].header.copy())
        for chip_id in xrange(1, 4):
            data = fits_data['data_%d' % chip_id].data

            if fits_data['uncertainty_%d' % chip_id].header['uncertainty_type'] == 'stddev':
                uncertainty = StdDevUncertainty(fits_data['uncertainty_%d' % chip_id].data)

            mask = fits_data['mask_%d' % chip_id].data.astype(bool)

            chips.append(NDData(data, uncertainty=uncertainty, mask=mask))

        return cls(chips, meta)


    def __init__(self, chips, meta):
        self.chips = chips
        self.meta = meta

        for chip in self.chips:
            chip.meta = meta


    def __add__(self, other):
        if not isinstance(other, GMOSCCDImage):
            raise TypeError
