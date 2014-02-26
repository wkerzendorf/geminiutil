"""Basic reduction class Prepare

Typical usage:
Prep = Prepare()   # can give bias_subslice, data_subslice, clip, gains, ...
infits = fits.open(fitsfile)
outfits = Prep(infits)

***TODO***
Ensure the out fits headers contain everything required to reproduce the result
"""

from geminiutil.gmos.util import prepare_frame#, prepare_slices

import logging
logger = logging.getLogger(__name__)
from sqlalchemy.orm import object_session

class GMOSPrepareFrame(object):  # will be base when we know what
    file_prefix = 'prep'

    def __init__(self, bias_subslice=[slice(None), slice(1,11)],
            data_subslice=None, bias_clip_sigma=3.,
            gain=None, read_noise=None, combine=True,
            overscan_std_threshold=3.):

        self.bias_subslice = bias_subslice
        self.data_subslice = data_subslice
        self.gain = gain
        self.bias_clip_sigma = bias_clip_sigma
        self.read_noise = read_noise
        self.combine = combine
        self.overscan_std_threshold = overscan_std_threshold

    def __call__(self, gmos_raw_object):
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
            fits file name to write out the cut_image to



        """


        fits_data = gmos_raw_object.fits.fits_data

        assert len(fits_data) - 1 == 3

        # subtract overscan, get useful part of detector, correct for gain,
        # and set read noise
        if self.read_noise is None:
            read_noise = [detector.readout_noise for detector in
                          gmos_raw_object.instrument_setup.detectors]
        else:
            read_noise = self.read_noise

        if self.gain is None:
            gain = [detector.gain for detector in
                    gmos_raw_object.instrument_setup.detectors]
        else:
            gain = self.gain


        fits_file = prepare_frame.prepare(fits_data,
                                    bias_subslice=self.bias_subslice,
                                    data_subslice=self.data_subslice,
                                    clip=self.bias_clip_sigma,
                                    gain=gain, read_noise=read_noise, combine=self.combine,
                                    overscan_std_threshold=self.overscan_std_threshold)
        # give each extension a name.  May later add error/mask extensions
        for i, extension in enumerate(fits_file):
            if i > 0:
                extension.name = 'chip{0:d}.data'.format(i)

        return fits_file


