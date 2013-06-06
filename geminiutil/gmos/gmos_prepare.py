"""Basic reduction class Prepare

Typical usage:
Prep = Prepare()   # can give bias_subslice, data_subslice, clip, gains, ...
infits = fits.open(fitsfile)
outfits = Prep(infits)

***TODO***
Ensure the out fits headers contain everything required to reproduce the result
"""

from .basic import prepare


class Prepare(object):
    def __init__(self, bias_subslice=[slice(None), slice(1,11)],
                 data_subslice=[slice(None), slice(-1)], clip=3.,
                 gain=None, read_noise=None, combine=True,
                 overscan_std_threshold=3.):
        """Class that extracts and combines bias- and gain-corrected extensions

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
