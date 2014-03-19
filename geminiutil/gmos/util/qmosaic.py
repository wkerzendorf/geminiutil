import astropy.io.fits as fits
import numpy as np


def qmosaic(red, green, blue):
    """Combine three Gemini/GMOS extensions into a single frame

    Offsets are based on CRPIX in the headers, rounded to the nearest integer
    pixel.  As a result, the WCS is correct for the middle chip, but only
    approximate for the left and right one.

    Parameters
    ----------
    red, green, blue : `~astropy.fits.ImageHDU` extensions
        fits extensions holding the CCD images from the left, middle, and right

    Returns
    -------
    `~astropy.fits.ImageHDU` extension with the combined image

    Notes
    -----
    At present, this does not properly update the headers.  E.g., BIAS, CCDSEC,
    etc., will all reflect the values from the middle chip.
    """
    gap_red_green = np.round(red.header['CRPIX1'] - red.shape[1] -
                             green.header['CRPIX1'])
    gap_green_blue = np.round(green.header['CRPIX1'] - green.shape[1] -
                              blue.header['CRPIX1'])
    combo = np.hstack((red.data,
                       np.zeros((red.shape[0], gap_red_green)),
                       green.data,
                       np.zeros((green.shape[0], gap_green_blue)),
                       blue.data))
    im = fits.ImageHDU(combo, green.header)
    im.header['CRPIX1'] += red.shape[1]+gap_red_green
    im.header.pop('CCDNAME')
    im.header.pop('CCDSIZE')
    return im


def rgb2mosaic(rgb):
    """Combine a Gemini/GMOS fits file with 3 extensions into a single frame

    Offsets are based on CRPIX in the headers, rounded to the nearest integer
    pixel.  As a result, the WCS is correct for the middle chip, but only
    approximate for the left and right one.

    Parameters
    ----------
    rgb : `~astropy.fits.HDUList`
        fits file with 3 extensions holding the CCD images from the
        left, middle, and right

    Returns
    -------
    `~astropy.fits.HDUList` with a single extension holding the combined image

    Notes
    -----
    At present, this does not properly update the headers.  E.g., BIAS, CCDSEC,
    etc., will all reflect the values from the middle chip.
    """
    return fits.HDUList([rgb[0], qmosaic(rgb[1],rgb[2],rgb[3])])
