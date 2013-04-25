import astropy.io.fits as fits
import numpy as np
from numpy.polynomial.polynomial import polyfit as polyfit

"""Basic reduction classes SubtractOverscan, CorrectGain, CombineHalves.

Typical usage:
SO = SubtractOverscan()   # can give bias_slice, data_slice
CG = CorrectGain()        # can force gains other than those of headers
infits = fits.open(fitsfile)
outfits = CorrectGain(SubtractOverscan(infits))

***TODO*** 
Add CombineHalves into CorrectGain?  Add all three together?
set up classes so they can used even as functions with default values?
Ensure the out fits headers contain everything required to reproduce the result
"""

class SubtractOverscan(object):
    def __init__(self, bias_slice=[slice(None), slice(1,11)], 
                 data_slice=[slice(None), slice(-1)], clip=3.):
        """Class that extracts bias-corrected exposed parts of GMOS extensions.

        Bias is determined from selected regions of overscan, clipping outliers

        Parameters
        ----------
        bias_slice: list of 2 slices or 2x2 int array
            slices in Y, X within the overscan region given in the 
            extension header 'BIASSEC' to use for determining the bias level
            (default: all Y, 1:11 in X, relative to far end of the overscan) 
        data_slice: list of 2 slices or 2x2 int array
            good parts of the array, relative to header 'DATASEC'
            (default: all Y, 0:-1 in X, i.e., all but first column read-out)
        clip: float
            number of standard deviations from the mean within which overscan
            values have to be to be included.  I.e., bias substracted is mean of
            overscan[where[abs(overscan-overscan.mean())<clip*overscan.std()]]

        Returns
        -------
        Class instance which can be applied to input fits files, using call
        method (which calls function correct_overscan with parameters set here)
        """
        self.bias_slice = self._interpret_slice(bias_slice)
        self.data_slice = self._interpret_slice(data_slice)
        self.clip = clip

    def _interpret_slice(self, in_slice):
        """Turn lists or arrays into list of slices."""
        if in_slice is None: return [slice(None), slice(None)]
        assert len(in_slice) == 2 # need two dimensions
        out_slice = []
        for slc in in_slice:
            if isinstance(slc, slice):
                out_slice += [slc]
            else:
                out_slice += [slice(slc[0] if slc[0] else None, 
                                    slc[1] if slc[1] else None)]
        return out_slice

    def __call__(self, im):
        return correct_overscan(im, self.bias_slice, self.data_slice, self.clip)

class CorrectGain(object):
    def __init__(self, gain=None, ron=None, error_opt=1):
        """Class that corrects images for gain and estimates uncertainties.
        
        Parameters
        ----------
        gain: sequence of floats, optional
            e-/ADU for each amplifiers (default: read from header)
        ron: sequence of floats, optional
            read-out noise (e-) for each amplifiers (default: read from header)
        error_opt: integer, optional
            whether to estimate uncertainties, and how to store them (default 1)
         None: do not estimate uncertainties
            1: store uncertainties as an additional plane in each extension
            2: store uncertainties as extensions with names En

        Notes
        -----
        Uncertainties are estimated using sqrt(abs(counts) + ron*ron),
        where counts=ADU*gain
        """
        self.gain = gain
        self.ron = ron
        self.error_opt = 0 if error_opt is None else error_opt

    def __call__(self, im):
        return correct_gain(im, gain=self.gain, readout_noise=self.ron,
                            error_opt=self.error_opt)

class CombineHalves(object):
    """Combine detector halves into single frames.

    ***TODO*** Maybe better just as a function? Combine with gain?
    """
    def __call__(self, im):
        return combine_halves(im)

def correct_overscan(amplifier, isleftamp, bias_slice=[slice(None), slice(1,11)],
                 data_slice=[slice(None), slice(-1)], sigma_clip=3.):
    """Extract bias-corrected, exposed parts of raw GMOS fits file

    Bias is determined from selected regions of overscan, clipping outliers

    Parameters
    ----------
    amplifier: `fits.ImageHDU` object

    bias_slice: list of 2 slices
        slices in Y, X within the overscan region given in the 
        extension header 'BIASSEC' to use for determining the bias level
        (default: all Y, 1:11 in X, relative to far end of the overscan) 
    data_slice: list of 2 slices
        good parts of the array, relative to header 'DATASEC'
        (default: all Y, 0:-1 in X, i.e., all but first column read-out)
    clip: float
        number of standard deviations from the mean within which overscan
        values have to be to be included.  I.e., bias substracted is mean of
        overscan[where[abs(overscan-overscan.mean())<clip*overscan.std()]]

    Returns
    -------
    output fits structure
        Same extensions as the input one, with bias-corrected exposed parts
     
    Note
    ----
    The output fits structure will have the following new headers:
    DATAUSED: section of initial data array extracted
    BIASUSED: section of initial data array used to determine bias level
    BIAS:     bias value subtracted (ADU)
    BIASSTD:  standard deviation of values going in to bias value (ADU)
              (should be similar to read noise; this is not checked here)

    The following headers are updated
    CRPIX1:   X reference pixel for the world coordinate system
    CRPIX1:   Y reference pixel for the world coordinate system
    CCDSEC:   the part of the chip represented by the data section
    DETSEC:   the part of the overall array represented by the data section

    ***TODO*** update history in primary extension
    """


    #    isleftamp = 'left' in amp.header['AMPNAME']

    bias_slice = sub_slices(sec2slice(amplifier.header['BIASSEC']),
                            reverse_yslice(bias_slice, doreverse=isleftamp))
    data_slice = sub_slices(sec2slice(amplifier.header['DATASEC']),
                            reverse_yslice(data_slice, doreverse=isleftamp))

    overscan = amplifier.data[bias_slice]

    if sigma_clip is None:
        clipped = overscan
    else:
        clipped = overscan[np.abs(overscan-overscan.mean()) <
                                    sigma_clip*overscan.std()]

    bias_estimate, bias_std = clipped.mean(), clipped.std()
    outamp = fits.ImageHDU(amplifier.data[data_slice]-bias_estimate,
                           amplifier.header)
    outamp.header['CRPIX1'] += data_slice[1].start
    outamp.header['CRPIX2'] += data_slice[0].start
    binning = np.fromstring(amplifier.header['CCDSUM'], sep=' ', dtype=np.int)
    outamp.header['CCDSEC'] = slice2sec(sub_slices(
        sec2slice(amplifier.header['CCDSEC']),
        reverse_yslice(data_slice, doreverse=isleftamp), factors=binning))
    outamp.header['DETSEC'] = slice2sec(sub_slices(
        sec2slice(amplifier.header['DETSEC']),
        reverse_yslice(data_slice, doreverse=isleftamp), factors=binning))
    outamp.header.pop('DATASEC') # all of image is now DATA
    outamp.header.pop('BIASSEC') # no BIAS section left
    outamp.header['DATAUSED'] = slice2sec(data_slice)
    outamp.header['BIASUSED'] = slice2sec(bias_slice)
    outamp.header['BIAS'] = bias_estimate
    outamp.header['BIASSTD'] = bias_std

    return outamp

def correct_gain(amplifier, gain=None, readout_noise=None):
    """Correct image for gain and possibly estimate uncertainties.
    
    Parameters
    ----------
    amplifier: `~astropy.io.fits.ImageHDU`
    gain: sequence of floats, optional
       e-/ADU for each amplifiers (default: read from header)
    ron: sequence of floats, optional
       read-out noise (e-) for each amplifiers (default: read from header)
    error_opt: integer, optional
       whether to estimate uncertainties, and how to store them (default 1)

    Notes
    -----
    Uncertainties are estimated using sqrt(abs(counts) + ron*ron),
    where counts=ADU*gain

    In the output fits structure, the following headers are updated/created:
    GAIN:     set to 1 in output
    GAINUSED: gain values were multiplied with
    RDNOISE:  read-out noise to be used for estimating uncertainties
    """


    if gain is None:
        gain = amplifier.header['GAIN']

    if readout_noise is None:
        readout_noise = amplifier.header['RDNOISE']

    corrected_amplifier_data = amplifier.data*gain

    corrected_amplifier = fits.ImageHDU(corrected_amplifier_data, amplifier.header.copy())
    corrected_amplifier.header['GAIN'] = 1.
    corrected_amplifier.header['GAINUSED'] = gain
    corrected_amplifier.header['RDNOISE'] = readout_noise

    return corrected_amplifier


def create_uncertainty_frame(amplifier, readout_noise=None, bias_uncertainty=None):
    """
    Create an standard deviation uncertainty frame for GMOS

    Parameters
    ----------

    amplifier: ~astropy.io.fits.ImageHDU

    readout_noise: ~float
        readout_noise in e-/ADU !!! TODO really in e-/ADU or e-?? !!!!
        if None will be read from the header 'RDNOISE'
    bias_uncertainty: ~float
        uncertainty from the overscan subtract
        if None will be read from the header 'BIASSTD'
    """

    if readout_noise is None:
        readout_noise = amplifier.header['RDNOISE']

    if bias_uncertainty is None:
        bias_uncertainty = amplifier.header['BIASSTD']

    uncertainty_data = np.sqrt(np.abs(amplifier.data) + readout_noise**2 + bias_uncertainty**2)

    uncertainty = fits.ImageHDU(uncertainty_data, amplifier.header.copy())
    uncertainty.header['uncertainty_type'] = 'stddev'

    return uncertainty



def combine_halves(im):
    outlist = [im[0]]
    if 'E' in im[2].name:
        pairlist = zip(im[1::4], im[3::4])+zip(im[2::4], im[4::4])
        # re-order??
    else:
        pairlist = zip(im[1::2], im[2::2])
    for amp in pairlist:
        left, right = ((amp[0], amp[1]) 
                       if 'left' in amp[0].header['AMPNAME']
                       else (amp[1], amp[0]))
        if left.data.ndim == 2: # no error planes
            combo = np.hstack([left.data, right.data])
        else:
            combo = np.dstack([left.data, right.data]).copy()
        chip = fits.ImageHDU(combo, left.header)
        # remove amplifier specific keys (maybe want to keep??)
        for key in ('AMPNAME', 'GAIN', 'RDNOISE', 'GAINUSED', 
                    'CCDSEC', 'DATAUSED', 'BIASUSED', 
                    'FRMNAME', 'FRAMEID'):
            chip.header.remove(key)
        ldet = sec2slice(left.header['DETSEC'])
        rdet = sec2slice(right.header['DETSEC'])
        chip.header['DETSEC'] = slice2sec((ldet[0], slice(ldet[1].start,
                                                          rdet[1].stop)))
        outlist += [chip]
    return fits.HDUList(outlist)


def reverse_yslice(in_slice, doreverse=True):
    """Reverse second dimension of in_slice (if doreverse=True)"""
    return [in_slice[0], 
            reverse_slice(in_slice[1]) if doreverse else in_slice[1]]

def reverse_slice(in_slice):
    """Reverse a slice, such that a[out]==a[::-1][in_slice][::-1]"""
    return slice(None if in_slice.stop is None else -in_slice.stop, 
                 None if in_slice.start is None else -in_slice.start)

def sub_slices(slices1, subslices, factors=None):
    """Get subslices of slices, such that a[out]==a[slice1][subslices]
    Both slices1 and subslices are lists of slices."""
    if factors is None: factors = [1]*len(slices1)
    return [sub_slice(slice1, subslice, factor) 
            for slice1,subslice,factor in zip(slices1, subslices, factors)]

def sub_slice(slice1, subslice, factor=1):
    """Get a subslice of a slice, such that a[out]==a[slice1][subslices]"""
    return slice(slice_part(slice1, subslice.start, factor, start=True),
                 slice_part(slice1, subslice.stop, factor, start=False))

def slice_part(slice1, pos, factor, start):
    if pos is None:
        return slice1.start if start else slice1.stop
    if pos >= 0:
        if slice1.start is None:
            return pos
        else:
            return slice1.start+pos*factor
    else:
        if slice1.stop is None:
            return pos
        else:
            return slice1.stop+pos*factor

def sec2slice(sec):
    """Convert a section in IRAF format to a list of slices.

    E.g., [1:10,5:20] becomes [slice(4,19),slice(0,9)]
    """
    secnums = np.fromstring(sec.strip('[]').replace(':',','),np.int,sep=',')-1
    return [slice(secnums[2],secnums[3]+1), slice(secnums[0],secnums[1]+1)]

def slice2sec(slice):
    """Convert (list of) python slice(s) to a string in IRAF format

    E.g., [slice(4,19),slice(0,9)] becomes[1:10,5:20]
    """
    sec = '['
    for i in range(len(slice)-1,-1,-1):
        sec = sec + repr(slice[i].start+1) + ':' + repr(slice[i].stop)
        if i>0: sec = sec + ','
    return sec + ']'
