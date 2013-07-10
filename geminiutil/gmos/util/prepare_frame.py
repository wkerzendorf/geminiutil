import astropy.io.fits as fits
import astropy.stats.funcs as stats
import numpy as np
import warnings

warnings.filterwarnings('ignore',
                        message='.+ a HIERARCH card will be created.')


def prepare(image, bias_subslice=[slice(None), slice(1,11)],
            data_subslice=[slice(None), slice(-1)], clip=3.,
            gain=None, read_noise=None, combine=True,
            overscan_std_threshold=3.):
    """Extract and combine bias- and gain-corrected extensions.

    Bias is determined from selected regions of overscan, clipping outliers.
    Then, active region have bias subtracted, and are multiplied by gains.
    If needed, detector halves are combined.

    Parameters
    ----------
    image: fits.HDUList
        input GMOS image with extensions for each amplifier
    bias_subslice: list of 2 slices
        slices in Y, X within the overscan region given in the
        extension header 'BIASSEC' to use for determining the bias level
        (default: all Y, 1:11 in X, relative to far end of the overscan;
         None: no subslice, i.e., use header)
    data_subslice: list of 2 slices
        good parts of the array, relative to header 'DATASEC'
        (default: all Y, 0:-1 in X, i.e., all but first column read out;
         None: no subslice, i.e., use header)
    clip: ~float
        number of standard deviations from the median within which overscan
        values have to be to be included.
    gain: sequence of ~float, optional
        e-/ADU for each amplifiers (default: header 'GAIN')
    read_noise: sequence of ~float, optional
        read-out noise (e-) for each amplifier (default: header 'RDNOISE')
    combine: ~bool, optional
        Combine halves if detectors were read out through two amplifiers
    overscan_std_threshold: ~float, optional
        Warn if standard deviation of overscan is larger than
        this value times the read noise (default: 3)

    Returns
    -------
    ~fits.HDUList
    """
    outlist = [image[0]]
    gain = gain or [None]*(len(image)-1)
    read_noise = read_noise or [None]*(len(image)-1)
    for amp, ampgain, ampron in zip(image[1:], gain, read_noise):
        isleftamp = 'left' in amp.header['AMPNAME']
        bias_slice = slice_slices(sec2slice(amp.header['BIASSEC']),
                                  adjust_subslices(bias_subslice,
                                                   reverse1=isleftamp))
        data_slice = slice_slices(sec2slice(amp.header['DATASEC']),
                                  adjust_subslices(data_subslice,
                                                   reverse1=isleftamp))
        bias_corrected = correct_overscan(amp, bias_slice, data_slice, clip)
        gain_corrected = correct_gain(bias_corrected, ampgain)
        if ampron:
            gain_corrected.header['RDNOISE'] = ampron
        if overscan_std_threshold:
            if(gain_corrected.header['BIASSTD'] *
               gain_corrected.header['APPLIED GAIN'] >
               overscan_std_threshold * gain_corrected.header['RDNOISE']):
                warnings.warn('Overscan standard deviation of {} '
                              '(gain={}; ron={}; amp={})'.format(
                                  gain_corrected.header['BIASSTD'],
                                  gain_corrected.header['APPLIED GAIN'],
                                  gain_corrected.header['RDNOISE'],
                                  gain_corrected.header['AMPNAME']))
        outlist += [gain_corrected]

    if combine and image[0].header['NAMPS'] == 2:
        outlist = [outlist[0]] + [combine_halves(amp1, amp2) for amp1, amp2
                                  in zip(outlist[1::2], outlist[2::2])]
    return fits.HDUList(outlist)


def correct_overscan(amplifier, bias_slice=None, data_slice=None, clip=3.):
    """Extract bias-corrected, exposed parts of raw GMOS fits file

    Bias is determined from selected regions of overscan, clipping outliers

    Parameters
    ----------
    amplifier: ~fits.ImageHDU
    bias_slice: list of 2 slices
        slices in Y, X to use for determining the bias level
        (default=None: use 'BIASSEC' from header)
    data_slice: list of 2 slices
        slices in Y, X that give good parts of the array
        (default=None: use 'DATASEC' from header)
    clip: ~float
        number of standard deviations from the median within which overscan
        values have to be to be included.

    Returns
    -------
    ~fits.ImageHDU
        bias-corrected exposed parts

    Note
    ----
    The output fits structure will have the following new headers:
    APPLIED DATASEC: section of initial data array extracted
    APPLIED BIASSEC: section of initial data array used to determine bias level
    BIAS:     bias value subtracted (ADU)
    BIASSTD:  standard deviation of values going in to bias value (ADU)
              (should be similar to read noise; this is not checked here)
    BIASUNC:  formal uncertainty in bias (BIASSTD/sqrt(N), in ADU)

    The following headers are updated
    CRPIX1:   X reference pixel for the world coordinate system
    CRPIX1:   Y reference pixel for the world coordinate system
    CCDSEC:   the part of the chip represented by the data section
    DETSEC:   the part of the overall array represented by the data section
    """

    if bias_slice is None:
        bias_slice = sec2slice(amplifier.header['BIASSEC'])

    if data_slice is None:
        data_slice = sec2slice(amplifier.header['DATASEC'])
        data_subslice = None
    else:
        # ensure we get numbers only, rather than relative indices,
        # otherwise cannot store it as DATASEC header later
        data_slice = [slice(*s.indices(shape))
                      for s, shape in zip(data_slice, amplifier.data.shape)]
        # get equivalent subslice of 'DATASEC' for updating CCDSEC, DETSEC
        data_subslice = get_subslices(sec2slice(amplifier.header['DATASEC']),
                                      data_slice)

    overscan = amplifier.data[bias_slice]

    if clip:
        clipped = stats.sigma_clip(overscan, clip, 1)
    else:
        clipped = overscan

    bias_estimate, bias_std = clipped.mean(), clipped.std()
    bias_unc = bias_std/np.sqrt(clipped.size)

    outamp = fits.ImageHDU(amplifier.data[data_slice]-bias_estimate,
                           amplifier.header)

    # add new headers with information on what was done
    outamp.header['APPLIED DATASEC'] \
        = slice2sec(data_slice), 'Section considered good'
    outamp.header['APPLIED BIASSEC'] \
        = slice2sec(bias_slice), 'Section used to estimate bias'
    outamp.header['BIAS'] = bias_estimate, 'ADU'
    outamp.header['BIASSTD'] = bias_std, 'ADU'
    outamp.header['BIASUNC'] = bias_unc, 'ADU'
    # adjust remove/existing headers as appropriate
    outamp.header.remove('BIASSEC')  # no BIAS section left
    outamp.header.remove('DATASEC')  # all of image is now DATA
    outamp.header['CRPIX1'] += data_slice[1].start
    outamp.header['CRPIX2'] += data_slice[0].start
    if data_subslice:
        binning = np.fromstring(amplifier.header['CCDSUM'], sep=' ',
                                dtype=np.int)
        ccd_subslice = adjust_subslices(data_subslice, binning)
        outamp.header['CCDSEC'] = slice2sec(slice_slices(
            sec2slice(amplifier.header['CCDSEC']), ccd_subslice))
        outamp.header['DETSEC'] = slice2sec(slice_slices(
            sec2slice(amplifier.header['DETSEC']), ccd_subslice))

    return outamp


def correct_gain(amplifier, gain=None):
    """Correct image for gain

    Parameters
    ----------
    amplifier: ~fits.ImageHDU
    gain: ~float, optional
       e-/ADU; default (~None): read from amplifier.header['GAIN']

    Returns
    -------
    ~fits.ImageHDU

    Notes
    -----
    In the output fits image, the following headers are updated/created:
    GAIN:      set to 1
    APPLIED GAIN:  gain values were multiplied with
    """

    if gain is None:
        gain = amplifier.header['GAIN']

    corrected = fits.ImageHDU(amplifier.data*gain, amplifier.header)
    corrected.header['GAIN'] = 1., 'Unity since data corrected for gain'
    corrected.header['APPLIED GAIN'] = gain, 'Gain multiplied with'

    return corrected


def combine_halves(amp1, amp2):
    """Combine two detector halves read out through 2 amplifiers into one.

    Parameters
    ----------
    amp1, amp2: ~fits.ImageHDU

    Returns
    -------
    ~fits.ImageHDU

    Notes
    -----
    Updates headers:
      CCDSEC, DETSEC, DATATYP
    Renames amplifier specific one by prefixing HIERARCH LEFT, HIERARCH RIGHT:
      AMPNAME, FRMNAME, FRAMEID, APPLIED DATASEC, APPLIED BIASSEC
      BIAS, BIASSTD, BIASUNC, GAIN, APPLIED GAIN, RDNOISE
    Also creates HIERARCH LEFT DATASEC and HIERARCH RIGHT DATASEC
    """
    left, right = ((amp1, amp2) if 'left' in amp1.header['AMPNAME']
                   else (amp2, amp1))
    assert 'left' in left.header['AMPNAME'] and \
        'right' in right.header['AMPNAME']

    chip = fits.ImageHDU(np.hstack([left.data, right.data]), left.header)

    # rename amplifier specific keys to left/right pairs
    for key in ('AMPNAME', 'FRMNAME', 'FRAMEID', 'GAIN', 'RDNOISE',
                'APPLIED BIASSEC', 'BIAS', 'BIASSTD', 'BIASUNC',
                'CCDSEC', 'DETSEC', 'APPLIED DATASEC', 'APPLIED GAIN'):
        if key in ('CCDSEC', 'DETSEC'):
            ls, rs = sec2slice(left.header[key]), sec2slice(right.header[key])
            if ls[0] == rs[0] and ls[1].stop == rs[1].start:
                chip.header[key] = slice2sec([ls[0],
                                              slice(ls[1].start, rs[1].stop)])
                continue
        chip.header['LEFT {}'.format(key)] = chip.header.pop(key)
        chip.header['RIGHT {}'.format(key)] = right.header[key]
    # and write the new data sections
    ls, rs = left.data.shape, right.data.shape
    chip.header['LEFT DATASEC'] \
        = slice2sec([slice(0, ls[0]), slice(0, ls[1])]), 'Part from left amp.'
    chip.header['RIGHT DATASEC'] \
        = slice2sec([slice(0, rs[0]),
                     slice(ls[1], ls[1]+rs[1])]), 'Part from right amp.'
    return chip


def adjust_subslices(subslices, factors=None, reverse1=False):
    """Scale subslices and possibly reverse horizontal direction."""
    if subslices is not None:
        if factors is not None:
            subslices = [multiply_slice(subslice, factor)
                         for subslice, factor in zip(subslices, factors)]
        if reverse1:
            subslices = [subslices[0], reverse_slice(subslices[1])]
    return subslices


def multiply_slice(in_slice, factor):
    """Multiply a slice with a given factor (ignoring step)."""
    return slice(None if in_slice.start is None else factor*in_slice.start,
                 None if in_slice.stop is None else factor*in_slice.stop)


def reverse_slice(in_slice):
    """Reverse a slice, such that a[out]==a[::-1][in_slice][::-1]"""
    return slice(None if in_slice.stop is None else -in_slice.stop,
                 None if in_slice.start is None else -in_slice.start)


def get_subslices(slices1, slices2):
    """Get subslices of slices1 such that a[slices1][out]==a[slices2]"""
    if slices2 is None:
        return slices1
    else:
        return [get_subslice(slice1, slice2)
                for slice1, slice2 in zip(slices1, slices2)]


def get_subslice(slice1, slice2):
    assert slice1.start <= slice2.start and slice1.stop >= slice2.stop
    return slice(None if slice2.start == slice1.start
                 else slice2.start-slice1.start,
                 None if slice2.stop == slice1.stop
                 else slice2.stop-slice1.stop)


def slice_slices(slices1, subslices):
    """Get subslices of slices, such that a[out]==a[slice1][subslices]
    Both slices1 and subslices are lists of slices."""
    if subslices is None:
        return slices1
    else:
        return [slice_slice(slice1, subslice)
                for slice1, subslice in zip(slices1, subslices)]


def slice_slice(slice1, subslice):
    """Get a subslice of a slice, such that a[out]==a[slice1][subslice]

    Parameters
    ----------
    slice1: ~slice
        Must have numbers for all start, stop (e.g., taken from a section)
    subslice: ~slice
        Can have None for parts that stay the same

    Returns
    -------
    ~slice, containing subslice of slice1
    """
    return slice(slice_part(slice1, subslice.start, start=True),
                 slice_part(slice1, subslice.stop, start=False))


def slice_part(slice1, pos, start):
    """Helper routine for ~sub_slice"""
    if pos is None:
        return slice1.start if start else slice1.stop
    if pos >= 0:
        return slice1.start+pos if slice1.start else pos
    else:
        return slice1.stop+pos if slice1.stop else pos


def sec2slice(sec):
    """Convert a section in IRAF format to a list of two slices.

    E.g., [1:10,5:20] becomes [slice(4,19),slice(0,9)]
    """
    secnums = np.fromstring(sec.strip('[]').replace(':',','),np.int,sep=',')-1
    return [slice(secnums[2],secnums[3]+1), slice(secnums[0],secnums[1]+1)]


def slice2sec(slices):
    """Convert (list of) python slice(s) to a string in IRAF format

    E.g., [slice(4,19),slice(0,9)] becomes[1:10,5:20]

    Note: slices have to contain numbers for start, stop, and step is ignored
    """
    sec = ','.join(['{0}:{1}'.format(slc.start+1, slc.stop)
                    for slc in reversed(slices)])
    return '[' + sec + ']'


def header_value(chip, key):
    """Get header value, combining information from two amplifiers if needed.

    For a given key, if header exists, return it.
    If not, retrieve 'LEFT '+key, 'RIGHT '+key, and return 1-dimensional array
    with left and right parts set following LEFT DATASEC, RIGHT DATASEC regions

    Parameters
    ----------
    chip: ~fits.ImageHDU
        holding fits headers for chip
    key: ~string
        header requested

    Returns
    -------
    1-dimensional array of header content (1 element) or left+right array
    """
    if key in chip.header:
        return np.array([chip.header[key]])

    values = leftright_header(chip, key)
    slices = [sec2slice(datasec) for datasec in
              leftright_header(chip, 'DATASEC')]
    assert slices[0][0] == slices[1][0] and slices[0][1] != slices[1][1]
    out = np.zeros(chip.data.shape[1])
    out[slices[0][1]], out[slices[1][1]] = values
    return out


def leftright_header(chip, key):
    """Get tuple with chip.header['LEFT '+key] & chip.header['RIGHT '+key]."""
    return (chip.header['LEFT {}'.format(key)],
            chip.header['RIGHT {}'.format(key)])


def multiext_data(im):
    """Return multiple extensions as (presumably) 3-dimensional array."""
    return np.array([ext.data for ext in im[1:]])


def multiext_header_value(im, key):
    """Return headers from multiple extensions as 3-dimensional array."""
    return np.array([header_value(ext, key)[np.newaxis,:] for ext in im[1:]])


def x_coords(chip, key='CCDSEC'):
    assert key in ['CCDSEC', 'DETSEC']
    slice_x = sec2slice(chip.header[key])[1]
    step = np.fromstring(chip.header['CCDSUM'], sep=' ', dtype=np.int)[1]
    return np.arange(slice_x.start-0.5+step/2., slice_x.stop-0.5+step/2., step)


def multiext_x_coords(im, key='CCDSEC'):
    out = np.array([x_coords(ext, key)[np.newaxis,:] for ext in im[1:]])
    if all(np.all(out[0] == t) for t in out[1:]):
        out = out[0].reshape(1, 1, -1)
    return out



