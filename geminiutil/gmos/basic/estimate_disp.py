# -*- coding: utf-8 -*-
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.special import erf
from astropy.table import Table


def smoothed_slit(x, cx, wx, sx):
    """Return a slit profile, for given slit width and instrumental broadening.

    Treats a slit image as a convolution of a box car and a Gaussian.

    Parameters
    ----------
    x: float array
        coordinates on detector
    cx, wc, sx: float (array)
        position of slit, slit width, and sigma of gaussian (in pixel units)
    """
    return (erf((x-cx+wx/2.)/(sx*np.sqrt(2.))) -
            erf((x-cx-wx/2.)/(sx*np.sqrt(2.))))/2.


def fake_arc(grid, wavelengths, amplitudes, width, sigma):
    """Make a fake arc spectrum

    Parameters
    ----------
    grid: float array
        wavelength grid
    wavelengths, amplitudes: float array
        wavelengths and amplitudes of fake lines to be added
        (use amplitudes=None if no line strengths are required)
    width, sigma : float
        slit width and sigma of other instrumental broadening

    Returns
    -------
    Table with columns 'w' holding the input grid and 'f' holding the fake arc
    """
    within_grid = np.logical_and(wavelengths > grid.min(),
                                 wavelengths < grid.max())
    g2d = smoothed_slit(grid, wavelengths[within_grid][:, np.newaxis],
                        width, sigma)
    if amplitudes is not None:
        g2d *= amplitudes[within_grid]

    return Table([grid, g2d.sum(axis=0)], names=['w', 'f'])


def fit_peak1d(x, y, peakloc, full=False):
    """Interpolate z = a+b*x+c*x^2 and return peak location.

    Constructs interpolating parabolae around given peak locations, using
    itself and its two neighbour.

    Parameters
    ----------
    x, y: float array
        x coordinates, and corresponding values
    peakloc: int or int array
        location of the peak(s) in the array (should hold highest or lowest
        value in among its neighbours)
    full: bool
        whether to return just location, or also peak value and the parameters
        of the interpolating parabolas (default: False)

    Returns
    -------
    x_peak: float or float array (same shape as peakloc)
        locations of the peak, as determined from the interpolating parabola
    y_peak: float or float array
        peak value expected from interpolating parabolas (if full=True)
    [a, b, c]: float list/array
        parameters of the interpolating parabolas (if full=True)
    """
    x0, y0 = x[peakloc], y[peakloc]
    xm, ym = x[peakloc-1]-x0, y[peakloc-1]-y0
    xp, yp = x[peakloc+1]-x0, y[peakloc+1]-y0
    #parabola y = ax^2+bx+c
    a = y0
    b = (yp*xm**2-ym*xp**2)/(xp*xm**2-xm*xp**2)
    c = (ym*xp-yp*xm)/(xp*xm**2-xm*xp**2)
    # location of maximum -b/2a, peak counts
    offset = -b/2./c
    x_peak, y_peak = x0+offset, a+b*offset+c*offset**2
    if full:
        return x_peak, y_peak, [a, b, c]
    else:
        return x_peak


def fit_peak2d(x, y, z, peakloc, full=False):
    """Interpolate z = a+b*x+c*y+d*x^2+e*xy+f*y^2 and return peak location.

    Constructs interpolating parabolae around given peak locations, using
    the six points given by itself, its direct neighbours, and the closest
    corner pixel.

    Parameters
    ----------
    x, y, z: float array
        x, y coordinates, and corresponding values
    peakloc: int list/array of dimension (2, n)
        location of the peak(s) in the array (should hold highest or lowest
        value in among their neighbours
    full: bool
        whether to return just location, or also peak value and the parameters
        of the interpolating parabolas (default: False)

    Returns
    -------
    x_peak, y_peak: float list/array (same shape as peakloc)
        locations of the peak, as determined from the interpolating parabola
    z_peak: float list/array
        peak value expected from interpolating parabolas (if full=True)
    [a, b, c, d, e, f]: float list/array
        parameters of the interpolating parabolas (if full=True)
    """
    xpx, zpx, [a, b, d] = fit_peak1d(x, z[:, peakloc[1]], peakloc[0], True)
    ypy, zpy, [a, c, f] = fit_peak1d(y, z[peakloc[0], :], peakloc[1], True)
    ix = peakloc[0] + np.where(-b/2./d < 0., -1, 1)
    iy = peakloc[1] + np.where(-c/2./f < 0., -1, 1)
    xx, yx = x[ix]-x[peakloc[0]], y[iy]-y[peakloc[1]]
    e = (z[ix,iy] - a - xx*(b+d*xx) - yx*(c+f*yx))/(xx*yx)
    denom = 4*d*f-e*e
    offset_x = (c*e-2.*b*f)/denom
    offset_y = (b*e-2.*c*d)/denom
    z_peak = a + offset_x*(b+offset_x*d+offset_y*e) + offset_y*(c+offset_y*f)
    x_peak = x[peakloc[0]] + offset_x
    y_peak = y[peakloc[1]] + offset_y
    if full:
        return x_peak, y_peak, z_peak, [a, b, c, d, e, f]
    else:
        return x_peak, y_peak


def estimate_disp(arc, wref, xgrid, dispgrid, comparison=None,
                  line_catalog=None, slit_width=7., sigma=1., full=False):
    """Estimate location of reference wavelength and the dispersion of an arc

    Iterates over a set of location and dispersion estimates, matching the
    model to the input at each, and determines the location of the best match
    by interpolation among the point with the best chi2 and its neighbours.

    Parameters
    ----------
    arc : float array or Table
        Arc spectrum; if given as table, should hold column 'f'; if a column
        'x' is present, those will be used; otherwise x=arange(len(arc))
    wref : float
        Reference wavelength for which location in arc is to be found
    xgrid : array
        Trial values within which reference wavelength will lie
    dispgrid : array
        Trial dispersions within which the real dispersion will lie
    comparison : Table or None
        if given, a Table with the model to match the arc to, containing
        columns 'w' and 'f'
    line_catalog : array or None
        wavelength catalog from which a model arc will be constructed
    slit_width, sigma : float
        slit width and sigma of other instrumental broadening, for
        constructing a fake arc model from the line catalog.
    full : bool
        whether to return diagnostic output

    Returns
    -------
    xrefbest, dwbest : float
        estimates of location of reference wavelength, dispersion
    shift : Table
        Match quality for grid.  Contains columns 'xref', 'chi2', 'p0', 'p1',
        with the trial positions, match qualities, and polynomial fit between
        input and model (input = p0 + p1*model); 'chi2', 'p0', and 'p1' have
        dimension equal to the number of dispersion positions tried, with
        the corresponding dispersions stored in shift.meta['dw']
        Only returned if Full=True
    model : Table
        Fake arc constructed from the wavelength catalog, with columns
        'w' and 'f'.  Only returned if Full=True, and comparison=None
    """
    try:
        data = arc['f']
    except:
        data = arc

    try:
        x = arc['x']
    except:
        x = np.arange(len(arc))

    if line_catalog is None:
        model = comparison['w', 'f']
    else:
        assert comparison is None
        xguess = xgrid.mean()
        xsize = xgrid.max()-xgrid.min()
        dispguess = dispgrid.mean()
        model = fake_arc(wref + dispguess*(np.arange(x.min()-xsize*1.2,
                                                     x.max()+xsize*1.2,
                                                     .5) - xguess),
                         line_catalog, None,
                         slit_width*abs(dispguess), sigma*abs(dispguess))
    model.sort('w')

    shift = Table([xgrid] +
                  [np.zeros((xgrid.size, dispgrid.size))]*3,
                  names=('xref','chi2','p0','p1'))
    shift.meta['dw'] = dispgrid
    for j, dw in enumerate(shift.meta['dw']):
        for i, xref in enumerate(shift['xref']):
            pfit,extra = poly.polyfit(data, np.interp(wref+dw*(x-xref),
                                                      model['w'], model['f']),
                                      1, full=True)
            shift['p0'][i,j], shift['p1'][i,j] = pfit
            shift['chi2'][i,j] = np.sqrt(extra[0])
    imin, jmin = np.unravel_index(shift['chi2'].data.argmin(),
                                  shift['chi2'].data.shape)
    imin = max(min(imin, shift['chi2'].data.shape[0]-2), 1)
    jmin = max(min(jmin, shift['chi2'].data.shape[1]-2), 1)

    xrefbest, dwbest = fit_peak2d(shift['xref'], shift.meta['dw'],
                                  shift['chi2'], (imin, jmin))
    if full:
        if line_catalog is not None:
            return xrefbest, dwbest, shift, model
        else:
            return xrefbest, dwbest, shift
    else:
        return xrefbest, dwbest
