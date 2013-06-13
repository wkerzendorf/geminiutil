"""Wavelength calibration.

Typical usage
-------------
Find peaks in a 1-dimensional arc spectrum
    import wavecal
    pos, peak = wavecal.search(arc, min_curvature=10., sigma=1.)

Then make a line table from those (or do it directly)
    from wavecal import LineTable
    linetab = LineTable(pos, peak)
    linetab = LineTable.fromsearch(arc, min_curvature=10., sigma=1.)

And calibrate it
    solguess = [1024., -0.458, 0.]    # xref-wavelength, dispersion, quadratic
    linetab.calibrate(solguess,
                      [(3.,1.5e4), (0.5,2e3), (0.13,0.), (0.06,0.)],
                      lab_wavelengths, refwave=4300., doplot=True)

GMOS example
------------
Assumes arctab is Table with extracted spectrum in ['x'], ['f'][:,(0,1,2)]
    lines = []
    for i, min_curvature in enumerate([10.,10.,7.]):
        lines.append(wavecal.LineTable.fromsearch(arctab['f'][:,i],
                                                  x=arctab['x'],
                                                  min_curvature=min_curvature,
                                                  sigma=1.))
    chip1solguess = [1019., -0.458, 0.]  # x(refwave), dispersion, quadratic
    lines[1].calibrate(chip1solguess,
                       [(3.,1.5e4), (0.5,2e3), (0.13,0.), (0.06,0.)],
                       lincat['w'], refwave=arctab.meta['grwlen']*10,
                       doplot=doplot)
    linesall = wavecal.ThreeChipLineTable.fromlist(lines, refchip=1)
    linesall.calibrate([lines[1].meta['par'][0], 2087.3, 0.54] +
                       list(lines[1].meta['par'][1:]) + [0.],
                       [(1.5,1000), (0.5,200), (0.13,0.), (0.04,0.)],
                       lincat['w'], doplot=doplot)
"""

from __future__ import print_function, division
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from astropy.table import Table
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import leastsq


class LineTable(Table):
    """List of lines identified in an arc spectrum and associated fits.

    Initialization
    --------------
    lintab = LineTable(pos, peak)
        pos, peak can be found using search(arc, ...)
    lintab = LineTable.fromsearch(...)
        initialize from search call directly

    Arc line identification
    -----------------------
    lintab.calibrate(...)
        needs an initial guess that is close, and a laboratory line list.
        See detailed help of calibrate.

    Application
    -----------
    w = lintab.fit(x)

    Fitting function
    ----------------
    For LineTable, the fitting function (_fitfunc) used is

    w = meta['refwave'] + \
        par[1]*(x-par[0])+par[2]*(x-par[0])**2+ ... +par[n]*(x-par[0])**n

    where n is set by the length of the array/list of guesses given when
    the identification is done.

    Save and restore
    ----------------
    lintab.write(...)
        uses Table writer; hdf5 keeps meta-data
    lintab = LineTable.read(...)
        uses Table reader, restores appropriate class from meta['class']
    """

    def __init__(self, *args, **kwargs):
        """Initialise line table using variables.  Details depend on class;
        see cls.init
        """
        table = kwargs.pop('table', None)
        if table is None:
            self.init(*args, **kwargs)
        else:
            assert len(args) == 0
            self._initfromtable(table)
            if kwargs:
                self.meta.update(kwargs)

    @classmethod
    def fromsearch(cls, *args, **kwargs):
        """Initialize LineTable with result of line search.
        See `search` for usage."""
        return cls(*search(*args, **kwargs))

    def init(self, pos, peak, **kwargs):
        """Initialize LineTable with positions and peak fluxes.
        Any other named keywords get added to meta.  For LineTable, this can
        be used to set the reference wavelength, by passing refwave=...
        (maybe more logically done as part of calibrate)
        """
        self._init(pos, peak, names=('x', 'peak'))

    def _init(self, *args, **kwargs):
        """Line table initializer shared between all LineTable classes."""
        names = kwargs.pop('names', None)
        assert names is not None
        assert len(names) == len(args)
        dtypes = kwargs.pop('dtypes', (np.float,)*len(args))
        assert len(dtypes) == len(args)
        Table.__init__(self,
                       args + (np.zeros((args[0].size,)),)*3,
                       names=names+('fit', 'id', 'resid'),
                       dtypes=dtypes+(np.float,)*3,
                       masked=True)
        self.meta['class'] = self.__class__.__name__
        self.meta['id.mask_value'] = self['id'].fill_value
        self.mask['id'] = self.mask['resid'] = True
        self.meta['indep'] = names[:-1]
        if kwargs:
            self.meta.update(kwargs)

    @classmethod
    def read(cls, *args, **kwargs):
        """Read line table from a file using Table reader.
        Expects to find the LineTable class in .meta['class']
        """
        table = Table.read(*args, **kwargs)
        # if table holds recognized line table class, use it
        if 'class' in table.meta:
            class_lookup = {'LineTable': LineTable,
                            'ThreeChipLineTable': ThreeChipLineTable}
            if table.meta['class'] in class_lookup:
                return class_lookup[table.meta['class']](table=table)
        #otherwise, just hope for the best
        return cls(table=table)

    def _initfromtable(self, table):
        """Initializer for table, used by classmethod read"""
        assert {'x', 'peak'}.issubset(table.columns)
        Table.__init__(self, table, masked=True)
        self.mask['id'] = self.mask['resid'] = (self['id'] ==
                                                self.meta['id.mask_value'])

    def match(self, catalog, mindist, minpeak=0.):
        """Find closest matches in catalog for each data['fit'].
        set data['ok'] and set data['id'] to catalog wavelength
        for those within mindist and with data['peak']>minpeak.
        Typically only called from within method calibrate."""
        ix = np.searchsorted(catalog, self['fit'])
        ix = np.clip(ix, 1, len(catalog)-1)
        ix -= self['fit']-catalog[ix-1] < catalog[ix]-self['fit']
        self['id'] = catalog[ix]
        self['resid'] = self['fit']-catalog[ix]
        self.mask['resid'] = False
        self.mask['id'] = np.logical_or(np.abs(self['resid']) > mindist,
                                        self['peak'] < minpeak)

    def fit(self, *args):
        """Return fitted wavelengths for input positions:
        x position for LineTable; chip,x for ThreeChipLineTable.
        If no arguments are given, returns wavelengths for all lines in table.
        """
        if len(args) == 0:
            args = (self[indep] for indep in self.meta['indep'])
        return self._fitfunc(self.meta['par'], *args)

    def estimate(self, *args):
        """Return wavelengths estimated using guess for solution."""
        if len(args) == 0:
            args = (self[indep] for indep in self.meta['indep'])
        return self._fitfunc(self.meta['guess'], *args)

    def setguess(self, coeff):
        """Set initial guess for solution.  Typically called from calibrate."""
        self.meta['guess'] = coeff
        self.setfit(coeff)

    def setfit(self, coeff):
        """Set fit parameters and use them to calculate wavelengths and
        residuals for all lines.  Typically called from calibrate."""
        self.meta['par'] = coeff
        self['fit'] = self.fit()
        self['resid'] = self['id']-self['fit']

    def rms(self):
        """Return root-mean-square of residuals of current solution."""
        return np.sqrt(np.mean(self['resid']**2))

    def _fitfunc(self, par, x):
        """LineTable fit function
        meta['refwave'] + \
            par[1]*(x-par[0])+par[2]*(x-par[0])**2+... +par[n]*(x-par[0])**n
        """
        px = Polynomial(np.hstack((self.meta['refwave'],par[1:])))
        return px(x-par[0])

    def _residual(self, par, wave, *args):
        """Residual function used by leastsq fit routine."""
        return wave-self._fitfunc(par,*args)

    def dofit(self):
        """Determine new best fit given current set of lines."""
        pfit,cov = leastsq(self._residual, self.meta['par'],
                           args=tuple(self[indep]
                                      for indep in('id',)+self.meta['indep']))
        self.setfit(pfit)

    def calibrate(self, guess, minlist, catalog, doplot=False, **kwargs):
        """Calibrate line lines given initial guesses

        Parameters
        ----------
        guess: ~float array or list
            Initial guesses for solution.  For LineTable, the length determines
            the polynomial degree that is fitted.
        minlist: ~list of ~tuple, with maximum distance and minimum strength
            Calibration is attempted iteratively, to allow e.g. ever weaker
            lines, but with ever more stringent constraint on the maximum
            deviation.  E.g., [(3.,1.5e4), (0.5,2e3), (0.13,0.), (0.06,0.)],
            would run a first attempt finding lines within 3 units of the line
            catalog and with peak>1.5e4, a second attempt within 0.5 units and
            peak>2e3, etc.
        catalog: ~float array
            list of laboratory arc line wavelengths
        doplot: ~bool, optional
            whether to show residuals (default: False)

        Any further keywords get stored in the meta data.  These can include
        (f not set at initialisation):

        refwave: ~float       # for LineTable, ThreeChipLineTable
            reference wavelength for which X position is sought
        refchip: ~int         # for ThreeChipLineTable
            reference chip
        """
        self.setguess(guess)
        if kwargs:
            self.meta.update(kwargs)
        if doplot:
            cols = iter('rbgmcybbbbbbbbbbbbbbbbbbb')
            import matplotlib.pylab as plt
        for mindist, minpeak in minlist:
            self.match(catalog, mindist, minpeak)
            self.dofit()
            if doplot:
                c = cols.next()
                plt.scatter(self['id'],
                            self['id']-self.estimate(),c=c)
                plt.plot(self['fit'],
                         self['fit']-self.estimate(), c=c)

    def write(self, *args, **kwargs):
        """Write to disk, using Table writer.  Can be restored with .read"""
        return Table.write(self.filled(), *args, **kwargs)


class ThreeChipLineTable(LineTable):
    """Lines in an arc spectrum over three chips and associated fits.

    Note
    ----
    The fitting function uses two polynomials, one in chip number to get an
    effective pixel position relative to the pixel at the reference wavelength,
    and another in this effective pixel number to get the wavelength:

    y = x - par[0] + par[1]*(chip-refchip) + par[2]*(chip-refchip)**2
    w = refwave + par[3]*y + par[4]*y**2 + ... + par[n]*y**(n-2)

    where refchip and refwave are taken from the meta data (they can
    be set in the call to calibrate)

    Above, the polynomial on chip number is quadratic (where the linear
    coefficient gives the average offset in pixels between chips, and the
    quadratic parameter accounts for differences in offset).  The degree of
    the polynomial on x position is determined by the number of guesses
    parameters given.
    """
    def init(self, chip, pos, peak, **kwargs):
        """Initialize ThreeChipLineTable with chips, positions, peak fluxes

        Parameters
        ----------
        chip, pos, peak: ~array
            chip number for arc line, fitted line position, line strength

        One can pass on additional parameters for storage in the meta-data
        here.  For ThreeChipLineTable, these would be
        refwave: ~float
            reference wavelength for which X position is sought
        refchip: ~int
            reference chip number
        """
        self._init(chip, pos, peak, names=('chip', 'x', 'peak'), **kwargs)

    @classmethod
    def fromlist(cls, tablist, **kwargs):
        """Initialize from a list of LineTables.

        Parameters
        ----------
        tablist: ~list of LineTable
            chip number assumed to be equal to list index

        One can pass on additional parameters for storage in the meta-data
        here.  For ThreeChipLineTable, these would be

        refwave: ~float
            reference wavelength for which X position is sought
        refchip: ~int
            reference chip number
        """
        return cls(np.hstack((np.ones(len(tab),dtype=np.int)*i
                              for i,tab in enumerate(tablist))),
                   np.hstack((tab['x'].data.data for tab in tablist)),
                   np.hstack((tab['peak'].data.data for tab in tablist)),
                   **kwargs)

    def _fitfunc(self, p, chip, x):
        """ThreeChipLineTable fit function:
        y = x - par[0] + par[1]*(chip-refchip) + par[2]*(chip-refchip)**2
        w = refwave + par[3]*y + par[4]*y**2 + ... + par[n]*y**(n-2)
        where refchip, refwave are taken from the Table meta data.
        """
        pchip = Polynomial([-p[0], p[1], p[2]])
        px = Polynomial(np.hstack((self.meta['refwave'], p[3:])))
        return px(x+pchip(chip-self.meta['refchip']))


def search(arc, min_counts=None, min_curvature=None, sigma=None, x=None):
    """Find peak locations and sizes in spectral arc, above some threshold.

    Peaks are found by looking for local maxima, possibly after
    smoothing the spectrum, and interpolating between the high point
    and its two surrounding ones using a parabola.  Weak peaks can be
    rejected either based on an absolute threshold or by having small
    curvature (i.e., the second-order term in the parabola).

    Parameters
    ----------
    arc: ~array
        One-dimensional arc spectrum
    min_counts: ~float, optional
        Absolute threshold for inclusion of line in line list (default: None)
    min_curvature: ~float, optional
        Threshold for curvature (quadrature term in parabola; default: None)
    sigma: ~float, optional
        Width of Gaussian (in pixels) with which to smooth arc before
        attempting to locate maximum (default: None; recommended: 1.)
    x: ~array, optional
        Positions associated with array.  Useful for arrays with holes,
        or to keep track of binning, etc.  If None (default), indices of arc
        are used (i.e., x is assumed to run from 0 to len(arc)-1).
    """
    assert arc.ndim == 1  # only 1-dimensional spectra

    if x is None:
        x = np.arange(len(arc))

    if sigma is not None:
        arc = gaussian_filter1d(arc, sigma, mode='nearest')

    peak = np.logical_and(arc[1:-1] > arc[:-2], arc[1:-1] > arc[2:])
    if min_counts is not None:
        peak = np.logical_and(peak, arc[1:-1] > min_counts)

    peakloc = np.where(peak)[0]+1
    x0, y0 = x[peakloc], arc[peakloc]
    xm, ym = x[peakloc-1]-x0, arc[peakloc-1]-y0
    xp, yp = x[peakloc+1]-x0, arc[peakloc+1]-y0
    # parabola y = ax^2+bx+c
    b = (yp*xm**2 - ym*xp**2) / (xp*xm**2 - xm*xp**2)
    a = (ym*xp - yp*xm) / (xp*xm**2 - xm*xp**2)
    # location of maximum -b/2a, peak counts
    offset = -b/2./a
    x, y = x0+offset, y0+b*offset+a*offset**2
    if min_curvature is None:
        return x, y
    else:
        ok = -a > min_curvature
        return x[ok], y[ok]


def read(*args, **kwargs):
    """Read line table and initialize class from meta['class']."""
    return LineTable.read(*args, **kwargs)
