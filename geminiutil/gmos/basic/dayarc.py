"""Helper functions for wavelength calibration from day-time long-slit arcs.

Sample usage
------------
arctab = dayarc.get_arcs(fits.open(arcfitsname))
lines = dayarc.find_lines(arctab, min_curvature=[5.,3.,2.])
chip1guess = [xrefguess, dwguess, 0.]
chippar = [2087.3, 1.5]
linesall = dayarc.calibrate(
                lines, arctab.meta['grwlen']*10, chip1guess, lincat['w'],
                minlist1=[(3.,1e3), (1.,3e2), (0.26,0.), (0.1,0.)],
                chippar=chippar,
                minlist2=[(2.,0.), (0.26,0.), (0.1,0.)],
                doplot=True)
arctab['w'] = linesall.fit(np.arange(3)[np.newaxis,:], arctab['x'])

Notes
-----
Above, arcfitsname would be a "prepared" file, for which bias subtraction and
gain correction has been done, and a relevant region extracted, e.g., by

PrepDayArc = gmos_prepare.Prepare(data_subslice=[slice(1150,1250),slice(-1)])
PrepDayArc(fits.open('some-raw-fits-file')).writeto(arcfitsname)
"""

from __future__ import print_function, division

import numpy as np

#import extract_psf.extract as extract
from astropy.table import Table

from prepare import multiext_data, multiext_header_value, multiext_x_coords
from wavecal import LineTable, ThreeChipLineTable


def get_arcs(image, skypol=2, clip=3, fferr=0.015):
    """Fit arc image and create table with one-dimensional arc spectra.

    Parameters
    ----------
    image: ~fits.HDUList
        prepared image set with one extension for each chip
    skypol: ~int, optional
        degree with which image is fit along spatial direction (default: 2)
    clip: ~float, optional
        maximum sigma with which points are allowed to deviate from fit
        (default: 3)
    fferr: ~float, optional
        Flat field uncertainty, used for estimating uncertainties
        (default: 0.015)

    Returns
    -------
    Table with columns 'x' (positions on exposed part of CCD) and
    extracted fluxes 'f' (with dimension equal to number of extensions)
    """
    data = multiext_data(image)
    ron = multiext_header_value(image, 'RDNOISE')
    arc, chi2, _, ntbadl, ntbadh = extract.fitsky(data, ron=ron, skypol=skypol,
                                                  clip=clip, fferr=0.015,
                                                  ibadlimit=7)
    x = multiext_x_coords(image, 'CCDSEC').squeeze()
    arctab = Table([x.T, arc.T, chi2.T], names=('x', 'f', 'chi2'))
    arctab.meta['grwlen'] = image[0].header['GRWLEN']
    return arctab


def find_lines(arctab, min_curvature=[10.,10.,7.], sigma=1.):
    """Find lines in one-dimensional arc spectra.

    Parameters
    ----------
    arctab: ~Table
        with 'x' and 'f', usually made with ~get_arcs
    min_curvature: ~list of ~float
        curvature thresholds for each extension; see ~wavecal.search
    sigma: ~float
        Gaussian smoothing kernal applied before searching for lines;
        see ~wavecal.search
    """
    if hasattr(min_curvature, '__len__'):
        assert len(min_curvature) == arctab['f'].shape[1]
    else:
        min_curvature = [min_curvature]*arctab['f'].shape[1]
    lines = []
    for i, min_curv in enumerate(min_curvature):
        lines.append(LineTable.fromsearch(arctab['f'][:,i], x=arctab['x'][:,i],
                                          min_curvature=min_curv,
                                          sigma=sigma))
    return lines


def calibrate(lines, refwave, chip1guess, catalog,
              minlist1=[(3.,1.5e4), (1.,2e3), (0.26,0.), (0.1,0.)],
              chippar=[2087.3, 0.54], extrapar=[0.],
              minlist2=None,
              doplot=False):
    if minlist2 is None:
        minlist2 = minlist1
    if refwave:
        lines[1].meta['refwave'] = refwave
    else:
        refwave = lines[1].meta['refwave']
    lines[1].calibrate(chip1guess, minlist1, catalog, doplot)
    linesall = ThreeChipLineTable.fromlist(lines, refchip=1, refwave=refwave)
    linesall.calibrate(np.hstack((lines[1].meta['par'][:1], chippar,
                                  lines[1].meta['par'][1:], extrapar)),
                       minlist2, catalog, doplot=doplot)
    return linesall
