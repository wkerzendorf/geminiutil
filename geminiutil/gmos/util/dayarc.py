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

import extract_psf.extract as extract

from astropy.table import Table
import astropy.units as u

from prepare_frame import multiext_data, multiext_header_value, multiext_x_coords
from wavecal import LineTable, ThreeChipLineTable
import estimate_disp


class GMOSDayArc(object):
    def __init__(self, line_catalog,
                 xref_grid=np.arange(-40, 40),
                 frac_disp_grid=np.linspace(0.97, 1.03, 31),
                 min_curvature=[5.,3.,2.],
                 minlist1=[(3.,1e3), (1.,3e2), (0.26,0.), (0.1,0.)],
                 minlist2=[(2.,0.), (0.26,0.), (0.1,0.)]):
        try:  # is this a quantity with a unit of length?
            line_catalog.to(u.m)
        except:  # convert to Å
            self.line_catalog = u.Quantity(line_catalog, u.Angstrom)
        else:
            self.line_catalog = line_catalog
        self.xref_grid = xref_grid
        self.frac_disp_grid = frac_disp_grid
        self.min_curvature = min_curvature
        self.minlist1 = minlist1
        self.minlist2 = minlist2
        # for if one really needs to twiddle
        self.get_arcs_skypol = 2
        self.get_arcs_clip = 3.

    def __call__(self, prepared_arc, doplot=False):
        """Find wavelength solution.
        prepared_arc should have a limited range in the spatial direction."""
        instrument_setup = prepared_arc.raw_fits.instrument_setup
        wref = (instrument_setup.grating_central_wavelength_value *
                instrument_setup.grating_central_wavelength_unit
                ).to(self.line_catalog.unit).value
        # instrument_setup.detectors[1].naxis1 includes overscan
        xref_guess = (prepared_arc.fits.fits_data[2].shape[-1] *
                      instrument_setup.x_binning // 2)
        disp_guess = -(instrument_setup.spectral_pixel_scale /
                       (instrument_setup.x_binning / u.pix)
                       ).to(self.line_catalog.unit).value
        # slit width not yet in data base for long slits
        slit_width = prepared_arc.fits.fits_data[0].header['maskname']
        slit_width = float(slit_width.replace('arcsec','')) * u.arcsec
        # get resolution in Å, and then in pixels
        slit_width = wref / instrument_setup.calculate_resolution(slit_width)
        slit_width /= abs(disp_guess)
        slit_sigma = abs(disp_guess)  # additional instrumental broadening

        # extract spectrum from image
        arctab = get_arcs(prepared_arc.fits.fits_data,
                                 skypol=self.get_arcs_skypol,
                                 clip=self.get_arcs_clip)
        # estimate dispersion and position of reference wavelength
        xref_estimate, disp_estimate, shift, fake \
            = estimate_disp(Table([arctab['x'][:,1], arctab['f'][:,1]]),
                            wref=wref,
                            xref_grid=xref_guess+self.xref_grid,
                            disp_grid=disp_guess*self.frac_disp_grid,
                            line_catalog=self.line_catalog.value,
                            slit_width=slit_width, sigma=slit_sigma,
                            full=True)
        # find lines
        lines = find_lines(arctab, min_curvature=self.min_curvature)
        # calibrate
        chip1guess = [xref_estimate, disp_estimate, 0.]
        chippar = [2087.3, 1.5]
        linesall = calibrate(
            lines, wref, chip1guess, self.line_catalog.value,
            minlist1=self.minlist1, chippar=chippar, minlist2=self.minlist2,
            doplot=doplot)
        linesall.meta['xref_estimate'] = xref_estimate
        linesall.meta['disp_estimate'] = disp_estimate
        arctab['w'] = linesall.fit(np.arange(3)[np.newaxis,:],
                                   arctab['x'])
        return arctab, linesall, shift, fake


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
