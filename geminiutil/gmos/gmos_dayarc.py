# -*- coding: utf-8 -*-
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.units import Quantity

from .basic import dayarc
from .basic.estimate_disp import estimate_disp


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
            self.line_catalog = Quantity(line_catalog, u.Angstrom)
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
        arctab = dayarc.get_arcs(prepared_arc.fits.fits_data,
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
        lines = dayarc.find_lines(arctab, min_curvature=self.min_curvature)
        # calibrate
        chip1guess = [xref_estimate, disp_estimate, 0.]
        chippar = [2087.3, 1.5]
        linesall = dayarc.calibrate(
            lines, wref, chip1guess, self.line_catalog.value,
            minlist1=self.minlist1, chippar=chippar, minlist2=self.minlist2,
            doplot=doplot)
        linesall.meta['xref_estimate'] = xref_estimate
        linesall.meta['disp_estimate'] = disp_estimate
        arctab['w'] = linesall.fit(np.arange(3)[np.newaxis,:],
                                   arctab['x'])
        return arctab, linesall, shift, fake
