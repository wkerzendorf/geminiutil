import numpy as np
from astropy.table import Table

import extract_psf.extract as extract


def extract_spectrum(gmos_slice, model_errors=1, ff_noise=0.03):

    slitno = gmos_slice.list_id
    raw_science = gmos_slice.science_set.science
    prep_science = raw_science.prepared_fits
    instrument_setup = raw_science.instrument_setup
    x_binning = instrument_setup.x_binning
    y_binning = instrument_setup.y_binning
    prep_fits = prep_science.fits.fits_data

    sci_slice = slice(*(np.ceil(edge).astype(int)
                        for edge in (gmos_slice.lower_edge,
                                     gmos_slice.upper_edge)))

    scidata = np.array([prep_fits[chip].data[sci_slice, :]
                        for chip in range(1,4)])
    read_noise = np.array([amp.header['RDNOISE']
                           for amp in prep_fits[1:]]).reshape(-1,1,1)

    tracepos = np.array([scidata.shape[1]/2. -
                         raw_science.mask.table[slitno]['specpos_y'] /
                         y_binning])

    e_source = scidata
    psf = extract.PSF(form=0, guesses=[[0., 0.], 8., (2.5, '@')])
    for i in range(model_errors+1):
        # if we want to use model errors, first get approximate model
        # (iteratively, if model_errors>1)
        error_estimate = extract.ccd_noise(e_source, read_noise, ff_noise)
        (out,eout,back,chi2,test,
         ntdiscard,ntbadl,ntbadh,
         nproblems) = extract.extract(scidata, tracepos, psf,
                                      e=error_estimate, skypol=0,
                                      ibadlimit=2, squeeze_dims=False,
                                      itesttype=101 if i < model_errors
                                      else 103)
        # set error source to test frame
        e_source = test

    # make table of extracted spectra
    # Note: out.shape=(nstar,norder,nwav) but need
    #       (nwav,nstar,norder) for nstar>1 -> transpose(2,0,1)
    #       (nwav,norder) for nstar=1 -> transpose(1,0)
    x_offset = 0
    ndisp = prep_fits[1].data.shape[-1]
    x = x_offset + x_binning*np.arange(0.5, ndisp) - 0.5
    scitab = Table([np.array(x, dtype=np.float32)] +
                   [a.transpose(2,0,1) for a in [out, eout, back, chi2]],
                   names=['x', 'f', 'e', 's', 'chi2'])

    # fits0_header = prep_sci_fits[0].header
    # for hdr in fits0_header:
    #     if hdr not in ('SIMPLE', 'BITPIX', 'EXTEND', ''):
    #         scitab.meta[hdr.lower()] = fits0_header[hdr]
    scitab.meta['nstars'] = len(tracepos)
    scitab.meta['tracepos'] = tracepos
    scitab.meta['psfform'] = psf.form
    scitab.meta['psfnpoldisp'] = psf.npoldisp
    scitab.meta['psfpar'] = psf.par
    scitab.meta['psferr'] = psf.err
    scitab.meta['psfchi2'] = psf.chi2
    scitab.meta['ndiscard'] = ntdiscard
    scitab.meta['nbad'] = [ntbadl, ntbadh]
    scitab.meta['nproblems'] = nproblems

    return scitab, test
