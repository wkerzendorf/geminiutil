***************************
GMOS Wavelength Calibration
***************************

Quick & Automatic Wavelength Calibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





Manual wavelength calibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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