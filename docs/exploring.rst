*********************
Exploring the Dataset
*********************



This guide is written with a GMOS MOS dataset in mind. But components should also work for future parts of this framework.

To load the database
>>> from geminiutil import gmos, base
>>> proj = gmos.GMOSMOSProject('sqlite:///mcsnr.db3')
>>>

Each FITS file has an entry in the database as ~FITSFile and we can query the FITS files:

>>> test_fits = proj.session.query(base.FITSFile).first()
>>> test_fits
<geminiutil.base.gemini_alchemy.FITSFile at 0x11611d10>
>>> test_fits.path, test_fits.fname
(u'/media/data1/mcsnr/gmos_data/raw', u'GS20121011S0254_BIAS.fits')

All FITS files that are actually recorded with the detector have an entry at ~GMOSMOSRawFITS. This object links back to
the FITSFile:

>>> test_fits.raw_fits
<gmos fits="GS20121011S0254_BIAS.fits" class="daycal" type="bias" object="bias">
>>> test_fits.mjd
56211.447782407406

All raw fits entries are linked to a type and class. We can find all types and classes

>>> proj.observation_classes
(u'daycal', u'acq', u'science', u'partnercal', u'acqcal', u'progcal')
>>> proj.observation_types
(u'bias', u'flat', u'object', u'arc')

We can also query for raw fits files that belong to a particular type

>>> proj.session.query(gmos.GMOSMOSRawFITS).join(base.ObservationClass).filter(base.ObservationClass.name=='science').first()
<gmos fits="S20121010S0086.fits" class="science" type="object" object="j0103.5-7247">

There are short cuts to do those queries:

>>> proj.science
[<gmos fits="S20121010S0086.fits" class="science" type="object" object="j0103.5-7247">,
 <gmos fits="S20121010S0089.fits" class="science" type="object" object="j0103.5-7247">,
 <gmos fits="S20121010S0092.fits" class="science" type="object" object="j0103.5-7247">,
 <gmos fits="S20121010S0095.fits" class="science" type="object" object="j0103.5-7247">,
...
]
>>> proj.flat
[<gmos fits="S20120917S0215.fits" class="daycal" type="flat" object="gcalflat">,
 <gmos fits="S20120917S0218.fits" class="daycal" type="flat" object="gcalflat">,
 <gmos fits="S20121009S0038.fits" class="daycal" type="flat" object="gcalflat">,
 ...
]

If you would like to get the query item back to do more filtering on your query that is also possible:

>>> test_fits = proj.science_query.filter_by(gmos.GMOSMOSRawFITS.mjd > 56215.).first()
>>> test_fits
<gmos fits="S20121103S0147.fits" class="science" type="object" object="j0459.9-7008">
>>> test_fits.mjd
56234.21007291666

Instrument Setup
^^^^^^^^^^^^^^^^

Each GMOSMOSRawFITS is also linked to an InstrumentSetup that knows about how GMOS was setup

>>> isetup = test_fits.instrument_setup
>>> isetup
<GMOS MOS Instrument Setup Filter1 open Filter2 open Grating B600+_G5323 Tilt 56.82 central wave=430.00 nm>
>>>  isetup.instrument
<Gemini Instrument gmos-s>



