**************************************
Initializing Database for `geminiutil`
**************************************


Initializing the database is an important step in using `geminiutil`. Each project is contained in a database that keeps
track of files and operations.

Initializing Database
^^^^^^^^^^^^^^^^^^^^^


>>> import geminiutil.base as base
>>> from geminiutil.gmos import GMOSMOSProject
>>> proj = GMOSMOSProject('sqlite:///mcsnr.db3')

First time, initialize to read in GMOS filter/grating information

>>> proj.initialize_database()

Adding Data to the DB
^^^^^^^^^^^^^^^^^^^^^

>>> proj.add_directory('/media/data1/mcsnr/gmos_data/raw', file_filter='S*S*.fits')

One can add different directories as well:

>>> proj.add_directory('/media/data1/mcsnr/gmos_data/mdf_dir')


Linking Masks
^^^^^^^^^^^^^

This step will link the observations to specific masks (only in the database). Longslit exposures have a link to masks as
well, however, these do not have corresponding files. Currently we only support the "xxarcsec" format for longslit mask names

>>> proj.link_masks()


Linking Arcs
^^^^^^^^^^^^^

This step will link the Longslit arcs to a GMOSLongSlitArc object that contains information about the arc used. In addition,
these GMOSLongSlitArcs have methods for calibration.

>>> proj.link_long_slit_arcs()



Grouping Science Sets
^^^^^^^^^^^^^^^^^^^^^

Science sets are a way to group each individual science observation with its calibration data. For this to work the function
must be given a dictionary linking science instrument setups to longslit arc instrument setups
first look at the science instrument_setups

>>> proj.science_instrument_setups
[<GMOS MOS Instrument Setup ID=11 Filter1 open Filter2 open Grating B600+_G5323 Tilt 56.82 central wave=430.00 nm>,
 <GMOS MOS Instrument Setup ID=12 Filter1 open Filter2 open Grating B600+_G5323 Tilt 56.05 central wave=470.00 nm>,
 <GMOS MOS Instrument Setup ID=13 Filter1 OG515_G0330 Filter2 open Grating R400+_G5325 Tilt 55.67 central wave=735.00 nm>,
 <GMOS MOS Instrument Setup ID=14 Filter1 OG515_G0330 Filter2 open Grating R400+_G5325 Tilt 55.15 central wave=775.00 nm>]

>>> proj.longslit_arcs_instrument_setups
Out[3]:
[<GMOS MOS Instrument Setup ID=17 Filter1 open Filter2 open Grating B600+_G5323 Tilt 56.43 central wave=450.00 nm>,
 <GMOS MOS Instrument Setup ID=25 Filter1 open Filter2 open Grating B600+_G5323 Tilt 57.20 central wave=410.00 nm>,
 <GMOS MOS Instrument Setup ID=26 Filter1 open Filter2 open Grating B600+_G5323 Tilt 55.67 central wave=490.00 nm>,
 <GMOS MOS Instrument Setup ID=16 Filter1 OG515_G0330 Filter2 open Grating R400+_G5325 Tilt 55.41 central wave=755.00 nm>,
 <GMOS MOS Instrument Setup ID=27 Filter1 OG515_G0330 Filter2 open Grating R400+_G5325 Tilt 56.05 central wave=705.00 nm>,
 <GMOS MOS Instrument Setup ID=30 Filter1 OG515_G0330 Filter2 open Grating R400+_G5325 Tilt 0.00 central wave=755.00 nm>,
 <GMOS MOS Instrument Setup ID=28 Filter1 OG515_G0330 Filter2 open Grating R400+_G5325 Tilt 55.92 central wave=715.00 nm>,
 <GMOS MOS Instrument Setup ID=29 Filter1 OG515_G0330 Filter2 open Grating R400+_G5325 Tilt 54.90 central wave=795.00 nm>]

>>> proj.link_science_sets({11:17, 12:17, 13:16, 14:16}
)

