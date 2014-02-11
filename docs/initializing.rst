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

>>> proj.initialize()

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

Grouping Science Sets
^^^^^^^^^^^^^^^^^^^^^

Science sets are a way to group each individual science observation with its calibration data

>>> proj.link_science_sets()

