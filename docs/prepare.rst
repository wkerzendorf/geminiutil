************************
Preparing the FITS files
************************

Overview
^^^^^^^^

Preparing does overscan, blah blah more detail needed


Actual preparation
^^^^^^^^^^^^^^^^^^

1. The GMOSMOSFITSFile now have a `prepare` function that will automatically call the `GMOSPrepareFrameSet` class with
 sensible defaults (or you can give it a prepare-object). This will just return a fits file and make no changes to the Database

2. THE `GMOSMOSFITSFile` also has a `prepare_to_database`-function that will check if a prepared file already exists in the database.
if not it will automatically create one in the directory that is listed by the user. If a prepared file exists, one can override
with `force=True` and the old one will be deleted and new one created.
