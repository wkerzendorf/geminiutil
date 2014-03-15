************************************
Extracting point sources from Slices
************************************


We again start out with having a project open called proj::

    >>> proj = GMOSMOSProject(dbname, work_dir=work_dir, echo=False)


Calculating Slice Geometries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Sources in the Database
^^^^^^^^^^^^^^^^^^^^^^^

For this to work the database has to be constructed and has linked the science sets together. Each slice is
associated with a number of MOSPointSources. These objects hold information about the priority and the slit
position for the current slice. Having a slice object it is easy to access the associated MOS Point Sources::

    >>> slice = proj.science_sets[0].slices[0]
    >>> mos_ps = slice.mos_point_sources[0]
    >>> mos_ps
    <MosPointSouce id=25036361 priority=0 slit position=0.0000>
    >>> mos_ps.slices
    [<GMOS MOS Slice lower_edge=741.76 upper_edge=755.46)>]