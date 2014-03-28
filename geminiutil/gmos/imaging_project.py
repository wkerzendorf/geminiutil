from geminiutil.gmos.gmos_project import GMOSProject


from geminiutil.gmos.alchemy.imaging_raw import GMOSImagingRawFITS


class GMOSImagingProject(GMOSProject):
    """
    TBD
    """
    def __init__(self, database_string, work_dir, echo=False):
        super(GMOSImagingProject, self).__init__(database_string, work_dir,
                                             GMOSImagingRawFITS, echo=echo)


    def add_fits_file(self, fname):
        """
        Add FITS file to database

        Parameters
        ----------

        fname: str
            FITS file name/path to add to database
        """

        GMOSImagingRawFITS.verify_fits_class(fname)
        GMOSImagingRawFITS.from_fits_file(fname, self.session)