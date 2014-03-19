from geminiutil.gmos.gmos_project import GMOSProject


from geminiutil.gmos.alchemy.imaging_raw import GMOSImagingRawFITS


class GMOSImagingProject(GMOSProject):
    """
    TBD
    """
    def __init__(self, database_string, work_dir, echo=False):
        super(GMOSImagingProject, self).__init__(database_string, work_dir,
                                             GMOSImagingRawFITS, echo=echo)

