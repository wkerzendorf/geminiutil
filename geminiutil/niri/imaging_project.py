from geminiutil.gmos.gmos_project import GMOSProject


from geminiutil.base.base_project import BaseProject


class NIRIImagingProject(GMOSProject):
    """
    TBD
    """
    def __init__(self, database_string, work_dir, echo=False):
        super(GMOSImagingProject, self).__init__(database_string, work_dir,
                                             GMOSImagingRawFITS, echo=echo)

