from geminiutil.gmos.gmos_project import GMOSProject
from geminiutil.gmos.gmos_alchemy import GMOSMOSRawFITS
class GMOSImagingProject(GMOSProject):
    """
    TBD
    """
    def __init__(self, database_string, work_dir, echo=False):
        super(GMOSImagingProject, self).__init__(database_string, work_dir,
                                             GMOSMOSRawFITS, echo=echo)

