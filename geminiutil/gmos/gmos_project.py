from ..base import BaseProject
from .gmos_alchemy import GMOSMOSRawFITS


class GMOSMOSProject(BaseProject):

    def __init__(self, database_string):
        super(GMOSMOSProject, self).__init__(database_string, GMOSMOSProject)


    def classify_added_fits(self, current_fits):
        pass

    #
