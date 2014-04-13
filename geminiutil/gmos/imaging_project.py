from geminiutil.gmos.gmos_project import GMOSProject

import geminiutil
from geminiutil.gmos.alchemy.imaging_raw import GMOSImagingRawFITS
from geminiutil.base.base_project import BaseProject

from geminiutil.gmos.util.read_standard_fields import read_standard_star_db

import os

default_configuration_dir = os.path.join(geminiutil.__path__[0],
                                         'data', 'gmos')


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
        return GMOSImagingRawFITS.from_fits_file(fname, self.session)

    def initialize_standard_stars(self):
        standard_star_db_fname = os.path.join(default_configuration_dir,
                                              'gmosnlandolt.dat')
        standard_star_db = read_standard_star_db(standard_star_db_fname)

        for line in standard_star_db:
            print line
        1/0
            
    def initialize_database(self, configuration_dir=None):
        self.initialize_standard_stars()
        super(GMOSImagingProject, self).initialize_database(
            configuration_dir=configuration_dir)




