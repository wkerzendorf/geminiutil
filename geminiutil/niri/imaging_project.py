import logging

from geminiutil.gmos.gmos_project import GMOSProject

from geminiutil.base.base_project import BaseProject
from geminiutil.niri.alchemy.imaging_alchemy import NIRIImagingRawFits, \
    NIRIClassifyError

logger = logging.getLogger(__name__)

class NIRIImagingProject(GMOSProject):
    """
    TBD
    """
    def __init__(self, database_string, work_dir, echo=False):
        super(NIRIImagingProject, self).__init__(database_string, work_dir,
                                             NIRIImagingRawFits, echo=echo)

    def add_fits_file(self, fname):
        """
        Add FITS file to database

        Parameters
        ----------

        fname: str
            FITS file name/path to add to database
        """

        NIRIImagingRawFits.verify_fits_class(fname)
        NIRIImagingRawFits.from_fits_file(fname, self.session)

    def classify_raw_fits(self, current_fits):
        pass
