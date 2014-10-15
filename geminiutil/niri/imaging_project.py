import logging

from geminiutil.gmos.gmos_project import GMOSProject

from geminiutil.base.base_project import BaseProject
from geminiutil.niri.alchemy.imaging_alchemy import NIRIImagingRawFITS, \
    NIRIClassifyError

from geminiutil.base.alchemy.file_alchemy import FITSClassifyError

logger = logging.getLogger(__name__)

class NIRIImagingProject(GMOSProject):
    """
    TBD
    """
    def __init__(self, database_string, work_dir, echo=False):
        super(NIRIImagingProject, self).__init__(database_string, work_dir,
                                             NIRIImagingRawFITS, echo=echo)

    def add_fits_file(self, fname):
        """
        Add FITS file to database

        Parameters
        ----------

        fname: str
            FITS file name/path to add to database
        """
        try:
            NIRIImagingRawFITS.verify_fits_class(fname)
        except FITSClassifyError:
            pass
        else:
            NIRIImagingRawFITS.from_fits_file(fname, self.session)

