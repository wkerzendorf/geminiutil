from .. import base
from ..base import BaseProject
from .gmos_alchemy import GMOSMOSRawFITS
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


class GMOSMOSProject(BaseProject):

    def __init__(self, database_string):
        super(GMOSMOSProject, self).__init__(database_string, GMOSMOSRawFITS)


    def classify_added_fits(self, current_fits):
        fits_object = self.add_gmos_raw_fits(current_fits)
        if fits_object is not None:
            return fits_object
        #self.add_gmos_mask(self)

    def add_gmos_raw_fits(self, fits_file):
        required_categories = [base.Object, base.Program, base.ObservationBlock, base.ObservationClass, base.ObservationType,
                               base.Instrument]

        required_keywords = [item.category_keyword for item in required_categories] + ['date-obs']

        if not all([keyword in fits_file.header for keyword in required_keywords]):
            logger.info("%s is not a normal raw gmos fits file", fits_file.fname)
            return

        object = base.Object.from_fits_object(fits_file)
        program = base.Program.from_fits_object(fits_file)
        observation_block = base.ObservationBlock.from_fits_object(fits_file)
        observation_class = base.ObservationClass.from_fits_object(fits_file)
        observation_type = base.ObservationType.from_fits_object(fits_file)
        instrument = base.Instrument.from_fits_object(fits_file)


        date_obs_str = '%sT%s' % (fits_file.header['date-obs'], fits_file.header['time-obs'])
        date_obs = datetime.strptime(date_obs_str, '%Y-%m-%dT%H:%M:%S.%f')

        gmos_raw = self.raw_fits_class(date_obs, instrument.id, observation_block.id, observation_class.id,
                                  observation_type.id, object.id, None, exclude=False)

        gmos_raw.id = fits_file.id

        self.session.add(gmos_raw)
        self.session.commit()




def add_gmos_mask(self, fits_file):
        necessary_keywords = ['GEMPRGID', 'OBSTYPE', 'ODFNAME']
        if not all([keyword in fits_file.header for keyword in necessary_keywords]):
            print "%s is not a normal gemini mask file" % fits_file.fname
            return




