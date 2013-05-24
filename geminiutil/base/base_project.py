from .gemini_alchemy import Base, FITSFile, Instrument, Program, ObservationBlock, ObservationClass, ObservationType

from sqlalchemy import engine, create_engine
from sqlalchemy.orm import sessionmaker, backref, relationship, object_session
from datetime import datetime
from glob import glob
import os
import logging


logger = logging.getLogger(__name__)


def get_category(session, category_str, category):
    if session.query(category).filter_by(name=category_str).count() == 0:
        current_category = category(category_str)
        session.add(current_category)
        session.commit()
    else:
        current_category = session.query(category).filter_by(name=category_str).one()

    return current_category

class BaseProject(object):

    def __init__(self, database_string, raw_fits_class, echo=False):
        self.metadata = Base.metadata
        self.engine = create_engine(database_string, echo=echo)
        self.metadata.bind = self.engine
        self.metadata.create_all()
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.conn = self.session.bind.connect()

        self.raw_fits_class = raw_fits_class



    def add_directory(self, directory, file_filter='*.fits'):
        for fname in sorted(glob(os.path.join(directory, file_filter))):
            logger.info('Adding %s to project', fname)
            current_fits = FITSFile.from_fits_file(fname)
            self.session.add(current_fits)
            self.session.commit()
            self.classify_added_fits(current_fits)

        self.session.commit()

    def classify_added_fits(self, current_fits):
        fits_object = self.add_gemini_raw_fits(current_fits)
        return fits_object



    def add_gemini_raw_fits(self, fits_file):
        necessary_keywords = ['instrume', 'object', 'obstype', 'obsclass', 'gemprgid', 'obsid', 'date-obs', 'time-obs']

        #print "Working on %s" % fits_file.fname

        instrument_str = fits_file.header['instrume'].lower().strip()
        current_instrument = get_category(self.session, instrument_str, Instrument)

        program_str = fits_file.header['gemprgid'].lower().strip()
        current_program = get_category(self.session, program_str, Program)

        observation_str = fits_file.header['obsid'].lower().strip()
        current_observation_block = get_category(self.session, observation_str, ObservationBlock)

        if current_observation_block.program_id is None:
            current_observation_block.program_id = current_program.id
            self.session.commit()

        obstype_str = fits_file.header['obstype'].lower().strip()
        current_observation_type = get_category(self.session, obstype_str, ObservationType)

        obsclass_str = fits_file.header['obsclass'].lower().strip()
        current_observation_class = get_category(self.session, obsclass_str, ObservationClass)

        date_obs_str = '%sT%s' % (fits_file.header['date-obs'], fits_file.header['time-obs'])
        date_obs = datetime.strptime(date_obs_str, '%Y-%m-%dT%H:%M:%S.%f')

        current_raw_fits = raw_fits_class(date_obs, current_instrument.id, current_observation_block.id, current_observation_class.id, current_observation_type.id)
        current_raw_fits.id = fits_file.id

        self.session.add(current_raw_fits)

        self.session.commit()
        return current_raw_fits

















