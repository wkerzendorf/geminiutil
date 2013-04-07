from .gemini_alchemy import Base, FitsFile, GeminiRawFITS, Instrument, Program, Observation

from sqlalchemy import engine, create_engine
from sqlalchemy.orm import sessionmaker, backref, relationship

from glob import glob
import os

class BaseProject(object):

    def __init__(self, database_string):
        self.metadata = Base.metadata
        self.engine = create_engine(database_string)
        self.metadata.bind = self.engine
        self.metadata.create_all()
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()


    def add_directory(self, directory, file_filter='*.fits'):
        for fname in glob(os.path.join(directory, file_filter)):
            current_fits = FitsFile.from_fits_file(fname)
            self.session.add(current_fits)

        self.session.commit()

    def scan_raw_gemini(self):
        necesarry_keywords = ['instrume', 'object', 'obstype', 'obsclass', 'gemprgid', 'obsid', 'date-obs', 'time-obs']
        for fits_file in self.session.query(FitsFile):
            if self.session.query(GeminiRawFITS).filter_by(id=fits_file.id) != []:
                continue
            if not all([keyword in fits_file.header for keyword in necesarry_keywords]):
                print "%s is not a normal raw gemini fits file"
                continue

            instrument_str = fits_file.header['instrume'].lower().strip()
            program_str = fits_file.header['gemprgid'].lower().strip()
            obsclass_str = fits_file.header['obsclass'].lower().strip()

            if self.session.query(Instrument).filter_by(name=instrument_str).count() == 0:
                current_instrument = Instrument(instrument_str)
                self.session.add(current_instrument)
                self.session.commit()
            else:
                current_instrument = self.session.query(Instrument).filter_by(name=instrument_str).one()

            if self.session.query(Program).filter_by(name=program_str).count() == 0:
                current_program = Program(program_str)
                self.session.add(current_program)
                self.session.commit()
            else:
                current_program = self.session.query(Program).filter_by(name=program_str)








