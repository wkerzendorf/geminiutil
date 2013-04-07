
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy import Column, ForeignKey

#sqlalchemy types
from sqlalchemy import String, Integer, Float, DateTime

from glob import glob
import numpy as np
from astropy.io import fits
import abc
import os

Base = declarative_base()


class FitsFile(object):
    #__metaclass__ = abc.ABCMeta



    @property
    def full_path(self):
        return os.path.join(self.path, self.fname)

    @property
    def fits_data(self):
        return fits.open(self.full_path)

    @property
    def header(self):
        return fits.getheader(self.full_path)

    @property
    def data(self):
        return fits.getdata(self.full_path)


    @property
    def shape(self):
        return self.header['naxis1'], self.header['naxis2']


class Program(Base):
    __tablename__ = 'program'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

class ObservationBlock(Base):
    __tablename__ = 'observation_block'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, description):
        self.description = description



class Instrument(Base):
    __tablename__ = 'instrument'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

class GenericGeminiFITS(Base, FitsFile):
    __tablename__ = 'eso_raw_fits'

    id = Column(Integer, primary_key=True)
    program_id = Column(Integer, ForeignKey('program.id'))
    observation_block_id = Column(Integer, ForeignKey('observation_block.id'))
    instrument_id = Column(Integer, ForeignKey('instrument.id'))
    date_obs = Column(DateTime)
    fname = Column(String)
    path = Column(String)

    observation_block = relationship(ObservationBlock)


    def __init__(self, program_id, observation_block_id, template_id, instrument_id, date_obs, fname, path):
        self.program_id = program_id
        self.observation_block_id = observation_block_id
        self.instrument_id = instrument_id
        self.date_obs = date_obs
        self.fname = fname
        self.path = path

    def __repr__(self):
        try:
            template_name = self.template.name
        except:
            template_name = 'N/A'

        try:
            obs_block_id = self.observation_block.id
        except:
            obs_block_id = 'N/A'


        return "<Gemini FITS Obsblock %s Date-Obs %s>" % (obs_block_id, self.date_obs)



