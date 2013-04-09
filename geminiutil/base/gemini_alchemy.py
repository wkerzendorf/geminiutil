
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy import Column, ForeignKey

from astropy.io import fits


#sqlalchemy types
from sqlalchemy import String, Integer, Float, DateTime, Boolean

from glob import glob
import numpy as np
from astropy.io import fits
import abc
import os
import hashlib

def hashfile(afile, hasher, blocksize=65536):
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.hexdigest()

Base = declarative_base()


class FitsFile(Base):
    __tablename__ = 'base_fits'

    id = Column(Integer, primary_key=True)
    fname = Column(String)
    path = Column(String)
    size = Column(Integer)
    extensions = Column(Integer)
    md5 = Column(String)
    scanned = Column(Boolean)

    @classmethod
    def from_fits_file(cls, fname):
        extensions = len(fits.open(fname))
        filename = os.path.basename(fname)
        path = os.path.abspath(os.path.dirname(fname))
        filesize = os.path.getsize(fname)
        md5_hash = hashfile(file(fname, 'rb'), hashlib.md5())

        return cls(filename, path, filesize, md5_hash, extensions)


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

    def __init__(self, fname, path, size, md5, extensions):
        self.fname = fname
        self.path = path
        self.size = size
        self.md5 = md5
        self.extensions = extensions


class Program(Base):
    __tablename__ = 'program'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Program %s>' % self.name


class object(Base):
    __tablename__ = 'object'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    ra = Column(Float)
    dec = Column(Float)
    description = Column(String)

    def __init__(self, name, ra, dec, description=None):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.description = description

    def __repr__(self):
        return '<Object %s>' % self.name


class ObservationBlock(Base):
    __tablename__ = 'observation_block'

    id = Column(Integer, primary_key=True)
    program_id = Column(Integer, ForeignKey('program.id'))
    name = Column(String)
    description = Column(String)

    program = relationship('Program')

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Observation Block %s>' % self.name


class ObservationType(Base):
    __tablename__ = 'observation_type'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Observation Type %s>' % self.name


class ObservationClass(Base):
    __tablename__ = 'observation_class'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Observation Class %s>' % self.name


class Instrument(Base):
    __tablename__ = 'instrument'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return "<Gemini Instrument %s>" % self.name

"""
class GeminiRawFITS(object):
    __tablename__ = 'raw_fits_file'







    def __init__(self, date_obs, instrument_id, observation_block_id, observation_class_id, observation_type_id, exclude=False):
        self.date_obs = date_obs
        self.instrument_id = instrument_id
        self.observation_block_id = observation_block_id
        self.observation_class_id = observation_class_id
        self.observation_type_id = observation_type_id
        self.exclude = exclude

    def __repr__(self):
        return "<raw FITS - Instrument %s Block %s Class %s Type %s>" % (self.instrument.name, self.observation_block.name,
                                    self.observation_class.name, self.observation_type.name)

"""





