
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy import Column, ForeignKey

from astropy.io import fits


#sqlalchemy types
from sqlalchemy import String, Integer, Float, DateTime

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

class Observation(Base):
    __tablename__ = 'observation'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)


    def __init__(self, name, description=None):
        self.name = name
        self.description = description


class ObservationType(Base):
    __tablename__ = 'observation_type'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.description = description

class ObservationClass(Base):
    __tablename__ = 'observation_class'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.description = description





class Instrument(Base):
    __tablename__ = 'instrument'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

def GeminiRawFITS(Base):
    __tablename__ = 'raw_fits_file'




