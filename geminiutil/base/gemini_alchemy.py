from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy import event
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, relationship, backref, object_session
from sqlalchemy import Column, ForeignKey

from astropy.io import fits
import logging
import os
import hashlib

logger = logging.getLogger(__name__)

#sqlalchemy types




from glob import glob
import numpy as np

def hashfile(afile, hasher, blocksize=65536):
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.hexdigest()

Base = declarative_base()





class CategoryBaseClass(object):
    category_keyword = None

    @classmethod
    def from_fits_object(cls, fits_object):
        session = object_session(fits_object)
        category_keyword_value = fits_object.header[cls.category_keyword].lower().strip()
        if session.query(cls).filter_by(name=category_keyword_value).count() == 0:
            logger.info('%s %s mentioned in %s not found in database - adding', cls.__name__ , category_keyword_value, fits_object.fname)
            new_category_object = cls(category_keyword_value)
            session.add(new_category_object)
            session.commit()
            return new_category_object
        else:
            category_object = session.query(cls).filter_by(name=category_keyword_value).one()
            logger.debug('%s %s mentioned in %s found in database - returning existing object %s', cls.__name__ , category_keyword_value, fits_object.fname, category_object)
            return category_object


class AbstractFileTable(Base):
    __abstract__ = True

    id = Column(Integer, primary_key=True)
    fname = Column(String)
    path = Column(String)
    size = Column(Integer)
    md5 = Column(String)

    work_dir = ''

    @classmethod
    def from_file(cls, full_fname, use_abspath=False):
        fname = os.path.basename(full_fname)
        if use_abspath:
            logger.warn('Using absolute paths is now discouraged - all files should be located in one workdirectory')
            path = os.path.abspath(os.path.dirname(full_fname))

        filesize = os.path.getsize(full_fname)
        md5_hash = hashfile(file(full_fname, 'rb'), hashlib.md5())

        return cls(fname=fname, path=path, size=filesize, md5=md5_hash)



class FITSFile(AbstractFileTable):
    __tablename__ = 'fits_file'

    extensions = Column(Integer)



    @classmethod
    def from_fits_file(cls, full_fname, use_abspath=False):
        extensions = len(fits.open(full_fname))

        fits_obj = cls.from_file(full_fname, use_abspath=use_abspath)
        fits_obj.extensions = extensions
        return fits_obj


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

    @property
    def children(self):
        session = object_session(self)
        return session.query(FITSFile).join(Operations, Operations.output_fits_id==FITSFile.id).filter(Operations.input_fits_id==self.id).all()


    @property
    def parents(self):
        session = object_session(self)
        return session.query(FITSFile).join(Operations, Operations.input_fits_id==FITSFile.id).filter(Operations.output_fits_id==self.id).all()


    def __init__(self, fname, path, size, md5, extensions):
        self.fname = fname
        self.path = path
        self.size = size
        self.md5 = md5
        self.extensions = extensions

    def __repr__(self):
        return "<FITS file ID {0:d} @ {1}>".format(self.id, self.full_path)

# standard decorator style
@event.listens_for(FITSFile, 'before_delete')
def receive_before_delete(mapper, connection, target):
    logger.info('Deleting FITS file {0}'.format(target.full_path))
    os.remove(target.full_path)


class Program(Base, CategoryBaseClass):
    __tablename__ = 'program'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    category_keyword = 'gemprgid'

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Program %s>' % self.name


class Object(Base, CategoryBaseClass):
    __tablename__ = 'object'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    ra = Column(Float)
    dec = Column(Float)
    description = Column(String)

    category_keyword = 'object'


    def __init__(self, name, ra=None, dec=None, description=None):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.description = description

    def __repr__(self):
        return '<Object %s>' % self.name


class ObservationBlock(Base, CategoryBaseClass):
    __tablename__ = 'observation_block'

    id = Column(Integer, primary_key=True)
    program_id = Column(Integer, ForeignKey('program.id'))
    name = Column(String)
    description = Column(String)

    program = relationship('Program')

    category_keyword = 'obsid'

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Observation Block %s>' % self.name


class ObservationType(Base, CategoryBaseClass):
    __tablename__ = 'observation_type'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    category_keyword = 'obstype'

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Observation Type %s>' % self.name


class ObservationClass(Base, CategoryBaseClass):
    __tablename__ = 'observation_class'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    category_keyword = 'obsclass'

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Gemini Observation Class %s>' % self.name


class Instrument(Base, CategoryBaseClass):
    __tablename__ = 'instrument'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    category_keyword = 'instrume'

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return "<Gemini Instrument %s>" % self.name


class Operations(Base):
    __tablename__ = 'operations'

    id = Column(Integer, primary_key=True)
    input_fits_id = Column(Integer, ForeignKey('fits_file.id'))
    operations_type_id = Column(Integer)
    operations_id = Column(Integer)
    output_fits_id = Column(Integer, ForeignKey('fits_file.id'))


    @classmethod
    def from_fits_objects(cls, input_fits_object, output_fits_object, operations_type_id, operations_id):
        current_operation = cls()
        current_operation.input_fits_object = input_fits_object
        current_operation.output_fits_object = output_fits_object
        current_operation.input_fits_id = input_fits_object.id
        current_operation.session = object_session(input_fits_object)
        current_operation.operations_id = operations_id
        current_operation.operations_type_id = operations_type_id
        return current_operation

    def commit(self):
        self.session.add(self.output_fits_object)
        self.session.commit()
        self.output_fits_id = self.output_fits_object.id
        self.session.add(self)
        self.session.commit()

class WaveCalType(Base):
    __tablename__ = 'wave_cal_type'

    #initialize with 0=guess, 1=arc, 2=sky


    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String, default=None)






