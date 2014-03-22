from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy import event
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, relationship, backref, object_session
from sqlalchemy import Column, ForeignKey

from astropy.io import fits
import logging
import os


logger = logging.getLogger(__name__)

#



from geminiutil.base.alchemy.base import Base



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


class PointSource(Base):
    """
    Table describing a point source in the database. Currently only supports
    name, ra, dec
    """

    __tablename__ = 'point_sources'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    ra = Column(Float)
    dec = Column(Float)


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






