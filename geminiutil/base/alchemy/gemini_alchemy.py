from sqlalchemy import String, Integer, Float

from sqlalchemy.ext.declarative import declared_attr

from sqlalchemy.orm import relationship, object_session
from sqlalchemy import Column, ForeignKey

import logging

logger = logging.getLogger(__name__)

from geminiutil.base.alchemy.base import Base

from geminiutil.base.alchemy.category_alchemy import ObservationType, \
    ObservationClass, ObservationBlock, Instrument, Object, Program

from geminiutil.base.alchemy.file_alchemy import FITSFile

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



class AbstractGeminiRawFITS(Base):
    __abstract__ = True


    @declared_attr
    def id(cls):
        return Column(Integer, ForeignKey('fits_file.id'), primary_key=True)


    @declared_attr
    def instrument_id(cls):
        return Column(Integer, ForeignKey('instrument.id'))

    @declared_attr
    def observation_block_id(cls):
        return Column(Integer, ForeignKey('observation_block.id'))

    @declared_attr
    def observation_class_id(cls):
        return Column(Integer, ForeignKey('observation_class.id'))

    @declared_attr
    def observation_type_id(cls):
        return Column(Integer, ForeignKey('observation_type.id'))

    @declared_attr
    def object_id(cls):
        return Column(Integer, ForeignKey('object.id'))


    ###### RELATIONSHIPS #######
    @declared_attr
    def fits(cls):
        return relationship(FITSFile, uselist=False)

    @declared_attr
    def instrument(cls):
        return relationship(Instrument, uselist=False)

    @declared_attr
    def observation_block(cls):
        return relationship(ObservationBlock, uselist=False)

    @declared_attr
    def observation_class(cls):
        return relationship(ObservationClass, uselist=False)

    @declared_attr
    def observation_type(cls):
        return relationship(ObservationType, uselist=False)

    @declared_attr
    def object(cls):
        return relationship(Object, uselist=False)
