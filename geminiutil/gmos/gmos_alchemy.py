from ..base import GeminiRawFITS, Base

from sqlalchemy import Column, ForeignKey

from sqlalchemy.orm import relationship, backref

#sqlalchemy types
from sqlalchemy import String, Integer, Float, DateTime, Boolean

class GMOSMDF(Base):
    __tablename__ = 'mdf'

    id = Column(Integer, ForeignKey('base_fits.id'))
    name = Column(String)
    program_id = Column(Integer, ForeignKey('program.id'))


    def __init__(self, name, program_id):
        self.name = name
        self.program_id = program_id


class GMOSMOSRawFITS(GeminiRawFITS, Base):
    __tablename__ = 'gmos_mos_raw_fits'


    id = Column(Integer, ForeignKey('base_fits.id'), primary_key=True)
    date_obs = Column(DateTime)
    instrument_id = Column(Integer, ForeignKey('instrument.id'))
    observation_block_id = Column(Integer, ForeignKey('observation_block.id'))
    observation_class_id = Column(Integer, ForeignKey('observation_class.id'))
    observation_type_id = Column(Integer, ForeignKey('observation_type.id'))
    exclude = Column(Boolean)
    mask_id = Column(Integer, ForeignKey('mdf.id'))


    fits = relationship(FitsFile, uselist=False, backref='raw_fits')
    instrument = relationship(Instrument, uselist=False, backref='raw_fits')
    observation_block = relationship(ObservationBlock, uselist=False, backref='raw_fits')
    observation_class = relationship(ObservationClass, uselist=False, backref='raw_fits')
    observation_type = relationship(ObservationType, uselist=False, backref='raw_fits')

    def __init(self, date_obs, instrument_id, observation_block_id, observation_class_id, observation_type_id, mask_id, object_id, exclude=False,):

        super(GeminiRawFITS, self).__init(self, date_obs, instrument_id, observation_block_id, observation_class_id, observation_type_id, exclude=False)
        self.mask_id = mask_id
        self.object_id = object_id

