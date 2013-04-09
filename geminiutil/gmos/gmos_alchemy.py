from ..base import Base, FITSFile, Instrument, ObservationType, ObservationClass, ObservationBlock, Object

from .. import base

from sqlalchemy import Column, ForeignKey

from sqlalchemy.orm import relationship, backref, object_session


from astropy.utils import misc

#sqlalchemy types
from sqlalchemy import String, Integer, Float, DateTime, Boolean

class GMOSMask(Base):
    __tablename__ = 'gmos_mask'

    id = Column(Integer, ForeignKey('fits_file.id'), primary_key=True)
    name = Column(String)
    program_id = Column(Integer, ForeignKey('program.id'))

    fits = relationship(base.FITSFile)

    @misc.lazyproperty
    def table(self):
        return self.fits.data

    @classmethod
    def from_fits_object(cls, fits_object):
        session = object_session(fits_object)
        mask_name = fits_object.header['DATALAB'].lower().strip()
        mask_program = session.query(base.Program).filter_by(name=fits_object.header['GEMPRGID'].lower().strip()).one()
        mask_object = cls(mask_name, mask_program.id)
        mask_object.id = fits_object.id
        return mask_object



    def __init__(self, name, program_id):
        self.name = name
        self.program_id = program_id


class GMOSMOSRawFITS(Base):
    __tablename__ = 'gmos_mos_raw_fits'


    id = Column(Integer, ForeignKey('fits_file.id'), primary_key=True)
    date_obs = Column(DateTime)
    instrument_id = Column(Integer, ForeignKey('instrument.id'))
    observation_block_id = Column(Integer, ForeignKey('observation_block.id'))
    observation_class_id = Column(Integer, ForeignKey('observation_class.id'))
    observation_type_id = Column(Integer, ForeignKey('observation_type.id'))
    object_id = Column(Integer, ForeignKey('object.id'))
    mask_id = Column(Integer, ForeignKey('gmos_mask.id'))
    exclude = Column(Boolean)


    fits = relationship(FITSFile, uselist=False, backref='raw_fits')
    instrument = relationship(Instrument, uselist=False, backref='raw_fits')
    observation_block = relationship(ObservationBlock, uselist=False, backref='raw_fits')
    observation_class = relationship(ObservationClass, uselist=False, backref='raw_fits')
    observation_type = relationship(ObservationType, uselist=False, backref='raw_fits')
    object = relationship(base.Object, uselist=False, backref='raw_fits')
    mask = relationship(GMOSMask, uselist=False, backref='raw_fits')

    def __init__(self, date_obs, instrument_id, observation_block_id, observation_class_id, observation_type_id, object_id, mask_id=None, exclude=False):
        self.date_obs = date_obs
        self.instrument_id = instrument_id
        self.observation_block_id = observation_block_id
        self.observation_class_id = observation_class_id
        self.observation_type_id = observation_type_id
        self.exclude = exclude

        self.mask_id = mask_id
        self.object_id = object_id

    def __repr__(self):
        return '<gmos fits="%s" class="%s" type="%s" object="%s">' % (self.fits.fname, self.observation_class.name,
                                                                  self.observation_type.name, self.object.name)

