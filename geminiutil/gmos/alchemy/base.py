from geminiutil.base.alchemy import gemini_alchemy


from sqlalchemy import Column, Table, ForeignKey
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy.orm import relationship, object_session

from sqlalchemy.ext.declarative import declared_attr

from astropy import time
from astropy.utils import misc

import os


import logging

logger = logging.getLogger(__name__)

class GMOSDatabaseDuplicate(Exception):
    pass

class GMOSNotPreparedError(Exception):
    #raised if a prepared fits is needed but none found for certain tasks
    pass

class GMOSMask(gemini_alchemy.Base):
    __tablename__ = 'gmos_mask'


    # GMOSMasks don't have the fits_id as the primary key anymore.
    #The reason for this change is that for a longslit exposure there exists no mask exposure.

    id = Column(Integer, primary_key=True)
    fits_id = Column(Integer, ForeignKey('fits_file.id'), default=None)
    name = Column(String)
    program_id = Column(Integer, ForeignKey('program.id'))

    fits = relationship(gemini_alchemy.FITSFile)

    @misc.lazyproperty
    def table(self):
        return self.fits.data

    @classmethod
    def from_fits_object(cls, fits_object):
        session = object_session(fits_object)
        mask_name = fits_object.header['DATALAB'].lower().strip()
        mask_program = session.query(gemini_alchemy.Program).filter_by(name=fits_object.header['GEMPRGID'].lower().strip()).one()
        mask_object = cls(mask_name, mask_program.id)
        mask_object.fits_id = fits_object.id
        return mask_object



    def __init__(self, name, program_id):
        self.name = name
        self.program_id = program_id

class AbstractGMOSRawFITS(gemini_alchemy.Base):
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

    @declared_attr
    def mask_id(cls):
        return Column(Integer, ForeignKey('gmos_mask.id'))

    @declared_attr
    def instrument_setup_id(cls):
        return Column(Integer, ForeignKey('gmos_mos_instrument_setup.id'))

    exclude = Column(Boolean)
    mjd = Column(Float)

    @declared_attr
    def fits(cls):
        return relationship(gemini_alchemy.FITSFile, uselist=False)

    @declared_attr
    def instrument(cls):
        return relationship(gemini_alchemy.Instrument, uselist=False)

    @declared_attr
    def observation_block(cls):
        return relationship(gemini_alchemy.ObservationBlock, uselist=False)

    @declared_attr
    def observation_class(cls):
        return relationship(gemini_alchemy.ObservationClass, uselist=False)

    @declared_attr
    def observation_type(cls):
        return relationship(gemini_alchemy.ObservationType, uselist=False)

    @declared_attr
    def object(cls):
        return relationship(gemini_alchemy.Object, uselist=False)

    @declared_attr
    def mask(cls):
        return relationship(GMOSMask, uselist=False)

    @declared_attr
    def instrument_setup(cls):
        return relationship('GMOSMOSInstrumentSetup')


    @property
    def associated_query(self):
        session = object_session(self)
        return session.query(self.__class__).filter_by(observation_block_id=self.observation_block_id)

    @property
    def associated(self):
        return self.associated_query.all()

    @property
    def date_obs(self):
        return time.Time(self.mjd, scale='utc', format='mjd')

    def __init__(self, mjd, instrument_id, observation_block_id, observation_class_id, observation_type_id,
                 object_id, mask_id=None, instrument_setup_id=None, exclude=False):
        self.mjd = mjd
        self.instrument_id = instrument_id
        self.observation_block_id = observation_block_id
        self.observation_class_id = observation_class_id
        self.observation_type_id = observation_type_id
        self.exclude = exclude

        self.mask_id = mask_id
        self.object_id = object_id
        self.instrument_setup_id = instrument_setup_id

    def __repr__(self):
        return '<gmos id ={0:d} fits="{1}" class="{2}" type="{3}" object="{4}">'.format(self.id,
                                                                                        self.fits.fname,
                                                                                        self.observation_class.name,
                                                                                        self.observation_type.name,
                                                                                        self.object.name)



