from geminiutil.base import gemini_alchemy


from sqlalchemy import Column, Table, ForeignKey
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy.orm import relationship, object_session


from astropy import time

import os


import logging

logger = logging.getLogger(__name__)

class GMOSDatabaseDuplicate(Exception):
    pass

class GMOSNotPreparedError(Exception):
    #raised if a prepared fits is needed but none found for certain tasks
    pass


class AbstractGMOSRawFITS(gemini_alchemy.Base):
    __abstract__ = True

    id = Column(Integer, ForeignKey('fits_file.id'), primary_key=True)
    mjd = Column(Float)
    instrument_id = Column(Integer, ForeignKey('instrument.id'))
    observation_block_id = Column(Integer, ForeignKey('observation_block.id'))
    observation_class_id = Column(Integer, ForeignKey('observation_class.id'))
    observation_type_id = Column(Integer, ForeignKey('observation_type.id'))
    object_id = Column(Integer, ForeignKey('object.id'))
    mask_id = Column(Integer, ForeignKey('gmos_mask.id'))


    instrument_setup_id = Column(Integer, ForeignKey('gmos_mos_instrument_setup.id'))

    exclude = Column(Boolean)


    fits = relationship(gemini_alchemy.FITSFile, uselist=False, backref='raw_fits')
    instrument = relationship(gemini_alchemy.Instrument, uselist=False, backref='raw_fits')
    observation_block = relationship(gemini_alchemy.ObservationBlock, uselist=False, backref='raw_fits')
    observation_class = relationship(gemini_alchemy.ObservationClass, uselist=False, backref='raw_fits')
    observation_type = relationship(gemini_alchemy.ObservationType, uselist=False, backref='raw_fits')
    object = relationship(gemini_alchemy.Object, uselist=False, backref='raw_fits')
    mask = relationship(gemini_alchemy.GMOSMask, uselist=False, backref='raw_fits')

    instrument_setup = relationship('GMOSMOSInstrumentSetup', backref='raw_fits')


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


