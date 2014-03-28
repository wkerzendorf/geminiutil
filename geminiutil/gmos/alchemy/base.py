from geminiutil.base.alchemy.base import Base
from geminiutil.base.alchemy import gemini_alchemy
from geminiutil.base.alchemy.file_alchemy import FITSFile, \
    AbstractCalibrationFileTable


from sqlalchemy import Column, Table, ForeignKey
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy.orm import relationship, object_session

from sqlalchemy.ext.declarative import declared_attr

from astropy import time
from astropy.utils import misc
from astropy import units as u

import os

import numpy as np

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

    fits = relationship(FITSFile)

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

class AbstractGMOSRawFITS(gemini_alchemy.AbstractGeminiRawFITS):
    """
    abstract class for gmos raw fits frames with mask relations and so on

    more info to be added
    """
    __abstract__ = True

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
        return relationship(FITSFile, uselist=False)

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

    def __repr__(self):
        return '<gmos id ={0:d} fits="{1}" class="{2}" type="{3}" object="{4}">'.format(self.id,
                                                                                        self.fits.fname,
                                                                                        self.observation_class.name,
                                                                                        self.observation_type.name,
                                                                                        self.object.name)



class GMOSFilter(AbstractCalibrationFileTable):
    __tablename__ = 'gmos_filters'

    name = Column(String)
    wavelength_start_value = Column(Float)
    wavelength_start_unit = Column(String)

    wavelength_end_value = Column(Float)
    wavelength_end_unit = Column(String)


    @classmethod
    def from_keyword(cls, filter_keyword, session):
        if filter_keyword.startswith('open'):
            return session.query(GMOSFilter).filter_by(name='open').one()
        else:
            return session.query(GMOSFilter).filter_by(name=filter_keyword).one()


    @property
    def wavelength_start(self):
        return u.Quantity(self.wavelength_start_value, self.wavelength_start_unit)

    @property
    def wavelength_end(self):
        return u.Quantity(self.wavelength_end_value, self.wavelength_end_unit)

    @misc.lazyproperty
    def wavelength(self):
        return np.loadtxt(self.full_path, usecols=(0,))

    @misc.lazyproperty
    def flux(self):
        return np.loadtxt(self.full_path, usecols=(1,))


    def __repr__(self):
        return "<GMOS Filter %s>" % self.name

    def __str__(self):
        return self.name


class GMOSGrating(Base):
    __tablename__ = 'gmos_gratings'

    id = Column(Integer, primary_key=True)

    name = Column(String)

    ruling_density_value = Column(Float)
    ruling_density_unit = u.Unit('1/mm') # lines/mm

    blaze_wavelength_value = Column(Float)
    blaze_wavelength_unit = u.Unit('nm')

    R = Column(Float)

    coverage_value = Column(Float)
    coverage_unit = u.Unit('nm')

    wavelength_start_value = Column(Float)
    wavelength_start_unit = u.Unit('nm')

    wavelength_end_value = Column(Float)
    wavelength_end_unit = u.Unit('nm')

    wavelength_offset_value = Column(Float)
    wavelength_offset_unit = u.Unit('nm')

    y_offset_value = Column(Float)
    y_offset_unit = u.Unit('pix')

    @classmethod
    def from_keyword(cls, grating_keyword, session):
        if grating_keyword.lower() == 'mirror':
            return session.query(GMOSGrating).filter_by(name='mirror').one()
        else:
            return session.query(GMOSGrating).filter_by(
                name=grating_keyword).one()



    def __getattr__(self, item):
        if item in ['ruling_density', 'blaze_wavelength', 'coverage', 'wavelength_start', 'wavelength_end',
                    'wavelength_offset', 'y_offset']:
            item_value = getattr(self, '%s_value' % item)
            item_unit = getattr(self, '%s_unit' % item)
            return u.Quantity(item_value, item_unit)
        else:
            raise AttributeError('%s has no attribute %s' % (self.__class__.__name__, item))

    def __repr__(self):
        return "<GMOS Grating %s>" % self.name

    def __str__(self):
        return self.name
