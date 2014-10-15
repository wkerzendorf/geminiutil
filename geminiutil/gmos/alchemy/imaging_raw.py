from sqlalchemy import ForeignKey, Column
from sqlalchemy import Integer, Float, func, Table
from sqlalchemy.orm import relationship



from geminiutil.base.alchemy.base import Base

from geminiutil.base.alchemy.file_alchemy import FITSClassifyError
from geminiutil.base.alchemy.gemini_alchemy import AbstractGeminiRawFITS
from geminiutil.base.alchemy.file_alchemy import DataFile, TemporaryFITSFile
from geminiutil.base.alchemy.category_alchemy import Instrument, Object, \
    Program, ObservationBlock, ObservationType, ObservationClass
from geminiutil.gmos.alchemy.base import GMOSFilter, GMOSGrating
from geminiutil.gmos.util import qmosaic, prepare_frame

import logging

from astropy import units as u
from astropy.io import fits
from astropy import time

logger = logging.getLogger(__name__)

class GMOSImagingFITSClassifyError(FITSClassifyError):
    pass

class GMOSImagingRawFITS(AbstractGeminiRawFITS):
    __tablename__ = 'gmos_imaging_raw_fits'


    categories = [Object, Program,
                  ObservationBlock,
                  ObservationClass,
                  ObservationType,
                  Instrument]


    id = Column(Integer, ForeignKey('temporary_fits_files.id'), primary_key=True)

    instrument_setup_id = Column(Integer,
                                 ForeignKey('gmos_imaging_instrument_setup.id'))

    fits = relationship('TemporaryFITSFile', uselist=False)
    instrument_setup = relationship('GMOSImagingInstrumentSetup')
    mjd = Column(Float)

    @property
    def datetime(self):
        return time.Time(self.mjd, format='mjd', scale='utc').datetime

    def __repr__(self):
        return '<gmos id ={0:d} fits="{1}" class="{2}" type="{3}" object="{4}">'\
            .format(self.id, self.fits.fname, self.observation_class.name,
                    self.observation_type.name, self.object.name)

    def reduce(self):
        prepared_frame = prepare_frame.prepare(self.fits.fits_data)
        return qmosaic.qmosaic(*prepared_frame[1:4])

    @classmethod
    def from_fits_file(cls, fname, session):
        """
        generate a

        """
        data_file_object = DataFile.from_file(fname)
        fits_object = TemporaryFITSFile()
        fits_object.data_file = data_file_object

        session.add(fits_object)
        session.commit()

        gmos_raw_image = cls()
        gmos_raw_image.fits = fits_object
        gmos_raw_image.object = Object.from_fits_object(fits_object)
        gmos_raw_image.program = Program.from_fits_object(fits_object)
        gmos_raw_image.observation_block = ObservationBlock.from_fits_object(fits_object)
        gmos_raw_image.observation_class = ObservationClass.from_fits_object(fits_object)
        gmos_raw_image.observation_type = ObservationType.from_fits_object(fits_object)
        gmos_raw_image.instrument = Instrument.from_fits_object(fits_object)
        gmos_raw_image.instrument_setup = GMOSImagingInstrumentSetup.from_fits_file(fname, session)

        date_obs_str = '{0}T{1}'.format(fits_object.header['date-obs'],
                                        fits_object.header['time-obs'])

        gmos_raw_image.mjd = time.Time(date_obs_str, scale='utc').mjd

        session.add(gmos_raw_image)
        session.commit()

        return gmos_raw_image


class GMOSImagingInstrumentSetup(Base):
    __tablename__ = 'gmos_imaging_instrument_setup'

    id = Column(Integer, primary_key=True)

    filter1_id = Column(Integer, ForeignKey('gmos_filters.id'))
    filter2_id = Column(Integer, ForeignKey('gmos_filters.id'))

    grating_id = Column(Integer, ForeignKey('gmos_gratings.id'))

    grating_slit_wavelength_value = Column(Float)
    grating_slit_wavelength_unit = u.Unit('nm')

    @property
    def grating_slit_wavelength(self):
        return self.grating_slit_wavelength_value *\
               self.grating_slit_wavelength_unit

    grating_central_wavelength_value = Column(Float)
    grating_central_wavelength_unit = u.Unit('nm')

    @property
    def grating_central_wavelength(self):
        return self.grating_central_wavelength_value *\
               self.grating_central_wavelength_unit

    grating_tilt_value = Column(Float)
    grating_tilt_unit = u.Unit('degree')

    grating_order = Column(Integer)

    instrument_id = Column(Integer, ForeignKey('instrument.id'))


    filter1 = relationship(GMOSFilter, primaryjoin=(GMOSFilter.id==filter1_id),
                                uselist=False)

    filter2 = relationship(GMOSFilter, primaryjoin=(GMOSFilter.id==filter2_id),
                                uselist=False)

    grating = relationship(GMOSGrating)

    instrument = relationship(Instrument)



    @classmethod
    def from_fits_file(cls, fits_file, session,
                       tilt_equivalency_threshold=0.001,
                       wavelength_equivalency_threshold=0.0001):
        header = fits.getheader(fits_file)

        filter1 = GMOSFilter.from_keyword(header['filter1'], session)
        filter2 = GMOSFilter.from_keyword(header['filter2'], session)

        grating = GMOSGrating.from_keyword(header['grating'], session)
        instrument = Instrument.from_keyword(header['instrume'], session)

        grating_central_wavelength = header['centwave']
        grating_slit_wavelength = header['grwlen']

        grating_tilt = header['grtilt']
        grating_order = header['grorder']


        #Checking if the same instrument setup already exists
        instrument_setup_query = session.query(cls).filter(
            cls.filter1_id==filter1.id, cls.filter2_id==filter2.id,
            cls.grating_id==grating.id, cls.instrument_id==instrument.id,
            func.abs(cls.grating_central_wavelength_value -
                     grating_central_wavelength)
            < wavelength_equivalency_threshold,
            func.abs(cls.grating_slit_wavelength_value -
                      grating_slit_wavelength) <
            wavelength_equivalency_threshold,
            func.abs(cls.grating_tilt_value -
                      grating_tilt) < tilt_equivalency_threshold)

        instrument_setup_object = instrument_setup_query.first()

        if instrument_setup_object is None:
            instrument_setup_object = cls(
                filter1_id=filter1.id,
                filter2_id=filter2.id,
                grating_id=grating.id,
                grating_central_wavelength_value=grating_central_wavelength,
                grating_slit_wavelength_value=grating_slit_wavelength,
                grating_tilt_value=grating_tilt,
                grating_order=grating_order,
                instrument_id=instrument.id)

            session.add(instrument_setup_object)
            session.commit()

        return instrument_setup_object




    def __getattr__(self, item):
        if item in ['grating_slit_wavelength', 'grating_central_wavelength', 'grating_tilt']:
            item_value = getattr(self, '%s_value' % item)
            item_unit = getattr(self, '%s_unit' % item)
            return u.Quantity(item_value, item_unit)
        else:
            return self.__getattribute__(item)


    def __repr__(self):
        return "<GMOS MOS Instrument Setup ID=%s Filter1 %s Filter2 %s Grating %s Tilt %.2f central wave=%.2f %s>" % \
                (self.id, self.filter1, self.filter2, self.grating, self.grating_tilt_value, self.grating_central_wavelength_value,
                self.grating_central_wavelength_unit)

"""
collection2raw_image = Table(
    'gmos_imaging_collection2raw_image', Base.metadata,
    Column('gmos_imaging_science_set_id', Integer,
           ForeignKey('gmos_imaging_science_set')),
    Column('gmos_imaging_raw_fits_id', Integer,
           ForeignKey('gmos_imaging_raw_fits.id')))

class GMOSImagingCollection(Base):
    __tablename__ = 'gmos_imaging_science_sets'

    id = Column(Integer, primary_key=True)
    raw_fits = relationship(GMOSImagingRawFITS, secondary=science_set2raw_image,
                            backref='science_set')



"""


