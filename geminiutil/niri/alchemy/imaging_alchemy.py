import logging

from astropy.io import fits

from sqlalchemy import Column, ForeignKey
from sqlalchemy import Integer
from sqlalchemy.orm import relationship


from geminiutil import base
from geminiutil.base.alchemy.file_alchemy import DataFile
from geminiutil.util import convert_camel_case2underscore
from geminiutil.base.alchemy.gemini_alchemy import AbstractGeminiRawFITS

from geminiutil.base.alchemy.base import Base

logger = logging.getLogger(__name__)

from astropy import time

class NIRIClassifyError(ValueError):
    pass


class NIRIImagingRawFits(AbstractGeminiRawFITS):
    __tablename__ = 'niri_imaging_raw_fits'

    categories = [base.Object, base.Program, base.ObservationBlock,
                            base.ObservationClass, base.ObservationType,
                            base.Instrument]


    id = Column(Integer, ForeignKey('temporary_fits_files.id'), primary_key=True)

    fits = relationship('TemporaryFITSFile', uselist=False)

    def __repr__(self):
        return '<niri id ={0:d} fits="{1}" class="{2}" type="{3}" object="{4}">'\
            .format(self.id, self.fits.fname, self.observation_class.name,
                    self.observation_type.name, self.object.name)





    @classmethod
    def verify_fits_class(cls, fname):
        """
        Class method to check if a FITS file matches the given class of FITS files
        This is often done by requiring a numner of different keywords

        Parameters
        ----------

        fname: str
            FITS filename
        """
        required_keywords = [item.category_keyword for item in cls.categories] + ['date-obs']
        fits_header = fits.getheader(fname)
        if not all([keyword in fits_header
                    for keyword in required_keywords]):
            raise NIRIClassifyError(
                "{0} is not a {1} fits file".format(fname, cls.__name__))



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

        niri_raw_image = cls()
        niri_raw_image.fits = fits_object
        niri_raw_image.object = base.Object.from_fits_object(fits_object)
        niri_raw_image.program = base.Program.from_fits_object(fits_object)
        niri_raw_image.observation_block = base.ObservationBlock.from_fits_object(fits_object)
        niri_raw_image.observation_class = base.ObservationClass.from_fits_object(fits_object)
        niri_raw_image.observation_type = base.ObservationType.from_fits_object(fits_object)
        niri_raw_image.instrument = base.Instrument.from_fits_object(fits_object)


        date_obs_str = '{0}T{1}'.format(fits_object.header['date-obs'],
                                        fits_object.header['time-obs'])

        niri_raw_image.mjd = time.Time(date_obs_str, scale='utc').mjd

        session.add(niri_raw_image)
        session.commit()

        return niri_raw_image


class TemporaryFITSFile(Base):
    __tablename__ = 'temporary_fits_files'

    id = Column(Integer, primary_key=True)
    data_file_id = Column(Integer, ForeignKey('data_files.id'))
    extensions = Column(Integer)

    data_file = relationship('DataFile', uselist=False)

    @property
    def fname(self):
        return self.data_file.fname


    @property
    def fits_data(self):
        return fits.open(self.data_file.full_path, ignore_missing_end=True)

    @property
    def header(self):
        return fits.getheader(self.data_file.full_path, ignore_missing_end=True)

    @property
    def data(self):
        return fits.getdata(self.data_file.full_path, ignore_missing_end=True)


    @property
    def shape(self):
        return self.header['naxis1'], self.header['naxis2']


    def __repr__(self):
        return "<FITS file ID {0:d} @ {1}>".format(self.id,
                                                   self.data_file.full_path)