import logging

from astropy.io import fits

from sqlalchemy.orm import relationship


from geminiutil import base
from geminiutil.util import convert_camel_case2underscore



logger = logging.getLogger(__name__)


class NIRIClassifyError(ValueError):
    pass

class NIRIImagingRawFits(Base):

        __tablename__ = 'niri_imaging_raw_fits'

        categories = [base.Object, base.Program, base.ObservationBlock,
                                base.ObservationClass, base.ObservationType,
                                base.Instrument]

        @classmethod
        def initialize_table(cls):
            for category in cls.categories:
                relationship_name = convert_camel_case2underscore(category.__name__)
                cls.__setattr__(relationship_name, relationship(category,
                                                                uselist=False))

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
                logger.debug("{0} is not a {1} fits file".format(fname, cls.__name__))
                raise NIRIClassifyError



        @classmethod
        def from_fits_file(cls, fname, session):
            """
            generate a

            """
            fits_object = base.FITSFile.from_file(fname)
            session.add(fits_object)
            session.commit()
            return cls.from_fits_object(fits_object)

        @classmethod
        def from_fits_object(cls, fits_object):

            object = base.Object.from_fits_object(fits_file)
            program = base.Program.from_fits_object(fits_file)
            observation_block = base.ObservationBlock.from_fits_object(fits_file)
            observation_class = base.ObservationClass.from_fits_object(fits_file)
            observation_type = base.ObservationType.from_fits_object(fits_file)
            instrument = base.Instrument.from_fits_object(fits_file)

        if len(fits_file.fits_data) == 4 or len(fits_file.fits_data) == 7:
            instrument_setup_id = GMOSMOSInstrumentSetup.from_fits_object(
                fits_file).id
        else:
            logger.warn('Unusual fits data with %d HDUs '
                        '(expecting either 4 or 7)', len(fits_file.fits_data))
            instrument_setup_id = None

        date_obs_str = '%sT%s' % (fits_file.header['date-obs'],
                                  fits_file.header['time-obs'])
        mjd = time.Time(date_obs_str, scale='utc').mjd

        gmos_raw = self.raw_fits_class(
            mjd=mjd, instrument_id=instrument.id,
            observation_block_id=observation_block.id,
            observation_class_id=observation_class.id,
            observation_type_id=observation_type.id, object_id=object.id,
            instrument_setup_id=instrument_setup_id)

        gmos_raw.id = fits_file.id

        self.session.add(gmos_raw)
        self.session.commit()

        return gmos_raw