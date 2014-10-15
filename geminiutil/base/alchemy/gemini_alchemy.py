import logging

from sqlalchemy import String, Integer
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship, object_session
from sqlalchemy import Column, ForeignKey

logger = logging.getLogger(__name__)

from geminiutil.base.alchemy.base import Base

from geminiutil.base.alchemy.category_alchemy import ObservationType, \
    ObservationClass, ObservationBlock, Instrument, Object

from geminiutil.base.alchemy.file_alchemy import FITSFile, FITSClassifyError


from astropy.io import fits


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

class AbstractGeminiRawFITS(Base):
    __abstract__ = True


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
            raise FITSClassifyError(
                "{0} is not a {1} fits file".format(fname, cls.__name__))


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

