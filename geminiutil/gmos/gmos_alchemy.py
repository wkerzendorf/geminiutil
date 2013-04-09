from ..base import GeminiRawFITS, Base

from sqlalchemy import Column, ForeignKey

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

    mdf_id = Column(Integer, ForeignKey('mdf.id'))


class FITS2MDF(Base):
    __tablename__ = 'fits2mdf'


    id = Column(Integer, ForeignKey('raw_fits_file.id'), primary_key=True)
    mdf_id = Column(Integer, ForeignKey('mdf.id'))

