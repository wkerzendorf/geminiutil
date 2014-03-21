from sqlalchemy import Column, Table, ForeignKey
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy.orm import relationship

from geminiutil.base.alchemy import gemini_alchemy
from geminiutil.base.alchemy.file_alchemy import AbstractFileTable

# association table
mos_point_source2mos_slice = Table('mos_point_source2point_source',
                                         gemini_alchemy.Base.metadata,
                                         Column('mos_point_source_id', Integer, ForeignKey('mos_point_sources.id')),
                                         Column('mos_slice_id', Integer, ForeignKey('gmos_mos_slices.id')))

class MOSPointSource(gemini_alchemy.Base):
    """
    This table describes a point source in the term of a slice in a mask. It contains a slit position and a priority,

    for the many-to-many relationship between the basic PointSource Table and the GMOSMOSSlice table

    It stores what point source is found on which slit and where.
    """

    __tablename__ = "mos_point_sources"

    id = Column(Integer, ForeignKey('point_sources.id'), primary_key=True)
    slit_position = Column(Float)
    priority = Column(Integer)

    point_source = relationship(gemini_alchemy.PointSource, backref='mos_point_source')

    slices = relationship('GMOSMOSSlice', secondary=mos_point_source2mos_slice, backref='mos_point_sources')


    def __repr__(self):
        return '<MosPointSouce id={0} priority={1} slit position={2:.4f}>'.format(self.id, self.priority, self.slit_position)





mos_spectrum2mos_slice = Table('mos_spectrum2mos_point_source',
                                      gemini_alchemy.Base.metadata,
                                         Column('mos_spectrum_id',
                                                Integer,
                                                ForeignKey('mos_spectra.id')),
                                         Column('mos_slice_id',
                                                Integer,
                                                ForeignKey('gmos_mos_slices.id')
                                         ))

class MOSSpectrum(AbstractFileTable):
    """
    This table holds Extracted Spectra from MOS Slices
    """

    __tablename__ = "mos_spectra"

    id = Column(Integer, primary_key=True)
    mos_point_source_id = Column(Integer, ForeignKey('mos_point_sources.id'))

    slices = relationship('GMOSMOSSlice', secondary=mos_spectrum2mos_slice,
                          backref='mos_spectra')

    mos_point_source = relationship(MOSPointSource, uselist=False,
                                    backref='mos_spectra')



