
from sqlalchemy import Column, Table, ForeignKey
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy.orm import relationship

from geminiutil.base import gemini_alchemy

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
