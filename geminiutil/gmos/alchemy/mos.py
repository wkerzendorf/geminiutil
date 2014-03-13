
from sqlalchemy import Column, Table, ForeignKey
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy.orm import relationship

from geminiutil import base


class MOSPointSource(base.Base):
    """
    This Table establishes a link with the point sources and the slices. It serves as an association table
    for the many-to-many relationship between the basic PointSource Table and the GMOSMOSSlice table

    It stores what point source is found on which slit and where.
    """

    __tablename__ = "mos_object"

    id = Column(Integer, primary_key=True)
    point_source_id = Column(Integer, ForeignKey('point_source.id'))
    slice_id = Column(Integer, ForeignKey('gmos_mos_slices.id'))
    slit_pos_y = Column(Float)

    point_source = relationship(base.PointSource, backref='mos_point_source')

