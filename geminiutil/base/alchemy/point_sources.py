from sqlalchemy import Column, Integer, String, Float, ForeignKey
from sqlalchemy.orm import relationship
from geminiutil.base.alchemy.base import Base


class PointSource(Base):
    """
    Table describing a point source in the database. Currently only supports
    name, ra, dec
    """

    __tablename__ = 'point_sources'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    ra = Column(Float)
    dec = Column(Float)

    magnitudes_id = Column(Integer, ForeignKey('point_source_magnitudes.id'))

    magnitudes = relationship('PointSourceMagnitude')


class PointSourceMagnitude(Base):

    __tablename__ = 'point_source_magnitudes'

    id = Column(Integer, primary_key=True)
    magnitude = Column(Float)
    point_source_id = Column(Integer)
    band_id = Column(Integer, ForeignKey('magnitude_bands.id'))


class MagnitudeBand(Base):
    __tablename__ = 'magnitude_bands'

    id = Column(Integer, primary_key=True)
    short_name = Column(String)
    long_name = Column(String)

class Field(Base):
    __tablename__ = 'fields'

    id = Column(Integer, primary_key=True)
    name = Column(String)


    @classmethod
    def from_name(cls, name, session):
        """
        Get Field object from name
        """

        name = name.lower()

        if session.query(cls).filter_by(name=name).count()>0:
            return session.query(cls).filter_by(name=name).one()
        else:
            new_field = cls(name=name)
            session.add(new_field)
            return new_field