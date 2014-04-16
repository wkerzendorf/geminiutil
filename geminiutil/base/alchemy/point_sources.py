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

    field_id = Column(Integer, ForeignKey('fields.id'))

    field = relationship('Field', uselist=False, backref='point_sources')

    def __repr__(self):
        return "<Point Source {0}>".format(self.name)



class PointSourceMagnitude(Base):

    __tablename__ = 'point_source_magnitudes'

    id = Column(Integer, primary_key=True)
    magnitude = Column(Float)
    point_source_id = Column(Integer, ForeignKey('point_sources.id'))
    band_id = Column(Integer, ForeignKey('photometry_bands.id'))

    point_source = relationship(PointSource, uselist=False, backref='magnitudes')
    band = relationship('PhotometryBand', uselist=False)

    def __repr__(self):
        return "<Magnitude {0}={1:2f}>".format(self.band.name, self.magnitude)


class PhotometryBand(Base):
    __tablename__ = 'photometry_bands'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    long_name = Column(String)


    @classmethod
    def from_name(cls, name, session):
        if session.query(cls).filter_by(name=name).count()>0:
            return session.query(cls).filter_by(name=name).one()
        else:
            new_photometry_band = cls(name=name)
            session.add(new_photometry_band)
            session.commit()
            return new_photometry_band

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
            session.commit()
            return new_field

    def __repr__(self):
        return "<Field {0}>".format(self.name)