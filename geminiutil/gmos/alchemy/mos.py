
import logging

from sqlalchemy import Column, Table, ForeignKey
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy.orm import relationship

from geminiutil.base.alchemy import gemini_alchemy
from geminiutil.base.alchemy.file_alchemy import AbstractFileTable, DataPathMixin, DataFile

import os

from astropy import table
from astropy import units as u

logger = logging.getLogger(__name__)

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





class MOSSpectrum(gemini_alchemy.Base, DataPathMixin):
    """
    This table holds Extracted Spectra from MOS Slices
    """

    __tablename__ = "mos_spectra"

    id = Column(Integer, ForeignKey('data_files.id'), primary_key=True)
    mos_point_source_id = Column(Integer, ForeignKey('mos_point_sources.id'))
    slice_id = Column(Integer, ForeignKey('gmos_mos_slices.id'))

    slice = relationship('GMOSMOSSlice', uselist=False,
                          backref='mos_spectra')

    mos_point_source = relationship(MOSPointSource, uselist=False,
                                    backref='mos_spectra')

    data_file = relationship(DataFile, uselist=False)

    hdf5_path = 'spectrum'

    @classmethod
    def from_table(cls, table, slice, mos_point_source, relative_path='spectra',
                   fname=None):
        """
        class method to add a new spectrum given a table from the extracted
        spectrum.

        Parameters
        ----------

        table: astropy.table.Table

        mos_point_source: MosPointSource

        relative_path: str
            relative path to work dir (default='spectra')

        fname: str
            filename of the extracted spectrum (default None)
            if None it automatically generates
            <raw_fits_fname>_slc<slc_id>_src<src_id>

        """

        if not os.path.exists(os.path.join(cls.work_dir, relative_path)):
            logger.warn('Relative dir {0} does not exist for adding spectrum '
                        '- creating'.format(relative_path))
            os.mkdir(os.path.join(cls.work_dir, relative_path))

        if fname is None:
            fname = '{0}_slc{1}_src{2}.h5'.format(
                slice.science_set.science.fits.fname.replace('.fits', ''),
                slice.id, mos_point_source.id)

        full_fname = os.path.join(cls.work_dir, relative_path, fname)

        table.write(full_fname, path=cls.hdf5_path, format='hdf5',
                    overwrite=True)

        data_file = DataFile.from_file(os.path.join(relative_path, fname))

        spectrum = cls(mos_point_source_id=mos_point_source.id,
                       slice_id=slice.id)

        spectrum.data_file = data_file

        return data_file, spectrum

    @property
    def table(self):
        return table.Table.read(self.data_file.full_path, path=self.hdf5_path,
                                format='hdf5')

    @property
    def wavelength(self):
        return self.table['wave'].T.flatten() * u.Angstrom

    @property
    def flux(self):
        return self.table['src'].T.flatten()
