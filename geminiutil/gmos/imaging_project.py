from geminiutil.gmos.gmos_project import GMOSProject

import geminiutil
from geminiutil.gmos.alchemy.imaging_raw import GMOSImagingRawFITS
from geminiutil.base.base_project import BaseProject

from geminiutil.gmos.util.read_standard_fields import read_standard_star_db

from geminiutil.base.alchemy.point_sources import Field, PointSourceMagnitude,\
    PointSource, PhotometryBand

import os

default_configuration_dir = os.path.join(geminiutil.__path__[0],
                                         'data', 'gmos')


class GMOSImagingProject(GMOSProject):
    """
    TBD
    """

    def __init__(self, database_string, work_dir, echo=False):
        super(GMOSImagingProject, self).__init__(database_string, work_dir,
                                             GMOSImagingRawFITS, echo=echo)





    def add_fits_file(self, fname):
        """
        Add FITS file to database

        Parameters
        ----------

        fname: str
            FITS file name/path to add to database
        """

        GMOSImagingRawFITS.verify_fits_class(fname)
        return GMOSImagingRawFITS.from_fits_file(fname, self.session)

    def initialize_standard_stars(self):
        standard_star_db_fname = os.path.join(default_configuration_dir,
                                              'standards', 'gmosnlandolt.dat')
        standard_star_db = read_standard_star_db(standard_star_db_fname)

        for id, line in standard_star_db.iterrows():
            stds_data = line
            name = stds_data.pop('name')
            field_name = stds_data.pop('fields')
            ra = stds_data.pop('ra')
            dec = stds_data.pop('dec')

            field_object = Field.from_name(field_name, self.session)
            point_source_object = PointSource(name=name, ra=ra, dec=dec)
            self.session.add(point_source_object)

            for band in stds_data.keys():
                band_object = PhotometryBand.from_name(band, self.session)
                magnitude_object = PointSourceMagnitude(magnitude=stds_data[band])
                magnitude_object.band = band_object
                point_source_object.magnitudes.append(magnitude_object)


            point_source_object.field = field_object

            self.session.commit()


            
    def initialize_database(self, configuration_dir=None):
        self.initialize_standard_stars()
        super(GMOSImagingProject, self).initialize_database(
            configuration_dir=configuration_dir)




