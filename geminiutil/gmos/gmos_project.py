from .. import base
from ..base import BaseProject, ObservationClass, ObservationType
import geminiutil
# following not yet used
# from geminiutil.base.gemini_alchemy import AbstractFileTable
# from .gmos_alchemy import GMOSDetector
from .gmos_alchemy import (GMOSMOSRawFITS, GMOSMask, GMOSFilter,
                           GMOSGrating, GMOSMOSInstrumentSetup,
                           GMOSMOSScienceSet, GMOSArcLamp,
                           GMOSLongSlitArc)
import logging
from sqlalchemy import func

from astropy import time

import numpy as np
import re
import os


logger = logging.getLogger(__name__)

default_configuration_dir = os.path.join(geminiutil.__path__[0],
                                         'data', 'gmos')


class GMOSDBError(ValueError):
    pass


class GMOSMOSProject(BaseProject):
    """GMOS multi-object spectroscopy project.

    Example
    -------
    Start/restart a project with a given database that holds information
    about all observations

    >>> from geminiutil.gmos.gmos_project import GMOSMOSProject
    >>> proj = GMOSMOSProject('sqlite:///mcsnr.db3')

    First time, initialize to read in GMOS filter/grating information

    >>> proj.initialize()

    Add observations to the database

    >>> proj.add_directory('/raw/mhvk/gemini/mcsnr', file_filter='S*S*.fits')

    Now, can get sets of relevant files.  E.g.,

    >>> print proj.observation_classes
    (u'daycal', u'acq', u'science', u'partnercal', u'acqcal', u'progcal')

    >>> proj.observation_types
    (u'flat', u'object', u'arc')

    >>> science_list = proj.science
    >>> daycal_list = proj.daycal
    """

    def __init__(self, database_string, work_dir, echo=False):
        super(GMOSMOSProject, self).__init__(database_string, work_dir,
                                             GMOSMOSRawFITS, echo=echo)

    @property
    def science_sets(self):
        return self.session.query(GMOSMOSScienceSet).all()

    #showing the different setups for science exposures
    @property
    def science_instrument_setups(self):
        return (self.session.query(GMOSMOSInstrumentSetup)
                .join(GMOSMOSRawFITS)
                .join(ObservationClass)
                .filter(ObservationClass.name == 'science')
                .all())

    #showing the different setups for longslit_arcs
    @property
    def longslit_arcs_instrument_setups(self):
        return (self.session.query(GMOSMOSInstrumentSetup)
                .join(GMOSMOSRawFITS).join(ObservationType).join(GMOSMask)
                .filter(ObservationType.name == 'arc',
                        GMOSMask.name.like('%arcsec'))
                .all())

    def __getattr__(self, item):
        if(item.endswith('_query') and
           item.replace('_query', '') in self.observation_types):
            return (self.session.query(GMOSMOSRawFITS)
                    .join(ObservationType)
                    .filter(ObservationType.name ==
                            item.replace('_query', '')))
        elif item in self.observation_types:
            return self.__getattr__(item+'_query').all()

        elif (item.endswith('_query') and
              item.replace('_query', '') in self.observation_classes):
            return (self.session.query(GMOSMOSRawFITS)
                    .join(ObservationClass)
                    .filter(ObservationClass.name ==
                            item.replace('_query', '')))

        elif item in self.observation_classes:
            return self.__getattr__(item+'_query').all()
        else:
            return self.__getattribute__(item)

    def classify_added_fits(self, current_fits):
        """
        Classify added FITS objects

        Parameters
        ----------

        current_fits: FITSFile
            FITSFile object to classify
        """

        fits_object = self.add_gmos_raw_fits(current_fits)
        if fits_object is not None:
            return fits_object
        else:
            fits_object = self.add_gmos_mask(current_fits)

        if fits_object is None:
            logger.warning('Could not classify fits file %s',
                           current_fits.fname)
        return fits_object
        #self.add_gmos_mask(self)

    def add_gmos_raw_fits(self, fits_file):
        required_categories = [base.Object, base.Program,
                               base.ObservationBlock, base.ObservationClass,
                               base.ObservationType, base.Instrument]

        required_keywords = [item.category_keyword
                             for item in required_categories] + ['date-obs']

        if not all([keyword in fits_file.header
                    for keyword in required_keywords]):
            logger.debug("%s is not a normal raw gmos fits file",
                         fits_file.fname)
            return

        object = base.Object.from_fits_object(fits_file)
        # not yet used: program = base.Program.from_fits_object(fits_file)
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

    def add_gmos_mask(self, fits_object):
        required_keywords = ['GEMPRGID', 'OBSTYPE', 'ODFNAME']
        if(not all([keyword in fits_object.header
                    for keyword in required_keywords]) and
           fits_object.header['OBSTYPE'].lower().strip() != 'mask'):
            logger.debug("%s is not a normal gemini mask file" %
                         fits_object.fname)
            return None

        gmos_mask = GMOSMask.from_fits_object(fits_object)

        self.session.add(gmos_mask)
        self.session.commit()
        logger.info('Added GMOS MDF %s', fits_object.fname)
        return gmos_mask

    def link_masks(self):
        """For each MOS observation, link it to the corresponding mask file."""
        for gmos_raw in self.session.query(GMOSMOSRawFITS):
            slit_mask_pattern = re.compile('G[SN]\d{4}.+')
            longslit_mask_pattern = re.compile('(\d+)\.(\d+)arcsec')
            if gmos_raw.mask_id is not None:
                logger.debug('Mask is already set for %s - moving on',
                             gmos_raw.fits.fname)
                continue
            mask_name = gmos_raw.fits.header['maskname']

            if longslit_mask_pattern.match(mask_name) is not None:

                if self.session.query(GMOSMask).filter_by(
                        name=mask_name.strip().lower()).count() == 0:
                    program_id = self.session.query(base.Program).filter_by(
                        name=(gmos_raw.fits.header['GEMPRGID']
                              .lower().strip())).one().id
                    longslit_placeholder_mask = GMOSMask(mask_name, program_id)
                    self.session.add(longslit_placeholder_mask)
                    self.session.commit()

            elif slit_mask_pattern.match(mask_name) is None:
                logger.warn('%s (in %s) doesn\'t seem to be a valid maskname',
                            mask_name, gmos_raw.fits.fname)
                continue

            masks_found = self.session.query(GMOSMask).filter_by(
                name=mask_name.strip().lower()).count()
            if masks_found == 0:
                logger.critical('Mask %s is required by %s but does not exist'
                                'in database', mask_name, gmos_raw.fits.fname)
                continue
            elif masks_found > 1:
                logger.warn('Mask %s is duplicate in the database - '
                            ' please check', mask_name)
                continue
            else:
                logger.info('Linking %s with mask %s',
                            gmos_raw.fits.fname, mask_name)
                mask = self.session.query(GMOSMask).filter_by(
                    name=mask_name.strip().lower()).one()
                gmos_raw.mask_id = mask.id

        self.session.commit()

    def link_longslit_arcs(self):
        logger.info('Linking Longslit Arcs')

        for longslit_arc in self.arc_query.join(GMOSMask).filter(
                GMOSMask.name.like('%arcsec')).all():
            current_arc_lamp = self.session.query(GMOSArcLamp).filter_by(
                name=longslit_arc.fits.header['GCALLAMP'].lower()).one()
            gmos_longslit_arc = GMOSLongSlitArc(
                id=longslit_arc.id, arc_lamp_id=current_arc_lamp.id)
            self.session.add(gmos_longslit_arc)

        self.session.commit()

    def link_science_sets(self):
        """
        Link individual science observations to their calibration data

        Parameters
        ----------

        science_instrument2longslit_instrument: ~np.recarray
            a dictionary mapping science_frame instrument_setup_id to
            longslit_arc_instrument_setup_id
        """

        science_frames = (self.session.query(GMOSMOSRawFITS)
                          .join(ObservationType).join(ObservationClass)
                          .filter(ObservationClass.name == 'science',
                                  ObservationType.name == 'object')
                          .all())
        for science_frame in science_frames:
            flat = (self.session.query(GMOSMOSRawFITS)
                    .join(ObservationType)
                    .filter(ObservationType.name == 'flat',
                            GMOSMOSRawFITS.mask_id == science_frame.mask_id,
                            GMOSMOSRawFITS.observation_block_id ==
                            science_frame.observation_block_id)
                    .order_by(func.abs(GMOSMOSRawFITS.mjd -
                                       science_frame.mjd))
                    .all())

            #### WRONG - just picking one flat for now ###
            flat = flat[0]

            mask_arc = (self.session.query(GMOSMOSRawFITS)
                        .join(ObservationType).join(ObservationClass)
                        .filter(
                            ObservationType.name == 'arc',
                            GMOSMOSRawFITS.mask_id == science_frame.mask_id,
                            GMOSMOSRawFITS.instrument_setup_id ==
                            science_frame.instrument_setup_id)
                        .order_by(func.abs(GMOSMOSRawFITS.mjd -
                                           science_frame.mjd))
                        .all())

            if len(mask_arc) != 1:
                logger.warn('Science Frame {0} has more than one mask arc:\n'
                            '{1} - selecting closest arc'
                            .format(science_frame,
                                    '\n'.join(map(str, mask_arc))))
            mask_arc = mask_arc[0]

            self.session.add(GMOSMOSScienceSet(id=science_frame.id,
                                               flat_id=flat.id,
                                               mask_arc_id=mask_arc.id,))
            logger.info('Link Science Frame {0} with:\nFlat: {1}\n'
                        'Mask Arc: {2}\n'
                        .format(science_frame, flat, mask_arc))
        self.session.commit()

    def link_database(self):
        """
            Linking the raw files into meaningful constructs (i.e. masks, ...)
        """
        self.link_masks()
        self.link_longslit_arcs()
        self.link_science_sets()

    def initialize_database(self, configuration_dir=None):
        """Read in GMOS filter/grating information, for matching to headers."""
        if configuration_dir is None:
            configuration_dir = default_configuration_dir

        self._initialize_gmos_filters(configuration_dir)
        self._initialize_gmos_gratings(configuration_dir)
        self._initialize_gmos_arcs(configuration_dir)

    def _initialize_gmos_arcs(self, configuration_dir=None):
        logger.info('Reading Arc information')
        cuar_arc = GMOSArcLamp(name='cuar', fname='CuAr.dat',
                               path=os.path.join('gmos', 'arcs'))
        self.session.add(cuar_arc)
        self.session.commit()

    def _initialize_gmos_filters(self, configuration_dir):

        logger.info('Reading Filter information')

        gmos_filters = np.recfromtxt(
            os.path.join(configuration_dir, 'GMOSfilters.dat'),
            names=['name', 'wave_start', 'wave_end', 'fname'])
        for line in gmos_filters:
            new_filter = GMOSFilter(name=line['name'],
                                    wavelength_start_value=line['wave_start'],
                                    wavelength_start_unit='nm',
                                    wavelength_end_value=line['wave_end'],
                                    wavelength_end_unit='nm',
                                    fname=line['fname'],
                                    path=os.path.join('gmos', 'filter_data'))
            self.session.add(new_filter)

        open_filter = GMOSFilter(name='open', wavelength_start_value=0,
                                 wavelength_start_unit='nm',
                                 wavelength_end_value=np.inf,
                                 wavelength_end_unit='nm',
                                 fname=None, path=None)
        self.session.add(open_filter)

    def _initialize_gmos_gratings(self, configuration_dir):
        gmos_gratings = np.recfromtxt(
            os.path.join(configuration_dir, 'GMOSgratings.dat'),
            names=['name', 'ruling_density', 'blaze_wave', 'R', 'coverage',
                   'wave_start', 'wave_end', 'wave_offset', 'y_offset'])
        logger.info('Reading grating information')
        for line in gmos_gratings:
            new_grating = GMOSGrating(
                name=line['name'],
                ruling_density_value=line['ruling_density'],
                blaze_wavelength_value=line['blaze_wave'],
                R=line['R'],
                coverage_value=line['coverage'],
                wavelength_start_value=line['wave_start'],
                wavelength_end_value=line['wave_end'],
                wavelength_offset_value=line['wave_offset'],
                y_offset_value=line['y_offset'])
            self.session.add(new_grating)

        mirror = GMOSGrating(name='mirror',
                             ruling_density_value=0.0,
                             blaze_wavelength_value=0.0, R=0.0,
                             coverage_value=np.inf,
                             wavelength_start_value=0.0,
                             wavelength_end_value=np.inf,
                             wavelength_offset_value=0.0,
                             y_offset_value=0.0)
        self.session.add(mirror)
        self.session.commit()
