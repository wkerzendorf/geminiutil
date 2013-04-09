from .. import base
from ..base import BaseProject
from .gmos_alchemy import GMOSMOSRawFITS, GMOSMask
import logging
from datetime import datetime

import re

logger = logging.getLogger(__name__)


class GMOSMOSProject(BaseProject):

    def __init__(self, database_string):
        super(GMOSMOSProject, self).__init__(database_string, GMOSMOSRawFITS)



    @property
    def science(self):
        return self.session.query(base.ObservationClass).filter_by(name='science').one().raw_fits


    def classify_added_fits(self, current_fits):
        fits_object = self.add_gmos_raw_fits(current_fits)
        if fits_object is not None:
            return fits_object
        else:
            fits_object = self.add_gmos_mask(current_fits)

        if fits_object is None:
            logger.warning('Could not classify fits file %s', current_fits.fname)
        return fits_object
        #self.add_gmos_mask(self)

    def add_gmos_raw_fits(self, fits_file):
        required_categories = [base.Object, base.Program, base.ObservationBlock, base.ObservationClass, base.ObservationType,
                               base.Instrument]

        required_keywords = [item.category_keyword for item in required_categories] + ['date-obs']

        if not all([keyword in fits_file.header for keyword in required_keywords]):
            logger.debug("%s is not a normal raw gmos fits file", fits_file.fname)
            return

        object = base.Object.from_fits_object(fits_file)
        program = base.Program.from_fits_object(fits_file)
        observation_block = base.ObservationBlock.from_fits_object(fits_file)
        observation_class = base.ObservationClass.from_fits_object(fits_file)
        observation_type = base.ObservationType.from_fits_object(fits_file)
        instrument = base.Instrument.from_fits_object(fits_file)


        date_obs_str = '%sT%s' % (fits_file.header['date-obs'], fits_file.header['time-obs'])
        date_obs = datetime.strptime(date_obs_str, '%Y-%m-%dT%H:%M:%S.%f')

        gmos_raw = self.raw_fits_class(date_obs, instrument.id, observation_block.id, observation_class.id, observation_type.id, object.id)

        gmos_raw.id = fits_file.id

        self.session.add(gmos_raw)
        self.session.commit()

        return gmos_raw




    def add_gmos_mask(self, fits_object):
        required_keywords = ['GEMPRGID', 'OBSTYPE', 'ODFNAME']
        if not all([keyword in fits_object.header for keyword in required_keywords]) and \
                        fits_object.header['OBSTYPE'].lower().strip() != 'mask':

            logger.debug("%s is not a normal gemini mask file" % fits_object.fname)
            return None

        gmos_mask = GMOSMask.from_fits_object(fits_object)

        self.session.add(gmos_mask)
        self.session.commit()
        logger.info('Added GMOS MDF %s', fits_object.fname)
        return gmos_mask

    def link_masks(self):
        for gmos_raw in self.session.query(GMOSMOSRawFITS):
            if gmos_raw.mask_id is not None:
                logger.debug('Mask is already set for %s - moving on', gmos_raw.fits.fname)
                continue
            mask_name = gmos_raw.fits.header['maskname']
            if re.match('G[SN]\d{4}.+', mask_name) is None:
                logger.warn('%s (in %s) doesn\'t seem to be a valid maskname', mask_name, gmos_raw.fits.fname)
                continue
            masks_found = self.session.query(GMOSMask).filter_by(name=mask_name.strip().lower()).count()
            if masks_found == 0:
                logger.critical('Mask %s is required by %s but does not exist in database', mask_name, gmos_raw.fits.fname)
                continue
            elif masks_found > 1:
                logger.warn('Mask %s is duplicate in the database - please check', mask_name)
                continue
            else:
                logger.info('Linking with %s with mask %s', gmos_raw.fits.fname, mask_name)
                mask = self.session.query(GMOSMask).filter_by(name=mask_name.strip().lower()).one()
                gmos_raw.mask_id = mask.id

        self.session.commit()









