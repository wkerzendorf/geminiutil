from geminiutil.base.gemini_alchemy import Base, FITSFile, Instrument, Program, \
    ObservationBlock, ObservationClass, ObservationType, AbstractFileTable

from geminiutil.base import gemini_alchemy as base_alchemy

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from datetime import datetime
from glob import glob
import os
import logging

import re

http_header_disposition_fitsgz = re.compile('^inline; filename=(.+\.fits)\.gz')
http_header_disposition_gz = re.compile('^inline; filename=(.+)\.gz')
http_header_disposition = re.compile('^inline; filename=(.+)\.gz')

logger = logging.getLogger(__name__)


def get_category(session, category_str, category):
    if session.query(category).filter_by(name=category_str).count() == 0:
        current_category = category(category_str)
        session.add(current_category)
        session.commit()
    else:
        current_category = session.query(category).filter_by(name=category_str).one()

    return current_category

class BaseProject(object):

    def __init__(self, database_string, work_dir, raw_fits_class=FITSFile, echo=False):
        self.metadata = Base.metadata
        self.engine = create_engine(database_string, echo=echo)
        self.metadata.bind = self.engine
        self.metadata.create_all()
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.conn = self.session.bind.connect()

        self.raw_fits_class = raw_fits_class

        if not os.path.exists(work_dir):
            raise ValueError('Working directory {0} does not exist'.format(work_dir))

        AbstractFileTable.work_dir = work_dir

        self.work_dir = work_dir


    @property
    def observation_classes(self):
        """Names of observation classes in the database (e.g., 'daycal')."""
        return [r.name for r in self.session.query(base_alchemy.ObservationClass.name)]

    @property
    def observation_types(self):
        """Names of observation types in the database (e.g., 'science')."""
        return [r.name for r in self.session.query(base_alchemy.ObservationType.name)]

    def __dir__(self):
        return self.__dict__.keys() + list(self.observation_classes) + list(self.observation_types)

    def download_raw_fits(self, fname, username, password, raw_directory='raw', chunk_size=1024):
        """
        Download files from a "cadcUrlList.txt" into the database.


        Parameters
        ----------

        fname: str
            name of the url list file

        username: str
            CADC username

        password: str
            CADC password

        raw_directory: str
            relative path to the raw fits directory (default='raw')

        chunk_size: int
            Number of bytes to read in one chunk when downloading (default=1024)
        """

        if not os.path.exists(os.path.join(self.work_dir, raw_directory)):
            print "warn"
            logger.warn('Raw directory {0} does not exist - creating'.format(os.path.join(self.work_dir, raw_directory)))
            os.mkdir(os.path.join(self.work_dir, raw_directory))


        try:
            import requests
        except ImportError:
            raise ImportError('The package requests is required for downloading files.')



        with open(fname) as cadc_fh:
            for line in cadc_fh:
                url_request = requests.get(line.strip(), auth=(username, password), stream=True)
                try:
                    fits_fname = http_header_disposition_fitsgz.match(url_request.headers['content-disposition']).groups()[0]

                except AttributeError:
                    try:
                        download_fname = http_header_disposition_gz.match(url_request.headers['content-disposition']).groups()[0]
                    except AttributeError:
                        download_fname = http_header_disposition.match(url_request.headers['content-disposition']).groups()[0]
                    finally:
                        file_type = 'other'
                else:
                    file_type = 'fits'


                if file_type == 'fits':
                    assert url_request.headers['content-encoding'] == 'gzip'
                    assert url_request.headers['content-type'] == 'application/fits'

                    #check if exists and move on if it does
                    if self.session.query(FITSFile).filter_by(fname=fits_fname).count() > 0:
                        current_fits = self.session.query(FITSFile).filter_by(fname=fits_fname).one()
                        assert url_request.headers['x-uncompressed-md5'] == current_fits.md5
                        logger.info('File {0} already exists - skip download'.format(fits_fname))
                        continue

                logger.info("Downloading file {0} from url {1}".format(fits_fname, line.strip('\r\n')))

                #writing to disk
                with open(os.path.join(self.work_dir, raw_directory, fits_fname), 'w') as local_fh:
                    for chunk in url_request.iter_content(chunk_size):
                        if chunk:
                            local_fh.write(chunk)

                if file_type == 'fits':
                    # Add to DB
                    current_fits = self.add_fits_file(os.path.join(raw_directory, fits_fname))

                    if url_request.headers['x-uncompressed-md5'] != current_fits.md5:
                        raise IOError('File {0} MD5 mismatch with downloaded version'.format(fits_fname))


    def add_directory(self, directory, file_filter='*.fits'):
        """
        Adding directory to the da
        """
        for fname in sorted(glob(os.path.join(directory, file_filter))):
            current_fits = self.add_fits_file(fname)
            self.classify_added_fits(current_fits)

    def add_fits_file(self, fname):
        """
        Add FITS file to database

        Parameters
        ----------

        fname: str
            FITS file name/path to add to database
        """

        logger.info('Adding {0} to project'.format(fname))
        current_fits = FITSFile.from_fits_file(fname)
        self.session.add(current_fits)
        self.session.commit()

        return current_fits


    def classify_added_fits(self, current_fits):
        fits_object = self.add_gemini_raw_fits(current_fits)
        return fits_object



    def add_gemini_raw_fits(self, fits_file):
        necessary_keywords = ['instrume', 'object', 'obstype', 'obsclass', 'gemprgid', 'obsid', 'date-obs', 'time-obs']

        #print "Working on %s" % fits_file.fname

        instrument_str = fits_file.header['instrume'].lower().strip()
        current_instrument = get_category(self.session, instrument_str, Instrument)

        program_str = fits_file.header['gemprgid'].lower().strip()
        current_program = get_category(self.session, program_str, Program)

        observation_str = fits_file.header['obsid'].lower().strip()
        current_observation_block = get_category(self.session, observation_str, ObservationBlock)

        if current_observation_block.program_id is None:
            current_observation_block.program_id = current_program.id
            self.session.commit()

        obstype_str = fits_file.header['obstype'].lower().strip()
        current_observation_type = get_category(self.session, obstype_str, ObservationType)

        obsclass_str = fits_file.header['obsclass'].lower().strip()
        current_observation_class = get_category(self.session, obsclass_str, ObservationClass)

        date_obs_str = '%sT%s' % (fits_file.header['date-obs'], fits_file.header['time-obs'])
        date_obs = datetime.strptime(date_obs_str, '%Y-%m-%dT%H:%M:%S.%f')

        current_raw_fits = raw_fits_class(date_obs, current_instrument.id, current_observation_block.id, current_observation_class.id, current_observation_type.id)
        current_raw_fits.id = fits_file.id

        self.session.add(current_raw_fits)

        self.session.commit()
        return current_raw_fits

















