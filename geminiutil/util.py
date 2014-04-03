import re
import os
import logging
import hashlib
import shutil

from glob import glob

from astropy.io import fits


logger = logging.getLogger(__name__)


http_header_disposition_gz = re.compile('^inline; filename=(.+)\.gz$')

def convert_camel_case2underscore(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

def hashfile(afile, hasher, blocksize=65536):
    """
    Create MD5 hash out of open file

    Parameters
    ----------

    afile: open file handle
    hasher: hash object
    blocksize: int
        number of bytes to be read simultaneously

    """
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.hexdigest()

def get_request(cadc_url, username, password):
    try:
        import requests
    except ImportError:
        raise ImportError('The package requests is required for downloading files.')

    url_request = requests.get(cadc_url, auth=(username, password), stream=True)
    try:
        download_fname = http_header_disposition_gz.match(url_request.headers['content-disposition']).groups()[0]
    except AttributeError:
        raise ValueError('Download not a gzip file')
    else:
        uncompressed_md5 = url_request.headers['x-uncompressed-md5']

    assert url_request.headers['content-encoding'] == 'gzip'

    return download_fname, uncompressed_md5, url_request


def download_cadc_url(cadc_url, username, password, chunk_size=1024):

    #check if exists and move on if it does

    download_fname, uncompressed_md5, url_request = get_request(cadc_url,
                                                                username,
                                                                password)

    logger.info('File {0} already exists - skip download'.format(download_fname))


    logger.info("Downloading file {0} from url {1}".format(download_fname,
                                                           cadc_url))


    #writing to disk
    with open(download_fname, 'w') as local_fh:
        for chunk in url_request.iter_content(chunk_size):
            if chunk:
                local_fh.write(chunk)

    if uncompressed_md5 == hashfile(download_fname, hashlib.md5()):
        raise IOError('File {0} MD5 mismatch with downloaded version'.format(download_fname))

    return download_fname



def download_cadc_list(fname, username, password, directory_dict=None, chunk_size=1024):
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

    if directory_dict is not None:
        all_files_full_path = []
        for path in directory_dict.values():
            all_files_full_path += glob(os.path.join(path, '*'))
        all_files = map(os.path.basename, all_files_full_path)
    else:
        all_files = []

    with open(fname) as cadc_fh:
        for line in cadc_fh:
            cadc_url = line.strip('\n\r')
            (download_fname, uncompressed_md5, url_request)\
                = get_request(cadc_url, username, password)

            if download_fname in all_files:
                i = all_files.index(download_fname)
                logger.info('File {0} already downloaded -- skipping'
                            .format(all_files_full_path[i]))
                continue

            download_fname = download_cadc_url(cadc_url, username, password)

            if directory_dict is not None:
                instrument = fits.getval(download_fname, 'instrume')
                shutil.move(download_fname, directory_dict[instrument])




