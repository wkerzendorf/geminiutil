from geminiutil.gmos.alchemy.base import AbstractGMOSRawFITS
import tempfile
import os
import logging

logger = logging.getLogger(__name__)



class GMOSImagingRawFITS(AbstractGMOSRawFITS):
    __tablename__ = 'gmos_imaging_raw_fits'

    def reduce(self):
        """
        Prepare, reduce, mosaic FITS File - currently using IRAF, to be replaced by our own constructs

        """
        reduce_dir = os.path.join(self.fits.work_dir, 'reduced')
        if not os.path.exists(reduce_dir):
            os.mkdir(reduce_dir)

        reduce_fname = os.path.join(reduce_dir, 'red-{0}'.format(self.fits.fname))

        if os.path.exists(reduce_fname):
            logger.warn('{0} exists - deleting')
            os.system('rm {0}'.format(reduce_fname))

        from pyraf import iraf
        prepare_temp_fname = tempfile.NamedTemporaryFile().name
        reduce_temp_fname = tempfile.NamedTemporaryFile().name

        iraf.gemini()
        iraf.gmos()

        iraf.gprepare(self.fits.full_path, rawpath='', outimag=prepare_temp_fname)
        iraf.gireduce(inimages=prepare_temp_fname, outimag=reduce_temp_fname, fl_over=True,
                      fl_trim=True, fl_bias=False, fl_dark=False, fl_qeco=False, fl_flat=False)

        iraf.gmosaic(inimages=reduce_temp_fname, outimages=reduce_fname)

        return reduce_fname
        #os.system('rm {0}'.format(temporary_file_name))