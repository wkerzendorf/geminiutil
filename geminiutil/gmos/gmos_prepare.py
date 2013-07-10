"""Basic reduction class Prepare

Typical usage:
Prep = Prepare()   # can give bias_subslice, data_subslice, clip, gains, ...
infits = fits.open(fitsfile)
outfits = Prep(infits)

***TODO***
Ensure the out fits headers contain everything required to reproduce the result
"""

from geminiutil.base import FITSFile
from geminiutil.gmos.util import prepare_data, prepare_slices
from geminiutil.gmos import GMOSMOSPrepared
import logging
logger = logging.getLogger(__name__)
from sqlalchemy.orm import object_session
import os

class GMOSPrepareFrame(object):  # will be base when we know what
    __tablename__ = 'gmos_prepare'

    file_prefix = 'prep'

    def __init__(self, bias_subslice=None,
                 data_subslice=None,
                 bias_clip_sigma=3.):
        self.bias_subslice = bias_subslice
        self.data_subslice = data_subslice
        self.bias_clip_sigma = bias_clip_sigma

    def __call__(self, gmos_raw_object, fname=None, destination_dir='.'):
        """
        Preparing the Image

        Parameters
        ----------

        gmos_raw_object :

        fname :
            output filename

        write_steps :
            write out the individual steps to the individual fits files

        write_cut_image :
            fits file name to write out the cut_image to



        """

        if fname is None:
            fname = '%s-%s' % (self.file_prefix, gmos_raw_object.fits.fname)

        full_path = os.path.join(destination_dir, fname)

        fits_data = gmos_raw_object.fits.fits_data

        assert len(fits_data) - 1 == 3

        # subtract overscan, get useful part of detector, correct for gain,
        # and set read noise
        read_noise = [detector.readout_noise for detector in
                      gmos_raw_object.instrument_setup.detectors]
        gain = [detector.gain for detector in
                gmos_raw_object.instrument_setup.detectors]
        fits_file = prepare_data.prepare(fits_data,
                                    bias_subslice=self.bias_subslice,
                                    data_subslice=self.data_subslice,
                                    clip=self.bias_clip_sigma,
                                    gain=gain, read_noise=read_noise)
        # give each extension a name.  May later add error/mask extensions
        for i, extension in enumerate(fits_file):
            if i > 0:
                extension.name = 'chip{0:d}.data'.format(i)

        fits_file.writeto(full_path, clobber=True)

        # read it back in and add to database
        fits_file = FITSFile.from_fits_file(full_path)

        session = object_session(gmos_raw_object)
        session.add(fits_file)
        session.commit()
        gmos_mos_prepared = GMOSMOSPrepared(id=fits_file.id,
                                            raw_fits_id=gmos_raw_object.id)
        session.add(gmos_mos_prepared)
        session.commit()
        return gmos_mos_prepared



class GMOSPrepareScienceSet(object):
    def __init__(self, prepare_function):
        self.prepare = prepare_function

    def __call__(self, science_set):

        for item in ['science', 'flat', 'mask_arc']:

            prepared_fits = self.prepare(getattr(science_set, item))

        slices = prepare_slices.calculate_slice_geometries(science_set)
        session = object_session(science_set)
        session.add_all(slices)
        session.commit()

