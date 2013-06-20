"""Basic reduction classes GMOSPrepare and GMOSCutSlits

***TODO***
Ensure the out fits headers contain everything required to reproduce the result
"""

from .. base import FITSFile
# import astropy.io.fits as fits
import logging
import os
from geminiutil.gmos.gmos_alchemy import GMOSMOSPrepared
from sqlalchemy.orm import object_session
from .basic import prepare, mask_cut

logger = logging.getLogger(__name__)


class GMOSPrepare(object):  # will be base when we know what
    __tablename__ = 'gmos_prepare'

    file_prefix = 'prep'

    def __init__(self, bias_subslice=None,
                 data_subslice=None,
                 bias_clip_sigma=3.):
        self.bias_subslice = bias_subslice
        self.data_subslice = data_subslice
        self.bias_clip_sigma = bias_clip_sigma

    def __call__(self, gmos_raw_object, fname=None, destination_dir='.',
                 write_steps=False, write_cut_image=None):
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
        fits_file = prepare.prepare(fits_data,
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


class GMOSPrepareScienceFrame(GMOSPrepare):
    def __call__(self, science_frame, **kwargs):
        prepare = super(GMOSPrepareScienceFrame, self).__call__
        if science_frame.raw_fits.prepared_fits is None:
            prepared_science = prepare(science_frame.raw_fits, **kwargs)

        if science_frame.flat.prepared_fits is None:
            prepared_flat = prepare(science_frame.flat, **kwargs)

        if science_frame.mask_arc.prepared_fits is None:
            prepared_mask_arc = prepare(science_frame.mask_arc, **kwargs)

        slices = mask_cut.calculate_slice_geometries(science_frame)
        session = object_session(science_frame)
        session.add_all(slices)
        session.commit()
        return (prepared_science, prepared_flat, prepared_mask_arc,
                science_frame.slices)


def prepare_science_frame(science_frame, destination_dir='.'):
    prepare = GMOSPrepare()

    if science_frame.raw_fits.prepared_fits is None:
        prepare(science_frame.raw_fits, destination_dir=destination_dir)

    if science_frame.flat.prepared_fits is None:
        prepare(science_frame.flat, destination_dir=destination_dir)

    if science_frame.mask_arc.prepared_fits is None:
        prepare(science_frame.mask_arc, destination_dir=destination_dir)

    slices = mask_cut.calculate_slice_geometries(science_frame)
    session = object_session(science_frame)
    session.add_all(slices)
    session.commit()
