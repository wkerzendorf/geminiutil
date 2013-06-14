"""Basic reduction classes GMOSPrepare and GMOSCutSlits

***TODO***
Ensure the out fits headers contain everything required to reproduce the result
"""

from .. base import FITSFile
import astropy.io.fits as fits
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

        final_hdu_list = [fits_data[0].copy()]

        # make sure we only have data for 3 fits files.
        # This does not yet work for GMOS-N
        assert len(fits_data) - 1 == 3

        for i in xrange(1, len(fits_data)):
            current_amplifier = fits_data[i]
            detector = gmos_raw_object.instrument_setup.detectors[i - 1]

            #####
            # Subtracting Overscan
            #####
            amplifier_data = prepare.correct_overscan(current_amplifier,
                                                      self.bias_subslice,
                                                      self.data_subslice,
                                                      self.bias_clip_sigma)
            #####
            #Correcting the amplifier data gain
            #####

            amplifier_data = prepare.correct_gain(amplifier_data,
                                                  gain=detector.gain)
            amplifier_data.name = 'chip%d.data' % i

            final_hdu_list += [amplifier_data]

            fits.HDUList(final_hdu_list).writeto(full_path, clobber=True)

        fits_file = FITSFile.from_fits_file(full_path)
        session = object_session(gmos_raw_object)
        session.add(fits_file)
        session.commit()
        gmos_mos_prepared = GMOSMOSPrepared(id=fits_file.id,
                                            raw_fits_id=gmos_raw_object.id)
        session.add(gmos_mos_prepared)
        session.commit()
        return gmos_mos_prepared



def prepare_science_frame(science_frame, destination_dir='.'):
    prepare = GMOSPrepare()

    if science_frame.raw_fits.prepared_fits is None:
        prepared_science = prepare(science_frame.raw_fits, destination_dir=destination_dir)

    if science_frame.flat.prepared_fits is None:
        prepared_flat = prepare(science_frame.flat, destination_dir=destination_dir)

    if science_frame.mask_arc.prepared_fits is None:
        prepared_mask_arc = prepare(science_frame.mask_arc, destination_dir=destination_dir)

    slices = mask_cut.calculate_slice_geometries(science_frame)
    session = object_session(science_frame)
    session.add_all(slices)
    session.commit()



