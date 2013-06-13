
import os
import yaml

from ..base import Base, FITSFile, Instrument, ObservationType, ObservationClass, ObservationBlock, Object

from .. import base

from scipy import interpolate
from sqlalchemy import Column, ForeignKey

from sqlalchemy.orm import relationship, backref, object_session
from sqlalchemy import func

from astropy.utils import misc
from astropy import units

from astropy import time

import numpy as np

detector_yaml_fname = os.path.join(os.path.dirname(__file__), 'data', 'gmos_detector_information.yml')
detector_information = yaml.load(file(detector_yaml_fname))

grating_eq, tilt = np.loadtxt(os.path.join(os.path.dirname(__file__), 'data', 'gratingeq.dat'), unpack=True)

# in the simplest case this is m*lambda = (1/ruling_density) * (sin(beta) + sin(alpha))
grating_eq_sorting = np.argsort(grating_eq)
grating_equation_interpolator = interpolate.interp1d(grating_eq[grating_eq_sorting], tilt[grating_eq_sorting])
import numpy as np

#sqlalchemy types
from sqlalchemy import String, Integer, Float, DateTime, Boolean

class GMOSMask(Base):
    __tablename__ = 'gmos_mask'

    id = Column(Integer, ForeignKey('fits_file.id'), primary_key=True)
    name = Column(String)
    program_id = Column(Integer, ForeignKey('program.id'))

    fits = relationship(base.FITSFile)

    @misc.lazyproperty
    def table(self):
        return self.fits.data

    @classmethod
    def from_fits_object(cls, fits_object):
        session = object_session(fits_object)
        mask_name = fits_object.header['DATALAB'].lower().strip()
        mask_program = session.query(base.Program).filter_by(name=fits_object.header['GEMPRGID'].lower().strip()).one()
        mask_object = cls(mask_name, mask_program.id)
        mask_object.id = fits_object.id
        return mask_object



    def __init__(self, name, program_id):
        self.name = name
        self.program_id = program_id


class GMOSDetector(Base):
    __tablename__ = 'gmos_detector'

    """
    Detector properties table
    """

    id = Column(Integer, primary_key=True)
    naxis1 = Column(Integer)
    naxis2 = Column(Integer)
    ccd_name = Column(String)
    readout_direction = Column(String)
    gain = Column(Float)
    readout_noise = Column(Float)
    x_binning = Column(Integer)
    y_binning = Column(Integer)
    frame_id = Column(Integer)
    instrument_id = Column(Integer, ForeignKey('instrument.id'))

    instrument = relationship(Instrument)

    @classmethod
    def from_fits_object(cls, fits_object, ccd_no):
        hdu = fits_object.fits_data[ccd_no]

        header = hdu.header
        session = object_session(fits_object)

        x_binning, y_binning = map(int, header['CCDSUM'].split())

        readout_direction = header['ampname'].split(',')[1].strip()

        instrument_id = Instrument.from_fits_object(fits_object).id

        detector_object = session.query(cls).filter(cls.naxis1==header['NAXIS1'], cls.naxis2==header['NAXIS2'],
                                  cls.ccd_name==header['CCDNAME'], cls.readout_direction==readout_direction,
                                  (func.abs(cls.gain - header['GAIN']) / header['GAIN']) < 0.0001,
                                  (func.abs(cls.readout_noise - header['RDNOISE']) / header['RDNOISE']) < 0.0001,
                                  cls.x_binning==x_binning, cls.y_binning==y_binning,
                                  cls.frame_id==int(header['FRAMEID']),
                                  cls.instrument_id==instrument_id).all()
        if detector_object == []:
            detector_object = cls(header['NAXIS1'], header['NAXIS2'], header['CCDNAME'], readout_direction, header['GAIN'],
                       header['RDNOISE'], x_binning, y_binning, header['FRAMEID'], instrument_id)
            session.add(detector_object)
            session.commit()
            return detector_object
        elif len(detector_object) == 1:
            return detector_object[0]
        else:
            raise ValueError('Found more than one detectors')

    def __init__(self, naxis1, naxis2, ccd_name, readout_direction, gain, read_noise, x_binning, y_binning, frame_id, instrument_id):
        self.naxis1 = naxis1
        self.naxis2 = naxis2
        self.ccd_name = ccd_name
        self.readout_direction = readout_direction
        self.gain = gain
        self.readout_noise = read_noise
        self.x_binning = x_binning
        self.y_binning = y_binning
        self.frame_id = frame_id
        self.instrument_id = instrument_id


    @misc.lazyproperty
    def pixel_scale(self):
        if self.ccd_name.startswith('EEV'):
            return detector_information['pixel_scale'][self.instrument.name]['eev'] * units.Unit('arcsec/pixel')
        else:
            raise NotImplemented('CCD %s not implemented yet' % self.ccd_name)

    @misc.lazyproperty
    def spectral_cutoff(self):
        if self.ccd_name.startswith('EEV'):
            return detector_information['spectral_cutoff']['eev'] * units.Unit('nm')
        else:
            raise NotImplemented('CCD %s not implemented yet' % self.ccd_name)

    def __repr__(self):
        return "<detector id=%d ccdname=%s xbin=%d ybin=%d gain=%.2f>" % (self.id, self.ccd_name, self.x_binning,
        self.y_binning, self.gain)





class GMOSFilter(Base):
    __tablename__ = 'gmos_filters'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    fname = Column(String)
    path = Column(String)
    wavelength_start_value = Column(Float)
    wavelength_start_unit = Column(String)

    wavelength_end_value = Column(Float)
    wavelength_end_unit = Column(String)


    @property
    def full_path(self):
        return os.path.join(self.path, self.fname)

    @property
    def wavelength_start(self):
        return units.Quantity(self.wavelength_start_value, self.wavelength_start_unit)

    @property
    def wavelength_end(self):
        return units.Quantity(self.wavelength_end_value, self.wavelength_end_unit)

    @misc.lazyproperty
    def wavelength(self):
        return np.loadtxt(self.full_path, usecols=(0,))

    @misc.lazyproperty
    def flux(self):
        return np.loadtxt(self.full_path, usecols=(1,))


    def __repr__(self):
        return "<GMOS Filter %s>" % self.name

    def __str__(self):
        return self.name


class GMOSGrating(Base):
    __tablename__ = 'gmos_gratings'

    id = Column(Integer, primary_key=True)

    name = Column(String)

    ruling_density_value = Column(Float)
    ruling_density_unit = units.Unit('1/mm') # lines/mm

    blaze_wavelength_value = Column(Float)
    blaze_wavelength_unit = units.Unit('nm')

    R = Column(Float)

    coverage_value = Column(Float)
    coverage_unit = units.Unit('nm')

    wavelength_start_value = Column(Float)
    wavelength_start_unit = units.Unit('nm')

    wavelength_end_value = Column(Float)
    wavelength_end_unit = units.Unit('nm')

    wavelength_offset_value = Column(Float)
    wavelength_offset_unit = units.Unit('nm')

    y_offset_value = Column(Float)
    y_offset_unit = units.Unit('pix')


    def __getattr__(self, item):
        if item in ['ruling_density', 'blaze_wavelength', 'coverage', 'wavelength_start', 'wavelength_end',
                    'wavelength_offset', 'y_offset']:
            item_value = getattr(self, '%s_value' % item)
            item_unit = getattr(self, '%s_unit' % item)
            return units.Quantity(item_value, item_unit)
        else:
            raise AttributeError('%s has no attribute %s' % (self.__class__.__name__, item))

    def __repr__(self):
        return "<GMOS Grating %s>" % self.name

    def __str__(self):
        return self.name

class GMOSMOSInstrumentSetup2Detector(Base):
    __tablename__ = 'gmos_mos_instrument_setup2gmosdetector'

    id = Column(Integer, primary_key=True)
    instrument_setup_id = Column(Integer, ForeignKey('gmos_mos_instrument_setup.id'))
    detector_id = Column(Integer, ForeignKey('gmos_detector.id'))
    fits_extension_id = Column(Integer)



class GMOSMOSInstrumentSetup(Base):
    __tablename__ = 'gmos_mos_instrument_setup'

    id = Column(Integer, primary_key=True)

    filter1_id = Column(Integer, ForeignKey('gmos_filters.id'))
    filter2_id = Column(Integer, ForeignKey('gmos_filters.id'))

    grating_id = Column(Integer, ForeignKey('gmos_gratings.id'))

    grating_slit_wavelength_value = Column(Float)
    grating_slit_wavelength_unit = units.Unit('nm')

    grating_central_wavelength_value = Column(Float)
    grating_central_wavelength_unit = units.Unit('nm')

    grating_tilt_value = Column(Float)
    grating_tilt_unit = units.Unit('degree')

    grating_order = Column(Integer)

    instrument_id = Column(Integer, ForeignKey('instrument.id'))

    #instrument_setup2detector_id = Column(Integer, ForeignKey('gmos_mos_instrument_setup2gmosdetector.id'))

    filter1 = relationship(GMOSFilter, primaryjoin=(GMOSFilter.id==filter1_id),
                                uselist=False)

    filter2 = relationship(GMOSFilter, primaryjoin=(GMOSFilter.id==filter2_id),
                                uselist=False)

    grating = relationship(GMOSGrating)

    instrument = relationship(base.Instrument)


    @classmethod
    def from_fits_object(cls, fits_object, equivalency_threshold = 0.0001):
        session = object_session(fits_object)
        header = fits_object.fits_data[0].header
        filter1 = header['filter1']
        filter2 = header['filter2']

        if filter1.startswith('open'):
            filter1_id = session.query(GMOSFilter).filter_by(name='open').one().id
        else:
            filter1_id = session.query(GMOSFilter).filter_by(name=filter1).one().id

        if filter2.startswith('open'):
            filter2_id = session.query(GMOSFilter).filter_by(name='open').one().id
        else:
            filter2_id = session.query(GMOSFilter).filter_by(name=filter2).one().id

        grating = header['grating']
        if grating.lower() == 'mirror':
            grating_id = session.query(GMOSGrating).filter_by(name='mirror').one().id
        else:
            grating_id = session.query(GMOSGrating).filter_by(name=header['grating']).one().id

        grating_central_wavelength = header['centwave']
        grating_slit_wavelength = header['grwlen']

        grating_tilt = header['grtilt']
        grating_order = header['grorder']

        instrument_id = base.Instrument.from_fits_object(fits_object).id

        instrument_setup_object = session.query(cls).filter(cls.filter1_id==filter1_id, cls.filter2_id==filter2_id,
            cls.grating_id==grating_id, cls.instrument_id==instrument_id,
            (func.abs(cls.grating_central_wavelength_value - grating_central_wavelength)
                                                    / grating_central_wavelength) < equivalency_threshold,
            (func.abs(cls.grating_slit_wavelength_value - grating_slit_wavelength)
                                                    / grating_slit_wavelength) < equivalency_threshold,
            (func.abs(cls.grating_tilt_value - grating_tilt)
                                                    / grating_tilt) < equivalency_threshold).all()
        if instrument_setup_object == []:
            instrument_setup_object = cls(filter1_id, filter2_id, grating_id, grating_central_wavelength,
                                          grating_slit_wavelength, grating_tilt, grating_order, instrument_id)

            session.add(instrument_setup_object)
            session.commit()

            for fits_extension_id in xrange(1, len(fits_object.fits_data)):
                current_detector_id = GMOSDetector.from_fits_object(fits_object, fits_extension_id).id
                session.add(GMOSMOSInstrumentSetup2Detector(instrument_setup_id=instrument_setup_object.id,
                                                            fits_extension_id=fits_extension_id,
                                                            detector_id=current_detector_id))
            session.commit()



            return instrument_setup_object

        elif len(instrument_setup_object) == 1:
            return instrument_setup_object[0]

        else:
            raise ValueError('More than one Instrument setup with the same setup found: %s' % instrument_setup_object)


    @property
    def detectors(self):
        session = object_session(self)
        detectors = session.query(GMOSDetector).join(GMOSMOSInstrumentSetup2Detector).\
            filter(GMOSMOSInstrumentSetup2Detector.instrument_setup_id==self.id).all()
        return detectors

    def __init__(self, filter1_id, filter2_id, grating_id, grating_central_wavelength_value,
                 grating_slit_wavelength_value, grating_tilt_value,
                 grating_order, instrument_id):
        self.filter1_id = filter1_id
        self.filter2_id = filter2_id
        self.grating_id = grating_id
        self.grating_central_wavelength_value = grating_central_wavelength_value
        self.grating_slit_wavelength_value = grating_slit_wavelength_value
        self.grating_tilt_value = grating_tilt_value
        self.grating_order = grating_order
        self.instrument_id = instrument_id

    def __getattr__(self, item):
        if item in ['grating_slit_wavelength', 'grating_central_wavelength', 'grating_tilt']:
            item_value = getattr(self, '%s_value' % item)
            item_unit = getattr(self, '%s_unit' % item)
            return units.Quantity(item_value, item_unit)
        else:
            return self.__getattribute__(item)

    @misc.lazyproperty
    def x_binning(self):
        x_binnings = np.array([detector.x_binning for detector in self.detectors])
        assert np.all(x_binnings == x_binnings[0])
        return x_binnings[0]

    @misc.lazyproperty
    def y_binning(self):
        y_binnings = np.array([detector.y_binning for detector in self.detectors])
        assert np.all(y_binnings == y_binnings[0])
        return y_binnings[0]

    @misc.lazyproperty
    def anamorphic_factor(self):
        return np.sin((self.calculated_grating_tilt + 50 * units.degree).to('rad').value) / \
               np.sin(self.calculated_grating_tilt.to('rad').value)


    @misc.lazyproperty
    def grating_equation_coefficient(self):
        return (self.grating.ruling_density * self.grating_central_wavelength).to(1).value


    @misc.lazyproperty
    def calculated_grating_tilt(self):
        return grating_equation_interpolator(self.grating_equation_coefficient) * units.degree


    @misc.lazyproperty
    def wavelength_start(self):
        wavelength_start_value = np.max([item.to('nm').value for item in [self.filter1.wavelength_start,
                                                                          self.filter2.wavelength_start,
                                                                          self.grating.wavelength_start]])
        return wavelength_start_value * units.Unit('nm')

    @misc.lazyproperty
    def wavelength_end(self):
        wavelength_end_value = np.min([item.to('nm').value for item in [self.filter1.wavelength_end,
                                                                        self.filter2.wavelength_end,
                                                                        self.detectors[-1].spectral_cutoff,
                                                                        self.grating.wavelength_end]])
        return wavelength_end_value * units.Unit('nm')

    @misc.lazyproperty
    def y_offset(self):
        return self.grating.y_offset / self.y_binning

    @misc.lazyproperty
    def y_distortion_coefficients(self):
        return detector_information['y-distortion'][self.instrument.name]

    @misc.lazyproperty
    def arcsec_per_mm(self):
        return detector_information['arcsecpermm'] * units.Unit('arcsec/mm')

    @misc.lazyproperty
    def x_pix_per_mm(self):
        return self.arcsec_per_mm / self.x_scale

    @misc.lazyproperty
    def y_pix_per_mm(self):
        return self.arcsec_per_mm / self.y_scale



    @misc.lazyproperty
    def chip_gap(self):
        if self.x_binning == 1:
            return detector_information['chip_gap']['unbinned']
        elif self.x_binning > 1:
            return np.int(np.round(detector_information['chip_gap']['binned'] / self.x_binning))

    def calculate_resolution(self, slit_width):
        """
            Calculate resolution

            Parameters
            ----------

            slit_width : `~astropy.units.Quantity`
                angle
        """
        return self.grating_equation_coefficient / (slit_width.to('rad').value * 81.0 *
                                             np.sin(self.calculated_grating_tilt.to('rad').value))

    @misc.lazyproperty
    def x_scale(self):
        return self.detectors[1].pixel_scale * self.x_binning

    @misc.lazyproperty
    def y_scale(self):
        return self.detectors[1].pixel_scale * self.y_binning


    @misc.lazyproperty
    def spectral_pixel_scale(self):
        xscale = self.x_scale.to('rad/pix').value
        spectral_pixel_scale_value = self.anamorphic_factor * xscale * self.grating_central_wavelength.to('nm').value * \
            81.0 * np.sin(self.calculated_grating_tilt.to('rad').value) / self.grating_equation_coefficient
        return spectral_pixel_scale_value * units.Unit('nm/pix')



    def __repr__(self):
        return "<GMOS MOS Instrument Setup Filter1 %s Filter2 %s Grating %s Tilt %.2f central wave=%.2f %s>" % \
                (self.filter1, self.filter2, self.grating, self.grating_tilt_value, self.grating_central_wavelength_value,
                self.grating_central_wavelength_unit)


class GMOSMOSRawFITS(Base):
    __tablename__ = 'gmos_mos_raw_fits'


    id = Column(Integer, ForeignKey('fits_file.id'), primary_key=True)
    mjd = Column(Float)
    instrument_id = Column(Integer, ForeignKey('instrument.id'))
    observation_block_id = Column(Integer, ForeignKey('observation_block.id'))
    observation_class_id = Column(Integer, ForeignKey('observation_class.id'))
    observation_type_id = Column(Integer, ForeignKey('observation_type.id'))
    object_id = Column(Integer, ForeignKey('object.id'))
    mask_id = Column(Integer, ForeignKey('gmos_mask.id'))


    instrument_setup_id = Column(Integer, ForeignKey('gmos_mos_instrument_setup.id'))

    exclude = Column(Boolean)


    fits = relationship(FITSFile, uselist=False, backref='raw_fits')
    instrument = relationship(Instrument, uselist=False, backref='raw_fits')
    observation_block = relationship(ObservationBlock, uselist=False, backref='raw_fits')
    observation_class = relationship(ObservationClass, uselist=False, backref='raw_fits')
    observation_type = relationship(ObservationType, uselist=False, backref='raw_fits')
    object = relationship(base.Object, uselist=False, backref='raw_fits')
    mask = relationship(GMOSMask, uselist=False, backref='raw_fits')

    instrument_setup = relationship(GMOSMOSInstrumentSetup)


    @property
    def associated_query(self):
        session = object_session(self)
        return session.query(GMOSMOSRawFITS).filter_by(observation_block_id=self.observation_block_id)

    @property
    def associated(self):
        return self.associated_query.all()

    @property
    def date_obs(self):
        return time.Time(self.mjd, scale='utc', format='mjd')

    def __init__(self, mjd, instrument_id, observation_block_id, observation_class_id, observation_type_id,
                 object_id, mask_id=None, instrument_setup_id=None, exclude=False):
        self.mjd = mjd
        self.instrument_id = instrument_id
        self.observation_block_id = observation_block_id
        self.observation_class_id = observation_class_id
        self.observation_type_id = observation_type_id
        self.exclude = exclude

        self.mask_id = mask_id
        self.object_id = object_id
        self.instrument_setup_id = instrument_setup_id
    def __repr__(self):
        return '<gmos fits="%s" class="%s" type="%s" object="%s">' % (self.fits.fname, self.observation_class.name,
                                                                  self.observation_type.name, self.object.name)


class GMOSMOSPrepared(Base):
    __tablename__ = 'gmos_mos_prepared'

    id = Column(Integer, ForeignKey('fits_file.id'), primary_key=True)
    raw_fits_id = Column(Integer, ForeignKey('gmos_mos_raw_fits.id'))

    raw_fits = relationship(GMOSMOSRawFITS, backref='prepared_fits')
    fits = relationship(FITSFile, uselist=False)
    #prepare_param_id = Column(Integer)


class GMOSMOSScience(Base):
    __tablename__ = 'gmos_mos_science'

    id = Column(Integer, ForeignKey('gmos_mos_raw_fits.id'), primary_key=True)
    flat_id = Column(Integer, ForeignKey('gmos_mos_raw_fits.id'))
    mask_arc_id = Column(Integer, ForeignKey('gmos_mos_raw_fits.id'))


    flat = relationship(GMOSMOSRawFITS, primaryjoin=(GMOSMOSRawFITS.id==flat_id),
                        backref=backref('flat2science', uselist=False))
    mask_arc = relationship(GMOSMOSRawFITS, primaryjoin=(GMOSMOSRawFITS.id==mask_arc_id),
                            backref=backref('mask2science', uselist=False))










