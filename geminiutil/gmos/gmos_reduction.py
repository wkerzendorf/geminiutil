import numpy as np
from astropy.io import fits
from .. base import FITSFile, Base, Operations

import os
from sqlalchemy import String, Integer, Float, DateTime, Boolean
from sqlalchemy import Column, ForeignKey

class GMOSMOSChipCombine(Base):
    __tablename__ = 'gmos_mos_chip_combine'

    operation_prefix = 'c'
    id = Column(Integer, primary_key=True)

    def __call__(self, fits_object, fname=None, reduction_dir='.'):
        new_data = combine_gmos_chips(fits_object)
        header = fits_object.header
        header[self.__tablename__] = 'DONE'
        if fname is None:
            fname = self.operation_prefix + fits_object.fname
        full_path = os.path.join(reduction_dir, fname)
        fits.PrimaryHDU(data=new_data, header=header).writeto(full_path)

        out_fits_object = FITSFile.from_fits_file(full_path)
        operations_object = Operations.from_fits_objects(input_fits_object=fits_object, out_fits_object=out_fits_object, operations_type_id=0, operations_id=self.id)
        return operations_object


def combine_gmos_chips(fits_object):
    chip_data = []
    new_data_shape = [-1, 0]
    column_shapes = []
    for i in xrange(3):
        current_data = fits_object.fits_data[i+1].data

        if new_data_shape[0] == -1:
            new_data_shape[0] = current_data.shape[0]
        else:
            assert new_data_shape[0] == current_data.shape[0]

        column_shapes.append(current_data.shape[1])
        chip_data.append(current_data)

    new_data_shape[1] = np.sum(column_shapes)
    new_data = np.zeros(new_data_shape, dtype=current_data.dtype)

    chip_column_start = 0
    chip_column_end = column_shapes[0]
    for i in xrange(3):
        new_data[:, chip_column_start:chip_column_end] = fits_object.fits_data[i+1].data
        if i < 2:
            chip_column_start += column_shapes[i]
            chip_column_end += column_shapes[i+1]
    return new_data





