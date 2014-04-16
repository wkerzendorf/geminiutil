import re

import numpy as np
import pandas as pd

from astropy import units as u

sexagesimal_pattern = re.compile('([+-]?\d+):(\d+):(\d+\.?\d+?)')
def convert_sexagesimal(sexagesimal_string):
    sexagesimal = map(float, sexagesimal_pattern.match(sexagesimal_string).groups())
    return sexagesimal[0] + sexagesimal[1] / 60. + sexagesimal[2] / 3600.


def read_standard_star_db(filename):
    column_names=('name', 'ra', 'dec', 'u', 'g', 'r', 'i', 'z')
    standard_fields_raw = np.recfromtxt(filename,
                                     converters={1:convert_sexagesimal,
                                                 2:convert_sexagesimal},
                                     usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                                     names=column_names)
    standard_fields = pd.DataFrame(standard_fields_raw)

    standard_fields['ra'] = u.hourangle.to(u.degree, standard_fields['ra'])

    standard_fields['fields'] = [re.split('[-_+]', item)[0] for item in standard_fields['name']]
    return standard_fields