"""
Compare the behavior of the accepatance rate with the data dimension.
"""

from study_abc_noise.model.conversion_reaction import *
from study_abc_noise.model import ConversionReactionModelVars as ModelVars
import pyabc
import numpy as np


mv = ModelVars()

# generate data
y = mv.generate_data()
y = y['y'].flatten()

pdf_max = normal_dty(y, y, mv.noise_std * np.ones(mv.n_t))

