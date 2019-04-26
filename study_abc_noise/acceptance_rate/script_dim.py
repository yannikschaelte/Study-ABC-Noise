"""
Compare the behavior of the accepatance rate with the data dimension.
"""

from study_abc_noise.model.conversion_reaction import *
from study_abc_noise.model import ConversionReactionModelVars as ModelVars
import pyabc
import numpy as np


list_n_t = np.array([5, 10, 15, 20, 25, 30])
n_r = 1


for n_t in list_n_t:
    for i_r in range(n_r):
        mv = ModelVars(n_t=n_t)
        # generate data
        y = mv.generate_data()
        y = y['y'].flatten()

pdf_max = normal_dty(y, y, mv.noise_std * np.ones(mv.n_t))

