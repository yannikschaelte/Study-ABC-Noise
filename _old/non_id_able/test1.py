"""
This is for testing noisy data with a noise-less model.
"""

import pyabc
from pyabc.visualization import plot_kde_2d
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import tempfile
from models import *

# VARIABLES

db_path = "sqlite:///db1.db"

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model_random,
                   parameter_priors=prior,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   sampler=sampler,
                   eps=eps)

abc.new(db_path, y_obs)
h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations, min_acceptance_rate=min_acceptance_rate)
h = pyabc.History(db_path)

# PLOT
#viz("test1", h)
