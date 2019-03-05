"""
This is for testing noisy data with a noisy model with unknown noise.
"""

import pyabc
from pyabc.visualization import plot_kde_2d
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import tempfile
from models_noise import *

# VARIABLES

db_path = "sqlite:///db_noise1.db"

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model_random_noise,
                   parameter_priors=prior_noise,
                   distance_function=distance_l2,
                   population_size=pop_size,
                   transitions=transition,
                   sampler=sampler,
                   eps=eps)

#abc.new(db_path, y_obs)
#h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations, min_acceptance_rate=min_acceptance_rate)
h = pyabc.History(db_path)

# PLOT
viz_noise("test_noise1", h, True)
