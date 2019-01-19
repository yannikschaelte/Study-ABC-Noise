"""
This is for testing model_random_uvar (variance from p),
estimating mean and variance.
"""

import pyabc
from models import *

# VARIABLEs

db_path = "sqlite:///db4.db"

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model_random_uvar,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)

# PLOT

visualize_uvar("test4", h)
visualize_uvar_animated("test4", h)
