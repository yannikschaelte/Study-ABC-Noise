"""
This is for testing the model2 (mean+variance), estimating
mean and variance.
"""


import pyabc
from models import *


# VARIABLEs

db_path = "sqlite:///db5.db"

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model2_random_uvar,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps)

abc.new(db_path, get_y_meas2())

h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)

# PLOT

visualize_uvar("test5", h)
visualize_uvar_animated("test5", h)
