"""
This is for testing a noisy model.
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

db_path = "sqlite:///db2.db"

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance_l2,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps)

#abc.new(db_path, y_obs)
#h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
h = pyabc.History(db_path)
# PLOT

visualize("test2", h)
