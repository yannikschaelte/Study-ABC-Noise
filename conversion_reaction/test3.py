"""
This is for testing the stochastic acceptor.
"""

import pyabc
import numpy as np
import scipy as sp
from models import *
from scipy import stats

# VARIABLES

db_path = "sqlite:///db3.db"

# ACCEPTOR

distr = stats.multivariate_normal([0, 0], noise**2 * np.asarray([[1, 0],[0, 1]]))
acceptor = pyabc.StochasticAcceptor(distr)

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=pyabc.distance_functions.NoDistance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.epsilon.NoEpsilon)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)

# PLOT

visualize("test3", h)
