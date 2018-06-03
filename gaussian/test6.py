"""
This is for testing the stochastic acceptor
with known variance.
"""


import pyabc
from models import *
from scipy import stats

# VARIABLEs

db_path = "sqlite:///db6.db"
distr = stats.multivariate_normal([0], [noise**2])
acceptor = pyabc.StochasticAcceptor(distr=distr)

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=500,
                   transitions=transition,
                   eps=eps,
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=1)

# PLOT

visualize("test6", h)
