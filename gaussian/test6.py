"""
This is for testing the stochastic acceptor
with known variance.

Attention:
    Set models.noise = 0.01 to get nicer results.
"""


import pyabc
from models import *
from scipy import stats

# VARIABLEs

db_path = "sqlite:///db6.db"
distr = stats.multivariate_normal([0], [noise**2])
nr_pops = 1
acceptor = pyabc.StochasticAcceptor(distribution=distr, nr_populations=nr_pops)
acceptor.max_temp = 100

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=100,
                   transitions=transition,
                   eps=eps,
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=nr_pops)

# PLOT

visualize("test6", h)
