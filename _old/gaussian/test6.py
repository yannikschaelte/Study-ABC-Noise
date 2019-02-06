"""
This is for testing the stochastic acceptor
with known variance.

Attention:
    Set models.noise = 0.01 to get nicer results.
"""


import pyabc
from models import *
from scipy import stats

import logging
logger = logging.getLogger("Acceptor")
logger.setLevel(logging.DEBUG)

# VARIABLEs
db_path = "sqlite:///db6.db"
nr_pops = 8
distance = pyabc.distance.NormalKernel(mean=[0.0], cov=[noise**2])
acceptor = pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_decay])

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=1, max_nr_populations=1)
print(h.get_all_populations())
# PLOT

visualize("test6_0", h)

distance = pyabc.distance.NormalKernel(mean=[0], cov=[noise**2])
acceptor = pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_decay])

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=1, max_nr_populations=max_nr_populations)
print(h.get_all_populations())
# PLOT

visualize("test6_1", h)

distance = pyabc.distance.NormalKernel(mean=[0], cov=[noise**2])
acceptor = pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_daly])

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=1, max_nr_populations=max_nr_populations)
print(h.get_all_populations())
# PLOT

visualize("test6_2", h)

distance = pyabc.distance.NormalKernel(mean=[0], cov=[noise**2])
acceptor = pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_exponential_decay])

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=1, max_nr_populations=max_nr_populations)
print(h.get_all_populations())
# PLOT

visualize("test6_3", h)

