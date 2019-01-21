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
nr_pops = 8
def pdf(x_0, x):
    return stats.multivariate_normal.pdf(np.array(list(x_0.values())) - np.array(list(x.values())), mean=[0], cov=[noise**2])
acceptor = pyabc.StochasticAcceptor(pdf=pdf, use_target_acceptance_rate=False)

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=pyabc.NoDistance(),
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=1, max_nr_populations=1)
print(h.get_all_populations())
# PLOT

visualize("test6_0", h)


acceptor = pyabc.StochasticAcceptor(pdf=pdf, use_target_acceptance_rate=False)

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=pyabc.NoDistance(),
                   population_size=pop_size,
                   transitions=transition,
                   eps=pyabc.NoEpsilon(),
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=1, max_nr_populations=max_nr_populations)
print(h.get_all_populations())
# PLOT

visualize("test6_1", h)
