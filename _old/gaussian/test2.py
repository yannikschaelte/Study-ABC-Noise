import pyabc
from models import *
from scipy import stats

# VARIABLEs

db_path = "sqlite:///db2.db"

# PERFORM ABC ANALYSIS

abc = pyabc.ABCSMC(models=model_random,
                   parameter_priors=prior,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)

# PLOT

visualize("test2", h, True)
visualize_animated("test2", h, True)

# AND WITH STOCHASTIC ACCEPTOR

distr = stats.multivariate_normal([0], [noise**2])
nr_pops = 8
acceptor = pyabc.StochasticAcceptor(distribution=distr, nr_populations=nr_pops)
acceptor.max_temp = 100

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior_uvar,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps,
                   acceptor=acceptor)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=nr_pops)

visualize("test2_stoch", h)
visualize_animated("test2_stoch", h)
