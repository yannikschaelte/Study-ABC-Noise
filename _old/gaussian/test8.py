import pyabc
from models import *

import logging
logger = logging.getLogger("Acceptor")
logger.setLevel(logging.DEBUG)

db_path = "sqlite:///db8.db"
distance = pyabc.distance.NormalKernel(mean=[0], cov=[noise**2])
acceptor = pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_ess, pyabc.acceptor.scheme_decay])

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

#visualize("test6_3", h)

