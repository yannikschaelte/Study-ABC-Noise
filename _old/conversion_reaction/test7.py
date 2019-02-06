"""
This is for testing the stochastic acceptor.
"""

import pyabc
import numpy as np
import scipy as sp
from models import *
from scipy import stats
import numpy as np


# VARIABLES

db_path = "sqlite:///db7.db"

# ACCEPTOR
#max_nr_populations = 40

def pdf(x_0, x):
    v_0 = np.array(list(x_0.values()))
    v = np.array(list(x.values()))
    mean = np.zeros(n_timepoints)
    cov = noise**2 * np.eye(n_timepoints)
    return stats.multivariate_normal.pdf(v_0 - v, mean=mean, cov=cov)

c = pickle.load(open("best_found_pdf_" + str(noise) + "_" + str(n_timepoints) + ".dat", "rb"))
c = np.log(c)


for acceptor, label in [
        (pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_acceptance_rate, pyabc.acceptor.scheme_decay], pdf_max_method=pyabc.acceptor.pdf_max_take_max_found), "adaptive"),
        (pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_decay], pdf_max_method=pyabc.acceptor.pdf_max_take_max_found), "decay"),
        (pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_exponential_decay], pdf_max_method=pyabc.acceptor.pdf_max_take_max_found), "exp_decay"),
        (pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_daly], pdf_max_method=pyabc.acceptor.pdf_max_take_max_found), "daly"),
        (pyabc.StochasticAcceptor(temp_schemes = [pyabc.acceptor.scheme_acceptance_rate, pyabc.acceptor.scheme_exponential_decay, pyabc.acceptor.scheme_decay, pyabc.acceptor.scheme_daly], pdf_max_method=pyabc.acceptor.pdf_max_take_max_found), "all")]:

    distance = pyabc.distance.IndependentNormalKernel(mean=np.zeros(n_timepoints), var = noise**2 * np.ones(n_timepoints), pdf_max = c)

    abc = pyabc.ABCSMC(models=model,
                       parameter_priors=prior,
                       distance_function=distance,
                       population_size=pop_size,
                       transitions=transition,
                       eps=pyabc.NoEpsilon(),
                       acceptor=acceptor,
                       sampler=sampler)

    abc.new(db_path, y_obs)
    h = abc.run(minimum_epsilon=1, max_nr_populations=max_nr_populations, min_acceptance_rate=min_acceptance_rate)
    h = pyabc.History(db_path)

    # PLOT
    #viz("test5_" + label, h)