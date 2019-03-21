"""
This is for testing noisy data with a noisy model with unknown noise.
"""

import pyabc
from pyabc.visualization import plot_kde_2d
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import tempfile
from models_noise import *

# VARIABLES

db_path = "sqlite:///db_noise2.db"

# PERFORM ABC ANALYSIS


def compute_var(par):
    return np.ones(n_timepoints) * par['noise']**2


distance = pyabc.distance.IndependentNormalKernel(mean=np.zeros(n_timepoints), var=compute_var, pdf_max=36.86)

acceptor = pyabc.StochasticAcceptor(temp_schemes =[pyabc.acceptor.scheme_acceptance_rate, pyabc.acceptor.scheme_decay], pdf_max_method=pyabc.acceptor.pdf_max_take_max_found)

abc = pyabc.ABCSMC(models=model_random_noise,
                   parameter_priors=prior_noise,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   sampler=sampler,
                   acceptor=acceptor,
                   eps=pyabc.NoEpsilon())


#abc.new(db_path, y_obs)
#h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations, min_acceptance_rate=min_acceptance_rate)
h = pyabc.History(db_path)

# PLOT
viz_noise("test_noise2", h, True)
