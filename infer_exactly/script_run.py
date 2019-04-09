"""
Illustrate the problem arising when not accounting for measurement noise.
Possible test cases: A deterministic model (-> ML estimator weighted by the
prior), and a simple SSA model (-> overestimation of parameters promoting
noise).

This means using a non-noisy model together with a deterministic acceptor,
and comparing that to the true analysis performed using either a noisy model,
or a stochastic acceptor.
"""

import pyabc
import datetime

from ..model import ConversionReactionModel
from ..settings import create_sampler

time = datetime.now().strftime("%Y%m%d_%H%M%S")
db_file = "sqlite:///db_" + time + ".db"

model = ConversionReactionModel()
y_obs = model.call_noisy(model.p_true)

abc = pyabc.ABCSMC(model = model.call,
                   parameter_priors = model.get_prior(),
                   distance_function = model.get_distance(),
                   population_size = model.get_pop_size(),
                   transitions = model.get_transition(),
                   eps = model.get_eps(),
                   acceptor = pyabc.UniformAcceptor(),
                   sampler = create_sampler())
abc.new(db_file, y_obs, gt_par=model.p_true)
abc.run(minimum_epsilon=model.eps_min,
        max_nr_populations=model.n_pop,
        min_acceptance_rate=model.min_acc_rate)
