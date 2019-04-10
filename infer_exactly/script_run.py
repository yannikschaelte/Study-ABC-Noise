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
import sys
sys.path.insert(0, '..')
from model import ConversionReactionModel
from util import create_sampler, get_timestamp


db_file = "sqlite:///db_" + get_timestamp() + ".db"
model = ConversionReactionModel()

# create data
y_obs = model.call_noisy(model.p_true)

# run with deterministic model and deterministic acceptor
# run with noisy model or stochastic acceptor
for m, a in zip([model.call, model.call_noisy, model.call],
                   [pyabc.UniformAcceptor(), pyabc.UniformAcceptor(),
                    pyabc.StochasticAcceptor(
                        temp_schemes=[pyabc.acceptor.scheme_acceptance_rate,
                                      pyabc.acceptor.scheme_decay])]):
    if isinstance(a, pyabc.acceptor.StochasticAcceptor):
        d = model.get_kernel()
        eps = pyabc.NoEpsilon()
    else:
        d = model.get_distance()
        eps = model.get_eps()
    abc = pyabc.ABCSMC(models = m,
        parameter_priors = model.get_prior(),
        distance_function = d,
        population_size = model.pop_size,
        transitions = model.get_transition(),
        eps = eps,
        acceptor = a,
        sampler = create_sampler())
    abc.new(db_file, y_obs, gt_par=model.p_true)
    abc.run(minimum_epsilon=model.eps_min,
            max_nr_populations=model.n_pop,
            min_acceptance_rate= model.min_acc_rate)
