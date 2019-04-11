"""
Compare all the different approaches of stochastic acceptance kernels
(Dalyy, Acceptance rate, Decay, ESS). Compare this to a setting with
a deterministic acceptance step but a noisy model output.

Note: It is problematic how to best choose termination criteria (minimal
acceptance rate, iteration, target epsilon variance).
"""


from study_abc_noise.model import (
    Gaussian1DModelVars, ConversionReactionModelVars)
from study_abc_noise.vars import (
    ModelVars, AnalysisVars, Task)
import pyabc
import numpy as np

n_rep = 10


# create model vars

list_model_vars = []
list_model_vars.append(Gaussian1DModelVars())
for n_t in [4, 7, 10, 13, 16]:
    model = ConversionReactionModelVars()
    model.n_t = n_t
    model.ts = np.linspace(0, 30, n_t)
    list_model_vars.append(model)

# create analysis vars

list_analysis_vars = []
list_analysis_vars.append(
    AnalysisVars(
        get_acceptor=lambda: pyabc.UniformAcceptor(), id_="uniform_acceptor"))
list_temp_schemes = [
    ([pyabc.acceptor.scheme_acceptance_rate], 'acc'),
    ([pyabc.acceptor.scheme_daly], 'daly'),
    ([pyabc.acceptor.scheme_decay], 'decay'),
    ([pyabc.acceptor.scheme_exponential_decay], 'exp_decay'),
    ([pyabc.acceptor.scheme_ess], 'ess'),
    ([pyabc.acceptor.scheme_acceptance_rate,
      pyabc.acceptor.scheme_decay], 'acc+dec'),
    ([pyabc.acceptor.scheme_acceptance_rate,
      pyabc.acceptor.scheme_ess,
      pyabc.acceptor.scheme_decay], 'acc+ess+dec'),
]
for temp_schemes, descr in list_temp_schemes:
    # note: need to bind input to lambda function, otherwise python silly
    get_acceptor = lambda temp_schemes=temp_schemes: \
        pyabc.StochasticAcceptor(temp_schemes=temp_schemes)
    analysis_vars = AnalysisVars(
        get_acceptor=get_acceptor,
        id_=f"stochastic_acceptor_{descr}")
    list_analysis_vars.append(analysis_vars)
list_analysis_vars.append(AnalysisVars(
    get_acceptor=lambda: pyabc.StochasticAcceptor(temp_schemes=[pyabc.acceptor.scheme_acceptance_rate]),
    n_pop=1,
    id_=f"stochastic_acceptor_rejection"))

# create tasks
tasks = []
for model_vars in list_model_vars:
    for analysis_vars in list_analysis_vars:
        for i_rep in range(n_rep):
            task = Task.from_vars(analysis_vars, model_vars, i_rep)
            tasks.append(task)

# run
for task in tasks:
    task.execute()
