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
from study_abc_noise.model import ConversionReactionModelVars
from study_abc_noise.vars import AnalysisVars, Task


# create analysis settings
list_analysis_vars = []
for acceptor, eps, id_ in [ 
        (pyabc.StochasticAcceptor(pdf_norm_method=pyabc.pdf_norm_from_kernel),
         pyabc.Temperature(),
         "stochastic_acceptor"),
        (pyabc.StochasticAcceptor(pdf_norm_method=pyabc.pdf_norm_from_kernel),
         pyabc.Temperature(initial_temperature=10),
         "stochastic_acceptor_init_temp_10")]:
    list_analysis_vars.append(
        AnalysisVars(
            get_acceptor=lambda acceptor=acceptor: acceptor,
            get_eps=lambda eps=eps: eps, id_=id_,
            n_acc=100))

# create tasks
tasks = []
for n_t in [2, 4, 6, 8, 10, 12, 15, 20]:
    for analysis_vars in list_analysis_vars:
        mv = ConversionReactionModelVars(n_t=n_t)
        tasks.append(Task.from_vars(analysis_vars, mv))

# run
for task in tasks:
    task.execute()
