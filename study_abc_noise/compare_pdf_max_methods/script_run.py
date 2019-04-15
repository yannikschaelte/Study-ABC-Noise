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
import cloudpickle as pickle
import os
import matplotlib.pyplot as plt
from study_abc_noise.model import ConversionReactionModelVars
from study_abc_noise.util import create_sampler, get_timestamp
from study_abc_noise.vars import AnalysisVars, Task
from study_abc_noise.optimize import get_and_store_optimal_kernel_value


n_rep = 10

mv = ConversionReactionModelVars()
y_obs = Task.get_data(mv, 0)
optimal_pdf_max = get_and_store_optimal_kernel_value(mv, y_obs, 0)
pdf_maxs = [None, optimal_pdf_max, None]
pdf_max_methods = [pyabc.acceptor.pdf_max_take_from_kernel,
                   pyabc.acceptor.pdf_max_take_from_kernel,
                   pyabc.acceptor.pdf_max_take_max_found]
# create analysis settings
list_model_vars = []
list_analysis_vars = []
ids = ["no_pdf_max_info", "optimal_pdf_max", "adaptive_pdf_max"]
for pdf_max, pdf_max_method, id_ in zip(pdf_maxs, pdf_max_methods, ids):
    acceptor = pyabc.StochasticAcceptor(
        temp_schemes=[
            pyabc.acceptor.scheme_acceptance_rate,
            pyabc.acceptor.scheme_decay],
        pdf_max_method=pdf_max_method)
    list_analysis_vars.append(
        AnalysisVars(
            get_acceptor=lambda acceptor=acceptor: acceptor, id_=id_))
    mv = ConversionReactionModelVars()
    mv.pdf_max = pdf_max
    list_model_vars.append(mv)

# create tasks
tasks = []
for model_vars, analysis_vars in zip(list_model_vars, list_analysis_vars):
    tasks.append(Task.from_vars(analysis_vars, model_vars, 0))

# run
for task in tasks:
    task.execute()

# save pdf_maxs
pdf_maxs = tasks[2].acceptor.pdf_maxs
plt.plot(pdf_maxs.keys(), pdf_maxs.values())
plt.savefig("pdf_maxs.png")

