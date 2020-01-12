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
from study_abc_noise.model import MRNATranscriptionModelVars as ModelVars
from study_abc_noise.util import create_sampler, get_timestamp
from study_abc_noise.vars import AnalysisVars, Task, get_data
from study_abc_noise.optimize import get_and_store_optimal_kernel_value


n_rep = 10
n_t = 10
n_acc = 100
noise_success_probability = 0.7

mv = ModelVars(n_acc=n_acc, n_t=n_t, noise_success_probability=noise_success_probability)
y_obs = get_data(mv, 0)
# optimal_pdf_max = get_and_store_optimal_kernel_value(mv, y_obs, 0)
# print(optimal_pdf_max)
pdf_maxs = [None, None]
pdf_max_methods = [pyabc.acceptor.pdf_max_take_max_found,
                pyabc.acceptor.pdf_max_take_from_kernel]
# create analysis settings
list_model_vars = []
list_analysis_vars = []
ids = ["no_pdf_max_info", "adaptive_pdf_max"]
for pdf_max, pdf_max_method, id_ in zip(pdf_maxs, pdf_max_methods, ids):
    acceptor = pyabc.StochasticAcceptor(
        temp_schemes=[
            pyabc.acceptor.scheme_acceptance_rate,
            pyabc.acceptor.scheme_exponential_decay],
        pdf_max_method=pdf_max_method)
    list_analysis_vars.append(
        AnalysisVars(
            get_acceptor=lambda acceptor=acceptor: acceptor, id_=id_))
    mv = ModelVars(n_acc=n_acc, n_t=n_t, noise_success_probability=noise_success_probability)
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
pdf_maxs = tasks[1].acceptor.pdf_maxs
with open("pdf_maxs.dat", 'wb') as f:
    pickle.dump(pdf_maxs, f)
