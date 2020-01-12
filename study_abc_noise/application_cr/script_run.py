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
from study_abc_noise.model import ConversionReactionModelVars as ModelVars
from study_abc_noise.vars import AnalysisVars, Task


arr_n_t = [3, 10, 30, 100, 300, 1000]
for n_t in arr_n_t:
    mv = ModelVars(n_t=n_t)
        acceptor = pyabc.StochasticAcceptor(
            

# create analysis settings
list_analysis_vars = []
for acceptor, id_ in [ 
        (pyabc.StochasticAcceptor(
            temp_schemes=[
                pyabc.acceptor.scheme_acceptance_rate,
                pyabc.acceptor.scheme_exponential_decay],
            pdf_max_method=pyabc.acceptor.pdf_max_take_max_found),
            "stochastic_acceptor_ada_c"),
        (pyabc.StochasticAcceptor(
            temp_schemes=[
                pyabc.acceptor.scheme_acceptance_rate,
                pyabc.acceptor.scheme_exponential_decay]),
            "stochastic_acceptor"),
        (pyabc.UniformAcceptor(), "deterministic"),
        (pyabc.UniformAcceptor(), "noisy_model")]:
    list_analysis_vars.append(
        AnalysisVars(
            get_acceptor=lambda acceptor=acceptor: acceptor, id_=id_))

# create tasks
tasks = []
for analysis_vars in list_analysis_vars:
    tasks.append(Task.from_vars(analysis_vars, mv))
# overwrite deterministic setting
tasks[2].model = mv.get_model()
tasks[2].eps_min = 0.0

# run
for task in tasks:
    task.execute()
