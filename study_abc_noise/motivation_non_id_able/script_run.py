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
from study_abc_noise.model import NonIdAblePrioredModelVars
from study_abc_noise.util import create_sampler, get_timestamp
from study_abc_noise.vars import AnalysisVars, Task


mv = NonIdAblePrioredModelVars()


# create analysis settings
list_analysis_vars = []
for acceptor, id_ in [
        (pyabc.UniformAcceptor(), "deterministic"),
        (pyabc.UniformAcceptor(), "noisy_model"),
        (pyabc.StochasticAcceptor(
            temp_schemes=[
                pyabc.acceptor.scheme_acceptance_rate,
                pyabc.acceptor.scheme_decay]), "stochastic_acceptor")]:
    list_analysis_vars.append(
        AnalysisVars(
            get_acceptor=lambda acceptor=acceptor: acceptor, id_=id_))

# create tasks
tasks = []
for analysis_vars in list_analysis_vars:
    tasks.append(Task.from_vars(analysis_vars, mv, 0))
# overwrite deterministic setting
tasks[0].model = mv.get_model()
tasks[0].eps_min = 0.0

# run
for task in tasks:
    task.execute()
