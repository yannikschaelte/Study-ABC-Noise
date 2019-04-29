import pyabc
import sys
from study_abc_noise.model import HodgkinHuxleyModelVars as ModelVars
from study_abc_noise.vars import AnalysisVars, Task


mv = ModelVars(n_acc=100, n_t=100)

# create analysis settings
list_analysis_vars = []
for acceptor, id_ in [
        (pyabc.StochasticAcceptor(
            temp_schemes=[
                pyabc.acceptor.scheme_acceptance_rate,
                pyabc.acceptor.scheme_exponential_decay],
            #pdf_max_method=pyabc.acceptor.pdf_max_take_max_found
            ),
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
tasks[1].model = mv.get_model()
tasks[1].eps_min = 0.0

# frun
for task in tasks:
    task.execute()
