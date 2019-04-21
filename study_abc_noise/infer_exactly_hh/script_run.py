import pyabc
import sys
from study_abc_noise.model import HodgkinHuxleyModelVars as ModelVars
from study_abc_noise.vars import AnalysisVars, Task


mv = ModelVars()

# create analysis settings
list_analysis_vars = []
for acceptor, id_ in [
        (pyabc.UniformAcceptor(), "deterministic"),
        (pyabc.UniformAcceptor(), "noisy_model"),
        (pyabc.StochasticAcceptor(
            temp_schemes=[
                pyabc.acceptor.scheme_acceptance_rate,
                pyabc.acceptor.scheme_exponential_decay]),
            "stochastic_acceptor")]:
