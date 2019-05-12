"""
Show that the approach can also be used to estimate noise parameters
like variances, which are often unknown. For the analysis, both a noisy
model output, and a stochastic acceptor can be employed to get correct
estimates.
"""

from study_abc_noise.model import ConversionReactionUVarModelVars as ModelVars
from study_abc_noise.vars import AnalysisVars, Task
import pyabc


mv = ModelVars(pdf_max=25.874)

uniform_acceptor = pyabc.UniformAcceptor()
stochastic_acceptor = pyabc.StochasticAcceptor(
    temp_schemes=[pyabc.acceptor.scheme_acceptance_rate,
                  pyabc.acceptor.scheme_decay],
    #pdf_max_method=pyabc.acceptor.pdf_max_take_max_found,
)
av_uniform = AnalysisVars(
    get_acceptor = lambda: uniform_acceptor,
    id_ = "noisy_model")
av_stochastic = AnalysisVars(
    get_acceptor = lambda: stochastic_acceptor,
    id_ = "stochastic_acceptor")

# create tasks
tasks = []
tasks.append(Task.from_vars(av_uniform, mv))
tasks.append(Task.from_vars(av_stochastic, mv))

# run
for task in tasks:
    task.execute()

