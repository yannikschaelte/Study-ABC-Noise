"""
Show that the approach can also be used to estimate noise parameters
like variances, which are often unknown. For the analysis, both a noisy
model output, and a stochastic acceptor can be employed to get correct
estimates.
"""

from study_abc_noise.model import ConversionReactionUVarModelVars as ModelVars
from study_abc_noise.vars import AnalysisVars, Task
import pyabc

mv = ModelVars(n_acc=100)#pdf_max=25.874)
pdf_max = mv.get_pdf_max()
print(pdf_max)

av_uniform = AnalysisVars(
    get_acceptor = lambda: pyabc.UniformAcceptor(),
    id_ = "noisy_model", n_pop=30)
av_stochastic = AnalysisVars(
    get_acceptor = lambda: pyabc.StochasticAcceptor(),#pdf_norm_method=pyabc.pdf_norm_from_kernel),
    get_eps = lambda: pyabc.Temperature(
        schemes=[pyabc.AcceptanceRateScheme(), pyabc.ExpDecayFixedRatioScheme()]),
    id_ = "stochastic_acceptor", n_pop = 10)

# create tasks
tasks = []
tasks.append(Task.from_vars(av_uniform, mv))
tasks.append(Task.from_vars(av_stochastic, mv))

for i in [0, 1]:
    print(tasks[i].acceptor, tasks[i].eps, tasks[i].distance)

# run
for task in tasks:
    task.execute()
