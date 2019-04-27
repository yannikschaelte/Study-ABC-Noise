import numpy as np
import pyabc
from study_abc_noise.model import ConversionReactionModelVars
from study_abc_noise.vars import AnalysisVars, Task


mv = ConversionReactionModelVars()

def l1(x, y):
    return np.sum(np.abs(x['y'] - y['y']) / mv.noise_std)

mv.get_distance = lambda: l1

analysis_vars = AnalysisVars(
    get_acceptor=lambda: pyabc.UniformAcceptor(),
    id_="deterministic_l1")

tasks = []
tasks.append(Task.from_vars(analysis_vars, mv))
tasks[0].model =mv.get_model()
tasks[0].eps_min = 0.0

# run
for task in tasks:
    task.execute()
