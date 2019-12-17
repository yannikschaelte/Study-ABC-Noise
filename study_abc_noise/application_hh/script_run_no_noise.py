from study_abc_noise.model import HodgkinHuxleyModelVars as ModelVars
from study_abc_noise.vars import AnalysisVars, Task
import pyabc


mv = ModelVars(noise_std=0.0)
av = AnalysisVars(get_acceptor=lambda: pyabc.UniformAcceptor(), id_="no_noise")

task = Task.from_vars(av, mv)
task.model = mv.get_model()
task.eps_min = 0.0

task.execute()
