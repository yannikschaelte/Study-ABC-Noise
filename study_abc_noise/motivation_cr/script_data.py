from study_abc_noise.read_pickle_file import read_pickle_file
from study_abc_noise.model import ConversionReaction1dModelVars as ModelVars
import matplotlib.pyplot as plt


data = read_pickle_file('data/cr_10_0.02__0.dat')
mv = ModelVars()
y = mv.get_model()(mv.p_true)

_, ax = plt.subplots()
ax.plot(mv.get_ts(), y['y'], '-', color='black', label='non-noisy model')
ax.plot(mv.get_ts(), data['y'], '*', color='black', label='noisy data')
ax.legend()
ax.set_xlabel("Time")
ax.set_ylabel("Model value")
plt.show()

