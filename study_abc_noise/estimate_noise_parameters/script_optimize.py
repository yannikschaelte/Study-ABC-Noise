import cloudpickle as pickle
from study_abc_noise.read_pickle_file import read_pickle_file
from study_abc_noise.model import ConversionReactionUVarModelVars as ModelVars
from study_abc_noise.optimize import get_optimal_kernel_value
import os
import numpy as np
import matplotlib.pyplot as plt
import pypesto.visualize


data_file = os.listdir('data')[0]
y_obs = read_pickle_file('data/' + data_file)
mv = ModelVars()
print(y_obs)

optimal_parameter, optimal_value = get_optimal_kernel_value(mv, y_obs)
print("optimal parameter: ", optimal_parameter)
print("optimal value: ", optimal_value)

distance_value = 2 * (optimal_value - 0.5 * (mv.n_t * np.log(2 * np.pi * mv.noise_std**2)))
print("optimal distance: ", distance_value)

pypesto.visualize.waterfall(result)
plt.show()
