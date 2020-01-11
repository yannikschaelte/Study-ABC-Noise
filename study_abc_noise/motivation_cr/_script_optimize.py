import cloudpickle as pickle
from study_abc_noise.read_pickle_file import read_pickle_file
from study_abc_noise.model import ConversionReactionModelVars
from study_abc_noise.optimize import get_optimal_kernel_value
import os
import numpy as np
import matplotlib.pyplot as plt
import pypesto.visualize


data_file = os.listdir('data')[0]
y_obs = read_pickle_file('data/' + data_file)
mv = ConversionReactionModelVars()
print(y_obs)

result = get_optimal_kernel_value(mv, y_obs)
optimal_parameter = result.optimize_result.get_for_key('x')[0]
optimal_value = result.optimize_result.get_for_key('fval')[0]
print("optimal parameter: ", optimal_parameter)
print("optimal value: ", optimal_value)

distance_value = 2 * (optimal_value - 0.5 * (mv.n_t * np.log(2 * np.pi * mv.noise_std**2)))
print("optimal distance: ", distance_value)

pypesto.visualize.waterfall(result)
plt.show()
