import cloudpickle as pickle
from study_abc_noise.read_pickle_file import read_pickle_file
from study_abc_noise.model import ConversionReactionModelVars
from study_abc_noise.optimize import get_optimal_kernel_value
import os
import numpy as np
import matplotlib.pyplot as plt
import pypesto.visualize
import cloudpickle as pickle


data_file = os.listdir('data')[0]
y_obs = read_pickle_file('data/' + data_file)
mv = ConversionReactionModelVars()
print(y_obs)

optimal_parameter, optimal_value = get_optimal_kernel_value(mv, y_obs)
with open("optimal_pdf_max.dat", 'wb') as f:
    pickle.dump(optimal_value, f)
print("optimal parameter: ", optimal_parameter)
print("optimal value: ", optimal_value)

distance_value = 2 * (optimal_value - 0.5 * (mv.n_t * np.log(2 * np.pi * mv.noise_std**2)))
print("optimal distance: ", distance_value)

#pypesto.visualize.waterfall(result)
#plt.show()
