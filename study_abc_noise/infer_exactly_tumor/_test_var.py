import pickle
import tumor2d
import numpy as np

noise_vector = pickle.load(open("noise_vector_v3.dat", 'rb'))

keys = ['growth_curve', 'extra_cellular_matrix_profile', 'proliferation_profile']  

_, _, data_var = tumor2d.load_default()
data_var[keys[1]] = data_var[keys[1]][:640]
data_var[keys[2]] = data_var[keys[2]][:345]

var = []
for key in keys:
    var.extend(data_var[key])
    print(data_var[key].size)
var = np.array(var)

total_var = np.sum(var * noise_vector)

print(total_var)

