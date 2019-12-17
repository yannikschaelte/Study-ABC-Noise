import pickle
import tumor2d
import numpy as np

noise_vector = pickle.load(open("noise_vector_v4.dat", 'rb'))
exact_data = pickle.load(open("exact_data_v4.dat", 'rb'))
noisy_data = pickle.load(open("noisy_data_v4.dat", 'rb'))
noise = {'growth_curve': 40, 'extra_cellular_matrix_profile': 0.15, 'proliferation_profile': 0.02}
keys = ['growth_curve', 'extra_cellular_matrix_profile', 'proliferation_profile']  

_, _, data_var = tumor2d.load_default()
data_var[keys[1]] = data_var[keys[1]][:640:10]
data_var[keys[2]] = data_var[keys[2]][:345:10]
print(data_var)
var = []
for key in keys:
    var.extend(data_var[key])
var = np.array(var)
weight = np.zeros_like(var)
weight[var != 0]  = 1 / var[var != 0]

weight = 1 / noise_vector**2
a = noise_vector
#a = np.sqrt(12) / 2 * noise_vector
eps = np.sum(a**2 * weight)

print(eps)

dist = 0
for key in keys:
    dist += sum((noisy_data[key] - exact_data[key])**2 / noise[key]**2)
print(dist)
