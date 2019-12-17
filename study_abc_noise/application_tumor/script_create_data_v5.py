import tumor2d
import numpy as np
import pandas as pd
import pickle

exact_data = tumor2d.simulate()
keys = ['growth_curve', 'extra_cellular_matrix_profile', 'proliferation_profile']
noise = {'growth_curve': 40, 'extra_cellular_matrix_profile': 0.15, 'proliferation_profile': 0.02}
noise_vector = []

# cut off superfluous radii
exact_data[keys[1]] = exact_data[keys[1]][:640]
exact_data[keys[2]] = exact_data[keys[2]][:345]

# thin out
exact_data[keys[0]] = exact_data[keys[0]][::2]
exact_data[keys[1]] = exact_data[keys[1]][::70]
exact_data[keys[2]] = exact_data[keys[2]][::37]

noisy_data = {}

for key in keys:
    n = len(exact_data[key])
    noisy_data[key] = exact_data[key] + noise[key] * np.random.randn(n)
    noise_vector.extend([noise[key]] * n)

noise_vector = np.array(noise_vector)

print(exact_data)
print(noisy_data)
print(noise_vector)

pickle.dump(exact_data, open("exact_data_v5.dat", "wb"))
pickle.dump(noisy_data, open("noisy_data_v5.dat", "wb"))
pickle.dump(noise_vector, open("noise_vector_v5.dat", "wb"))
