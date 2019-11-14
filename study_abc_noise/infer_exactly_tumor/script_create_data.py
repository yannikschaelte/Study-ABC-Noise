import tumor2d
import numpy as np
import pandas as pd
import pickle

exact_data = tumor2d.simulate()
keys = ['growth_curve', 'extra_cellular_matrix_profile', 'proliferation_profile']
noise = {'growth_curve': 5, 'extra_cellular_matrix_profile': 1e-2, 'proliferation_profile': 1e-2}
noise_vector = []

noisy_data = {}

for key in keys:
    n = len(exact_data[key])
    noisy_data[key] = exact_data[key] + noise[key] * np.random.randn(n)
    noise_vector.extend([noise[key]] * n)

noise_vector = np.array(noise_vector)

print(exact_data)
print(noisy_data)
print(noise_vector)

pickle.dump(exact_data, open("exact_data.dat", "wb"))
pickle.dump(noisy_data, open("noisy_data.dat", "wb"))
pickle.dump(noise_vector, open("noise_vector.dat", "wb"))
