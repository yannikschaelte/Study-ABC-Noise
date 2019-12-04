import tumor2d
import numpy as np
import pandas as pd
import pickle

exact_data = tumor2d.simulate()
keys = ['growth_curve', 'extra_cellular_matrix_profile', 'proliferation_profile']
noise = {'growth_curve': 30, 'extra_cellular_matrix_profile': 0.2, 'proliferation_profile': 0.03}
noise_vector = []

# cut off superfluous radii
for key in keys[1:]:
    exact_data[key] = exact_data[key][:600]

noisy_data = {}

for key in keys:
    n = len(exact_data[key])
    noisy_data[key] = exact_data[key] + noise[key] * np.random.randn(n)
    noise_vector.extend([noise[key]] * n)

noise_vector = np.array(noise_vector)

print(exact_data)
print(noisy_data)
print(noise_vector)

pickle.dump(exact_data, open("exact_data_v3.dat", "wb"))
pickle.dump(noisy_data, open("noisy_data_v3.dat", "wb"))
pickle.dump(noise_vector, open("noise_vector_v3.dat", "wb"))
