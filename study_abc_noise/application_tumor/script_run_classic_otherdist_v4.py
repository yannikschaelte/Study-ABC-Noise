import pyabc
import tumor2d
import pickle
import numpy as np

# for debugging
import logging
df_logger = logging.getLogger('Distance')
df_logger.setLevel(logging.DEBUG)
df_logger = logging.getLogger('Acceptor')
df_logger.setLevel(logging.DEBUG)
df_logger = logging.getLogger('Epsilon')
df_logger.setLevel(logging.DEBUG)

noisy_data = pickle.load(open("noisy_data_v4.dat", "rb"))
noise_vector = pickle.load(open("noise_vector_v4.dat", "rb"))
noise = {'growth_curve': 40, 'extra_cellular_matrix_profile': 0.15, 'proliferation_profile': 0.02}
keys = ['growth_curve', 'extra_cellular_matrix_profile', 'proliferation_profile']  

limits = dict(log_division_rate=(-3, -1),
              log_division_depth=(1, 3),
              log_initial_spheroid_radius=(0, 1.2),
              log_initial_quiescent_cell_fraction=(-5, 0),
              log_ecm_production_rate=(-5, 0),
              log_ecm_degradation_rate=(-5, 0),
              log_ecm_division_threshold=(-5, 0))

prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a)
                              for key, (a,b) in limits.items()})


_, _, data_var = tumor2d.load_default()
data_var[keys[1]] = data_var[keys[1]][:640][::10]
data_var[keys[2]] = data_var[keys[2]][:345][::10]

def model(p):
    sim = tumor2d.log_model(p)
    sim[keys[1]] = sim[keys[1]][:640][::10]
    sim[keys[2]] = sim[keys[2]][:345][::10]
    return sim

#distance = tumor2d.Tumor2DDistance(data_var)

def distance(x, x0):
    dist = 0.0
    for key in keys:
        dist += np.sum((x[key] - x0[key])**2 / noise[key]**2)
    return dist

sampler = pyabc.sampler.RedisEvalParallelSampler(host="icb-mona", port=8778)

abc = pyabc.ABCSMC(model, prior, distance, sampler=sampler,
                   population_size=500)
db_path="sqlite:///tumor2d_incorrect_otherdist_v4.db"
abc.new(db_path, noisy_data)
abc.run()
