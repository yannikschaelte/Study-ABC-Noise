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
              log_ecm_division_threshold=(-5, 0),
              std_growth_curve=(10, 100),
              std_extra_cellular_matrix_profile=(0.05, 0.5),
              std_proliferation_profile=(0.01, 0.03))

def model(p):
    p = {key: val for key, val in p.items() if 'std_' not in key}
    sim = tumor2d.log_model(p)
    sim[keys[1]] = sim[keys[1]][:640][::10]
    sim[keys[2]] = sim[keys[2]][:345][::10]
    return sim

def get_var(p):
    var_vector = []
    for key in keys:
        n = len(noisy_data[key])
        var_vector.extend([p['std_' + key]**2] * n)
    return np.array(var_vector)

prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a)
                              for key, (a,b) in limits.items()})

acceptor = pyabc.StochasticAcceptor()
temperature = pyabc.Temperature(schemes=[pyabc.AcceptanceRateScheme(), pyabc.ExpDecayFixedRatioScheme(alpha=0.6)])
kernel = pyabc.IndependentNormalKernel(keys=keys, var=get_var)

sampler = pyabc.sampler.RedisEvalParallelSampler(host="icb-mona", port=8774)
#sampler = pyabc.sampler.MulticoreEvalParallelSampler(daemon=False)

abc = pyabc.ABCSMC(model, prior, kernel, sampler=sampler,
                   acceptor=acceptor, eps=temperature, population_size=500)
db_path="sqlite:///tumor2d_stoch_acc_v4_uvar.db"
abc.new(db_path, noisy_data)
abc.run()
