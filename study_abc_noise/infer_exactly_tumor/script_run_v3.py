import pyabc
import tumor2d
import pickle

# for debugging
import logging
df_logger = logging.getLogger('Distance')
df_logger.setLevel(logging.DEBUG)
df_logger = logging.getLogger('Acceptor')
df_logger.setLevel(logging.DEBUG)
df_logger = logging.getLogger('Epsilon')
df_logger.setLevel(logging.DEBUG)

noisy_data = pickle.load(open("noisy_data_v3.dat", "rb"))
noise_vector = pickle.load(open("noise_vector_v3.dat", "rb"))
keys = ['growth_curve', 'extra_cellular_matrix_profile', 'proliferation_profile']  

limits = dict(log_division_rate=(-3, -1),
              log_division_depth=(1, 3),
              log_initial_spheroid_radius=(0, 1.2),
              log_initial_quiescent_cell_fraction=(-5, 0),
              log_ecm_production_rate=(-5, 0),
              log_ecm_degradation_rate=(-5, 0),
              log_ecm_division_threshold=(-5, 0))

def model(p):
    sim = tumor2d.log_model(p)
    sim[keys[1]] = sim[keys[1]][:640]
    sim[keys[2]] = sim[keys[2]][:345]
    return sim

prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a)
                              for key, (a,b) in limits.items()})

acceptor = pyabc.StochasticAcceptor()
temperature = pyabc.Temperature(schemes=[pyabc.AcceptanceRateScheme(), pyabc.ExpDecayFixedRatioScheme(alpha=0.6)])
kernel = pyabc.IndependentNormalKernel(keys=keys, var=noise_vector**2)

sampler = pyabc.sampler.RedisEvalParallelSampler(host="icb-mona", port=8775)
#sampler = pyabc.sampler.MulticoreEvalParallelSampler(daemon=False)

abc = pyabc.ABCSMC(model, prior, kernel, sampler=sampler,
                   acceptor=acceptor, eps=temperature, population_size=500)
db_path="sqlite:///tumor2d_stoch_acc_v3.db"
abc.new(db_path, noisy_data)
abc.run()
