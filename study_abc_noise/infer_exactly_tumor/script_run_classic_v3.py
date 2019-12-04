import pyabc
import tumor2d
import pickle

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

prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a)
                              for key, (a,b) in limits.items()})


_, _, data_var = tumor2d.load_default()
for key in keys[1:]:
    data_var[key] = data_var[key][:600]

distance = tumor2d.Tumor2DDistance(data_var)

sampler = pyabc.sampler.RedisEvalParallelSampler(host="icb-mona", port=8776)

abc = pyabc.ABCSMC(tumor2d.log_model, prior, distance, sampler=sampler,
                   population_size=500)
db_path="sqlite:///tumor2d_incorrect_v3.db"
abc.new(db_path, noisy_data)
abc.run()
