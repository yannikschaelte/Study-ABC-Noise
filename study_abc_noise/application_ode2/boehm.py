import os
os.environ["OMP_NUM_THREADS"] = "1"
import pyabc
import numpy as np
import pypesto
import logging
#%matplotlib inline

# for debugging
df_logger = logging.getLogger('Distance')
df_logger.setLevel(logging.DEBUG)
df_logger = logging.getLogger('Acceptor')
df_logger.setLevel(logging.DEBUG)
df_logger = logging.getLogger('Epsilon')
df_logger.setLevel(logging.DEBUG)

importer = pypesto.PetabImporter.from_folder("/home/icb/yannik.schaelte/benchmark-models/hackathon_contributions_new_data_format/Boehm_JProteomeRes2014")
problem = importer.petab_problem
objective = importer.create_objective()

#print(problem.get_dynamic_simulation_parameters())
#print(problem.get_optimization_parameters())
#df = problem.measurement_df
#print(df)
#print(df[['observableId', 'time', 'measurement']])
#print(objective(problem.x_nominal, return_dict=True)['rdatas'][0]['y'])

#problem.parameter_df

# data

mdf = problem.measurement_df
data = {}
keys = ['pSTAT5A_rel', 'pSTAT5B_rel', 'rSTAT5A_rel']
for key in keys:
    data[key] = np.array(mdf[mdf['observableId'] == key]['measurement'])
parameters = problem.get_optimization_parameters()
refval = {}
pdf = problem.parameter_df.reset_index()
for p in parameters:
    refval[p] = float(pdf[pdf['parameterId'] == p]['nominalValue'])

def model(p):
    if isinstance(p, (np.ndarray, list)):
        p_vector = p
    else:
        p_vector = np.zeros(len(parameters))
        for ip, p_id in enumerate(parameters):
            p_vector[ip] = p[p_id]

    rdatas =  objective(p_vector, return_dict=True)['rdatas']
    y = rdatas[0]['y']

    y_pSTAT5A_rel = np.array(y[:, 0])
    y_pSTAT5B_rel = np.array(y[:, 1])
    y_rSTAT5A_rel = np.array(y[:, 2])
    return {'pSTAT5A_rel': y_pSTAT5A_rel,
            'pSTAT5B_rel': y_pSTAT5B_rel,
            'rSTAT5A_rel': y_rSTAT5A_rel}

def get_var(p):
    sd_pSTAT5A_rel = p['sd_pSTAT5A_rel']
    sd_pSTAT5B_rel = p['sd_pSTAT5B_rel']
    sd_rSTAT5A_rel = p['sd_rSTAT5A_rel']
    var = []
    for key in keys:
        var.extend([(10**p[f'sd_{key}'])**2] * len(data[key]))
    return np.array(var)

# prior
limits = {}
for p in parameters:
    row = pdf[pdf['parameterId'] == p]
    limits[p] = (float(row['lowerBound']), float(row['upperBound']))
prior = pyabc.Distribution(**{p: pyabc.RV('uniform', a, b-a) for p, (a,b) in limits.items()})

class PDFNorm:

    def __init__(self):
        self.hit = False

    def __call__(
            self,
            prev_pdf_norm, get_weighted_distances, prev_temp, acceptance_rate, **kwargs):
        pdf_norm = pyabc.pdf_norm_max_found(prev_pdf_norm=prev_pdf_norm, get_weighted_distances=get_weighted_distances)
        print(" best: ", pdf_norm)
        if prev_temp is None or (acceptance_rate >= 0.1 and not self.hit):
            return pdf_norm
        self.hit = True
        temp = 0.5 * prev_temp
        offset = temp * np.log(20)
        used_norm = pdf_norm - offset
        used_norm = max(prev_pdf_norm, used_norm)
        print(" offsetted: ", pdf_norm - offset)
        return used_norm

# abc functions

distance = pyabc.IndependentNormalKernel(var=get_var, keys=keys)
acceptor = pyabc.StochasticAcceptor(pdf_norm_method=PDFNorm(), log_file="acceptor_log.log")
eps = pyabc.Temperature(schemes=[pyabc.AcceptanceRateScheme(target_rate=0.5, min_rate=1e-1),
                                 pyabc.ExpDecayFixedRatioScheme(alpha=0.5)],
                                 log_file="temperature_log.log")

#type(problem.x_nominal)
#print(model(problem.x_nominal))

#sampler = pyabc.sampler.SingleCoreSampler()
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=20)
#sampler = pyabc.sampler.RedisEvalParallelSampler(host="icb-sarah", port=8777)
abc = pyabc.ABCSMC(model, prior, distance, acceptor=acceptor, eps=eps, population_size=5000,
                   sampler=sampler)
abc.new("sqlite:///h_boehm4.db", data)
abc.run()

#h = pyabc.History("sqlite:///h_boehm.db")
#pyabc.visualization.plot_kde_matrix_highlevel(h, refval=refval, limits=limits)


