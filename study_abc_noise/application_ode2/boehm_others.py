#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import numpy as np
#import scipy as sp

#def model(y, t, k1, k2, k3, epoRA):
#    x1, x2, x3, x4 = y

#    dx1 = - k1 * x1 * epoRA
#    dx2 = - k2 * x2**2 + k1 * x1 * epoRA
#    dx3 = - k3 * x3 + 0.5 * k2 * x2**2
#    dx4 = k3 * x3

#    return [dx1, dx2, dx3, dx4]


# In[1]:

import os
os.environ["OMP_NUM_THREADS"] = "1"
import pyabc
import numpy as np
import pypesto
#get_ipython().run_line_magic('matplotlib', 'inline')

# for debugging
import logging
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
print(problem.get_optimization_parameters())
#df = problem.measurement_df
#print(df)
#print(df[['observableId', 'time', 'measurement']])
#print(objective(problem.x_nominal, return_dict=True)['rdatas'][0]['y'])

problem.parameter_df


# In[3]:


# data

mdf = problem.measurement_df
data = {}
keys = ['pSTAT5A_rel', 'pSTAT5B_rel', 'rSTAT5A_rel']
for key in keys:
    data[key] = np.array(mdf[mdf['observableId'] == key]['measurement'])
parameters = problem.get_optimization_parameters()
pdf = problem.parameter_df.reset_index()
refval = {}
for p in parameters:
    refval[p] = float(pdf[pdf['parameterId'] == p]['nominalValue'])

def model_noisy(p):
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
    return {'pSTAT5A_rel': y_pSTAT5A_rel + 10**p['sd_pSTAT5A_rel'] * np.random.randn(len(y_pSTAT5A_rel)),
            'pSTAT5B_rel': y_pSTAT5B_rel + 10**p['sd_pSTAT5B_rel'] * np.random.randn(len(y_pSTAT5B_rel)),
            'rSTAT5A_rel': y_rSTAT5A_rel + 10**p['sd_rSTAT5A_rel'] * np.random.randn(len(y_rSTAT5A_rel))}

# prior
limits = {}
for p in parameters:
    row = pdf[pdf['parameterId'] == p]
    limits[p] = (float(row['lowerBound']), float(row['upperBound']))
prior = pyabc.Distribution(**{p: pyabc.RV('uniform', a, b-a) for p, (a,b) in limits.items()})

sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=20)
# In[14]:


#sampler = pyabc.sampler.SingleCoreSampler()
def distance_l2(x, y):
    d = 0
    for key in keys:
        d += np.sum((x[key] - y[key])**2)
    return d
abc = pyabc.ABCSMC(model_noisy, prior, distance_l2, sampler=sampler, population_size=5000)
abc.new("sqlite:///h_boehm_noisymodel.db", data)
abc.run(min_acceptance_rate=1e-3)


# In[15]:


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

#sampler = pyabc.sampler.SingleCoreSampler()
def distance_l2(x, y):
    d = 0
    for key in keys:
        d += np.sum((x[key] - y[key])**2)
    return d
abc = pyabc.ABCSMC(model, prior, distance_l2, sampler=sampler, population_size=5000)
abc.new("sqlite:///h_boehm_incorrect.db", data)
abc.run(min_acceptance_rate=1e-3)


# In[ ]:




