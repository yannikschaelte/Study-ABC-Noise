"""
Simple gaussian model, mean (and variance) to be inferred.
"""

import numpy as np
import scipy as sp
import pickle
import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt


# VARIABLES

noise = 1
noise_model = np.random.randn

# prior
prior_lb = -5
prior_ub = 5
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)
                              for key in ['th0']})

# MODEL


def model(p):
    return {'y0': p['th0']}


def model_random(p):
    return {'y0': p['th0'} + noise * noise_model()}


# TRUE VARIABLES

th0_true = 0 # mean
th_true = {'th0': th0_true}
y_true = model(th_true)

# for unknown variance
th1_true = 1 # standard deviation
th_true_uvar = {'th0': th0_true, 'th1': th1_true}

# MESAURED DATA

_y_meas = None


def get_y_meas():
    global _y_meas
    if _y_meas is not None:
        return _y_meas

    y_meas_file = "y_meas.dat"
    try:
        y_meas = pickle.load(open(y_meas_file, 'rb'))
    except Exception:
        y_true = model(th_true)
        y_meas = {'y0': y_true['y0'] + noise * noise_model()}
        pickle.dump(y_meas, open(y_meas_file, 'wb'))

    _y_meas = y_meas
    return y_meas


# PYABC PARAMETERS
distance = pyabc.PNormDistance(p=2)
pop_size = 50
transition = pyabc.MultivariateNormalTransition()
eps = pyabc.MedianEpsilon()
max_nr_populations = 20
sampler = pyabc.sampler.SingleCoreSampler()
