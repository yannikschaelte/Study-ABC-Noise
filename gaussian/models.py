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
prior_lb = 0
prior_ub = 5
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)
                              for key in ['th0']})
prior_uvar = pyabc.Distribution(**{key: pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)
                                   for key in ['th0', 'th1']})

# MODEL


def model(p):
    return {'y0': p['th0']}


def model_random(p):
    return {'y0': p['th0'] + noise * noise_model()}


def model_random_uvar(p):
    return {'y0': p['th0'] + p['th1'] * noise_model()}


def model2(p):
    """
    Model2 returns a the sample mean and standard deviation.
    """
    xs = np.asarray([p['th0'] for _ in range(10)])
    return {'y0': np.mean(xs), 'y1': np.std(xs)}


def model2_random(p):
    xs = np.asarray([p['th0'] + noise * noise_model() for _ in range(10)])
    return {'y0': np.mean(xs), 'y1': np.std(xs)}


def model2_random_uvar(p):
    xs = np.asarray([p['th0'] + p['th1'] * noise_model() for _ in range(10)])
    return {'y0': np.mean(xs), 'y1': np.std(xs)}


# TRUE VARIABLES

th0_true = 2.5 # mean
th_true = {'th0': th0_true}
y_true = model(th_true)

# for unknown variance
th1_true = 1 # standard deviation
th_true_uvar = {'th0': th0_true, 'th1': th1_true}

# MEASURED DATA

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


def get_y_meas2():
    global _y_meas
    if _y_meas is not None:
        return _y_meas

    y_meas_file = "y_meas2.dat"
    try:
        y_meas = picle.load(open(y_meas_file, 'rb'))
    except Exception:
        y_meas = model2_random(th_true)
        pickle.dump(y_meas, open(y_meas_file, 'wb'))

    _y_meas = y_meas
    return y_meas


# PYABC PARAMETERS
distance = pyabc.PNormDistance(p=2)
pop_size = 50
transition = pyabc.MultivariateNormalTransition()
eps = pyabc.MedianEpsilon()
max_nr_populations = 8
sampler = pyabc.sampler.SingleCoreSampler()

# VISUALIZATION

def visualize(label, history):
    t = history.max_t

    df, w = history.get_distribution(m=0, t=t)
    ax = pyabc.visualization.plot_kde_1d(df, w,
                                          'th0',
                                          xmin=prior_lb, xmax=prior_ub, numx=300)
    ax.axvline(th0_true, color='k', linestyle='dashed')
    plt.savefig(label + "_kde_1d_" + str(t))
    plt.close()


def visualize_uvar(label, history):
    t = history.max_t

    df, w = history.get_distribution(m=0, t=t)
    ax = pyabc.visualization.plot_kde_matrix(df, w,
            limits={key: (prior_lb, prior_ub)
                    for key in ['th0', 'th1']})
    plt.savefig(label + "_kde_2d_" + str(t))
    plt.close()
