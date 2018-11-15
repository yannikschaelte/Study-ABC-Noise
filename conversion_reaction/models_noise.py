import pyabc
import pyabc.visualization
import copy
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.integrate as integrate
import os
from models import *


limits_noise = copy.deepcopy(limits)
limits_noise['noise'] = (0, 0.5)

prior_noise = pyabc.Distribution(**{key: pyabc.RV('uniform', bounds[0], bounds[1])
                                    for key, bounds in limits_noise.items()})


def model_random_noise(p):
    """
    Observations. Assume also noise variance unknown.
    """
    y = x(p)[1, :] + p['noise'] * noise_model(1, n_timepoints)
    return {'y': y.flatten()}


def distance_l2(x, y):
    return np.power( (x['y'] - y['y']), 2).sum()


th_true_noise = {'th0': th0_true, 'th1': th1_true, 'noise': noise}

def viz_noise(label, history):
    for t in range(1, history.max_t + 1):
        filename = label + "_kde_2d_" + str(t) + ".png"
        print(filename)
        if os.path.isfile(filename):
            continue
        df, w = history.get_distribution(m=0, t=t)
        axes = pyabc.visualization.plot_kde_matrix(
            df, w, numx=1000, numy=1000,
            limits=limits_noise, refval=th_true_noise)

        plt.savefig(filename)
        plt.close()

