import pyabc
import pyabc.visualization
import copy
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.integrate as integrate
import os
from models import *


limits['noise'] = (0, 0.2)

prior_noise = pyabc.Distribution(**{key: pyabc.RV('uniform', bounds[0], bounds[1])
                                    for key, bounds in limits.items()})


def model_random_noise(p):
    """
    Observations. Assume also noise variance unknown.
    """
    y = x(p)[1, :] + p['noise'] * noise_model(1, n_timepoints)
    return {'y': y.flatten()}


th_true_noise = {'th0': th0_true, 'th1': th1_true, 'noise': noise}


def pdf_true_noise(p):
    """
    Unscaled posterior density.
    """
    if type(p) is list:
        p = {'th0': p[0], 'th1': p[1], 'noise': p[2]}

    # prior value
    prior_val = prior_noise.pdf(p)

    # observed data
    y_bar = y_obs['y'].flatten()

    # data for p
    y = model(p)['y'].flatten()

    # likelihood value
    dim = len(y)
    sigma = p['noise'] * np.ones(dim)
    likelihood_val = normal_dty(y_bar, y, sigma)

    # posterior value
    unscaled_posterior = likelihood_val * prior_val

    return unscaled_posterior


# VISUALIZATION

def for_plot_pdf_true_noise(num_integral=False):
    for_plot_pdf_true_file = "for_plot_pdf_true_noise.dat"
    try:
        xs_0, ys_0, xs_1, ys_1, xs_noise, ys_noise = pickle.load(open(for_plot_pdf_true_file, 'rb'))
    except Exception as e:
        print(e)
        n_mesh = 200
        # th0
        def marginal_0(th0):
            return integrate.dblquad(lambda th1, noise: pdf_true_noise(
                {'th0': th0, 'th1': th1, 'noise': noise}),
                limits['th1'][0], limits['th1'][1],
                lambda x: limits['noise'][0], lambda x: limits['noise'][1])[0]
        if num_integral:
            integral_0 = 1.0
        else:
            integral_0 = integrate.quad(marginal_0, limits['th0'][0], limits['th0'][1])[0]
        xs_0 = np.linspace(limits['th0'][0], limits['th0'][1], n_mesh)
        ys_0 = []
        for x in xs_0:
            ys_0.append(marginal_0(x) / integral_0)
        if num_integral:
            sum_ = sum(ys_0) * (limits['th0'][1] - limits['th0'][0]) / (n_mesh - 1)
            ys_0 = [val / sum_ for val in ys_0]
        print(ys_0)

        # th1
        def marginal_1(th1):
            return integrate.dblquad(lambda th0, noise: pdf_true_noise(
                {'th0': th0, 'th1': th1, 'noise': noise}),
                limits['th0'][0], limits['th0'][1],
                lambda x: limits['noise'][0], lambda x: limits['noise'][1])[0]
        integral_1 = integrate.quad(marginal_1, limits['th1'][0], limits['th1'][1])[0]
        xs_1 = np.linspace(limits['th1'][0], limits['th1'][1], n_mesh)
        ys_1 = []
        for x in xs_1:
            ys_1.append(marginal_1(x) / integral_1)

        # noise
        def marginal_noise(noise):
            return integrate.dblquad(lambda th0, th1: pdf_true_noise(
                {'th0': th0, 'th1': th1, 'noise': noise}),
                limits['th0'][0], limits['th0'][1],
                lambda x: limits['th1'][0], lambda x: limits['th1'][1])[0]
        integral_noise = integrate.quad(marginal_noise, limits['noise'][0], limits['noise'][1])[0]
        xs_noise = np.linspace(limits['noise'][0], limits['noise'][1], n_mesh)
        ys_noise = []
        for x in xs_noise:
            ys_noise.append(marginal_noise(x) / integral_noise)
        pickle.dump((xs_0, ys_0, xs_1, ys_1, xs_noise, ys_noise), open(for_plot_pdf_true_file, 'wb'))
    
    return xs_0, ys_0, xs_1, ys_1, xs_noise, ys_noise


def viz_noise(label, history, show_true=True):
    # compute true posterior
    if show_true:
        xs_0, ys_0, xs_1, ys_1, xs_noise, ys_noise = for_plot_pdf_true_noise()

    # plot abc posteriors
    for t in range(1, history.max_t + 1):
        filename = label + "_kde_2d_" + str(t) + ".png"
        print(filename)
        if os.path.isfile(filename):
            continue
        df, w = history.get_distribution(m=0, t=t)
        axes = pyabc.visualization.plot_kde_matrix(
            df, w, #numx=1000, numy=1000,
            limits=limits, refval=th_true_noise)
        
        if show_true:
            axes[0, 0].plot(xs_0, ys_0, '-', color='k', alpha=0.75)
            axes[1, 1].plot(xs_1, ys_1, '-', color='k', alpha=0.75)
            axes[2, 2].plot(xs_noise, ys_noise, '-', color='k', alpha=0.75)

        plt.savefig(filename)
        plt.close()

