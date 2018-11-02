from gillespie import gillespie
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyabc
import pyabc.visualization
import pickle
import os

# number of replicates
n_r = 10
# true parameter
# note: translation rate = c>1 * mRNA decay rate
# 0.01, 0.12, 0.006, 0.0002
p_true = {'transcription': 1,
          'translation': 1,
          'mRNA decay': 0.1,
          'protein decay': 0.04}
# time
max_t = 2e2
n_timepoints = 20
timepoints = np.linspace(0, max_t, n_timepoints)
# reaction
pre = sp.array([[0, 0], [1, 0], [1, 0], [0, 1]], dtype=int)
post = sp.array([[1, 0], [1, 1], [0, 0], [0, 0]], dtype=int)
# initial state
x0 = sp.array([0, 0])


def model_raw(p):
    """
    Species:
    0. mRNA
    1. protein

    Reaction network:
    0. transcription:  0        --> mRNA
    1. translation:    mRNA     --> mRNA + protein
    2. mRNA decay:     mRNA     --> 0
    3. protein decay:  protein  --> 0
    """
    rates = np.array([p['transcription'],
                      p['translation'],
                      p['mRNA decay'],
                      p['protein decay']])
    ts = []
    xs = []
    for j in range(0, n_r):
        t, x = gillespie(x0, rates, pre, post, max_t)
        ts.append(t)
        xs.append(x)
    return {'ts': ts, 'xs': xs}


def _for_timepoints(y_raw):
    xs = np.zeros((n_r, n_timepoints, 2))
    for j in range(0, n_r):
        x_ind = sp.searchsorted(y_raw['ts'][j], timepoints)
        xs[j, ...] = y_raw['xs'][j][x_ind, ...]
    return {'t': timepoints, 'xs': xs}


def model(p):
    y_raw = model_raw(p)
    return _for_timepoints(y_raw)


def sumstat_p(y):
    """
    Assume only protein counts can be observed.
    """
    s = {}
    mean_p = np.mean(y['xs'][:, :, 1], axis=0).flatten()
    for j in range(0, n_timepoints):
        s['mean_p_' + str(j)] = mean_p[j]
    if n_r > 1:
        std_p = np.sqrt(np.var(y['xs'][:, :, 1], axis=0)).flatten()
        for j in range(0, n_timepoints):
            s['std_p_' + str(j)] = std_p[j]
    return s


def sumstat_mp(y):
    """
    Assume also mrna counts can be observed.
    """
    s = sumstat_p(y)
    mean_m = np.mean(y['xs'][:, :, 0], axis=0).flatten()
    for j in range(0, n_timepoints):
        s['mean_m_' + str(j)] = mean_m[j]
    if n_r > 1:
        std_m = np.sqrt(np.var(y['xs'][:, :, 1], axis=0)).flatten()
        for j in range(0, n_timepoints):
            s['std_m_' + str(j)] = std_m[j]
    return s


# observed data
_y_raw_obs = None
_y_obs = None


def get_y_raw_obs():
    global _y_raw_obs
    if _y_raw_obs is not None:
        return _y_raw_obs
    y_raw_file = "y_raw_" + str(n_r) + ".dat"
    try:
        y_raw_obs = pickle.load(open(y_raw_file, 'rb'))
    except Exception:
        y_raw_obs = model_raw(p_true)
        pickle.dump(y_raw_obs, open(y_raw_file, 'wb'))
    _y_raw_obs = y_raw_obs
    return _y_raw_obs


y_raw_obs = get_y_raw_obs()


def get_y_obs():
    global _y_obs
    if _y_obs is not None:
        return _y_obs
    y_file = "y_" + str(n_r) + ".dat"
    try:
        y_obs = pickle.load(open(y_file, 'rb'))
    except Exception:
        y_obs = _for_timepoints(y_raw_obs)
        pickle.dump(y_obs, open(y_file, 'wb'))
    _y_obs = y_obs
    return _y_obs


y_obs = get_y_obs()


# observed sumstats
sumstat_p_obs = sumstat_p(y_obs)
sumstat_mp_obs = sumstat_mp(y_obs)

# prior
limits = {key: (0, 3*rate) for key, rate in p_true.items()}
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', bounds[0], bounds[1])
                              for key, bounds in limits.items()})

# pyabc stuff
distance = pyabc.AdaptivePNormDistance(p=2)
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=16)
#sampler = pyabc.sampler.RedisEvalParallelSampler(host="wastl", port=8775)
max_nr_populations = 20
pop_size = 100


# visualize


def visualize_y_raw(y_raw):
    _, ax = plt.subplots()
    for j in range(0, n_r):
        ax.step(y_raw['ts'][j], y_raw['xs'][j][:, 0], color='b', alpha=0.2)
        ax.step(y_raw['ts'][j], y_raw['xs'][j][:, 1], color='r', alpha=0.2)
    
    custom_lines = [mpl.lines.Line2D([0], [0], color='b'),
                    mpl.lines.Line2D([0], [0], color='r')]
    ax.legend(custom_lines, ['mRNAs', 'Proteins'])

    ax.set_xlabel("Time [au]")
    ax.set_ylabel("Molecules")
    plt.savefig("obs_raw.png")


def visualize_y(y):
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    ax = axes[0]
    for j in range(0, n_r):
        ax.plot(y['t'], y['xs'][j, :, 0], color='b', alpha=0.2)
    mean_m = np.mean(y['xs'][:, :, 0], axis=0)
    ax.plot(y['t'], mean_m, 'b')
    std_m = np.sqrt(np.var(y['xs'][:, :, 0], axis=0))
    ax.fill_between(y['t'], mean_m - std_m , mean_m + std_m, color='b', alpha=0.1)
    ax.set_xlabel("Time [au]")
    ax.set_ylabel("mRNAs")

    ax = axes[1]
    for j in range(0, n_r):
        ax.plot(y['t'], y['xs'][j, :, 1], color='r', alpha=0.2)
    mean_p = np.mean(y['xs'][:, :, 1], axis=0)
    ax.plot(y['t'], mean_p, 'r')
    std_p = np.sqrt(np.var(y['xs'][:, :, 1], axis=0))
    ax.fill_between(y['t'], mean_p - std_p, mean_p + std_p, color='r', alpha=0.1)
    
    #custom_lines = [mpl.lines.Line2D([0], [0], color='b'),
    #                mpl.lines.Line2D([0], [0], color='r')]
    #ax.legend(custom_lines, ['mRNAs', 'Proteins'])
    
    ax.set_xlabel("Time [au]")
    ax.set_ylabel("Proteins")

    fig.tight_layout()

    plt.savefig("obs.png")


def visualize(label, history):
    for t in range(1, history.max_t + 1):
        df, w = history.get_distribution(m=0, t=t)
        print("t=" + str(t))
        ax = pyabc.visualization.plot_kde_matrix(
            df, w,
            limits=limits,
            refval=p_true)
        plt.savefig(label + "_kde_matrix_" + str(t) + ".png")
        plt.close()
