import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pickle
import scipy.integrate as integrate
import os
import logging

logger = logging.getLogger("Acceptor")
logger.setLevel(logging.DEBUG)

# VARIABLES

# noise variance
noise = 0.02
# gaussian error model
noise_model = np.random.randn

# prior
limits = {'th0': (0, 0.4), 'th1': (0, 0.4)}

prior = pyabc.Distribution(**{key: pyabc.RV('uniform', bounds[0], bounds[1])
                              for key, bounds in limits.items()})


# MODEL

# timepoints
n_timepoints = 15
timepoints = np.linspace(0, 30, n_timepoints)

# initial concentrations (normalized to 1) 
x0 = np.array([1, 0])


def f(x, t0, th0, th1):
    """
    Differential equation.
    """
    x0, x1 = x
    dx0 = - th0 * x0 + th1 * x1
    dx1 = th0 * x0 - th1 * x1
    return dx0, dx1


def x_numeric(p):
    """
    States.
    Solve ODE numerically.
    """
    th0 = p['th0']
    th1 = p['th1']

    def cur_f(t, x):
        return f(x, t, th0, th1)
    
    sol = sp.integrate.solve_ivp(fun=cur_f, 
                                 t_span=(min(timepoints), max(timepoints)), 
                                 y0=x0, 
                                 method='BDF', 
                                 t_eval=timepoints)
    return sol.y


def x(p):
    """
    States.
    Use analytic solution of ODE.
    Returns array of shape n_x * n_t.
    """
    th0 = p['th0']
    th1 = p['th1']
    sol = np.zeros((2, n_timepoints))
    x0 = np.array([[1], [0]])
    for ix, t in enumerate(timepoints):
        e = np.exp( - (th0 + th1) * t)
        A = 1 / (- th0 - th1) * np.array([[- th1 - th0 * e, - th1 + th1 * e],
                                          [- th0 + th0 * e, - th0 - th1 * e]])

        sol[:, ix] = np.dot(A, x0).flatten()
    return sol


def model(p):
    """
    Observations. Do not account for noise.
    Assume only species 1 is observable.
    """
    y = x(p)[1, :]
    return {'y': y.flatten()} 


def model_random(p):
    """
    Observations. Account for noise.
    """
    y = x(p)[1, :] + noise * noise_model(1, n_timepoints)
    return {'y': y.flatten()}


def distance_l2(x, y):
    """
    Simple l2 distance.
    """
    return np.power( (x['y'] - y['y']) / noise, 2 ).sum()


# TRUE VALUES

th0_true, th1_true = [0.06, 0.08]
th_true = {'th0': th0_true, 'th1': th1_true}
y_true = model(th_true)


# OBSERVED DATA
_y_obs = None


def get_y_obs():
    global _y_obs
    if _y_obs is None:
        y_obs_file = "y_obs_" + str(noise) + "_" + str(n_timepoints) + ".dat"
        try:
            y_obs = pickle.load(open(y_obs_file, 'rb'))
        except Exception:
            y_obs = model_random(th_true)
            pickle.dump(y_obs, open(y_obs_file, 'wb'))
        _y_obs = y_obs
    return _y_obs


y_obs = get_y_obs()


def normal_dty_1d(y_bar, y, sigma):
    dty = ( 1 / np.sqrt( 2 * np.pi * sigma**2 ) 
            * np.exp( - ( (y_bar - y) / sigma)**2 / 2) )
    return dty


def normal_dty(y_bar, y, sigma):
    """
    y_bar: size dim
        point at which to evaluate the density
    y, sigma: size dim
        For N(y, sigma).
    """
    dim = len(y_bar)
    dties = np.zeros(dim)
    for j in range(dim):
        dties[j] = normal_dty_1d(y_bar[j], y[j], sigma[j])
    dty = np.prod(dties)
    return dty


def pdf_true(p):
    """
    Unscaled posterior density.
    """
    if type(p) is list:
        p = {'th0': p[0], 'th1': p[1]}
    
    # prior value
    prior_val = prior.pdf(p)
    
    # observed data
    y_bar = y_obs['y'].flatten()

    # data for p
    y = model(p)['y'].flatten()

    # likelihood value
    dim = len(y)
    sigma = noise * np.ones(dim)
    likelihood_val = normal_dty(y_bar, y, sigma)

    # posterior value
    unscaled_posterior = likelihood_val * prior_val

    return unscaled_posterior


# VISUALIZATION

def for_plot_pdf_true():
    for_plot_pdf_true_file = "for_plot_pdf_true_" + str(noise) + "_" + str(n_timepoints) + ".dat"
    try:
        xs_0, ys_0, xs_1, ys_1, zs = pickle.load(open(for_plot_pdf_true_file, 'rb'))
    except Exception as e:
        print(e)
        n_mesh = 200
        # th0
        def marginal_0(th0):
            return integrate.quad(lambda th1: pdf_true({'th0': th0, 'th1': th1}),
                                  limits['th1'][0], limits['th1'][1])[0]
        integral_0 = integrate.quad(marginal_0, limits['th0'][0], limits['th0'][1])[0]
        xs_0 = np.linspace(limits['th0'][0], limits['th0'][1], n_mesh)
        ys_0 = []
        for x in xs_0:
            ys_0.append(marginal_0(x) / integral_0)
    
        # th1
        def marginal_1(th1):
            return integrate.quad(lambda th0: pdf_true({'th0': th0, 'th1': th1}),
                                  limits['th0'][0], limits['th0'][1])[0]
        integral_1 = integrate.quad(marginal_1, limits['th1'][0], limits['th1'][1])[0]
        xs_1 = np.linspace(limits['th1'][0], limits['th1'][1], n_mesh)
        ys_1 = []
        for x in xs_1:
            ys_1.append(marginal_1(x) / integral_1)

        # th0, th1
        zs = np.zeros((n_mesh, n_mesh))
        for i0, v0 in enumerate(xs_0):
            for i1, v1 in enumerate(xs_1):
                zs[i0, i1] = pdf_true({'th0': v0, 'th1': v1})

        pickle.dump((xs_0, ys_0, xs_1, ys_1, zs), open(for_plot_pdf_true_file, 'wb'))

    return xs_0, ys_0, xs_1, ys_1, zs


def viz(label, history, show_true=True):
    # compute true posterior
    xs_0, ys_0, xs_1, ys_1, zs = for_plot_pdf_true()

    # plot abc posteriors
    for t in range(1, history.max_t + 1):
        filename = label + "_kde_2d_" + str(t) + ".png"
        print(filename)
        if os.path.isfile(filename):
            continue
        df, w = history.get_distribution(m=0, t=t)
        axes = pyabc.visualization.plot_kde_matrix(
            df, w, numx=1000, numy=1000,
            limits=limits,
            refval=th_true)
    
        axes[0, 0].plot(xs_0, ys_0, '-', color='k', alpha=0.75)
        axes[1, 1].plot(xs_1, ys_1, '-', color='k', alpha=0.75)
        axes[1, 0].contour(xs_0, xs_1, zs.transpose(), cmap='Greys')
        plt.savefig(filename)
        plt.close()


def viz_system(x):
    _, ax = plt.subplots()
    ax.plot(timepoints, x[0, :], 'x-', color='C1', label="A")
    ax.plot(timepoints, x[1, :], 'x-', color='C2', label="B")
    plt.xlabel("Time [au]")
    plt.ylabel("Concentration [au]")
    plt.legend()
    plt.savefig("viz_system.png")

def viz_data(y, label):
    _, ax = plt.subplots()
    ax.plot(timepoints, y['y'], 'x-', color='C2', label="B")
    plt.xlabel("Time [au]")
    plt.ylabel("Concentration [au]")
    plt.legend()
    plt.savefig("viz_data_" + label + ".png")

def viz_both(y_true, y_obs, label):
    _, ax = plt.subplots()
    ax.plot(timepoints, y_true['y'], 'x-', color='C0', label = "ODE solution")
    ax.plot(timepoints, y_obs['y'], 'x', color='C2', label = "measured data")
    plt.legend()
    plt.savefig("viz_data_" + label + ".png")

def viz_fit(df, label):
    _, ax = plt.subplots()
    ax.plot(timepoints, y_true['y'], 'x-', color='C0', label="Noise-free data")
    ax.plot(timepoints, y_obs['y'], 'x-', color='C2', label="Observation")
    n_particles = len(df)
    sumstats = df['sumstat_y']
    #for j in range(0, n_particles):
    #    sumstat = sumstats.iloc[j]
    #    ax.plot(timepoints, sumstat, 'x-', color='C3', alpha=0.2)
    mean = np.mean(np.array(sumstats), axis=0)
    std = np.sqrt(np.var(np.array(sumstats), axis=0))
    #ax.fill_between(timepoints, mean - std, mean + std, color='C3', alpha=0.1)
    err = [std, std]
    print(err)
    ax.errorbar(timepoints, mean, yerr=np.array(err), capsize=5, color='C3', label="ABC posterior")
    #ax.plot(timepoints, mean, 'x-', color='C4')
    plt.legend()
    plt.xlabel("Time [au]")
    plt.ylabel("Concentration [au]")
    plt.savefig("viz_fit_" + label + ".png")

def viz_eps(list_h, list_label):
    list_eps = []
    for h in list_h:
        list_eps.append(np.array(h.get_all_populations()['epsilon']))
    _, ax = plt.subplots()
    for ix, eps_schedule in enumerate(list_eps):
        ax.plot(np.log(eps_schedule[1:]), 'x-', label=list_label[ix])
    plt.xlabel("Population index")
    plt.ylabel("Log(Epsilon)")
    plt.legend()
    plt.savefig("viz_eps.png")


def viz_samples(list_h, list_label):
    list_samples = []
    for h in list_h:
        list_samples.append(np.array(h.get_all_populations()['samples']))
    _, ax = plt.subplots()
    for ix, sample_schedule in enumerate(list_samples):
        ax.plot(np.log(sample_schedule[1:]), 'x-', label=list_label[ix])
    plt.xlabel("Population index")
    plt.ylabel("Log(#Samples)")
    plt.legend()
    plt.savefig("viz_samples.png")


# pyabc parameters
distance = distance_l2
pop_size = 500  # 500
transition = pyabc.MultivariateNormalTransition()
eps = pyabc.MedianEpsilon()
max_nr_populations = 20  # 20
min_acceptance_rate = 1e-6
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=40)
