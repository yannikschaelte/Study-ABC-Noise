import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pickle
import scipy.integrate as integrate


# VARIABLES

# noise variance
noise = 0.05
# gaussian error model
noise_model = np.random.randn

# prior
prior_lb = 0
prior_ub = 0.2
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)
                              for key in ['th0', 'th1']})


# MODEL

# timepoints
n_timepoints = 50
timepoints = np.arange(n_timepoints)

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
    # y = x(p)[1, :]
    y = x(p)[1, :]
    return {'y': y.flatten()} 


def model_random(p):
    """
    Observations. Account for noise.
    """
    #y = x(p)[1, :] + noise * noise_model(n_timepoints)
    y = x(p)[1, :] + noise * noise_model(1, n_timepoints)
    return {'y': y.flatten()}


def model_random_unknownnoise(p):
    """
    Observations when also noise std is a parameter.
    """
    #y = x(p)[1, :] + p['noise'] * noise_model(n_timepoints)
    y = x(p)[1, :] + p['noise'] * noise_model(1, n_timepoints)
    return {'y': y}


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
        y_obs_file = "y_obs.dat"
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
    for_plot_pdf_true_file = "for_plot_pdf_true.dat"
    try:
        xs_0, ys_0, xs_1, ys_1, zs = pickle.load(open(for_plot_pdf_true_file, 'rb'))
    except Exception as e:
        print(e)
        n_mesh = 200
        # th0
        def marginal_0(th0):
            return integrate.quad(lambda th1: pdf_true({'th0': th0, 'th1': th1}),
                                  prior_lb, prior_ub)[0]
        integral_0 = integrate.quad(marginal_0, prior_lb, prior_ub)[0]
        xs_0 = np.linspace(prior_lb, prior_ub, n_mesh)
        ys_0 = []
        for x in xs_0:
            ys_0.append(marginal_0(x) / integral_0)
    
        # th1
        def marginal_1(th1):
            return integrate.quad(lambda th0: pdf_true({'th0': th0, 'th1': th1}),
                                  prior_lb, prior_ub)[0]
        integral_1 = integrate.quad(marginal_1, prior_lb, prior_ub)[0]
        xs_1 = np.linspace(prior_lb, prior_ub, n_mesh)
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
        df, w = history.get_distribution(m=0, t=t)
        axes = pyabc.visualization.plot_kde_matrix(
            df, w, numx=1000, numy=1000,
            limits={key: (prior_lb, prior_ub)
                    for key in ['th0', 'th1']},
            refval=th_true)
    
        axes[0, 0].plot(xs_0, ys_0, '-', color='k', alpha=0.75)
        axes[1, 1].plot(xs_1, ys_1, '-', color='k', alpha=0.75)
        axes[1, 0].contour(xs_0, xs_1, zs.transpose(), cmap='Greys')
        plt.savefig(label + "_kde_2d_" + str(t))
        plt.close()


def viz_system(x):
    _, ax = plt.subplots()
    ax.plot(timepoints, x.transpose(), 'x-')
    plt.savefig("viz_system.png")


def viz_data(y, label):
    _, ax = plt.subplots()
    ax.plot(timepoints, y['y'], 'x-')
    plt.savefig("viz_data_" + label + ".png")

# pyabc parameters
distance = distance_l2
pop_size = 500  # 500
transition = pyabc.MultivariateNormalTransition()
eps = pyabc.MedianEpsilon()
max_nr_populations = 40  # 20
min_acceptance_rate = 1e-6
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=16)
