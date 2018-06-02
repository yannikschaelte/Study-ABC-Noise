import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pickle


# VARIABLES

noise = 0.5
noise_model = np.random.randn

# number of states
n_x = 4

# number of observables
n_y = 2

# (number of) timepoints
n_t = 50
t = np.linspace(0, 5, n_t)

# number of experiments (so far only 1 allowed)
# for n_e > 1, some loops and squeezes would be required
n_e = 1

# standard deviation and variance
std_dev = 0.05
variance = std_dev ** 2

# prior
prior_lb = -10
prior_ub = 5
prior = pyabc.Distribution(
        **{key: pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)
           for key in ['th0', 'th1', 'th2', 'th3']})

# pyabc parameters
pop_size = 50
transition = pyabc.MultivariateNormalTransition()
eps = pyabc.MedianEpsilon()
max_nr_populations = 10
# sampler = pyabc.sampler.RedisEvalParallelSampler(host="wastl", port=8765)
sampler = pyabc.sampler.SingleCoreSampler()

# true parameters
th_true = {'th0': 1.1770, 'th1': -2.3714, 'th2': -0.4827, 'th3': -5.5387}

# initial concentrations and measured data

_x0 = None
_y_meas = None


def get_x0():
    global _x0
    if _x0 is not None:
        return _x0

    x0_file = "x0.dat"
    try:
        x0 = pickle.load(open(x0_file, 'rb'))
    except Exception:
        x0 = np.exp(np.random.normal(0, 0.5, n_x))
        pickle.dump(x0, open(x0_file, 'wb'))
    
    _x0 = x0
    return x0


def get_y_meas():
    global _y_meas
    if _y_meas is not None:
        return _y_meas

    y_meas_file = "y_meas.dat"
    try:
        y_meas = pickle.load(open(y_meas_file, 'rb'))
    except Exception:
        y_true = model(th_true)
        y_meas = {'y0': y_true['y0'] + noise * noise_model(n_t),
                  'y3': y_true['y3'] + noise * noise_model(n_t)}
        pickle.dump(y_meas, open(y_meas_file, 'wb'))

    _y_meas = y_meas
    return y_meas


# MODEL


def f(t, x, th0, th1, th2, th3):
    x0, x1, x2, x3 = x
    dx0 = - th0*x0*x1 + th1*x2
    dx1 = - th0*x0*x1 + (th1+th2)*x2 - th3*x1*x3
    dx2 = th0*x0*x1 - (th1+th2)*x2 + th3*x1*x3
    dx3 = th2*x2 - th3*x1*x3
    
    return dx0, dx1, dx2, dx3


def x(p):
    th = np.exp([p['th0'], p['th1'], p['th2'], p['th3']])

    def cur_f(t, x):
        return f(t, x, th[0], th[1], th[2], th[3])

    sol = sp.integrate.solve_ivp(fun=cur_f, 
                                 t_span=(min(t), max(t)),
                                 y0=get_x0(),
                                 method='BDF',
                                 t_eval=t,
                                 rtol=1e-5,
                                 atol=1e-10)        
    
    return sol.y


def model(p):
    y = x(p)
    
    return {'y0': y[0, :], 'y3': y[3, :]}


def model_random(p):
    y = x(p)

    return {'y0': y[0, :] + noise * noise_model(n_t),
            'y3': y[3, :] + noise * noise_model(n_t)}

def normalize_sum_stats(x):
    x_flat = {}
    for key, value in x.items():
        for j in range(len(value)):
            x_flat[(key, j)] = value[j]
    
    return x_flat


def distance(x, y):
    x = normalize_sum_stats(x)
    y = normalize_sum_stats(y)

    return pow(
           sum(pow(abs(x[key]-y[key]), 2) 
               if key in x and key in y else 0
               for key in x),
           1/2)

# VISUALIZATION

def visualize(label, history):
    t = history.max_t

    df, w = history.get_distribution(m=0, t=t)
    pyabc.visualization.plot_kde_matrix(
            df, w, 
            limits={key: (prior_lb, prior_ub)
                    for key in ['th0', 'th1', 'th2', 'th3']},
            refval=th_true)
    plt.savefig(label + "_kde_matrix_" + str(t))
    plt.close()

