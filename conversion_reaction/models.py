import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pickle

# VARIABLES

noise = 0.2
noise_model = np.random.randn

# prior
prior_lb = 0
prior_ub = 1
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)
                              for key in ['th0', 'th1']})

# MODEL

n_timepoints = 11
timepoints = sp.arange(n_timepoints)
x0 = sp.array([1, 0])


def f(x, t0, th0, th1):
    x0, x1 = x
    dx0 = - th0 * x0 + th1 * x1
    dx1 = th0 * x0 - th1 * x1
    return dx0, dx1


def x(p):
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


def model(p):
    y = x(p)[1, :]
    return {'y': y} 


def model_random(p):
    y = x(p)[1, :] + noise * noise_model(n_timepoints)
    return {'y': y}


def model_random_unknownnoise(p):
    y = x(p)[1, :] + p['noise'] * noise_model(n_timepoints)
    return {'y': y}


def distance(x, y):
    return sp.absolute(x['y'] - y['y']).sum()


class ArrayPNormDistance(pyabc.PNormDistance):

    def __init__(self):
        super().__init__(p=1)

    def initialize(self, t, sample_from_prior):
        sum_stats = []
        for sum_stat in sample_from_prior:
            sum_stats.append(normalize_sum_stat(sum_stat))
        super().initialize(t, sum_stats)

    def __call__(self, t, x, y):
        x = normalize_sum_stat(x)
        y = normalize_sum_stat(y)
        return super().__call__(t, x, y)


def normalize_sum_stat(x):
    x_flat = {}
    for key, value in x.items():
        for j in range(len(value)):
            x_flat[(key, j)] = value[j]
    return x_flat


# TRUE VALUES

th0_true, th1_true = sp.exp([-2.5, -2])
th_true = {'th0': th0_true, 'th1': th1_true}
y_true = model(th_true)

# MEASURED DATA


def get_y_meas():
    y_meas_file = "y_meas.dat"
    try:
        y_meas = pickle.load(open(y_meas_file, 'rb'))
    except Exception:
        y_meas = {'y': y_true['y'] + noise * noise_model(n_timepoints)}
        pickle.dump(y_meas, open(y_meas_file, 'wb'))

    return y_meas

# VISUALIZATION

def visualize(label, history):
    t = history.max_t

    df, w = history.get_distribution(m=0, t=t)
    ax = pyabc.visualization.plot_kde_2d(df, w,
                                         'th0', 'th1',
                                         xmin=prior_lb, xmax=prior_ub, numx=300,
                                         ymin=prior_lb, ymax=prior_ub, numy=300)
    ax.scatter([th0_true], [th1_true],
               color='C1',
               label='$\Theta$ true = {:.3f}, {:.3f}'.format(th0_true, th1_true))
    ax.set_title("Posterior t={}".format(t))
    ax.legend()
    plt.savefig(label + "_kde_2d_" + str(t))
    plt.close()


# pyabc parameters
distance = ArrayPNormDistance()
pop_size = 50
transition = pyabc.MultivariateNormalTransition()
eps = pyabc.MedianEpsilon()
max_nr_populations = 10
