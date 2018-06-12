from gillespie import gillespie
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pyabc
import pyabc.visualization
import pickle
import os
import tempfile
import subprocess
tempdir = tempfile.mkdtemp()

# CONSTANTS

N_TEST_TIMES = 20


# MODELS

def noise_model(high, size):
    return np.random.randint(low=-high+1, high=high, size=size)  


class Model1:
    __name__ = 'Model1'
    x0 = sp.array([40, 3])
    pre = sp.array([[1, 1]], dtype=int)
    post = sp.array([0, 2])
    true_rate = {'r0': 2.3}
    limits = {'r0': (0, 10)}
    max_t = 0.1
    _obs = None

    def __init__(self, noise_range=1):
        self.noise_range = noise_range

    def extract_rates(self, par):
        return sp.array([par['r0']])

    def __call__(self, par):
        t, X = gillespie(self.x0,
                         self.extract_rates(par),
                         self.pre, self.post,
                         self.max_t)
        X = X + noise_model(high = self.noise_range, size=X.shape)       
        return {'t': t, 'X': X}
    
    def obs(self):
        if self._obs is not None:
            return self._obs

        obs_file = self.__name__ + "_" + str(self.noise_range) + ".dat"
        try:
            obs = pickle.load(open(obs_file, 'rb'))
        except Exception:
            obs = self(self.true_rate)
            pickle.dump(obs, open(obs_file, 'wb'))

        self._obs = obs
        return obs

    def visualize_obs(self, obs):
        ax = plt.gca()
        ax.step(obs['t'], obs['X'])
        ax.set_xlabel("Time")
        ax.set_ylabel("Concentration")
        plt.show()

    def visualize(self, label, history):
        t = history.max_t

        df, w = history.get_distribution(m=0, t=t)
        keys = list(self.true_rate.keys())
        if len(keys) == 1:
            ax = pyabc.visualization.plot_kde_1d(df, w,
                                                 keys[0],
                                                 xmin=0, xmax=10,
                                                 numx=300, refval=self.true_rate)
        else:
            ax = pyabc.visualization.plot_kde_matrix(df, w, limits=self.limits, refval=self.true_rate)

        plt.savefig(label + "_kde_" + str(t))
        plt.close()

    def visualize_animated(self, label, history):
        for t in range(history.n_populations):
            df, w = history.get_distribution(m=0, t=t)
            keys = list(self.true_rate.keys())
            if len(keys) == 1:
                pyabc.visualization.plot_kde_1d(df, w, keys[0], xmin=0, xmax=10,
                                                numx=300, refval=self.true_rate)
            else:
                pyabc.visualization.plot_kde_matrix(df, w, limits=self.limits, refval=self.true_rate)

            plt.title("Iteration " + str(t))
            plt.savefig(os.path.join(tempdir, f"{t:0>2}.png"))
            plt.close()
        subprocess.call("convert -delay 50 " + os.path.join(tempdir, "*.png") + " " + label + ".gif", shell=True)


class Model2(Model1):
    __name__ = 'Model2'
    pre = sp.array([[1, 0]], dtype=int)
    post = sp.array([[0, 1]])


class MRNAModel(Model1):
    """
    Species:
    0. mRNAs
    1. proteins

    Reaction network:
    0. transcription:      0        -->  mRNA
    1. translation:        mRNA     -->  mRNA + protein
    2. mRNA decay:         mRNA     -->  0
    3. protein decay:      protein  -->  0
    """
    __name__ = "MRNAModel"
    x0 = sp.array([0, 0])
    pre = sp.array([[0, 0], [1, 0], [1, 0], [0, 1]], dtype=int)
    post = sp.array([[1, 0], [1, 1], [0, 0], [0, 0]], dtype=int)
    true_rate = {'r0': 0.1, 'r1': 0.1, 'r2': 0.1, 'r3': 0.002}
    limits = {key: (0, 0.2) for key in true_rate}
    max_t = 1000

    def extract_rates(self, par):
        return sp.array([par['r0'], par['r1'], par['r2'], par['r3']])


class BirthDeathModel(Model1):
    x0 = sp.array([20])
    pre = sp.array([[0], [1]])
    post = sp.array([[1], [0]])
    true_rate = {'r0': 0, 'r1': 0.1}


# DISTANCES


def distance1(x, y):
    t_test_times = sp.linspace(0, Model1.max_t, N_TEST_TIMES)

    xt_ind = sp.searchsorted(x['t'], t_test_times) - 1
    yt_ind = sp.searchsorted(y['t'], t_test_times) - 1
    error = (sp.absolute(x['X'][:, 1][xt_ind]
                       - y['X'][:, 1][yt_ind]).sum()
             / t_test_times.size)
    return error


def distance_mrna(x, y):
    t_test_times = sp.linspace(0, MRNAModel.max_t, N_TEST_TIMES)

    xt_ind = sp.searchsorted(x['t'], t_test_times) - 1
    yt_ind = sp.searchsorted(y['t'], t_test_times) - 1
    error = (sp.absolute(x['X'][:, 1][xt_ind]
                       - y['X'][:, 1][yt_ind]).sum()
             / t_test_times.size)
    return error


# ABC STUFF


prior1 = pyabc.Distribution(**{key: pyabc.RV('uniform', a, b-a) for key, (a, b) in Model1.limits.items()})
prior_mrna = pyabc.Distribution(**{key: pyabc.RV('uniform', a, b-a) for key, (a, b) in MRNAModel.limits.items()})
pop_size = 100
max_nr_populations = 15
max_nr_populations_mrna = 15
