from ..vars import ModelVars
import numpy as np
import pyabc
import ssa
import matplotlib.pyplot as plt
from collections import OrderedDict


class DimerizationModelVars(ModelVars):

    def __init__(self,
                 n_acc: int = 1000,
                 noise_success_probability: float = 0.90,
                 pdf_max: float = None,
                 n_t : int = 32,
                 t_max: int = 100):
        super().__init__(
            p_true = OrderedDict([('k1', np.log(1.0)), ('k2', np.log(0.04)),
                                  ('k3', np.log(0.002)), ('k4', np.log(0.5))]),
            n_acc = n_acc,
            pdf_max=pdf_max)
        self.limits = OrderedDict([('k1', (-2, 2)), ('k2', (-3, 1)),
                                   ('k3', ( -5, -1)), ('k4', (-3, 1))])
        self.reactants = np.array([[1, 0, 0], [0, 1, 0], [2, 0, 0], [0, 1, 0]])
        self.products = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [2, 0, 0]])
        self.n_t = n_t
        self.t_max = t_max
        self.output = ssa.output.ArrayOutput(
            np.linspace(0, self.t_max, self.n_t))
        self.x0 = np.array([1e4, 0, 0])
        self.noise_success_probability = noise_success_probability

    def get_id(self):
        return f"dimerization_{self.n_t}_{self.t_max}_" \
               f"{self.noise_success_probability}"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', bounds[0], bounds[1])
               for key, bounds in self.limits.items()})

    def get_distance(self):
        def l2(x, y):
            return sum([np.sum(np.power((x[key] - y[key]), 2)) for key in ['S1', 'S2']])
        return l2

    def get_kernel(self):
        kernel = pyabc.distance.BinomialKernel(
            p=self.noise_success_probability,
            keys=['S1'],# 'S2'],
            ret_scale=pyabc.distance.SCALE_LOG,
            pdf_max=self.pdf_max)
        return kernel

    def get_model(self):
        def model(p):
            k = np.array([p['k1'], p['k2'], p['k3'], p['k4']])
            k = np.exp(k)
            ssa_model = ssa.Model(
                reactants=self.reactants, products=self.products,
                t_max=self.t_max, x0=self.x0, k=k,
                output=self.output)
            result = ssa_model.simulate(n_reps=1)
            #import matplotlib.pyplot as plt
            #plt.plot(result.list_xs)
            res = {'t': result.list_ts[0],
                   'S1': result.list_xs[0][:, 0].flatten()
                   }# 'S2': result.list_xs[0][:, 1].flatten()}
            print(res)
            return res

        return model

    def get_model_noisy(self):
        def model_noisy(p):
            y = self.get_model()(p)
            for key in ['S1']:# 'S2']:
                successes = np.random.binomial(
                    n=y[key].flatten().astype(dtype=int),
                    p=self.noise_success_probability)
                y[key] = successes
            return y
        return model_noisy
