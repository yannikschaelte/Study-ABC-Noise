from ..vars import ModelVars
import numpy as np
import pyabc
import ssa
import copy
import matplotlib.pyplot as plt
from collections import OrderedDict


class MRNATranscriptionModelVars(ModelVars):

    def __init__(self,
                 n_acc: int = 1000,
                 noise_success_probability: float = 0.90,
                 pdf_max: float = None,
                 n_t : int = 30,
                 t_max: int = 90,
                 noise_model='binom'):
        super().__init__(
            p_true = OrderedDict([('transcription', 10),
                                  ('decay', 0.1),]),
            n_acc = n_acc,
            pdf_max=pdf_max)
        self.limits = OrderedDict([('transcription', (0, 30)),
                                   ('decay', (0, 0.2))])
        self.reactants = np.array([[0], [1]])
        self.products = np.array([[1], [0]])
        self.n_t = n_t
        self.t_max = t_max
        self.output = ssa.output.ArrayOutput(
            np.linspace(0, self.t_max, self.n_t))
        self.x0 = np.array([0])
        self.noise_success_probability = noise_success_probability
        self.noise_model = noise_model

    def get_id(self):
        return f"mrna_transcription_{self.n_t}_{self.t_max}_" \
               f"{self.noise_success_probability}"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', bounds[0], bounds[1])
               for key, bounds in self.limits.items()})

    def get_distance(self):
        def l2(x, y):
            return np.sum(np.power( (x['mrna'] - y['mrna']), 2))
        return l2

    def get_kernel(self):
        if self.noise_model == 'binom':
            kernel = pyabc.distance.BinomialKernel(
                p=self.noise_success_probability,
                keys=['mrna'],
                ret_scale=pyabc.distance.SCALE_LOG,
                pdf_max=self.pdf_max)
        elif self.noise_model == 'poisson':
            kernel = pyabc.distane.PoissonKernel(
                keys=['mrna'],
                ret_scale=pyabc.distance.SCALE_LOG)
        return kernel

    def get_model(self):
        def model(p):
            # anything (e.g. pd.DataSeries) -> OrderedDict
            p = OrderedDict([('transcription', p.get('transcription', 10)),
                             ('decay', p.get('decay', 0.1))])
            k = np.array(list(p.values()))  # order should be correct
            ssa_model = ssa.Model(
                reactants=self.reactants, products=self.products,
                t_max=self.t_max, x0=self.x0, k=k,
                output=self.output)
            result = ssa_model.simulate(n_reps=1)
            return {'t': result.list_ts[0],
                    'mrna': result.list_xs[0].flatten()}
        return model

    def add_noise_to_data(self, y):
        successes = copy.deepcopy(y)
        if self.noise_model == 'binom':
            successes['mrna'] = np.random.binomial(
                n=y['mrna'].flatten().astype(dtype=int),
                p=self.noise_success_probability)
        elif self.noise_model == 'poisson':
            successes['mrna'] = np.random.poisson(
                lam=y['mrna'].flatten().astype(dtype=int))
        return successes

    def get_model_noisy(self):
        def model_noisy(p):
            y = self.get_model()(p)
            if self.noise_model == 'binom':
                successes = np.random.binomial(
                    n=y['mrna'].flatten().astype(dtype=int),
                    p=self.noise_success_probability)
            elif self.noise_model == 'poisson':
                successes = np.random.poisson(
                    lam=y['mrna'].flatten().astype(dtype=int))
            y['mrna'] = successes
            return y
        return model_noisy


class MRNATranscription1dModelVars(ModelVars):

    def __init__(self,
                 n_acc: int = 1000,
                 noise_success_probability: float = 0.90,
                 pdf_max: float = None,
                 n_t : int = 30,
                 t_max: int = 90):
        super().__init__(
            p_true = OrderedDict([('transcription', 10)]),
            n_acc = n_acc,
            pdf_max=pdf_max)
        self.limits = OrderedDict([('transcription', (0, 30))])
        self.reactants = np.array([[0], [1]])
        self.products = np.array([[1], [0]])
        self.n_t = n_t
        self.t_max = t_max
        self.output = ssa.output.ArrayOutput(
            np.linspace(0, self.t_max, self.n_t))
        self.x0 = np.array([0])
        self.noise_success_probability = noise_success_probability

    def get_id(self):
        return f"mrna_transcription_{self.n_t}_{self.t_max}_" \
               f"{self.noise_success_probability}"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', bounds[0], bounds[1])
               for key, bounds in self.limits.items()})

    def get_distance(self):
        def l2(x, y):
            return np.sum(np.power( (x['mrna'] - y['mrna']), 2))
        return l2

    def get_kernel(self):
        kernel = pyabc.distance.BinomialKernel(
            p=self.noise_success_probability,
            keys=['mrna'],
            ret_scale=pyabc.distance.SCALE_LOG,
            pdf_max=self.pdf_max)
        return kernel

    def get_model(self):
        def model(p):
            # anything (e.g. pd.DataSeries) -> OrderedDict
            p = OrderedDict([('transcription', p.get('transcription', 10)),
                             ('decay', p.get('decay', 0.1))])
            k = np.array(list(p.values()))  # order should be correct
            ssa_model = ssa.Model(
                reactants=self.reactants, products=self.products,
                t_max=self.t_max, x0=self.x0, k=k,
                output=self.output)
            result = ssa_model.simulate(n_reps=1)
            return {'t': result.list_ts[0],
                    'mrna': result.list_xs[0].flatten()}
        return model

    def add_noise_to_data(self, y):
        successes = copy.deepcopy(y)
        successes['mrna'] = np.random.binomial(
            n=y['mrna'].flatten().astype(dtype=int),
            p=self.noise_success_probability)
        return successes

    def get_model_noisy(self):
        def model_noisy(p):
            y = self.get_model()(p)
            successes = np.random.binomial(
                n=y['mrna'].flatten().astype(dtype=int),
                p=self.noise_success_probability)
            y['mrna'] = successes
            return y
        return model_noisy


