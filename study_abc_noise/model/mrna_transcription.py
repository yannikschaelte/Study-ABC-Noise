from ..vars import ModelVars
import numpy as np
import pyabc
import ssa
import matplotlib.pyplot as plt
from collections import OrderedDict


class MRNATranscriptionModelVars(ModelVars):

    def __init__(self):
        super().__init__(
            p_true = OrderedDict([('translation', 1),
                                  ('decay', 0.1),]))
        self.limits = OrderedDict([('transcription', (0, 3)),
                                   ('decay', (0, 0.2))])
        self.reactants = np.array([[0], [1]])
        self.products = np.array([[1], [0]])
        self.n_t = 100
        self.t_max = 100
        self.output = ssa.output.ArrayOutput(
                np.linspace(0, self.t_max, self.n_t))
        self.x0 = np.array([0])
        self.noise_success_probability = 0.95

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
            pdf_max=self.pdf_max)
        return kernel

    def get_model(self):
        def model(p):
            k = np.array(list(p.values()))  # order should be correct
            ssa_model = ssa.Model(
                reactants=self.reactants, products=self.products,
                t_max=self.t_max, x0=self.x0, k=k,
                output=self.output)
            result = ssa_model.simulate(n_reps=1)
            return {'t': result.list_ts[0],
                    'mrna': result.list_xs[0].flatten()}
        return model

    def get_model_noisy(self):
        def model_noisy(p):
            y = self.get_model()(p)
            successes = np.random.binomial(
                n=y['mrna'].flatten().astype(dtype=int),
                p=self.noise_success_probability)
            y['mrna'] = successes
            return y
        return model_noisy

    def generate_data(self):
        y = self.get_model_noisy()(self.p_true)
        return y
