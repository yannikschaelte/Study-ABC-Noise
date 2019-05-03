import ssa
from ..vars import ModelVars


class MjpModelVars(ModelVars):

    def __init__(self,
                 p_true, limits,
                 reactants, products,
                 n_t, t_max, x0, noise_success_probability,
                 n_acc: int = 1000, pdf_max: float = None):
        super().__init__(n_acc=n_acc, pdf_max=pdf_max, p_true=p_true)
        self.limits = limits
        self.reactants = reactants
        self.products = products
        self.n_t = n_t
        self.t_max = t_max
        self.x0 = x0
        selfoutput = ssa.output.ArrayOutput(
            np.linspace(0, self.t_max, self.n_t))
        self.noise_success_probabilty = noise_success_probability

    def get_id(self):
        return "MjpModelVars"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', bounds[0], bounds[1])
               for key, bounds in self.limits.items()})

    def get_distance(self):
        def l2(x, y):
            return np.sum([np.sum(np.power( (x[key] - y[key]), 2))
                           for key in self.limits])
        return l2

    def get_kernel(self):
        kernel = pyabc.distance.BinomialKernel(
            p=self.noise_success_probability,
            ret_scale=pyabc.distancel.RET_SCALE_LOG,
            pdf_max=self.pdf_max)
        return kernel

    def get_model_noisy(self):
        def model_noisy(p):
            y = self.get_model()(p)
            for key in y:
                successes = np.random.binomial(
                    n=y[key].flatten().astype(dtype=int),
                    p=self.noise_success_probability)
                y[key] = successes
            return y
        return model_noisy


class MjpConversion0ModelVars(MjpModelVars):

    def __init__(self):
        super().__init__(
            p_true = OrderedDict([('
