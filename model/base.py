import pyabc


class Model(ABC):

    def __init__(self):
        pass

    def get_p_true(self):
        return None

    def get_prior(self):
        raise NotImplementedError()

    def get_transitions(self):
        return pyabc.MultivariateNormalTransition()

    def get_distance(self):
        return pyabc.PNormDistance(p=2)

    def get_eps(self):
        return pyabc.MedianEpsilon()

    def get_pop_size(self):
        return 1000

    def get_n_pop(self):
        return 20

    def get_eps_min(self):
        return 0.0

    def get_id(self):
        raise NotImplementedError()

    def call(self):
        raise NotImplementedError()

    def call_noisy(self):
        raise NotImplementedError()

    def get_true_unnormalized_pdf_fun(y_obs):
        raise NotImplementedError()
