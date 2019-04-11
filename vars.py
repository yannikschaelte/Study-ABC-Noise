import pyabc
from abc import ABC
import os
import cloudpickle as pickle


class AnalysisVars(ABC):

    def __init__(
            self,
            acceptor,
            transition=None,
            eps=None,
            n_acc: int = 1000,
            n_pop: int = 20,
            eps_min: flat = 0.0,
            id_ = None):
        if acceptor is None:
            acceptor = pyabc.UniformAcceptor()
        self.acceptor = acceptor
        if transition is None:
            transition = pyabc.MultivariateNormalTransition()
        self.transition = transition
        if eps is None:
            if isinstance(self.acceptor, pyabc.StochasticAcceptor):
                eps = pyabc.NoEpsilon()
            else:
                eps = pyabc.MedianEpsilon()
        self.eps = eps
        self.n_acc = n_acc
        self.n_pop = n_pop
        self.eps_min = eps_min
        self.id = id_

class ModelVars(ABC):

    def __init__(
            self,
            n_acc: int = None,
            n_pop: int = None,
            p_true: dict):
        self.n_acc = n_acc
        self.n_pop = n_pop
        self.p_true = p_true

    def get_id(self):
        raise NotImplementedError()

    def get_prior(self):
        raise NotImplementedError()

    def get_distance(self):
        """
        Distance to use for deterministic acceptance.
        """
        return pyabc.PNormDistance(p=2)

    def get_kernel(self):
        """
        Kernel to use for stochsatic acceptance.
        """
        raise NotImplementedError()

    def get_model(self):
        raise NotImplementedError()

    def get_model_noisy(self):
        raise NotImplementedError()

    def generate_data(self):
        raise NotImplementedError()


class Task(ABC):

    def __init__(
            self,
            acceptor,
            transition,
            eps,
            distance,
            model,
            prior,
            sampler,
            n_acc,
            n_pop,
            eps_min,
            p_true,
            y_obs,
            analysis_id,
            model_id,
            i_rep: int = 0):
        self.acceptor = acceptor
        self.transition = transition
        self.eps = eps
        self.distance = distance
        self.model = model
        self.prior = prior
        self.sampler = sampler
        self.n_acc = n_acc
        self.n_pop = n_pop
        self.eps_min = eps_min
        self.p_true = p_true
        self.y_obs = y_obs
        self.analysis_id = analysis_id
        self.model_id = model_id
        self.i_rep = i_rep

    @staticmethod
    def from_vars(analysis_vars: AnalysisVars,
                  model_vars: ModelVars,
                  i_rep: int = 0):
        acceptor = analysis_vars.acceptor
        transition = analysis_vars.transition
        if isinstance(acceptor, pyabc.StochasticAcceptor):
            eps = pyabc.NoEpsilon()
            model = model_vars.get_model()
            distance = model_vars.get_kernel()
        else:
            eps = analysis_vars.eps
            model = model_vars.get_model_noisy()
            distance = model_vars.get_distance()
        prior = model_vars.get_prior()
        sampler = create_sampler()
        n_acc = analysis_vars.n_acc if model_vars.n_acc is None \
            else model_vars.n_acc
        n_pop = analysis_vars.n_pop if model_vars.n_pop is None \
            else model_vars.n_pop
        eps_min = analysis_vars.eps_min
        p_true = model_vars.p_true
        y_obs = get_data(model_vars, i_rep)
        analysis_id = analysis_vars.id
        model_id = model_vars.get_id()

        return Task(
            acceptor=acceptor, transition=transition, eps=eps,
            distance=distance, model=model, prior=prior,
            n_acc=n_acc, n_pop=n_pop, eps_min=eps_min, p_true=p_true,
            y_obs=y_obs,
            analysis_id=analysis_id, model_id=model_id,
            i_rep=i_rep)

    @staticmethod
    def get_data(model_vars, i_rep):
        data_id = f"{model_vars.get_id()}_{i_rep}"
        os.mkdir("data")
        filename = "data/" + data_id + ".dat"
        if os.path.isfile(filename):
            with open(filename, 'rb') as f:
                return pickle.load(f)
        y_obs = self.model_vars.generate_data()
        with open(filename, 'wb') as f:
            pickle.dump(y_obs, f)
        return y_obs

    def execute(self):
        y_obs = self.get_data()
        result_id = f"{self.model_id}_{self.analysis_id}_{self.i_rep}"
        db_file = f"db_{result_id}.db"
        
        abc = pyabc.ABCSMC(
            models = self.model,
            parameter_priors = self.prior,
            distance_function = self.distance,
            population_size = self.n_pop,
            transitions = self.transition,
            eps = self.eps,
            acceptor = self.acceptor,
            sampler = self.sampler)
        abc.new(db_file, self.y_obs, gt_par=self.p_true)
        abc.run(minimum_epsilon=self.eps_min)
