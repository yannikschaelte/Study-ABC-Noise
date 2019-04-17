import pyabc
from abc import ABC
import os
import cloudpickle as pickle
from .util import create_sampler


class AnalysisVars(ABC):

    def __init__(
            self,
            get_acceptor,
            get_transition=None,
            get_eps=None,
            n_acc: int = 1000,
            n_pop: int = 20,
            eps_min: float = 0.0,
            min_acc_rate: float = 0.0,
            id_ = None):
        if get_acceptor is None:
            get_acceptor = lambda: pyabc.UniformAcceptor()
        self.get_acceptor = get_acceptor
        if get_transition is None:
            get_transition = lambda: pyabc.MultivariateNormalTransition()
        self.get_transition = get_transition
        if get_eps is None:
            get_eps = lambda: pyabc.MedianEpsilon()
        self.get_eps = get_eps
        self.n_acc = n_acc
        self.n_pop = n_pop
        self.eps_min = eps_min
        self.min_acc_rate = min_acc_rate
        self.id = id_

class ModelVars(ABC):

    def __init__(
            self,
            p_true: dict,
            n_acc: int = None,
            n_pop: int = None,
            pdf_max: float = None):
        self.n_acc = n_acc
        self.n_pop = n_pop
        self.p_true = p_true
        self.pdf_max = pdf_max

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
            min_acc_rate,
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
        self.min_acc_rate = min_acc_rate
        self.p_true = p_true
        self.y_obs = y_obs
        self.analysis_id = analysis_id
        self.model_id = model_id
        self.i_rep = i_rep

    @staticmethod
    def from_vars(analysis_vars: AnalysisVars,
                  model_vars: ModelVars,
                  i_rep: int = 0):
        acceptor = analysis_vars.get_acceptor()
        transition = analysis_vars.get_transition()
        eps_min = analysis_vars.eps_min
        if isinstance(acceptor, pyabc.StochasticAcceptor):
            eps = pyabc.NoEpsilon()
            model = model_vars.get_model()
            distance = model_vars.get_kernel()
            eps_min = 1.0
        else:
            eps = analysis_vars.get_eps()
            model = model_vars.get_model_noisy()
            distance = model_vars.get_distance()
        prior = model_vars.get_prior()
        sampler = create_sampler()
        n_acc = analysis_vars.n_acc if model_vars.n_acc is None \
            else model_vars.n_acc
        n_pop = analysis_vars.n_pop if model_vars.n_pop is None \
            else model_vars.n_pop
        min_acc_rate = analysis_vars.min_acc_rate
        p_true = model_vars.p_true
        y_obs = Task.get_data(model_vars, i_rep)
        analysis_id = analysis_vars.id
        model_id = model_vars.get_id()

        return Task(
            acceptor=acceptor, transition=transition, eps=eps,
            distance=distance, model=model, prior=prior, sampler=sampler,
            n_acc=n_acc, n_pop=n_pop, eps_min=eps_min,
            min_acc_rate=min_acc_rate, p_true=p_true,
            y_obs=y_obs,
            analysis_id=analysis_id, model_id=model_id,
            i_rep=i_rep)

    @staticmethod
    def get_data(model_vars, i_rep):
        data_id = f"{model_vars.get_id()}__{i_rep}"
        if not os.path.exists("data"):
            os.mkdir("data")
        filename = "data/" + data_id + ".dat"
        if os.path.isfile(filename):
            with open(filename, 'rb') as f:
                y_obs = pickle.load(f)
                return y_obs
        y_obs = model_vars.generate_data()
        with open(filename, 'wb') as f:
            pickle.dump(y_obs, f)
        return y_obs

    def execute(self):
        result_id = f"{self.model_id}__{self.analysis_id}__{self.i_rep}"
        db_file = f"db_{result_id}.db"

        print("Result id: ", result_id)
        if os.path.isfile(db_file):
            print("Skipping since exists already.")
            return

        abc = pyabc.ABCSMC(
            models = self.model,
            parameter_priors = self.prior,
            distance_function = self.distance,
            population_size = self.n_acc,
            transitions = self.transition,
            eps = self.eps,
            acceptor = self.acceptor,
            sampler = self.sampler)
        abc.new("sqlite:///" + db_file, self.y_obs, gt_par=self.p_true)
        abc.run(minimum_epsilon=self.eps_min, min_acceptance_rate=self.min_acc_rate, max_nr_populations=self.n_pop)
