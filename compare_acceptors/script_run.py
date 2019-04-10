"""
Compare all the different approaches of stochastic acceptance kernels
(Dalyy, Acceptance rate, Decay, ESS). Compare this to a setting with
a deterministic acceptance step but a noisy model output.

Note: It is problematic how to best choose termination criteria (minimal
acceptance rate, iteration, target epsilon variance).
"""


from settings import n_rep, n_acc, n_pop
from models import Gaussian1DModel
import pyabc

# create models

models = []
models.append(Gaussian1DModel())
for n_t in [5, 10, 15, 20]:
    model = ConversionReactionModel()
    model.n_t = n_t
    model.ts = np.linspace(0, 30, n_t)
    models.append(model)

# create acceptors


# run

for model in models:
    for i_rep in range(n_rep):
        y_obs = model.generate_data()
        db_file = model.get_id() + f"_{n_rep}_{i_rep}_{n_acc}_{n_pop}" + get_timestamp()

        for acceptor in acceptors:
            id_ = f"{n_rep}_{n_acc}_{n_pop}"
            if isinstance(acceptor, pyabc.acceptor.StochasticAcceptor):
                call = model.call
                dist = model.get_kernel
                eps = pyabc.NoEpsilon()
            else:
                m = model.call_noisy
                dist = model.get_distance
                eps = model.get_eps()
            abc = pyabc.ABCSMC(
                models = m,
                parameter_priors = model.get_prior(),
                distance_function = dist,
                population_size = model.pop_size,
                transitions = model.get_transition(),
                eps = eps,
                acceptor = acceptor,
                sampler = create_sampler())
            abc.new(db_file, y_obs, gt_par=model.p_true)
            abc.run(minimum_epsilon=model.eps_min)
