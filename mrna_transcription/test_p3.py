from model import *
import pyabc

db_file= "sqlite:///db2.db"

abc = pyabc.ABCSMC(models=model_random,
                   parameter_priors=prior,
                   distance_function=pyabc.NoDistance(),
                   summary_statistics=sumstat,
                   population_size=pop_size,
                   eps=pyabc.NoEpsilon(),
                   acceptor=stochastic_acceptor,
                   sampler=sampler)
abc.new(db=db_file, observed_sum_stat=sumstat(y_obs))
h = abc.run(minimum_epsilon=1, max_nr_populations=max_nr_populations)
h = pyabc.History(db_file)
#h.id = 1
visualize("viz3_" + str(n_r), h)
