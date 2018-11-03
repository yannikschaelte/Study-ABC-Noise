from model import *
import pyabc

db_file = "sqlite:///db2_" + str(n_r) + ".db"
#db_file = "sqlite:///db2_" + str(n_r) + "_noweights.db"

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   summary_statistics=sumstat_mp,
                   population_size=pop_size,
                   sampler=sampler)
#abc.new(db=db_file, observed_sum_stat=sumstat_mp_obs)
#h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
h = pyabc.History(db_file)
h.id = 1
visualize("pic2_" + str(n_r), h)
#visualize("pic2_" + str(n_r) + "_noweights", h)
