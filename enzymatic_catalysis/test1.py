import pyabc
from models import *

db_path = "sqlite:///db1.db"

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps,
                   sampler=sampler)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)

visualize("test1", h)