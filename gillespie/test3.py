"""
This is for testing noise data with a noise model.
"""

import models
import pyabc


db_file = "sqlite:///db3.db"
obs = models.Model1(noise_range=2).obs()

abc = pyabc.ABCSMC(models=models.Model1(noise_range=2),
                   parameter_priors=models.prior1,
                   distance_function=models.distance1,
                   population_size=models.pop_size)

abc.new(db_file, obs)

h = abc.run(minimum_epsilon=0, max_nr_populations=models.max_nr_populations)

models.Model1().visualize("test3", h)
