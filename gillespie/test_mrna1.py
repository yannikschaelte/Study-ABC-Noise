"""
This is for testing non-noise data with a non-noise model.
"""

import models
import pyabc


db_file = "sqlite:///db_mrna1.db"
obs = models.MRNAModel(noise_range=1).obs()

abc = pyabc.ABCSMC(models=models.MRNAModel(noise_range=1),
                   parameter_priors=models.prior_mrna,
                   distance_function=models.distance_mrna,
                   population_size=models.pop_size)

abc.new(db_file, obs)

h = abc.run(minimum_epsilon=0, max_nr_populations=models.max_nr_populations_mrna)

models.MRNAModel().visualize("test_mrna1", h)
models.MRNAModel().visualize_animated("test_mrna1", h)
