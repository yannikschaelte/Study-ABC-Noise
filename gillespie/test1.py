import models
import pyabc


db_file = "sqlite:///db1.db"
obs = models.Model1().obs()

abc = pyabc.ABCSMC(models=models.Model1(noise_range=1),
                   parameter_priors=models.prior1,
                   distance_function=models.distance1,
                   population_size=models.pop_size)

abc.new(db_file, obs)

abc.run(minimum_epsilon=0, max_nr_populations=models.max_nr_populations)
