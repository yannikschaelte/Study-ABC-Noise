import models
import pyabc


db_file = "sqlite:///db1.db"

abc = pyabc.ABCSMC(model=models.Model1(),
                   prior=prior,
                   distance=distance,
                   population_size=pop_size)

abc.new(db_file, models.obs)
