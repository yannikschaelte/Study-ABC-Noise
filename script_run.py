from settings import n_rep, n_acc, n_pop
from models import Gaussian1DModel
import pyabc

# generate models
models = [Gaussian1DModel()]

# run

for model in models:
    for i_rep in range(n_rep):
        id_ = f"{n_rep}_{n_acc}_{n_pop}_{i_rep}" 
        abc = pyabc.ABCSMC()
        y_obs = model.generate_data()
        model.save_data(y_obs, id_)
        db_name = model.get_history_identifier(id_)
        abc.new()
        abc.run(n_pop)
