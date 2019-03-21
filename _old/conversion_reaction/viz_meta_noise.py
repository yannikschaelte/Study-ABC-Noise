import pyabc
from models import viz_eps, viz_samples


list_label = ["Noisy model", "Deterministic model", ]
list_h = [pyabc.History("sqlite:///db_noise1.db"), pyabc.History("sqlite:///db_noise2.db"), ]


viz_eps(list_h, list_label)
viz_samples(list_h, list_label)
