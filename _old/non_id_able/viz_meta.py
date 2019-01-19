import pyabc
from models import viz_eps, viz_samples


list_label = ["Noisy model", "Deterministic model", "Stochastic acceptor"]
list_h = [pyabc.History("sqlite:///db1.db"), pyabc.History("sqlite:///db2.db"), pyabc.History("sqlite:///db3.db")]


viz_eps(list_h, list_label)
viz_samples(list_h, list_label)
