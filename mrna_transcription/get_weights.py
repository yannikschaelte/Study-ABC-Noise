import pyabc
import numpy as np

db_file = "sqlite:///db2_1.db"
h = pyabc.History(db_file)
h.id = 1

for t in range(1, h.max_t + 1):
    d_weighted = h.get_weighted_distances(t)
    w = np.array(d_weighted['w'])
    n_eff = np.sum(w)**2 / np.sum(w**2)
    print(n_eff)
