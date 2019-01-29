import matplotlib.pyplot as plt
import pyabc
import numpy as np
import matplotlib
from models import *

ss = []
n_tests = 5
labels = ["acc + decay", "decay", "exp_decay", "daly", "all"]
for _ in range(n_tests):
    ss.append(np.zeros(max_nr_populations))

h = pyabc.History("sqlite:///db5.db")

for j in range(n_tests):
    h.id = j + 1
    s = np.asarray(h.get_all_populations()['samples'][1:])
    ss[j][:len(s)] = s

print(ss)
print([sum(ss[j]) for j in range(n_tests)])

for k in range(max_nr_populations):
    plt.bar(np.arange(n_tests), tuple([ss[j][k] for j in range(n_tests)]),
            bottom=tuple([sum(ss[j][l] for l in range(k)) for j in range(n_tests)]))
plt.xticks(np.arange(n_tests), labels)
plt.title("total required samples")
plt.xlabel("Method")
plt.ylabel("Samples")
plt.savefig("total_samples.png")
