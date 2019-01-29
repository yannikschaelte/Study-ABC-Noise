import matplotlib.pyplot as plt
import pyabc
import numpy as np
import matplotlib
from models import *

ss = []
n_tests = 10
labels = ["acc-dec", "dec", "exp", "daly", "all",
          "acc-dec-c-opt", "dec-c-opt", "exp-c-opt", "daly-c-opt", "all-c-opt"]

for _ in range(n_tests):
    ss.append(np.zeros(max_nr_populations))

h = pyabc.History("sqlite:///db5.db")
for i in range(n_tests / 2):
    h.id = i + 1
    s = np.asarray(h.get_all_populations()['samples'][1:])
    ss[i][:len(s)] = s

h = pyabc.History("sqlite:///db6.db")
for i in range(n_tests / 2):
    h.id = i + 1
    s = np.asarray(h.get_all_populations()['samples'][1:])
    ss[int(n_tests / 2 + i)][:len(s)] = s

print(ss)

for j in range(max_nr_populations):
    plt.bar(np.arange(n_tests),tuple([ss[p][j] for p in range(n_tests)]),
            bottom=(tuple([sum(ss[p][k] for k in range(j)) for p in range(n_tests)])))
plt.xticks(np.arange(n_tests), labels, rotation='vertical')
plt.title("Total required samples")
plt.ylabel("Samples")
plt.xlabel("Method")
plt.ylabel("Samples")
plt.savefig("total_samples.png")
