import matplotlib.pyplot as plt
import pyabc
import numpy as np
import matplotlib
from models import *

ss = []
n_tests = 6
labels = ["acc-dec", "dec", "exp", "daly", "all", "all-bestc"]
for _ in range(n_tests):
    ss.append(np.zeros(max_nr_populations))

h = pyabc.History("sqlite:///db5.db")
for i in range(n_tests - 1):
    h.id=i + 1
    s = np.asarray(h.get_all_populations()['samples'][1:])
    ss[i][:len(s)] = s

h = pyabc.History("sqlite:///db6.db")
s = np.asarray(h.get_all_populations()['samples'][1:])
ss[5][:len(s)] = s

print(ss)

for j in range(max_nr_populations):
    plt.bar(np.arange(n_tests),tuple([ss[p][j] for p in range(n_tests)]),
            bottom=(tuple([sum(ss[p][k] for k in range(j)) for p in range(n_tests)])))
plt.xticks(np.arange(n_tests), labels, rotation='vertical')
plt.title("Total required samples")
plt.ylabel("Samples")
plt.xlabel("Method")
plt.gcf().tight_layout()
plt.savefig("test6_bar")
