import matplotlib.pyplot as plt
import pyabc
import numpy as np
import matplotlib
from models import *

ss = []
n_tests = 5
labels = ["Adaptive", ""]
for _ in range(4):
    ss.append(np.zeros(max_nr_populations))

h = pyabc.History("sqlite:///db5.db")


h.id = 1
s2 = np.asarray(h.get_all_populations()['samples'][1:])

s1 = np.zeros_like(s2)
s3 = np.zeros_like(s2)
s4 = np.zeros_like(s2)
s5 = np.zeros_like(s2)

h.id = 1
s1[0] = h.get_all_populations()['samples'][1]

h.id = 3
tmp = np.asarray(h.get_all_populations()['samples'][1:])
s3[:len(tmp)] = tmp

h.id = 4
tmp = np.asarray(h.get_all_populations()['samples'][1:])
s4[:len(tmp)] = tmp

h = pyabc.History("sqlite:///db7.db")
h.id = 1
tmp = np.asarray(h.get_all_populations()['samples'][1:])
s5[:len(tmp)] = tmp

l = len(s2)
print(s1, s2, s3, s4, s5)
for j in range(l):
    plt.bar(np.arange(5),(s1[j], s2[j], s3[j], s4[j], s5[j]),
            bottom=(sum(s1[k] for k in range(j)), sum(s2[k] for k in range(j)), sum(s3[k] for k in range(j)), sum(s4[k] for k in range(j)), sum(s5[k] for k in range(j))))
plt.xticks(np.arange(5),("Rejection ABC", "Decay ABC-SMC", "Daly ABC-SMC", "Exponential Decay ABC-SMC", "Adaptive ABC-SMC"),rotation='vertical')
plt.title("Total required samples")
plt.ylabel("Samples")
plt.xlabel("Method")
plt.gcf().tight_layout()
plt.savefig("test6_bar")
