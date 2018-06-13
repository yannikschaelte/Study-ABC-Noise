import matplotlib.pyplot as plt
import pyabc
import numpy as np

h = pyabc.History("sqlite:///db6.db")
h.id = 1
s2 = np.asarray(h.get_all_populations()['samples'][1:])
s1 = np.zeros_like(s2)
h.id = 2
s1[0] = h.get_all_populations()['samples'][1]
l = len(s2)
print(s1, s2)
for j in range(l):
    plt.bar(np.arange(2),(s1[j], s2[j]),
            bottom=(sum(s1[k] for k in range(j)), sum(s2[k] for k in range(j))))
plt.xticks(np.arange(2),("Rejection ABC", "ABC-SMC"))
plt.ylabel("Samples")
plt.xlabel("Method")
plt.savefig("test6_bar")
