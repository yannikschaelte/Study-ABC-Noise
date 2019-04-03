import matplotlib.pyplot as plt
import pyabc
import numpy as np
import matplotlib

from models import *

histories = []
labels = []
for j  in range(1, 5):
    h = pyabc.History("sqlite:///db6.db")
    h.id = j
    histories.append(h)
labels.extend(["Rejection ABC-SMC", "Decay ABC-SMC", "Daly ABC-SMC", "Exp-Decay ABC-SMC"])

h = pyabc.History("sqlite:///db7.db")
h.id = 1
histories.append(h)
labels.append("Adaptive ABC-SMC")
h.id = 2
histories.append(h)
labels.append("Adaptive+Decay ABC-SMC")

h = pyabc.History("sqlite:///db8.db")
histories.append(h)
labels.append("ESS ABC-SMC")



pyabc.visualization.plot_histogram_1d(histories[1], x='th0', bins=30)
plt.savefig("hist.png")
pyabc.visualization.plot_sample_numbers(histories, labels, rotation=45)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(histories, labels, scale='log10')
plt.savefig("epsilons.png")
