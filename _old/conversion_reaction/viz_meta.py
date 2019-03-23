import pyabc
from models import viz_eps, viz_samples
import matplotlib.pyplot as plt


list_label = ["Noisy model", "Deterministic model", "Stochastic acceptor", "Adaptive Stochastic acceptor"]
list_h = [pyabc.History("sqlite:///db1.db"), pyabc.History("sqlite:///db2.db"), pyabc.History("sqlite:///db3.db"), pyabc.History("sqlite:///db5.db")]


pyabc.visualization.plot_confidence_intervals(list_h[0], confidences=[0.5, 0.75, 0.9, 0.95])
plt.savefig("confidences.png")
pyabc.visualization.plot_sample_numbers(list_h, list_label, rotation=45)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(list_h, list_label)
plt.savefig("epsilons.png")
pyabc.visualization.plot_histogram_1d(list_h[0], 'th0', bins=40)
plt.savefig("hist_th0.png")
pyabc.visualization.plot_histogram_1d(list_h[0], 'th1', bins=40)
plt.savefig("hist_th1.png")
pyabc.visualization.plot_histogram_2d(list_h[0], 'th0', 'th1', bins=40)
plt.savefig("hist_2d.png")
pyabc.visualization.plot_histogram_matrix(list_h[0], bins=40)
plt.savefig("hist_matrix.png")
