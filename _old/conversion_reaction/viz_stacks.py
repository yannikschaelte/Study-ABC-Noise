import matplotlib.pyplot as plt
import pyabc
import numpy as np
import matplotlib
from models import *

n_tests = 15
histories = []
for idb in [5,5,5]:
    for i in range(int(n_tests / 3)):
        history = pyabc.History("sqlite:///db6_0.02_15.db")
        history.id = i + 1
        histories.append(history)    

labels = ["acc-dec", "dec", "exp", "daly", "all",
          "acc-dec-c-opt", "dec-c-opt", "exp-c-opt", "daly-c-opt", "all-c-opt",
          "acc_dep_c-learned", "dec-c-learned", "exp-c-learned", "daly-c-learned", "all-c-learned"]


pyabc.visualization.plot_sample_numbers(histories, labels, rotation=45)
plt.savefig("total_samples.png")
