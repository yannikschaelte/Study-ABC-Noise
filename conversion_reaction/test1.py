import pyabc
from pyabc.visualization import plot_kde_2d
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import tempfile
from models import *

# VARIABLES

db_path = "sqlite:///db1.db"

# PERFORM ABC ANALYSIS

prior = pyabc.Distribution(th0=pyabc.RV('uniform', 0, 1),
                           th1=pyabc.RV('uniform', 0, 1))

distance = ArrayPNormDistance()
pop_size = 50
transition = pyabc.LocalTransition(k_fraction=.3)
eps = pyabc.MedianEpsilon(median_multiplier=.7)

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps)

abc.new(db_path, data_meas)

h = abc.run(minimum_epsilon=0, max_nr_populations=5)

# PLOT

t = h.max_t

ax = plot_kde_2d(*h.get_distribution(m=0, t=t),
                 'th0', 'th1',
                 xmin=0, xmax=1, numx=300,
                 ymin=0, ymax=1, numy=300)
ax.scatter([th0_true], [th1_true],
           color='C1',
           label='$\Theta$ true = {:.3f}, {:.3f}'.format(th0_true, th1_true))
ax.set_title("Posterior t={}".format(t))
ax.legend()
plt.savefig("test1")
