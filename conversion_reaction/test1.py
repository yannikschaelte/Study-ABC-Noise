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

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   population_size=pop_size,
                   transitions=transition,
                   eps=eps)

abc.new(db_path, get_y_meas())

h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)

# PLOT

visualize("test1", h)
