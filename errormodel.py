import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

x_axis = np.arange(-5, 5, 0.01)
plt.plot(x_axis, stats.norm.pdf(x_axis, 0, 1))
plt.savefig("normal_error")
