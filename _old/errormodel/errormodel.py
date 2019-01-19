import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
x_axis = np.arange(-5, 5, 0.01)
plt.plot(x_axis, stats.norm.pdf(x_axis, 0, 1))
plt.savefig("normal_error")

plt.figure()
a=np.sqrt(3)
plt.step([-5, -a, a, 5], [0, 1/(2*a), 0, 0], where='post')
plt.savefig("step_error")
