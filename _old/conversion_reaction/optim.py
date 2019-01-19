import pypesto
import pypesto.visualize
from models import *
import numpy as np
import matplotlib.pyplot as plt


def fun(x):
    th = {'th0': x[0], 'th1': x[1]}
    return np.power( (model(th)['y'] - y_obs['y']) / noise, 2 ).sum()


objective = pypesto.Objective(fun=fun)
problem = pypesto.Problem(objective=objective, lb=[0, 0], ub=[0.4, 0.4])
optimizer = pypesto.ScipyOptimizer()
result = pypesto.minimize(problem=problem, optimizer=optimizer, n_starts=20)

print(result.optimize_result.as_dataframe({'fval'}))

pypesto.visualize.waterfall(result)
plt.savefig("viz_optim.png")
