import pypesto
import pypesto.visualize
from models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pickle

def pdf(y_0, y):
    v_0 = np.array(list(y_0.values()))
    v = np.array(list(y.values()))
    mean = np.zeros(n_timepoints)
    cov = noise**2 * np.eye(n_timepoints)
    return stats.multivariate_normal.pdf(v_0 - v, mean=mean, cov=cov)


def obj(x):
    x = {'th0': x[0], 'th1': x[1]}
    y = model(x)
    return - pdf(y_obs, y)


lb = [limits['th0'][0], limits['th1'][0]]
ub = [limits['th0'][1], limits['th1'][1]]

objective = pypesto.Objective(fun=obj)
problem = pypesto.Problem(objective=objective, lb=lb, ub=ub)
optimizer = pypesto.ScipyOptimizer(options={'maxiter': 10000, 'maxfun': 10000, 'disp': False})
result = pypesto.minimize(problem=problem, optimizer=optimizer, n_starts=100)

pdf_at_y_obs = pdf(y_obs, y_obs)
pdf_at_th_true = pdf(y_obs, model(th_true))
best_found_pdf = - result.optimize_result.get_for_key('fval')[0]

pickle.dump(best_found_pdf, open("best_found_pdf_" + str(noise) + "_" + str(n_timepoints) + ".dat", "wb"))

print(result.optimize_result.as_dataframe(['fval', 'n_fval', 'n_grad']))

print("pdf_at_y_obs, pdf_at_th_true, best_found_pdf, pdf_at_y_obs / best_found_pdf \n",
        "{:.3E} {:.3E} {:.3E} {:.3E}".format(pdf_at_y_obs, pdf_at_th_true, best_found_pdf, pdf_at_y_obs / best_found_pdf))

pypesto.visualize.waterfall(result)
plt.savefig("viz_optim.png")
