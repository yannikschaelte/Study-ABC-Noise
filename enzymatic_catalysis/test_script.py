from models import *
import scipy as sp
import numpy as np

def cur_f(t, x):
    return f(t, x, th_true['th0'], th_true['th1'], th_true['th2'], th_true['th3'])

sol = sp.integrate.solve_ivp(fun=cur_f,
                             t_span=(min(t), max(t)),
                             y0=get_x0(),
                             t_eval=t)



print(sol)


