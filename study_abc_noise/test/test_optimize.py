from study_abc_noise.optimize import optimize
from study_abc_noise.model import ConversionReactionModelVars as ModelVars


mv = ModelVars()
y_obs = mv.generate_data()
lb = [mv.limits['p0'][0], mv.limits['p1'][0]]
ub = [mv.limits['p0'][1], mv.limits['p1'][1]]
def fun(p):
    p = {'p0': p[0], 'p1': p[1]}
    y = mv.get_model()(p)
    return - mv.get_kernel()(y, y_obs)
best_parameter, best_distance = optimize(fun, lb, ub)
print(best_parameter, best_distance)
print(mv.p_true, mv.get_kernel()(y_obs, y_obs))
