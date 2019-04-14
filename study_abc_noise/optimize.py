import pypesto


def get_optimal_kernel_value(model_vars, y_obs):
    keys = model_vars.p_true.keys()
    lb = [model_vars.limits[key][0] for key in keys]
    ub = [model_vars.limits[key][1] for key in keys]
    model = model_vars.get_model()
    kernel = model_vars.get_kernel()

    def obj_fun(p):
        p = {key: p[i] for i, key in enumerate(keys)}
        y = model(p)
        return - kernel(y, y_obs)

    return _get_optimal_kernel_value(obj_fun, lb, ub)


def _get_optimal_kernel_value(fun, lb, ub):
    objective = pypesto.Objective(fun=fun)
    problem = pypesto.Problem(objective=objective, lb=lb, ub=ub)
    optimizer = pypesto.ScipyOptimizer(options={
        'maxiter': 10000, 'maxfun': 10000, 'disp': False})
    result = pypesto.minimize(
        problem=problem, optimizer=optimizer, n_starts=100)
    return result
