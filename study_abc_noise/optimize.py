import pypesto


def optimize(fun, lb, ub):
    objective = pypesto.Objective(fun=fun)
    problem = pypesto.Problem(objective=objective, lb=lb, ub=ub)
    optimizer = pypesto.ScipyOptimizer(options={
        'maxiter': 10000, 'maxfun': 10000, 'disp': False})
    result = pypesto.minimize(
        problem=problem, optimizer=optimizer, n_starts=100)
    best_parameter = result.optimize_result.get_for_key('x')[0]
    best_distance = result.optimize_result.get_for_key('fval')[0]
    return best_parameter, best_distance
