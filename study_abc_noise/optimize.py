import pypesto
import os
import cloudpickle as pickle


def get_optimal_kernel_value(model_vars, y_obs):
    result = multistart_on_kernel(model_vars, y_obs)
    opt_p = result.optimize_result.get_for_key('x')[0]
    opt_fval = result.optimize_result.get_for_key('fval')[0]
    return opt_p, - opt_fval


def multistart_on_kernel(model_vars, y_obs, kernel):
    keys = model_vars.p_true.keys()
    lb = [model_vars.limits[key][0] for key in keys]
    ub = [model_vars.limits[key][1] for key in keys]
    model = model_vars.get_model()
    if kernel is None:
        kernel = model_vars.get_kernel()

    def obj_fun(p):
        p = {key: p[i] for i, key in enumerate(keys)}
        y = model(p)
        return - kernel(y, y_obs)

    return _multistart_on_kernel(obj_fun, lb, ub)


def _multistart_on_kernel(fun, lb, ub):
    objective = pypesto.Objective(fun=fun)
    problem = pypesto.Problem(objective=objective, lb=lb, ub=ub)
    optimizer = pypesto.ScipyOptimizer(options={
        'maxiter': 10000, 'maxfun': 10000, 'disp': False})
    result = pypesto.minimize(
        problem=problem, optimizer=optimizer, n_starts=100)
    return result


def get_and_store_optimal_kernel_value(model_vars, y_obs, i_rep): 
    pdf_max_id = f"pdf_max__{model_vars.get_id()}__{i_rep}"
    if not os.path.exists("data"):
        os.mkdir("data")
    filename = "data/" + pdf_max_id + ".dat"
    if os.path.isfile(filename):
        with open(filename, 'rb') as f:
            pdf_max = pickle.load(f)
            return pdf_max
    _, pdf_max = get_optimal_kernel_value(model_vars, y_obs)
    with open(filename, 'wb') as f:
        pickle.dump(pdf_max, f)
    return pdf_max
