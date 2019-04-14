from study_abc_noise.optimize import get_optimal_kernel_value
from study_abc_noise.model import ConversionReactionModelVars as ModelVars


mv = ModelVars()
y_obs = mv.generate_data()
best_parameter, best_distance = get_optimal_kernel_value(mv, y_obs)
print(best_parameter, best_distance)
print(mv.p_true, mv.get_kernel()(y_obs, y_obs))
