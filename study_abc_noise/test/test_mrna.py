from study_abc_noise.model import MRNATranscriptionModelVars as ModelVars

mv = ModelVars()
y_obs = mv.generate_data()
print(y_obs)
mv.get_model()(mv.p_true)
mv.get_model_noisy()(mv.p_true)
