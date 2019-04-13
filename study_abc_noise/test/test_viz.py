from study_abc_noise.model import ConversionReactionModelVars, NonIdAbleModelVars
import matplotlib.pyplot as plt


mv = ConversionReactionModelVars()
mv = NonIdAbleModelVars()
y = mv.generate_data()
mv.viz_x()
plt.show()
mv.viz_data_and_sim(y)
plt.show()
