import pyabc
import os
import matplotlib.pyplot as plt
from study_abc_noise.model import ConversionReactionModelVars


mv = ConversionReactionModelVars()


db_files = [f for f in os.listdir('.') if os.path.isfile(f) and "db_" in f]
print(f"Using db files {db_files}")

histories = []
labels = []
for db_file in db_files:
    id_ = db_file.split('__')[1]
    h = pyabc.History("sqlite:///" + db_file)
    h.id = 1
    histories.append(h)
    labels.append(id_)


gt_par = h.get_population(t=-1).get_list()[0].parameter


pyabc.visualization.plot_sample_numbers(histories, labels)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(histories, labels, scale="log10")
plt.savefig("epsilons.png")
for i in range(len(histories)):
    pyabc.visualization.plot_histogram_matrix(histories[i])
    plt.savefig("hist_" + str(i + 1) + ".png")
    df, w = histories[i].get_distribution()
    pyabc.visualization.plot_kde_matrix(df, w, refval=gt_par)  #, limits=mv.limits)
    plt.savefig("kde_" + str(i + 1) + ".png")
