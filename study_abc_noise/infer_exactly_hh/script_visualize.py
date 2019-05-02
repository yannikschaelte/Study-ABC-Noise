import pyabc
import os
import matplotlib.pyplot as plt
from study_abc_noise.model import HodgkinHuxleyModelVars as ModelVars


mv = ModelVars()


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


pyabc.visualization.plot_sample_numbers(histories, labels, rotation=45)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(histories, labels, scale="log10")
plt.savefig("epsilons.png")
for h, label in zip(histories, labels):
    pyabc.visualization.plot_histogram_matrix(h)
    plt.savefig(f"hist_{label}.png")
    for t in range(h.max_t + 1):
        df, w = h.get_distribution(t=t)
        pyabc.visualization.plot_kde_matrix(df, w, refval=gt_par, limits=mv.limits)
        plt.savefig(f"kde_{label}_{t}.png")
        plt.close()
pyabc.visualization.plot_effective_sample_sizes(histories, labels, rotation=45)
plt.savefig("ess.png")
