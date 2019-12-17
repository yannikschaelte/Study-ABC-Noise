import pyabc
import os
import matplotlib.pyplot as plt
from study_abc_noise.model import ConversionReactionUVarModelVars as ModelVars


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


pyabc.visualization.plot_sample_numbers(histories, labels)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(histories, labels, scale="log10")
plt.savefig("epsilons.png")
for h, label in zip(histories, labels):
    #for t in range(0, h.max_t + 1):
    #    pyabc.visualization.plot_histogram_matrix(h, t=t)
    #    plt.savefig("hist_" + label + ".png")
    _, axes = plt.subplots(1, 3)
    for (i, par) in enumerate(["p0", "p1", "std"]):
        for t in range(0, h.max_t + 1):
            pyabc.visualization.plot_kde_1d_highlevel(
                h, t=t, x=par, refval=gt_par, xmin=mv.limits[par][0], xmax=mv.limits[par][1],
                ax=axes[i], label=f"Iteration {t}")
    plt.legend()
    plt.gcf().set_size_inches((18, 6))
    plt.tight_layout()
    plt.savefig("kde_" + label +  ".png")
