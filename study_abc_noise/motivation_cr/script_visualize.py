import pyabc
import os
import matplotlib.pyplot as plt
from study_abc_noise.model import ConversionReaction1dModelVars as ModelVars


mv = ModelVars()


db_files = [f for f in os.listdir('.') if os.path.isfile(f) and ".db" in f]
print(f"Using db files {db_files}")

histories = []
labels = []
for db_file in db_files:
    id_ = db_file[2:-3]
    h = pyabc.History("sqlite:///" + db_file)
    h.id = 1
    histories.append(h)
    labels.append(id_)
print(f"Labels: {labels}")

gt_par = h.get_population(t=-1).get_list()[0].parameter


pyabc.visualization.plot_sample_numbers(histories, labels)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(histories, labels, scale="log10")
plt.savefig("epsilons.png")
for h, label in zip(histories, labels):
    for t in range(h.max_t + 1):
        pyabc.visualization.plot_histogram_1d(
            h, x='p0', t=t, xmin=mv.limits['p0'][0], xmax=mv.limits['p0'][1],
            refval=mv.p_true, bins=40)
        plt.savefig(f"hist_{label}_{t}.png")
        pyabc.visualization.plot_kde_1d_highlevel(
            h, x='p0', t=t, xmin=mv.limits['p0'][0], xmax=mv.limits['p0'][1],
            refval=mv.p_true)
        plt.savefig(f"kde_{label}_{t}.png")
        plt.close()
pyabc.visualization.plot_effective_sample_sizes(histories, labels)
plt.savefig("ess.png")
