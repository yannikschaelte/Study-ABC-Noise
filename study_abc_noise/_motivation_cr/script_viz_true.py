import pyabc
import os
import matplotlib.pyplot as plt
from study_abc_noise.model.conversion_reaction import *
from study_abc_noise.read_pickle_file import read_pickle_file


mv = ConversionReaction1dModelVars()

files = [f for f in os.listdir('data')]
print(files)
y_file = files[0]
y_obs = read_pickle_file("data/" + y_file)
print(y_obs)
posterior_scaled = get_posterior_scaled_1d(mv, y_obs)

xs = np.linspace(mv.limits['p0'][0], mv.limits['p0'][1], 200)
true_vals = [posterior_scaled([x]) for x in xs]

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


gt_par = h.get_population(t=0).get_list()[0].parameter


pyabc.visualization.plot_sample_numbers(histories, labels)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(histories, labels, scale="log10")
plt.savefig("epsilons.png")
for h, label in zip(histories, labels):
    for t in range(h.max_t + 1):
        pyabc.visualization.plot_histogram_1d(h, 'p0', t=t, xmin=mv.limits['p0'][0], xmax=mv.limits['p0'][1], bins=40)
        plt.plot(xs, true_vals, color='C4', label="True posterior")
        plt.legend()
        plt.savefig(f"hist_true_{label}_{t}.png")
    for t in range(h.max_t + 1):
        df, w = h.get_distribution(t=t)
        pyabc.visualization.plot_kde_1d(df, w, 'p0', xmin=mv.limits['p0'][0], xmax=mv.limits['p0'][1])
        plt.plot(xs, true_vals, color='C4', label="True posterior")
        plt.legend()
        plt.savefig(f"kde_true_{label}_{t}.png")
        plt.close()
pyabc.visualization.plot_effective_sample_sizes(histories, labels)
plt.savefig("ess.png")
