import pyabc
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from study_abc_noise.model.conversion_reaction import *
from study_abc_noise.read_pickle_file import read_pickle_file

mv = ConversionReaction1dModelVars()

dir_ = 'results_20190506_t15/'
db_files = [dir_ + f for f in ['db_cr_10_0.02__deterministic__0__0.db',
                               'db_cr_10_0.02__deterministic_l1__0__0.db',
                               'db_cr_10_0.02__noisy_model__0__0.db']]
labels = ['No noise, l2 distance',
          'No noise, l1 distance',
          'Noisy model']

histories = [pyabc.History("sqlite:///" + f) for f in db_files]


ax_epsilons = plt.subplot2grid((2, 3), (0, 0), colspan=3)
axes_h = [plt.subplot2grid((2, 3), (1, i)) for i in [0, 1, 2]]
color_truth = 'C0'
colors_h = ['C1', 'C2', 'C3']

# plot epsilons
pyabc.visualization.plot_epsilons(histories, labels, scale='log10', size=(9, 6), colors=colors_h, ax=ax_epsilons)


# find correct values
y_file = [f for f in os.listdir('data')][0]
y_obs = read_pickle_file("data/" + y_file)
posterior_scaled = get_posterior_scaled_1d(mv, y_obs)
xs = np.linspace(mv.limits['p0'][0], mv.limits['p0'][1], 200)
true_vals = [posterior_scaled([x]) for x in xs]

# plot histories
for i, (history, label, color) in enumerate(zip(histories, labels, colors_h)):
    pyabc.visualization.plot_histogram_1d(
        history, 'p0', xmin=mv.limits['p0'][0], xmax=mv.limits['p0'][1],
        size=(3, 2), bins=40, ax=axes_h[i], color=color)
    axes_h[i].plot(xs, true_vals, color=color_truth, label="True posterior")
    axes_h[i].set_title(label)
    axes_h[i].set_xlabel("Parameter $\\theta_0$")

axes_h[0].set_ylabel("Density")
# create legend for hists
lines = [Line2D([0], [0], color=color_truth)]
axes_h[0].legend(lines, ["True posterior"])

plt.gcf().set_size_inches((9, 6))
plt.gcf().tight_layout()
plt.savefig("motivation.png")

