import pyabc
import os
import matplotlib.pyplot as plt
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

# plot epsilons
ax_epsilons = pyabc.visualization.plot_epsilons(histories, labels, scale='log10', size=(9, 6))
plt.savefig("fig_epsilons.png")

# find correct values
y_file = [f for f in os.listdir('data')][0]
y_obs = read_pickle_file("data/" + y_file)
posterior_scaled = get_posterior_scaled_1d(mv, y_obs)
xs = np.linspace(mv.limits['p0'][0], mv.limits['p0'][1], 200)
true_vals = [posterior_scaled([x]) for x in xs]

# plot histories
axes_hist = []
for i, (history, label) in enumerate(zip(histories, labels)):
    ax = pyabc.visualization.plot_histogram_1d(
        history, 'p0', xmin=mv.limits['p0'][0], xmax=mv.limits['p0'][1],
        size=(3, 2), bins=40)
    plt.plot(xs, true_vals, color='C4', label="True posterior")
    axes_hist.append(ax)
    plt.savefig(f"fig_hist_{i}.png")
