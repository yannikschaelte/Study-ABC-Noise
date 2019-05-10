import pyabc
import os
import matplotlib.pyplot as plt
from study_abc_noise.model.conversion_reaction import *
from study_abc_noise.read_pickle_file import read_pickle_file


mv = ConversionReaction1dModelVars()

dir_ = 'results_20190510_t15/'
db_base = 'db_cr_10_0.02__'
db_end = '__0__0.db'
db_files = [dir_ + f for f in [db_base + 'deterministic' + db_end,
                               db_base + 'noisy_model' + db_end,
                               db_base + 'stochastic_acceptor' + db_end]]
labels = ['No noise, l2 distance',
          'Noisy model',
          'Modified acceptance']
ids = ['deterministic', 'noisy_model', 'stochastic_acceptor']

times = [3, 9, None]

# find correct values
y_file = [f for f in os.listdir('data')][0]
y_obs = read_pickle_file("data/" + y_file)
posterior_scaled = get_posterior_scaled_1d(mv, y_obs)
xs = np.linspace(mv.limits['p0'][0], mv.limits['p0'][1], 200)
true_vals = [posterior_scaled([x]) for x in xs]

color_truth = 'C0'

histories = [pyabc.History("sqlite:///" + f) for f in db_files]

for h, t, id_ in zip(histories, times, ids):
    ax = pyabc.visualization.plot_histogram_1d(
        h, 'p0', xmin=mv.limits['p0'][0], xmax=mv.limits['p0'][1],
        t=t,
        size=(3, 2), bins=40, color='C4')
    plt.plot(xs, true_vals, color=color_truth)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(None)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(f"hist_{id_}_{t}.png", format='png')
    plt.savefig(f"hist_{id_}_{t}.svg", format='svg')
