import pyabc
import os
import matplotlib.pyplot as plt


db_files = [f for f in os.listdir('.') if os.path.isfile(f) and "db_" in f]
if not db_files:
    raise ValueError("No database found.")
db_file = "sqlite:///" + db_files[-1]
print(f"Using db file {db_file}")


histories = []
labels = []
h = pyabc.History(db_file)
h.id = 1
histories.append(h)
labels.append("Deterministic model + uniform acceptor")

h = pyabc.History(db_file)
h.id = 2
histories.append(h)
labels.append("Noisy model + uniform acceptor")

h = pyabc.History(db_file)
h.id = 3
histories.append(h)
labels.append("Deterministic model + stochastic acceptor")


gt_par = h.get_population(t=-1).get_list()[0].parameter


pyabc.visualization.plot_sample_numbers(histories, labels)
plt.savefig("samples.png")
pyabc.visualization.plot_epsilons(histories, labels, scale="log10")
plt.savefig("epsilons.png")
for i in range(len(histories)):
    pyabc.visualization.plot_histogram_matrix(histories[i])
    plt.savefig("hist_" + str(i + 1) + ".png")
    df, w = histories[i].get_distribution()
    pyabc.visualization.plot_kde_matrix(df, w, refval=gt_par)
    plt.savefig("kde_" + str(i + 1) + ".png")
