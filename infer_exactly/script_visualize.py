import pyabc
import os

db_files = [f for f in os.listdir('.') if isfile(f) and "db_" in f]
if not db_files:
    raise ValueError("No database found.")
db_file = "sqlite:///" + db_files[-1]
print(f"Taking db file {db_file}")

histories = []
labels = []
h = pyabc.History(db_file)
h.id = 0
histories.append(h)
labels.append("Deterministic model + uniform acceptor")

h = pyabc.History(db_file)
h.id = 1
histories.append(h)
labels.append("Noisy model + uniform acceptor")

h = pyabc.History(db_file)
h.id = 2
histories.append(h)
labels.append("Deterministic model + stochastic acceptor")


pyabc.visualize.plot_sample_numbers(histories, labels)
