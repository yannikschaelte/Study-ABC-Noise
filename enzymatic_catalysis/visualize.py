import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
from models import *

def visualize_tmp(db_path, id):
    h = pyabc.History("sqlite:///" + db_path)
    h.id = id
    t = h.max_t
    print(t)
    df, w = h.get_distribution(m=0, t=t)
    pyabc.visualization.plot_kde_matrix(
            df, w,
            limits={key: (prior_lb, prior_ub)
                    for key in ['th0', 'th1', 'th2', 'th3']})
    plt.show()

visualize_tmp("db2.db", 3)
