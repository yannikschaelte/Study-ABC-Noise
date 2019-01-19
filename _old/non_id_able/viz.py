import pyabc
from models import *


for i in range(2, 3):
    db_path = "sqlite:///db" + str(i) + ".db"
    h = pyabc.History(db_path)
    viz("test" + str(i), h, show_true=False)
    #viz_fit(h.get_population_extended(), "test" + str(i))
