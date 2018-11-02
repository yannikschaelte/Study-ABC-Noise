import pyabc
from models import *


for i in range(1, 4):
    db_path = "sqlite:///db" + str(i) + ".db"
    h = pyabc.History(db_path)
    viz("test" + str(i), h)
    viz_fit(h.get_population_extended(), "test" + str(i))
