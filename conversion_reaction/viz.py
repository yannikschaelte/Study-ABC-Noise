import pyabc
from models import *


for i in range(1, 4):
    db_path = "sqlite:///db" + str(i) + ".db"
    h = pyabc.History(db_path)
    viz("test" + str(i), h)
