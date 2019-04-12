import cloudpickle as pickle
import sys

with open(sys.argv[1], 'rb') as f:
    var = pickle.load(f)

print(var)
