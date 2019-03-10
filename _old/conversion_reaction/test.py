import pyabc
import time
import pickle

start_time = time.time()
h = pyabc.History("sqlite:///db5.db")
print("open db: ", time.time() - start_time)
start_time = time.time()
print(h.get_population_extended())
print("get_weighted_distances: ", time.time() - start_time)
start_time = time.time()
h.get_weighted_sum_stats()
print("get_weighted_sum_stats: ", time.time() - start_time)
start_time = time.time()
pop = h.get_population()
print(len(pickle.dumps(pop)))
print("get_population: ", time.time() - start_time)
for p in pop._list:
    print(p.parameter)
print(len(pop._list))
