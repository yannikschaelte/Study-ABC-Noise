from study_abc_noise.read_pickle_file import read_pickle_file
import sys

var = read_pickle_file(sys.argv[1])
print(var)
