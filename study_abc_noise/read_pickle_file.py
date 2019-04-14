import cloudpickle as pickle


def read_pickle_file(filename):
    with open(filename, 'rb') as f:
        var = pickle.load(f)
    return var
