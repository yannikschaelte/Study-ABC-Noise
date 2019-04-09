import pyabc


def create_sampler():
    return pyabc.MulticoreEvalParallelSampler(n_procs=16)
