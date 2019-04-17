import pyabc
import logging
import datetime


logger = logging.getLogger("Acceptor")
logger.setLevel(logging.DEBUG)
logger = logging.getLogger("Distance")
logger.setLevel(logging.DEBUG)


def create_sampler():
    # return pyabc.sampler.SingleCoreSampler()
    return pyabc.sampler.MulticoreEvalParallelSampler(n_procs=40)


def get_timestamp():
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
