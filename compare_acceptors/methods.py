"""
Methods to use
--------------

classical Daly: epsp0 fixed,  eps0 = 1000, Imax = 5e6, terminate if kt <= 0.001 or eps1 = 1
improved Daly
acceptance rate
decay
something else
and combinations



for each model: specific 
"""

import pyabc


class AcceptorStruct(dict):

    def __init__(method, n_pop, description, noisy):
        super().__init__(method=method, n_pop=n_pop, description=description, noisy=noisy)

class AcceptorFactory.

    def __init__():
        self.n_acc = 6

    def __len__(self):
        return self.n_acc

    def get(self, i):
        n_pop = None
        if i == 0:
            method = pyabc.StochasticAcceptor(
                temp_schemes=[pyabc.acceptor.scheme_acceptance_rate])
            descr = "Acceptance rate"
        elif i == 1:
            method = pyabc.StochasticAcceptor(
                temp_schemes=[pyabc.acceptor.scheme_acceptance_rate, pyabc.acceptor.scheme_decay])
            descr = "Acceptance rate + Decay"
        elif i == 2:
            method = pyabc.StochasticAcceptor(
                temp_schemes=[pyabc.acceptor.scheme_daly])
            descr = "Daly"
        return AcceptorStruct(method=method,n_pop=n_pop, descr=descr)
    
