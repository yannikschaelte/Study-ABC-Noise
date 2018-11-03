import numpy as np
import scipy as sp


def h(x, pre, c):
    """
    Compute the reaction propensities.
    """
    return (x**pre).prod(1) * c


def gillespie(x, c, pre, post, max_t):
    """
    Gillespie simulation.

    Parameters
    ----------

    x: 1D array of size n_species
        The initial_numbers.

    c: 1D array of size n_reactions
        The reaction rates.

    pre: array of size n_reactions * n_species
        What is to be consumed.

    post: array of size n_reactions * n_species
        What is to be produced.

    max_t: int
        Simulate up to time max_t.

    Returns
    -------

    t, X: 1D array, 2D array
        t: The time points.
        X: The history of the species.
           ``X.shape == (t.size, x.size)``

    """

    t = 0
    t_store = [t]
    x_store = [x.copy()]
    S = post - pre

    while t < max_t:
        h_vec = h(x, pre, c)
        h0 = h_vec.sum()
        if np.isclose(h0, 0):
            break
        delta_t = sp.random.exponential(1 / h0)
        
        # no reaction can occur any more
        if not sp.isfinite(delta_t):
            t_store.append(max_t)
            x_store.append(x)
            break

        reaction = sp.random.choice(c.size, p=h_vec/h0)
        t = t + delta_t
        x = x + S[reaction]

        t_store.append(t)
        x_store.append(x)

    return sp.asarray(t_store), sp.asarray(x_store)
