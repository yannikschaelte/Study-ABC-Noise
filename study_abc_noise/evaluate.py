from scipy.stats import kstest


def perform_kstest(rvs, cdf):
    """
    Perform Kolmogorov-Smirnov test for goodness of fit for each parameter
    separately.
    
    Parameters
    ----------

    rvs: The samples. List of dicts, the keys being the parameter ids.
    cdf: Dict of callables, the keys being the parameter ids, and returning
    the cumulative distribution functions for each.
    """
    # call kstest in each dimension
    # TODO: Correct for multiple testing?
