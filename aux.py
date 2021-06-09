import numpy as np
from scipy import linalg


def isPD(X):
    """Function to check if a matrix is positive definite or not. if cholesky can be calculated it is treated as PD.

    Parameters
    ----------
    X : np.ndarray
        Input matrix

    Returns
    -------
    bool
        is positive definite or not
    """
    try:
        _ = np.linalg.cholesky(X)
        return True
    except np.linalg.LinAlgError:
        return False


def matrix_decomposition(matrix):
    """

    Parameters
    ----------
    matrix : np.ndarray
        [description]

    Returns
    -------
    [type]
        [description]
    """
    L, D, perm = linalg.ldl(matrix)
    D[D < 1e-8] = 0
    independent_variables = np.where(np.diag(D) != 0)[0]
    cholesky = L @ np.sqrt(D)
    cholesky = cholesky[:, independent_variables]

    return cholesky
