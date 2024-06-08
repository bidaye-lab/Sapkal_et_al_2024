import numpy as np

# registration
from pystackreg import StackReg

# parallel processing
from joblib import Parallel, delayed, parallel_backend


def get_tmats(arr, reg, n_bin=1):
    """Get transformation matrices for frames in array

    Apply one of the methods in pystackreg.StackReg to align frames in array
    to the mean of the array and return the transformation matrices for each frame.

    Optionally average frames in bins before registration (non-overlapping rolling mean)

    Parameters
    ----------
    arr : numpy.ndarray
        3D array with dimensions TYX
    reg : pystackreg.StackReg
        Method used for registration
    n_bin : int, optional
        Bin size for averaging before registration, by default 1

    Returns
    -------
    tmats : list of numpy.ndarray
        Transformation matrices for each frame
    """

    # bin time series
    n_y, n_x = arr.shape[-2:]
    arr_b = np.mean(np.reshape(arr, (-1, n_bin, n_y, n_x)), axis=1)

    # get mean image
    ref = np.mean(arr_b, axis=0)

    print(
        f"INFO getting transformation matrices using n_bin = {n_bin} and registration {reg}"
    )

    with parallel_backend("loky", n_jobs=-1):
        tmats_b = Parallel()(delayed(StackReg(reg).register)(ref, img) for img in arr_b)

    # non-binned tmats
    tmats = [i for i in tmats_b for _ in range(n_bin)]

    return tmats


def align(arr, tmats, reg):
    """Apply tranformation matrices to align array

    Uses the transformation matrices from get_tmats to align the frames in arr.

    Parameters
    ----------
    arr : numpy.ndarray
        3D array with dimensions TYX
    tmats : list of numpy.ndarray
        Transformation matrices for each frame
    reg : pysatckreg.StackReg
        Method used for registration

    Returns
    -------
    arr_a : numpy.ndarray
        Aligned 3D array with dimensions TYX
    """

    print(f"INFO aligning array using tranformation matrices and registration {reg}")

    with parallel_backend("loky", n_jobs=-1):
        arr_a = Parallel()(
            delayed(StackReg(reg).transform)(i, j) for i, j in zip(arr, tmats)
        )

    arr_a = np.array(arr_a)

    # replace 0 by mean
    arr_a[arr_a == 0] = np.mean(arr_a)

    return arr_a
