import numpy as np
from scipy.ndimage import gaussian_filter, percentile_filter

# imageJ ROIs
from read_roi import read_roi_zip

# images
from PIL import Image, ImageDraw
from skimage.draw import polygon, ellipse
from skimage.color import gray2rgb
from matplotlib import colors


def read_imagej_rois(p_roi, res_yx):
    """Convert ROIs exported from ImageJ to boolean array

    ImageJ ROIs can be exported as a zip file. This function
    reads the zip file and converts the ROIs to a boolean array
    with the shape (n_roi, n_y, n_x).

    Only polygon and oval ROIs are supported.

    Parameters
    ----------
    p_roi : path-like
        Path to the ROI zip file
    res_yx : tuple of int
        Y and X resolution of the image in which to place the ROIs

    Returns
    -------
    r : numpy.ndarray
        3D array with shape (n_roi, n_y, n_x) with boolean values
    """

    d_roi = read_roi_zip(p_roi)
    n_y, n_x = res_yx
    # all-false array with shape: (n_roi, n_y,  n_x)
    r = np.zeros((len(d_roi), n_y, n_x)).astype(bool)

    # set rois mask to true
    for i, v in enumerate(d_roi.values()):
        if v["type"] == "polygon":
            x, y = polygon(v["x"], v["y"])

        elif v["type"] == "oval":
            r_x = v["width"] / 2
            r_y = v["height"] / 2
            c_x = v["left"] + r_x
            c_y = v["top"] + r_y
            x, y = ellipse(c_x, c_y, r_x, r_y)

        else:
            print(
                f'WARNING skipping ROI {i+1}, because it has type {v["type"]} not implemented'
            )
            continue

        # ignore out of bounds
        m = (y < n_y) * (x < n_x)
        x, y = x[m], y[m]

        r[i, y, x] = True

    return r


def draw_rois(img, rois):
    """Draw ROIs on grayscale image

    Image can be, for example, the mean image. ROIs are extracted
    using the function read_imagej_rois.

    Parameters
    ----------
    img : numpy.ndarray
        Grayscale image with shape (n_y, n_x)
    rois : numpy.ndarray
        3D array with shape (n_roi, n_y, n_x) with boolean values

    Returns
    -------
    img : numpy.ndarray
        RGB image with ROIs drawn on top
    """

    img = gray2rgb(img).astype("uint8")

    for i, r in enumerate(rois):
        rgb = tuple([int(f * 255) for f in colors.to_rgb(f"C{i}")])
        img[r] = rgb

    img = Image.fromarray(img)
    draw = ImageDraw.Draw(img)

    for i, r in enumerate(rois):
        r_ = gaussian_filter(r.astype(float), sigma=1)
        x, y = np.unravel_index(r_.argmax(), r_.shape)
        draw.text((y, x), str(i + 1))

    img = np.array(img)

    return img


def get_mean_trace(rois, arr):
    """Get mean trace per ROI from array

    For each ROI in the array rois, extract the mean trace
    from the array arr.
    Since the ROIs are boolean arrays, the mean trace is
    the mean of all pixels within the ROI per frame.

    Parameters
    ----------
    rois : numpy.ndarray
        3D array with shape (n_roi, n_y, n_x) with boolean values
    arr : numpy.ndarray
        3D array with shape (n_frames, n_y, n_x)

    Returns
    -------
    arr_m : numpy.ndarray
        2D array with shape (n_roi, n_frames) with mean values per ROI
    """

    # array with shape: (n_rois, n_samples)
    arr_m = np.zeros((len(rois), len(arr)))

    for i, r in enumerate(rois):
        arr_m[i] = np.mean(arr[:, r], axis=1)

    return arr_m


def subtract_background(arr, sigma=0):
    """Subtract background ROI from all other ROIs

    This function assumes that the last ROI in the array
    is the background ROI.

    Parameters
    ----------
    arr : numpy.ndarray
        2D array with shape (n_roi, n_frames) with mean values per ROI
    sigma : int, optional
        Standard deviation of Gaussian kernel for smoothing last ROI, by default 0

    Returns
    -------
    arr : numpy.ndarray
        2D array with shape (n_roi - 1, n_frames) with background ROI subtracted
    """
    a = arr[:-1]  # all but last ROI
    b = arr[-1]  # last

    if sigma:
        b = gaussian_filter(
            arr[-1], sigma=sigma
        )  # this would smooth background ROI before subtraction
    arr = a - b  # subtract last

    return arr


def subtract_baseline(arr, percentile, size):
    """Subtract baseline using percentile filter

    This function calculates the baseline as the `perc`th percentile
    within a rolling window of size `size`. The baseline is subtracted
    as (x - x0) / x0.

    Parameters
    ----------
    arr : numpy.ndarray
        2D array with shape (n_roi, n_frames) with mean values per ROI
    percentile : float
        Percentile used to calculate baseline
    size : int
        Rolling window size

    Returns
    -------
    dxx : numpy.ndarray
        Baseline-subtracted array with shape (n_roi, n_frames)
    """

    x = arr
    x0 = percentile_filter(x, percentile=percentile, size=(1, size))
    dxx = (x - x0) / x

    return dxx
