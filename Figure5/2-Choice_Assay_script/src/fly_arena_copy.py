import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

import matplotlib.pyplot as plt
from matplotlib import colors

from scipy.io import loadmat
import h5py

from skvideo.io import vread, vwrite
from PIL import ImageDraw, Image


def load_track(matlab_file):
    """Load fly tracker data from matlab file

    Matlab file generated with
    https://github.com/jstaf/fly_tracker

    Parameters
    ----------
    matlab_file : path-like
        Path to matlab file

    Returns
    -------
    data : np.ndarray
        Array of shape (n_flies, n_frames, 2) with x and y coordinates
    """

    try:
        # matlab file pre 7.3
        m = loadmat(matlab_file, squeeze_me=True, struct_as_record=False)
        data = vars(m["trk"])["data"][:, :, [0, 1]]
    except NotImplementedError:
        # matlab file since 7.3
        with h5py.File(matlab_file, "r") as f:
            data = f["trk"]["data"][()].T[:, :, [0, 1]]

    return data


def load_separation(separation_file):
    """File with separation line between left and right side of chamber

    The line is defined by a series of points drawn for example with ImageJ.
    It does not need to span the whole chamber, but it need to cover the
    part of the y axis where the flies are walking.

    Parameters
    ----------
    separation_file : path-like
        Path to file with separation line

    Returns
    -------
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates
    """

    sep = np.loadtxt(separation_file).astype(int)

    return sep


def get_line(sep):
    """Get line separating left and right side of chamber based on separation points

    Line is defines a dict with ints for values and keys mapping y to x pixels

    Parameters
    ----------
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates

    Returns
    -------
    line : dict
        Dictionary with y as keys and x as values
    """

    # define line separating left and right
    line = dict()
    for p1, p2 in zip(sep, sep[1:]):
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        pxl = np.max([dx, dy])

        x = np.linspace(p1[0], p2[0], pxl + 1).astype(int)
        y = np.linspace(p1[1], p2[1], pxl + 1).astype(int)

        line = {**line, **{j: i for i, j in zip(x, y)}}

    return line


def get_sides(data, line):
    """Analyze fly trajectories and assign frames to left or right side of chamber

    Results are returned as a dictionary with fly as keys (0-indexed) and the
    following values:
    - x_pxl: x coordinates in pixels after removing nan
    - y_pxl: y coordinates in pixels after removing nan
    - left_mask: boolean mask with True for frames on left side of chamber
    - right_mask: boolean mask with True for frames on right side of chamber
    - nan_frames: number of frames with nan coordinates

    Parameters
    ----------
    data : np.ndarray
        Array of shape (n_flies, n_frames, 2) with x and y coordinates
    line : dict
        Dictionary with y as keys and x as values defining line separating left and right

    Returns
    -------
    results : dict
        Dictionary with fly as keys and dict with results as values
    """

    results = dict()

    for i_fly, xy in enumerate(data):
        # drop nan
        fnan = np.isnan(xy).any(axis=1)
        xy = xy[~fnan]

        # get x and y
        x = xy[:, 0]
        y = xy[:, 1]

        # masks for left and right of line
        bl = np.array([line[int(j)] >= i for i, j in zip(x, y)])
        br = ~bl

        results[i_fly] = {
            "x_pxl": x,
            "y_pxl": y,
            "left_mask": bl,
            "right_mask": br,
            "nan_frames": fnan.sum(),
        }

    return results

# def combined_plot_trajectory(p_png, sep, line, results, path=""):
    """Plot fly trajectory and separation line

    Parameters
    ----------
    p_png : path-like
        Example frame to plot ontop of
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates
    line : dict
        Dictionary with y as keys and x as values defining line separating left and right
    res : dict
        Dictionary with results for a single fly
    path : path-like, optional
        If not '', save plot to disk and close figure, by default ''
    """

    # load first frame
    img = Image.open(p_png)

    fig, ax = plt.subplots()

    # plot photo
    ax.imshow(img)

    # plot line defining points
    ax.scatter(sep[:, 0], sep[:, 1], zorder=99, color="k")

    # plot line definition
    y, x = line.keys(), line.values()
    ax.plot(x, y, color="C3", zorder=98)

    for fly, res in results.items():
        if fly in [0,1,5]:
            # plot trajectories
            cmap_paired = plt.cm.tab20.colors
            x, y = res["x_pxl"], res["y_pxl"]
            l, r = res["left_mask"], res["right_mask"]

            ax.scatter(x[l], y[l], marker=",", s=1, color=cmap_paired[2 * 0])
            ax.scatter(x[r], y[r], marker=",", s=1, color=cmap_paired[2 * 0 + 1])

    if path:
        fig.savefig(path)
        plt.close(fig)

def plot_trajectory(p_png, sep, line, res, path=""):
    """Plot fly trajectory and separation line

    Parameters
    ----------
    p_png : path-like
        Example frame to plot ontop of
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates
    line : dict
        Dictionary with y as keys and x as values defining line separating left and right
    res : dict
        Dictionary with results for a single fly
    path : path-like, optional
        If not '', save plot to disk and close figure, by default ''
    """

    # load first frame
    img = Image.open(p_png)

    fig, ax = plt.subplots()

    # plot photo
    ax.imshow(img)

    # plot line defining points
    ax.scatter(sep[:, 0], sep[:, 1], zorder=99, color="k")

    # plot line definition
    y, x = line.keys(), line.values()
    ax.plot(x, y, color="C3", zorder=98)

    # plot trajectories
    cmap_paired = plt.cm.tab20.colors
    x, y = res["x_pxl"], res["y_pxl"]
    l, r = res["left_mask"], res["right_mask"]

    ax.scatter(x[l], y[l], marker=",", s=1, color=cmap_paired[2 * 0])
    ax.scatter(x[r], y[r], marker=",", s=1, color=cmap_paired[2 * 0 + 1])

    if path:
        fig.savefig(path)
        plt.close(fig)


def summary_df(results):
    """Generate summary dataframe from results dict

       Dataframe has fly as index and the following columns:
        - n_frames_left: number of frames on left side of chamber
        - n_frames_right: number of frames on right side of chamber
        - ratio_frames_right_left: ratio of frames on right side of chamber to left side
        - n_frames_dropped: number of frames with nan coordinates


    Parameters
    ----------
    results : dict
        Dictionary with fly as keys and dict with results as values

    Returns
    -------
    df : pd.DataFrame
        Summary dataframe
    """

    df = pd.DataFrame(
        columns=[
            "n_frames_left",
            "n_frames_right",
            "ratio_frames_right_left",
            "n_frames_dropped",
        ]
    )
    df.index.name = "fly"

    for i, res in results.items():
        fly = i + 1

        ml, mr = res["left_mask"], res["right_mask"]

        nfl, nfr = ml.sum(), mr.sum()
        rfrl = nfr / nfl if nfl > 0 else 0
        nan = res["nan_frames"]

        
        df.loc[fly, :] = [
            nfl,
            nfr,
            rfrl,
            nan,
        ]

    return df
