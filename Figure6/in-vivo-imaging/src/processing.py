
import numpy as np
import pandas as pd

# signal processing
from scipy.ndimage import gaussian_filter
from scipy.stats import zscore
from scipy.ndimage import uniform_filter1d


#############
# tiff stacks

def split_channels(stack, n_z, n_ch):
    '''Split stack into 2 channels

    Assumes that stack is ordered as frames, channels, ypxl, xpxl (TCYX)
    This is equivalent to 
    - open TIF in ImageJ
    - Image > Hyperstacks > Stack to Hyperstack
        Order: xyczt (default)
        Channels: n_ch
        Slices: n_z

    Parameters
    ----------
    stack : numpy.ndarray
        4D array with dimensions TCYX
    n_z : int
        Number of z slices
    n_ch : int
        Number of channels

    Returns
    -------
    ch1 : numpy.ndarray
        4D array for 1st channel with dimensions TZYX
    ch2 : numpy.ndarray
        4D array for 2nd channel with dimensions TZYX
    '''

    # order is reversed in stack
    n_y, n_x = stack.shape[-2:]

    # reshape array channel, z, time, x, y
    arr = np.reshape(stack.T, (n_x, n_y, n_ch, n_z, -1), order='F')
    print(f'INFO reshaped transposed stack to {arr.shape}')

    # separate channels
    ch1, ch2 = arr[:, :, 0], arr[:, :, 1]
    ch1, ch2 = ch1.T, ch2.T
    print(f'INFO split into 2 channels and tranposed back: shape {ch1.shape}')

    return ch1, ch2

def maxproj_z(arr):
    '''Calculate maximum projection along z dimension

    Parameters
    ----------
    arr : numpy.ndarray
        4D array with dimensions TZYX

    Returns
    -------
    arr_z : numpy.ndarray
        3D array with dimensions TYX
    '''

    # collapse z: max procection
    arr_z = np.max(arr, axis=1)
    print(f'INFO max proj along z dim: shape {arr_z.shape}')
    
    return arr_z

def smooth_xy(arr, sigma):
    '''Smooth array along x and y with isotropic Gaussian filter

    Parameters
    ----------
    arr : numpy.ndarray
        3D array with dimensions TYX
    sigma : float
        Standard deviation of Gaussian kernel

    Returns
    -------
    arr_smth : numpy.ndarray
        Smoothed array
    '''
    
    # smooth only along x and y direction
    arr_smth = gaussian_filter(arr, sigma=(0, sigma, sigma))

    return arr_smth


#############
# stimulation
# TODO document once stim scrips is finalized

def get_onsets(arr, r_stim, r_rec, w=1, thresh=0.1, n_max=None):

    # convert to binary array based on thresh
    arr_bin = arr.copy()
    arr_bin[ arr_bin < thresh] = 0
    arr_bin[ arr_bin != 0 ] = 1
    
    # moving average and binary again
    arr_avg = uniform_filter1d(arr_bin, int(w*r_stim))
    arr_avg[ arr_avg != 0] = 1

    # indices for onset in arr_avg
    idx_avg = np.where(np.diff(arr_avg) == 1)[0]

    # onsets in arr
    idx = idx_avg + int(w*r_stim / 2)

    # convert from sample to s
    fs = idx / r_stim * r_rec
    f_max = len(arr) / r_stim * r_rec

    fs = np.round(fs, 0).astype(int)
    f_max = np.round(f_max, 0).astype(int)

    # limit number of stims
    fs = fs[:n_max]

    return fs, f_max

##########
# behavior
def resample(arr, sr, sr_new):
    '''Resample array to new sample rate

    Parameters
    ----------
    arr : numpy.ndarray
        1D array to be resampled
    sr : float
        Sample rate of input array
    sr_new : float
        New sample rate

    Returns
    -------
    y_new : numpy.ndarray
        Resampled array
    '''

    l = len(arr)
    x = np.linspace(0, 1, l)

    l_new = int(l / sr * sr_new)
    x_new = np.linspace(0, 1, l_new)

    y_new = np.interp(x_new, x, arr)

    return y_new

def resample_to_behavior(ca, beh, ball, f_ca, f_ball, f_beh, ca2=None):
    '''Resample imaging and ball velocity data to match behavior sample rate

    Returns dataframe with columns:
    - roi_{i}: fluorescence trace for ROI i
    - ball_{j}: ball velocity in x, y, or z direction
    - beh_{k}: behavior event k

    Optionally, normalize 1st channel by 2nd channel data (if ca2 is not None).


    Parameters
    ----------
    ca : numpy.ndarray
        2D array with shape (n_roi, n_frames) containing fluorescence traces
    beh : dict
        Mapping between behavior event name and data as returned by `load_behavior`
    ball : numpy.ndarray
        Ball velocity data in x, y, and z direction as returned by `load_ball`
    f_ca : float
        Sample rate of imaging data
    f_ball : float
        Sample rate of ball velocity data
    f_beh : float
        Sample rate of behavior data
    ca2 : numpy.array, optional
        Imaging data for 2nd channel. If given, calculate 1st/2nd channel ratio instead, by default None

    Returns
    -------
    df : pandas.DataFrame
        Dataframe with resampled data
    '''

    y = np.zeros_like(ca[0])
    y = resample(y, f_ca, f_beh)
    x = np.arange(len(y))

    df = pd.DataFrame(index=x)

    for i, y in enumerate(ca):
        y = resample(y, f_ca, f_beh)
        if ca2 is not None:
            y2 = resample(ca2[i], f_ca, f_beh)
            y /= y2
        df.loc[:, f'roi_{i+1}'] = y

    for j, y in zip('xyz', ball.T):
        y = resample(y, f_ball, f_beh)
        df.loc[:, f'ball_{j}'] = y[x]

    for k, v in beh.items():

        # boxcar
        y = np.zeros_like(x)

        arr = np.array(v, ndmin=2)
        if arr.any():
            for f_i, f_f, _ in arr:
                y[f_i:f_f+1] = 1

        df.loc[:, f'beh_{k}'] = y

        # delta 
        yi = np.zeros_like(x)
        yf = np.zeros_like(x)

        arr = np.array(v, ndmin=2)
        if arr.any():
            for f_i, f_f, _ in arr:
                yi[f_i:f_i+1] = 1
                yf[f_f:f_f+1] = 1

        df.loc[:, f'behi_{k}'] = yi
        df.loc[:, f'behf_{k}'] = yf

    return df

def align2events(df, beh, f, dt, use_z_roi=True):
    '''Align data to behavior events

    Selects data around behavior events and aligns to event onset.
    Returns new dataframe with aligned data.

    This will select all columns that start with `roi_` (or `z_roi_`) and `ball_` 
    and align +- `dt` to the behavior event `beh`. If trial start or end is within
    `dt` of the behavior event, the event is skipped.

    If `use_z_roi` is True, align z-scored traces, i.e. z_roi columns, otherwise align roi columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with resampled data as returned by `resample_to_behavior`
    beh : str
        Name of behavior column to align to
    f : float
        Sample rate of behavior data
    dt : float
        Time window in s around behavior event
    use_z_roi : bool, optional
        Use z-scored traces, by default True

    Returns
    -------
    data : pandas.DataFrame
        Data aligned to behavior events
    '''

    if use_z_roi:
        cols_roi = [ c for c in df.columns if c.startswith('z_roi_') ]
    else:
        cols_roi = [ c for c in df.columns if c.startswith('roi_') ]

    cols_ball = [ c for c in df.columns if c.startswith('ball_') ]
    cols = cols_roi + cols_ball
    
    dn = int(dt * f)
    s = np.arange(-dn, dn + 1)
    l = []
    n_good, n_bad = 0, 0

    for (fly, trial), d in df.groupby(['fly', 'trial']):
        ds = d.loc[:, beh]
        i_ons = np.flatnonzero(np.diff(ds) == 1) + 1
        idx_ons = ds.iloc[i_ons].index

        for idx in idx_ons:
            try: 
                for c in cols:
                    y = d.loc[idx - dn : idx + dn, c]
                    df_new = pd.DataFrame({
                        's': s,
                        't': s / f,
                        'fly': fly,
                        'trial': trial,
                        'y': y,
                        'match': c,
                        'onset': idx,
                    })
                    l.append(df_new)
                n_good += 1

            except ValueError:
                print(f'WARNING skipping fly {fly} trial {trial} event {idx}')
                n_bad += 1

    data = pd.concat(l)
    data.attrs['n_bad'] = n_bad
    data.attrs['n_good'] = n_good
    data.attrs['beh'] = beh
    data.attrs['n_fly'] = len(data.groupby(['fly']).groups.keys())
    data.attrs['n_trial'] = len(data.groupby(['fly', 'trial']).groups.keys())

    return data

def zscore_cols(df, col_start):
    '''Z-score columns in dataframe that start with `col_start`

    This will match `col_start` with `df.columns` and add new columns with prefix `z_`
    that contain the z-scored data.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with data
    col_start : str
        Prefix of columns to z-score

    Returns
    -------
    df : pandas.DataFrame
        Dataframe with z-scored data added
    '''

    cols = [ c for c in df.columns if c.startswith(col_start) ]
    cols_z = [ f'z_{c}' for c in cols ]

    df.loc[:, cols_z] = df.loc[:, cols].apply(zscore).values

    return df

def ca_kernel(tau_on, tau_off, f):
    '''Definition of the Calcium sensor kernel

    Parameters
    ----------
    tau_on : float
        Exponent for signal rise
    tau_off : float
        Expoenent for signal decay
    f : float
        Sample rate

    Returns
    -------
    y : numpy.ndarray
        1D array describing kernel
    '''

    x = np.arange(0, 10*tau_off, 1 / f)
    
    ca_trace = lambda a, b, c: - np.exp(-a / b) + np.exp(-a / c)
    y = ca_trace(x, tau_on, tau_off)

    # onset aligned to center of kernel
    y = np.pad(y, pad_width=(len(y), 0))

    # normalize area under curve to 1
    y = y / y.sum()

    return y

def convolute_ca_kernel(df, f, tau_on, tau_off):
    '''Convolute ball velocity and behavior data with Calcium kernel

    This will convolute all columns that start with `ball_` or `beh_`
    with the Calcium kernel and add new columns with prefix `conv_`.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with data
    f : float
        Sample rate of data
    tau_on : float
        Rise time in s of Calcium kernel
    tau_off : float
        Decay time in s of Calcium kernel

    Returns
    -------
    df : pandas.DataFrame
        Dataframe with convoluted data added
    '''

    kern = ca_kernel(tau_on, tau_off, f)

    cols = [ c for c in df.columns if c.startswith('ball_') or c.startswith('beh')]
    cols_new = [ f'conv_{c}' for c in cols ]
    df.loc[df.index, cols_new] = df.loc[df.index, cols].apply(np.convolve, axis=0, **{'v': kern, 'mode': 'same'}).values

    return df

def calculate_pearson(df, beh):
    '''Calculate Pearson correlation coefficients

    This will select all columns that start with `z_roi_`, `z_conv_ball_`, and `conv_{beh}_` 
    and calculate the Pearson correlation coefficient between them.

    Note that in the output dataframe, the prefix `z_`, `z_conv_` and `conv_` are removed.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with data
    beh : str
        String to match behavior colums `conv_{beh}_`, e.g. 'beh', 'behi', 'behf'

    Returns
    -------
    d : pandas.DataFrame
        Matrix with Pearson correlation coefficients
    '''


    c_roi = [ c for c in df.columns if c.startswith('z_roi_')]
    c_ball = [ c for c in df.columns if c.startswith('z_conv_ball_') ]
    c_beh = [ c for c in df.columns if c.startswith(f'conv_{beh}_')]
    cols = c_roi + c_ball + c_beh

    with np.errstate(divide='ignore', invalid='ignore'): # ignore division by 0 warning: TODO better check array before calling corrcoef
        c = np.corrcoef(df.loc[:, cols].T.values)
        # c = df.loc[:, cols].corr(method='pearson').values

    c[np.diag_indices(len(c))] = np.nan

    c_roi = [ c.replace('z_', '') for c in c_roi ] 
    c_ball = [ c.replace('z_conv_', '') for c in c_ball ]
    c_beh = [ c.replace(f'conv_{beh}_', '') for c in c_beh ]
    cols = c_roi + c_ball + c_beh

    d = pd.DataFrame(data=c, index=cols, columns=cols)

    return d

def calculate_ccf(df, dt, f, col1, col2, col2_):
    '''Calculate cross-correlation function (CCF)

    Match the start of column names with `col1`, `col2`, and `col2_` to define
    two sets of columns for which to calculate the CCF. The CCF is calculated
    for all combinations of columns in the two sets.
    The CCF is calculated for each fly, trial, and condition separately.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with time series data
    dt : float
        Time lag for CCF in s
    f : float
        Sample rate of data in Hz
    col1 : str
        Match with start of column name defines first set of columns
    col2 : str
        Match with start of column name defines second set of columns
    col2_ : str
        Match with start of column name also defines second set of columns

    Returns
    -------
    df_c : pandas.DataFrame
        Dataframe with CCF data
    '''
    n = 2* dt * f + 1
    t = np.linspace(-dt, dt, n)

    cols_1 = [ c for c in df.columns if c.startswith(col1)]
    cols_2 = [ c for c in df.columns if c.startswith(col2) or c.startswith(col2_)]

    l = []
    for k, df_k in df.groupby(['cond', 'fly', 'trial']):
        for c1 in cols_1:
            for c2 in cols_2:

                y1 = df_k.loc[:, c1].values
                y2 = df_k.loc[:, c2].values
                
                a = y1
                b = np.pad(y2, dt * f)
                c = np.correlate(a, b, mode='valid')

                # equivalent to matlab norm
                norm = np.sqrt( np.sum(a**2) * np.sum(b**2))
                with np.errstate(divide='ignore', invalid='ignore'): # ignore division by 0 warning: TODO better check array before calling corrcoef
                    c = c / norm

                # c = c / len(a)

                cond, fly, trial = k
                d = pd.DataFrame(data={
                    't': t,
                    'ccf': c,
                    'r': c1,
                    'b': c2,
                    'cond': cond,
                    'fly': fly,
                    'trial': trial,
                })

                l.append(d)

    df_c = pd.concat(l, ignore_index=True)
    
    return df_c

def select_event(df, col, f, dt):
    '''Select only data around behavior events

    Filter the data to contain only +- `dt` around behavior events in `col`.

    Return filtered dataframe, where other behavior columns are dropped.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with time series data
    col : str
        Column name around which to select data
    f : float
        Sample rate of data in Hz
    dt : float
        Time window in s around behavior event

    Returns
    -------
    df : pandas.DataFrame
        Filtered data
    '''

    dn = int(dt * f)

    l = []
    for k, d in df.groupby(['fly', 'trial']):
        ds = d.loc[:, col]
        i_ons = np.flatnonzero(np.diff(ds) == 1) + 1
        idx_ons = ds.iloc[i_ons].index
        for idx in idx_ons:
            try: 
                d_ = d.loc[idx - dn : idx + dn, :]
                l.append(d_)

            except ValueError:
                print(f'WARNING skipping fly {k[0]} trial {k[1]} event {idx}')

    df = pd.concat(l, ignore_index=True)

    # drop other behavior columns
    c_beh = [ c for c in df.columns if 'beh' in c]
    c_drop = [ c for c in c_beh if col.split('_')[-1] not in c]
    df = df.drop(columns=c_drop)

    return df