import numpy as np
from pathlib import Path

# plotting
import matplotlib.pylab as plt
from matplotlib.colors import CenteredNorm
import seaborn as sns

# video/image handling
from skvideo.io import vwrite
from PIL import Image

# local code
from src.processing import calculate_pearson, calculate_ccf, resample

###############
# visualization

def save_dual_movie(file, arr1, arr2, cmap='viridis', fps=30):
    '''Save arr1 and arr2 as side-by-side video.

    arr1 and arr2 are concatenated along the x-axis and therefore
    have to have the same height.

    Parameters
    ----------
    file : path-like
        Path to save file
    arr1 : numpy.ndarray
        Left 3D array of shape (n_frames, width, height)
    arr2 : numpy.ndarray
        Right 3D array of shape (n_frames, width, height)
    cmap : str, optional
        matplotlib colormap to use, by default 'viridis'
    fps : float, optional
        Frames rate for the movie, by default 30
    '''

    arr = np.concatenate([arr1, arr2], axis=-1) # along x dim

    cm = plt.get_cmap(cmap)
    arr = cm(arr.astype(int)) * 255

    print(f'INFO writing video to file {file}')
    vwrite(file, arr, inputdict={'-framerate': str(fps)})


def save_img(file, img):
    '''Save image to file

    Normalize data range to [0...255] and save as
    uint8 image file.

    Parameters
    ----------
    file : path-like
        Path to file on disk
    img : numpy.ndarray
        Image to be saved
    '''

    img = img.astype('float')
    img -= img.min()
    img /= img.max()
    img = (img * 255).astype('uint8')

    img = Image.fromarray(img)
    print(f'INFO saving normalized image to {file}')
    img.save(file)



def plot_aligned(df_al, ylims_roi=(None, None), ylims_ball=(None, None), path=''):
    "Plot averaged data after alignment to behavior event"

    fig, axarr = plt.subplots(nrows=2, figsize=(10, 10))

    ax = axarr[0]
    df = df_al.loc[ df_al.loc[:, 'match'].str.contains('roi') ]

    ax.axvline(0, ls=':', lw=.5, c='gray')
    sns.lineplot(ax=ax, data=df, x='t', y='y', hue='match', errorbar='se', palette='muted')
    ax.margins(x=0)
    ax.set_ylim(ylims_roi)


    ax.set_xlabel('time [s]')
    ax.set_ylabel('average zscored intensity')
    ax.set_title('''\
    Ca signal aligned to {}
    averaged over {} flies ({} trials)
    total of {} events ( {} skipped)
    '''.format(
        df.attrs['beh'],
        df.attrs['n_fly'],
        df.attrs['n_trial'],
        df.attrs['n_good'],
        df.attrs['n_bad'])
        )

    ax = axarr[1]
    df = df_al.loc[ df_al.loc[:, 'match'].str.contains('ball') ]

    ax.axvline(0, ls=':', lw=.5, c='gray')
    sns.lineplot(ax=ax, data=df, x='t', y='y', hue='match', errorbar='se', palette='muted')
    ax.margins(x=0)

    ax.set_xlabel('time [s]')
    ax.set_ylabel('average ball velocity')
    ax.set_ylim(ylims_ball)

    fig.tight_layout()
    if path:
        fig.savefig(path)
        plt.close(fig)


def plot_data(df, f, zroi=True, path=''):

    df = df.copy()

    c_roi = [ c for c in df.columns if c.startswith('roi_')]
    c_ball = [ c for c in df.columns if c.startswith('ball_')]
    c_beh = [ c for c in df.columns if c.startswith('beh_')]

    cols = c_roi + c_ball + c_beh

    # n = len(c_roi) + len(c_ball) + c_beh
    n = len(cols)
    fig, axarr = plt.subplots(nrows=n, figsize=(20, 2.5*n))
    x = df.index / f


    for ax, c in zip(axarr, cols):
        

        ax2 = ax.twinx()
        for i, c2 in enumerate(c_beh):
            y = df.loc[:, c2].values
            mask = y == 1
            offset = 0.1 * i
            y = y -  offset
            color = f'C{i}'
            ax2.scatter(x[mask], y[mask], color=color, marker='|', label=c2.replace('beh_', ''))
            ax2.axhline(1 - offset, color=color, lw=0.5, ls=':')
        ax2.set_ylim(-2, 1.1)
        ax2.set_yticklabels([])
        if c == c_roi[0]:
            ax2.legend()

        with np.errstate(divide='ignore', invalid='ignore'): # ignore division by 0 warning
            
            if c.startswith('ball') or c.startswith('beh_'):
                y = df.loc[:, c]
                y_c = df.loc[:, f'conv_{c}']
                y = y / y.max() * y_c.max()
                ax.plot(x, y, color='gray', ls='--', lw=1)
                ax.plot(x, y_c)

                if c.startswith('beh'):
                    y_c = df.loc[:, f'conv_{c}'.replace('beh', 'behi')]
                    y_c /= y_c.max()
                    ax.plot(x, y_c, color='C2', ls=':', lw=1)

                    y_c = df.loc[:, f'conv_{c}'.replace('beh', 'behf')]
                    y_c /= y_c.max()
                    ax.plot(x, y_c, color='C3', ls=':', lw=1)

            elif c.startswith('roi'):
                if zroi:
                    y = df.loc[:, f'z_{c}']
                else:
                    y = df.loc[:, c]
                ax.plot(x, y)
        

        ax.set_title(c)
        ax.set_xlabel('time [s]')
        ax.margins(x=0)

    fig.tight_layout()

    if path:
        fig.savefig(path)
        plt.close(fig)


def plot_corr_heatmap(df, beh, path=''):
    "Calculate and plot correlation heatmap for all ROIs, ball velocity, and behavior events"

    d = calculate_pearson(df, beh)

    fig, ax = plt.subplots(figsize=(10, 8))

    sns.heatmap(ax=ax, data=d, square=True, cmap='coolwarm', center=0)

    fig.tight_layout()
    if path:
        fig.savefig(path)
        plt.close(fig)

def plot_ccf(df, f, col1='z_roi_', col2='z_conv_ball_', col2_='conv_behi_', pool_fly=True, path=''):
    "Calculate and lot cross-correlation functions between (i) `col1` and (ii) `col2` + `col2_`"

    d = calculate_ccf(df, dt=2, f=f, col1=col1, col2=col2, col2_=col2_)
    d = d.dropna(axis=0)

    kw_args = {
        'x': 't',
        'y': 'ccf',
        'col': 'r', 
        'hue': 'b',
    }
    if not pool_fly:
        kw_args['hue'] = 'fly'
        kw_args['style'] = 'trial'
        kw_args['row'] = 'b'

    g = sns.relplot(data=d, kind='line', facet_kws={'sharey': False}, errorbar='se', **kw_args)
    g.set_axis_labels('time lag [s]', 'norm CCF')

    for ax in g.axes.flatten():
        ax.axvline(0, ls=':', lw=.5, c='gray', zorder=-1)
        ax.margins(x=0)
    
    fig = g.figure
    fig.tight_layout()
    if path:
        fig.savefig(path)
        plt.close(fig)

def plot_corrmap(arr1, arr2, df, b, f_ca, f_beh, cmap='bwr', path=''):
    "Calculate and plot pixel-wise correlation maps for two channels"
    
    y = df.loc[:, b].values

    fun = lambda x, y: np.corrcoef(resample(x, f_ca, f_beh), y)[0, 1]
    
    npy1 = Path(path).parent / f'corrmap_{b}_1.npy'
    npy2 = Path(path).parent / f'corrmap_{b}_2.npy'

    if npy1.exists() and npy2.exists():
        print(f'INFO found precalculated files, loading from disk')
        print(f'     {npy1}')
        print(f'     {npy2}')
        img1 = np.load(npy1)
        img2 = np.load(npy2)

    else:
        with np.errstate(divide='ignore', invalid='ignore'): # ignore division by 0 warning
            img1 = np.apply_along_axis(fun, 0, arr1, y)
            img2 = np.apply_along_axis(fun, 0, arr2, y)
        
        np.save(npy1, img1)
        np.save(npy2, img2)

    fig, axarr = plt.subplots(ncols=2, figsize=(14, 4))

    ax = axarr[0]
    im = ax.imshow(img1, cmap=cmap, norm=CenteredNorm())
    ax.set_title(f'{b} | ch 1')

    plt.colorbar(im, ax=ax)

    ax = axarr[1]
    ax.imshow(img2, cmap=cmap, norm=CenteredNorm())
    ax.set_title(f'{b} | ch 2')

    fig.tight_layout()
    if path:
        fig.savefig(path)
        plt.close(fig)