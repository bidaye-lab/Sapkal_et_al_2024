import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from src.file_handling import fname, load_tiff, write_tiff, save_params_json
from src.visualization import (
    plot_corr_heatmap,
    plot_ccf,
    plot_aligned,
    plot_corrmap,
    save_dual_movie,
    save_img,
)
from src.processing import align2events, select_event
from src.registration import align, get_tmats


def generate_plots(params):
    "Generate plots for data files in `params['p_out_all']`"

    # load parameters
    p_out_all = params["p_out_all"]
    f_beh = params["f_beh"]
    drop_rois = params["drop_rois"]
    dt_beh = params["dt_beh"]
    dt_align = params["dt_align"]
    s_align = params["s_align"]

    # save params
    save_params_json(params, p_out_all / "params_plots.json")

    # merge all trials and flies
    for p_df in p_out_all.glob("all_data_*.parquet"):
        # get method name from file name
        method = p_df.stem.split("_")[-1]

        # create plot folder for each method
        p_plot = p_out_all / method
        p_plot.mkdir(exist_ok=True)

        # read data from disk
        df = pd.read_parquet(p_df)
        df = df.fillna(0)  # TODO workaround because of missing behavior

        # remove ROIs
        df = df.drop(columns=[f"z_roi_{i}" for i in drop_rois])

        # pearson correlation heatmap (selection of columns: see src.processing.calculate_pearson)
        plot_corr_heatmap(df, beh="behi", path=p_plot / "heatmap.svg")

        # pearson heatmaps around behavior events
        # loop through all behavior
        cols = [c for c in df.columns if c.startswith("beh_")]
        for col in cols:
            # select df around behavoir
            d = select_event(df, col, f_beh, dt_beh)

            # generate plot
            plot_corr_heatmap(d, beh="behi", path=p_plot / f"heatmap_{col}.svg")

        # cross-correlation functions behavior/ball and ROIs (averaged)
        plot_ccf(df, f=f_beh, pool_fly=True, path=p_plot / "ccf.svg")

        # same (not averaged)
        plot_ccf(df, f=f_beh, pool_fly=False, path=p_plot / "ccf_indv.svg")

        # plot aligned data
        # smooth velocity
        df_ = df.copy()
        df_.loc[:, ["ball_x", "ball_y", "ball_z"]] = gaussian_filter(
            df_.loc[:, ["ball_x", "ball_y", "ball_z"]].values, (s_align * f_beh, 0)
        )

        # cycle through all behavoirs
        cols = [c for c in df_.columns if c.startswith("beh_")]
        for col in cols:
            # align to behavior
            df_al = align2events(df_, col, f_beh, dt_align)
            plot_aligned(df_al, path=p_plot / f"aligned_to_{col}.svg")


def spatial_corrmap(params, p_tifs):
    "Calculate and plot spatial correlation maps for ball velocities and behavior events"
    
    # load parameters
    p_out = params["p_out"]
    p_parent = params["parent_dir"]
    f_ca, f_beh = params["f_ca"], params["f_beh"]

    # define tif file for session of interest
    for p_tif in p_tifs:
        # define and create output folder
        p_plot = fname(p_tif, "", old_root=p_parent, new_root=p_out).parent / "corrmaps"
        p_plot.mkdir(exist_ok=True)
        print(p_plot)

        # load registered tif files
        ch1 = load_tiff(fname(p_tif, "ch1reg.tif", old_root=p_parent, new_root=p_out))
        ch2 = load_tiff(fname(p_tif, "ch2reg.tif", old_root=p_parent, new_root=p_out))

        # load preprocessed behavior data
        df = pd.read_parquet(
            fname(p_tif, f"data_ch1.parquet", old_root=p_parent, new_root=p_out)
        )  # either of the data files work, since we are only interested in the behavior and ball columns

        # loop through all conv behi and conv ball columns
        cols = [
            c
            for c in df.columns
            if c.startswith("conv_behi") or c.startswith("conv_ball")
        ]
        for col in cols:
            plot_corrmap(
                ch1,
                ch2,
                df,
                col,
                f_ca=f_ca,
                f_beh=f_beh,
                path=p_plot / f"{col}_1xy.svg",
            )


def spatial_corrmap_pooled(params, p_tifs, pool_name):
    "Calculate and plot pooled spatial correlation maps"

    # load parameters
    p_out = params["p_out"]
    p_parent = params["parent_dir"]
    reg = params["reg"]
    f_ca, f_beh = params["f_ca"], params["f_beh"]

    p_pool = p_out / pool_name
    p_plot = p_pool / "corrmaps"	
    p_plot.mkdir(exist_ok=True, parents=True)

    try:
        ch1 = load_tiff(p_pool / "ch1.tif")
        ch2 = load_tiff(p_pool / "ch2.tif")

    except FileNotFoundError:
        print('INFO Combining TIF files')

        # combine tifs
        l1, l2 = [], []
        for p_tif in p_tifs:
            # load uncorrected tif files
            ch1 = load_tiff(fname(p_tif, "ch1.tif", old_root=p_parent, new_root=p_out))
            ch2 = load_tiff(fname(p_tif, "ch2.tif", old_root=p_parent, new_root=p_out))
            l1.append(ch1)
            l2.append(ch2)

        # combine to single array
        ch1 = np.concatenate(l1, axis=0)
        ch2 = np.concatenate(l2, axis=0)

        # register
        print('INFO Motion correction')
        tmats = get_tmats(ch2, reg)
        ch1_a = align(ch1, tmats, reg)
        ch2_a = align(ch2, tmats, reg)

        # mean image
        ch1_am = np.mean(ch1_a, axis=0)
        ch2_am = np.mean(ch2_a, axis=0)

        # save to `pool_name` folder
        write_tiff(p_pool / "ch1.tif", ch1_a.astype("int16"))
        write_tiff(p_pool / "ch2.tif", ch2_a.astype("int16"))

        write_tiff(p_pool / "ch1reg.tif", ch1_a.astype("int16"))
        write_tiff(p_pool / "ch2reg.tif", ch2_a.astype("int16"))

        save_img(p_pool / "ch1mean.bmp", ch1_am)
        save_img(p_pool / "ch2mean.bmp", ch2_am)

        save_dual_movie(p_pool / "ch1ch2.mp4", ch1, ch2)
        save_dual_movie(p_pool / "ch1reg.mp4", ch1, ch1_a)
        save_dual_movie(p_pool / "ch2reg.mp4", ch2, ch2_a)

    # load preprocessed behavior data
    l = []
    for p_tif in p_tifs:
        p_prq = fname(p_tif, f"data_ch1.parquet", old_root=p_parent, new_root=p_out)
        df = pd.read_parquet(p_prq)  # either of the data works
        print(f'INFO Loading {p_prq}')
        l.append(df)
    df = pd.concat(l, ignore_index=True)

    # loop through all conv behi and conv ball columns
    cols = [
        c for c in df.columns if c.startswith("conv_behi") or c.startswith("conv_ball")
    ]
    for col in cols:
        print(f'INFO Plotting {col}')
        plot_corrmap(
            ch1, ch2, df, col, f_ca=f_ca, f_beh=f_beh, path=p_plot / f"{col}_1xy.svg"
        )
