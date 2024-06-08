# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: everything
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
from pystackreg import StackReg

from scripts.batch_processing import motion_correction_based_on_ch2, extract_traces, merge_imaging_and_behavior, merge_sessions
from scripts.batch_analysis import generate_plots, spatial_corrmap, spatial_corrmap_pooled

# %% [markdown]
# # User settings

# %%
params = {
    "n_ch": 2,  # number of channels
    "n_z": 15,  # number of z slices
    "xy_smth": 4,  # smoothing in xy before registration (sigma 2D Gaussian)
    "perc": 10,  # settings for background subtraction: percentile
    "winsize": 50,  # and window size (winsize / f_ca = window in seconds)
    "f_ca": 2,  # sample rate imaging
    "f_ball": 50,  # sample rate ball velocity
    "f_beh": 200,  # sample rate behavior
    "ca_tau_on":  0.13, # rise and decay time [s] of Ca kernel in seconds, used for covolving ball and behavior
    "ca_tau_off": 0.63, # taken from https://elifesciences.org/articles/23496#s4
    "beh_keys": [  # keys in behavior file to consider
        "pushing",
        "hindLegRub",
        "frontLegRub",
        "headGrooming",
        "abdomenGrooming",
        "PER",
        "midLeft+hindLegRub",
        "midRight+hindLegRub",
        # 'ObL1+R1'
    ],
    # transformation used for registration https://pystackreg.readthedocs.io/en/latest/readme.html#usage
    "reg": StackReg.SCALED_ROTATION,
    # path to folder
    "parent_dir": r"Y:\Nino_2P_for_Salil\for_Nico\stop1_filteredData_2\stop1-GCaMP6f-tdTomato_VNC_smth4",
    "overwrite": True,  # force overwriting files
}

# ensure `parent_dir` is Path object
params["parent_dir"] = Path(params["parent_dir"])

# selection rule for tif files
params["p_tifs"] = [
    *params["parent_dir"].glob("**/trials_to_register/*/trial*_00???.tif")
]

# set and create output folders
params["p_out"] = Path(r'C:\temp\imaging_pipeline\test_run')
params["p_out_all"] = params["p_out"] / "all_data"
params["p_out_all"].mkdir(exist_ok=True, parents=True)

# print info
print("INFO: Saving files to {}".format(params["p_out"]))
print("INFO: Found {} files:".format(len(params["p_tifs"])))
for p in params["p_tifs"]:
    print(p)

# %% [markdown]
# # Processing
# The processing part contains
# - motion correction
# - fluorescence trace extraction
# - merging behavior and imaging data
#
# All output is stored in the `params["p_out"]` folder.
# Here, the folder structure of all `params["p_tifs"]` files is regenerated relative to `params["parent_dir"]`.
# These folders contain all intermediate files of the processing pipeline, 
# as well as a `plots` folder with analysis for that specific TIF file.
#
# The output folder further contains an `all_data` folder,
# which contains analysis pooled over all TIF files.

# %% [markdown]
# ## Motion correction
# The motion correction is done with [pystackreg](https://pystackreg.readthedocs.io/en/latest/readme.html#summary), which is a python port of the TurboReg/StackReg extension in ImageJ.
#
# The processing steps are
# - maximum projection along z
# - smoothing with 2D Gaussian along xy (`params['xy_smth']`)
# - alignment based on channel 2 (`params['reg']`)
#
# Following files are generated for each TIF file, split by channel:
# - uncorrected maxproj data (`*ch1.tif` and `ch2.tif`)
# - motion-corrected maxproj data (`*ch1reg.tif` and `ch2reg.tif`)
# - mean maxproj image (`*ch1mean.bpm` and `*ch2mean.bpm`)
# - side-by-side movies (quick quality control)
#     - maxproj channel 1 and 2 (`*ch1ch2.mp4`)
#     - uncorrected and motion-corrected maxproj channel 1 (`*ch1reg.mp4`)
#     - uncorrected and motion-corrected maxproj channel 2 (`*ch2reg.mp4`)
#
#
# For more details, look into the `scripts.batch_processing.motion_correction` function.
#
#

# %%
motion_correction_based_on_ch2(params)

# %% [markdown]
# ## Trace extraction
# This part loads the motion-corrected TIF files saved to disk
# and extracts the fluorescence traces for ROIs defined in the ImageJ Roi.zip file.
#
# Check `ch1mean_rois.bmp` if ROIs are parsed correctly.
#
# Following files are generated for each TIF file, split by channel.
# The name indicates for which data the trace was extracted:
# - `ch1raw/ch2raw`: raw data
# - `r12`: ratio of ch1 to ch2
# - `dch1/dch2`: baseline subtracted from raw data
# - `dr12`: baseline subtracted from ratio
# - `ch1/ch2`: background (last) ROI subtracted from raw data
#

# %%
extract_traces(params)

# %% [markdown]
# ## Merge imaging and behavior data
# Imaging data, behavior data, and ball velocity data is recorded at different sampling rates.
# Here, they are resampled to the sampling rate of the behavior data and merged into one pandas.DataFrame.
#
# The ball velocity file is expected to be called `trial16.mat` for the TIF file `trial16_00001.tif`.
#
# The behavior data file is expteced to be called `*-actions.mat`

# %%
merge_imaging_and_behavior(params)

# %% [markdown]
# As a final processing step,
# the data from all TIF files is merged into one pandas.DataFrame and 
# saved to `params["all_data"] / all_data_{method}.parquet` (see _Trace extraction_).

# %%
merge_sessions(params)

# %% [markdown]
# # Analysis and plots
#
# This part analyzes the fluorescence traces and behavior data and gerenerates plots.
#
# The following plots are generated by `generate_plots` in the respective `params['p_out_all'] / {method}` folder
# - `heatmap.svg`: Pearson correlation coefficients between all combinations of z-scored traces, Ca-kernel-convolved and z-scored ball velocities, and Ca-kernel-convolved behavior onsets
# - `heatmap_{beh}.svg`: same as `heatmap.svg`, but with data filtered around behavior event
# - `ccf.svg`: cross-correlation function between (i) z-scored traces and (ii) Ca-kernel-convolved and z-scored ball velocities and Ca-kernel-convolved behavior onsets
# - `ccf_indv.svg`: same as `ccf.svg`, but not averaged over flies/trials
# - `aligned_to_{beh}.svg`: average z-scored traces and smoothed (not Ca-kernel-convolved!) ball velocities
#
#

# %%
# analysis-specific parameters
params.update({
    "drop_rois": [7, 8, 9], # drop these ROIs for all plots in `generate_plots`
    "dt_beh": 5, # +- time window in s around behavior events in `heatmap_{beh}.svg`
    "dt_align": 5, # +- time window in s to plot in `aligned_to_{beh}.svg`
    "s_align": .25 # STD of Gaussian to smooth velocity in `aligned_to_{beh}.svg`
})

# %%
generate_plots(params)

# %% [markdown]
# ## spatial correlation maps
# Spatial correlation maps show the pixel-wise correlation between
# (i) motion-corrected raw data and
# (ii) Ca-kernel-convolved ball velocities and Ca-kernel-convolved behavior onsets.
# Because these are calculated per TIF file, they are saved in the `params['p_out'] / path_to_tif / 'corrmaps'` folder
#
# The following plots are generated by `spatial_corrmap` 
# (see `plot` folder for individual TIF files)
# - `{col}_1xy.svg`: spatial correlation maps for `col` separated by channel 1 and 2
# - `corrmap_{col}_{ch}.npy`: not a plot, but correlation data file to be reused for plotting
#
# Because the spatial correlation maps take a long time to calculate,
# we can specifically choose the TIF files for which to calculate the maps.

# %%
# plot spatial correlation maps for subset of TIF files
p_tifs = params["p_tifs"][:1]
spatial_corrmap(params, p_tifs)

# %% [markdown]
# Here we calculate the spatial correlation map pooled over multiple TIF files (e.g. all trials of one fly).
#
# This requires motion-correcting all trials together.
# Since the plots contain data from multiple TIF files, they
# are saved in `params['p_out'] / pool_name / 'correlation_maps'`

# %%
# choose TIF files to pool
p_tifs = [ *params['parent_dir'].glob('female4/trials_to_register/*/*0000?.tif') ]
# define unique name for folder
pool_name = 'female4_all_trials'
# run 
spatial_corrmap_pooled(params, p_tifs, pool_name)

# %%
