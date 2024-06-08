Analysis of calcium imaging data recorded in the fruit fly

# How to use this repo
For more information on the structure of this repo, 
see this [template repo](https://github.com/bidaye-lab/template_data_pipelines).

## Analysis pipelines
The script files `scripts/*.py` are workflows to analyse imaging data.

|script file|use case|
|---|---|
|`full_pipeline.py`|start-to-finish analysis pipeline merging raw data, ImageJ ROIs, ball velocities, and behavior|

## old scripts
These old scripts need to become part of the new code structure.
- `scripts/old_notebooks/ca_imaging_stim.ipynb`
- ~~`scripts/old_notebooks/ca_imaging.ipynb`~~
- `scripts/old_notebooks/focused_heatmaps.ipynb`
- `scripts/old_notebooks/hacky_aligned_behavior.ipynb`
- `scripts/old_notebooks/pretty_plot.ipynb`


## Installation
```
# create conda environment with necessary dependencies
conda create -n imaging_analysis -f environment.yml

# get source code
git clone https://github.com/bidaye-lab/imaging_analysis

# install code as local local python module
cd imaging_analysis
pip install -e .

# create notebooks
python src/create_notebooks.py
```

## Update
```bash
# pull from github
cd imaging_analysis
git pull origin main

# recreate notebooks
conda activate imaging_analysis
python src/create_notebooks.py
```
