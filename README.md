# ERalpha_Steiner2025
Scripts relating to analysis of SMT and differential transcription data in Hunter and Steiner et al., 2025
Truncated source data is available in dat/ for trial runs of scripts
Full source data available at: 10.5281/zenodo.18635439

Requirements:
R 4.4.1 <br>
<br>
Python 3.12.8 <br>
PyMC 5.25.1 <br>
pytensor 2.23.0 <br>
jupyter 1.1.1 <br>
matplotlib 3.10.0 <br>
scikit-learn 1.6.1 <br>
scipy 1.15.0 <br>

<br>
For a full list of installed packages, see package.reqs for complete R session info and python environment details
<br>
model_fit_residence_time.ipynb - Jupyter notebook for fitting residence time data using PyMC. <br>
dat - input data for scripts <br>
scripts - R scripts which produce the indicated figure <br>


## Installation
For both R and Python installs, installation should be rapid (< 10 minutes total). 
### Jupyter Notebook
The packages for running the Jupyter notebooks are best installed using conda: <br>

```
conda env create -f environment.yml
conda activate pymc_env_hunter_steiner_2026
python -m ipykernel install --user --name pymc_env_hunter_steiner_2026
jupyter notebook
```
The jupyter notebook can be run end-to-end to remake all figures present in the manuscript. User data can be substituted in, provided it follows the same format as `dat/df_tracks.csv`. Each field-of-view will take approximately 1 minute to fit on a personal computer.

<br>
Navigate to the notebook path in the Jupyter hub window. Change the working path to the repository install location <br>

### R Scripts
The libraries for running R scripts are best installed by restoring the R env snapshot:

```
R
renv::restore(lockfile="renv.lock")
```
Each R script should run rapidly, generating each specified figure in under a minute.

