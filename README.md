# Fiol et al (2026)
This repo contains the notebooks to reproduce the figures of the Fiol et al. 2026 paper on sea ice short term predictability.

It is divided into four sub-folders: a personal Python library containing homemade functions, one with the scripts to compute and analyse the Spatial Probability Score (SPS, Goessling and Jung, 2018) applied to the sea ice edge, one with the scripts to compute and analyse the Continuous Rank Probability Score (CRPS, Hersbach, 2000; Candille et al., 2015; Leroux et al., 2022), and one with the scripts related to Lagrangian trajectories. Other scripts not organised in the above directories are also provided: 

- illustrations4article.ipynb: script to create illustrations for the article with perturbated experiments

- TUV_masks.ipynb: script that create the masks of our domain of interest

- figures_unperturbed_experiment.ipynb: script to create the illustrations made with the unperturbed experiment

- Analyse_distribution_defo.ipynb and Probability_maps_defo.ipynb: the first script looks at the hourly deformation distributions for two different rheologies in order to define a threshold to use to compute probability maps with the second script

All the scripts are Jupyter notebooks coded in Python. The libraries necessary to be able to run all the scripts are:

- numpy

- cartopy.crs

- xarray

- scipy.stats

- pandas (only for IABP_csv_to_NetCDF.ipynb)

- cftime (only for IABP_csv_to_NetCDF.ipynb)

- glob (only for IABP_csv_to_NetCDF.ipynb)

- csv (for writing files in preprocess_seeding.ipynb)

	For graphical display:

- matplotlib.pyplot, matplotlib.patches, matplotlib.gridspec, matplotlib.colors

- cmocean

- mpl_toolkits.axes_grid1.inset_locator

- cartopy.feature, cartopy.mpl.gridliner

Many scripts use functions defined in the homemade Python library: Python_library. For the scripts to work, make sure to install Python_library library and that Python is able to locate it: put *export PYTHONPATH=<path_to_Python_library\>:${PYTHONPATH}* in your .bashrc, .profile or equivalent.

Some scripts also used functions from ensdam or sitrack.

By L. Fiol with contributions from  S. Leroux, P. Rampal and J.-M. Brankart
