# About 

This directory contains intermediate files produced the analysis. The files located here are not results (e.g. figures and tables), but they are required to produce the files in the Result directory. 

After running the Run_all.sh script, this directory will contain two subdirectories, **Data_files** and **STS**. A **NLME** subdirectory is not a part of this directory, as I think it makes more sense to have the mlxtran-files (and consequently the Monolix output files) together with the Monolix model txt-files in the code-directory. 

### Data_files

This directory contains processed version of the SUC2-data files in the Data directory. More specifically, it contains: 
  
* **Data_monolix_SUC2.csv**, a filtered combination of the two experimentally identical data sets in the Data directory. As the data sets in the Data are not in readable format for Monolix, and some cells showed unrealistic behaviour, this filtered file is the one used for the NLME and the STS parameter estimations. 
* **Data_mean_SUC2.csv**, a per time point mean averaged version of the Data_monolix_SUC2.csv file. This file is used for generating start-guesses for the feedback cascade model (see method in paper for details). 

### STS

This directory contains the files (not plots and tables) produced by the STS parameter estimation. For a model, there is first a folder stating the name of the optimisation algorithm used. Within that, there are the following files: 

* **Cov_mat.csv**, covariance matrix of the logarithmed individual parameters. 
* **Mean_param.csv**, mean value of the logarithmed individual parameters.
* **Estimated_param.csv**, the individual parameter estimates.
* **Individual_fits.csv**, induvial fits for each cell.
* **Simulated_cells.csv**, simulated cells obtained by drawing random parameters from the estimated population parameters. 
