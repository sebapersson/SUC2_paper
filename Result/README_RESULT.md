This directory contains, after running the Run_all-script, all the results (tables and figures) presented in the paper. Results belonging to a specific part of the analysis is in the same subdirectory (description of each subdirectory follows below). 

### STS 

This directory contains, for a given model, the diagnostic STS plots (parameter values can be found in intermediate directory). These plots are: 

* IPRED, individual observations versus observed data. 
* qq-plot of all parameters. 
* Simulation plots for each state (plots obtained from typically simulating 10,000 cells) 

### Monolix_processed 

This directory contains, for a given model, the diagnostic plots and parameters values obtained from estimating parameters using Monolix. These plots and files are: 

* IPRED, individual observations versus observed data. 
* IWRES, individually weighted residuals. 
* qq-plot of all parameters. 
* Simulation plots for each state (plots obtained from typically simulating 10,000 cells) 
* Individual plots of randomly drawn cells. 
* qq-plot of the residuals. 
* A csv-file and LaTeX-table-file (txt) with the estimated population parameters. 

The subdirectories for a model can also contain special plots, like the deletion bar plot for the Snf1_feedback model. 

### Experimental_data 

This directory contains plots for the experimental (invertase assay) data. To keep a consistent graphical profile in the paper, the experimental data was plotted in R. 

### Data set 

Contains a plot of the data set, and an illustrative plot used for Fig. 2B in the paper. 
