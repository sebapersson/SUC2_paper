# About 

This directory contains all the code required to produce the result from the data in the data directory. Files performing a specific part of the analysis are located together in a subdirectory (description of each subdirectory below). For example, the run and process STS scripts are found in the STS subdirectory. Each script contains a description at the top. 

The Run_all.sh script calls, in the correct order, all scripts. Consequently, after running this script the results presented in the paper can be found in the result directory. Bellows follows a description of the Run_all-script. 

### Run_all.sh 

The run all script performs the following steps: 

1. Processes the data in the data-directory into a Monolix readable data set, and a mean-value data set, using the scripts in the **Process_plot_data** directory. Note, the data set is also plotted in this step. 
2.  Generates start guesses by a multiple shooting method using the mean value data. This is achieved using the script in the **Generate_start_guess* 
3.  Performs the STS optimisation and processes the result using scripts in the **STS** directory. 
4.  Performs the NLME parameter estimation using the models in the **Monolix_code** directory. Note the script will not call Monolix. Rather, if Monolix has not been run the script will instruct the reader on how to run the Monolix estimation. 
5. Processes the Monolix output using scripts in the **Monolix_code** directory. 
6. Plots the experimental data, and some supplementary material, using scripts in the **Other_plots** directory. 

It should be noted that the STS, but also simulating cells for NLME, requires several functions. These functions can be found in the **Julia_functions** directory. 
