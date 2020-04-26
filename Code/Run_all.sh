 #!/bin/bash

# This script will reproduce all the modeling results presented
# in the paper.
# Firstly it will process the raw data into csv tidy data and plot the data set.
# Secondly, it will perform the STS optimisation and provide 
# instructions for the Monolix optimisation.
# Thirdly, given estimated parameters this script will call the processing
# functions, which will produce model diagnostic plots and tables.
# Fourthly, other plots, like the experimental invertase data, is plotted. 

# Usage: This script will always carry out all the analysis. It should
# be run from the Code-directory. It will exit if this is not the case.
# NOTE, the user should in the run_julia-function provide the path to
# the Julia executeble. 
# NOTE, the STS-optimization and simulation of cells is quite time-consuming.
# Hence, it might take up to eight hours to run this script. 

# -------------------------------------------------------------------------
# Start of functions 
# -------------------------------------------------------------------------

# Function that will check if a directory exists. If the directory doesn't
# exist the function will create it. Note that input order matters 
# Args:
#     $1 Directory name 
check_if_dir_exists ()
{
    # Create directory if it doesn't exist
    if [ ! -d $1 ]; then
	mkdir $1
    fi
}


# Function that will check if the results from the monolix
# parameter estimation exists
# Args:
#    $1 Model name (name of the model)
check_monolix_files_exist ()
{
    if [ ! -d $1 ]; then
	echo "Error: The directory for"
	echo $1
	echo "does not exist. Set up the parameters estimation using the"
	echo "instructions provided in the How_run.md file found in the"
	echo "model directory in the monolix folder"
	exit 1
    fi

    # Check that charts data exist
    cd $1
    if [ ! -d ChartsData ];then
	echo "Error: Charts data does not exist for model" 
	echo $1
	echo "Run the estimation (or only create plots if there is a result" 
	echo "using export charts). Else run the analysis, checking the" 
	echo "how_to_run.md file in the monolix folder for a model"
	exit 1
    fi
    cd ..
}


# Function that will check the return code from a file, if the
# file does not return zero, the program will exit
# Args:
#    $1 Return code from last file run
#    $2 Short message about the error 
check_return_code ()
{
    if [ $? != 0 ];then
	echo "Error code != 0, error in:"
	echo $2
	exit 1
    fi
}


# Function that will run the julia executble
# Args:
#    the command line arguments passed to julia
run_julia ()
{
    ~/julia-1.3.1-linux-x86_64/julia-1.3.1/bin/julia "$@"
}

# -------------------------------------------------------------------------
# End of functions, start of run all script 
# -------------------------------------------------------------------------

# Check that the script is run from Shell directory
currentDir=${PWD##*/}
if [ ! $currentDir == "Code" ]; then
    >&2 echo "The script most be run from Code directory"
    exit 1
fi


## Process and plot the data
echo "Processing the data set into monolix format and"
echo "plotting the data set"
cd Process_plot_data/
Rscript ./Process_into_monolix_format.R 2> /dev/null
check_return_code $? "Process into monolix format"
Rscript ./Produce_into_mean_data.R 2> /dev/null
check_return_code $? "Produce mean data"
Rscript Plot_data.R 2> /dev/null
check_return_code $? "Plotting the data set"
cd ..
echo ""

# Generate the start guesses 
cd ./Generate_start_guess/
echo "Generating start guesses"
run_julia "Generate_start_guess.jl" "simple_feedback" "3" "100"
check_return_code $? "Generating start guess"
run_julia "Generate_start_guess.jl" "snf1_feedback" "3" "100"
check_return_code $? "Generating start guess"

# Run the STS optimisation
cd ../STS
echo "Running STS"
run_julia "Perform_STS.jl" "simple_feedback" "3"
check_return_code $? "Doing STS"

# Process the STS result 
echo "Processing the STS result"
Rscript Process_STS.R 2> /dev/null
check_return_code $? "Processing STS results"
echo ""
cd ..

# Do the NLME parameter estimation
echo "Checking that NLME-estimation has been performed"
# Simple feedback model 
cd Monolix_code/Simple_feedback/
check_monolix_files_exist "Simple_feedback"
cd ../
# Snf1 model
cd Snf1_feedback
check_monolix_files_exist "Snf1_feedback"
cd ../..

echo ""
echo "Check OK! Moving to processing NLME-result"
cd Process_monolix/
./Create_cov_mat.py
check_return_code $? "Creating covariance matrix"
echo "Simulating cells"
run_julia "Simulate_cells.jl" "simple_feedback" "10000"
check_return_code $? "Simulating cells"
run_julia "Simulate_cells.jl" "snf1_feedback" "10000"
check_return_code $? "Simulating cells"
run_julia "Simulate_cells.jl" "snf1_feedback_d" "10000"
check_return_code $? "Simulating cells"

echo "Processing the monolix result to graphs"
Rscript Process_monolix_result.R 2> /dev/null
check_return_code $? "Processing monolix result"

# Add the processing of the graphs in the paper 
cd ../Other_plots
Rscript Sigmoid_function.R 2> /dev/null
check_return_code $? "Plotting sigmoid"
Rscript Experimental_data.R 2> /dev/null
check_return_code $? "Experimental data"

exit 0 
