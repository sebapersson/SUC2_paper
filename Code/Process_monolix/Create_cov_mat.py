#!/usr/bin/env python
import numpy as np
import pandas as pd
import itertools
import os

# This file will process the correlation parameters that are outputted from
# the parameter estimation in Monolix into a covariance matrix that can
# be used to simulate the population behaviour. A model is added by
# adding the relative path (from this file) to the monolix result
# directory (where the populationParameters.txt file is).
# Output:
#    A csv-file saved in the monolix result directory which contains
#      the covariance matrix.


# Function that will find all indices that match a certain character in a string
# Args:
#    string, the string
#    character, the character
# Returns:
#    All the matching indices 
def find_ch(string, character):
    return [i for i, ltr in enumerate(string) if ltr == character]


# Function that will given a path to populationParameters.txt will calculate the 
# covariance matrix and write it to a csv-file. 
# Args:
#    dir_pop_param, the directory of the population parameters 
#    path_save, the path to where the file should be saved 
# Returns:
#    void 
def create_cor_mat(dir_pop, n_var):
    # Read the population parameters
    path_pop_param = dir_pop + "populationParameters.txt"
    data_pop = pd.read_csv(path_pop_param)
    data_pop_values = data_pop[["parameter", "value"]]
    # Remove the variances sigma 
    data_pop_values = data_pop_values.iloc[0:-n_var]
    
    # Get the number of parameters, use the fact that no parameter start on omega
    n_parameters = 0
    for parameter in data_pop_values["parameter"]:
        if parameter[0:5] == "omega":
            n_parameters += 1
    
    # Allocate the correlation matrix
    cov_mat = np.zeros((n_parameters, n_parameters))
    
    # Associate each parameter with a value
    dict_param = {}
    i, j = 0, 0
    while len(dict_param.keys()) < n_parameters:
        # Find the first underline
        parameter = data_pop_values.iloc[i, 0]
        pos_underline = parameter.rfind("_")
        key_dict = parameter[0:pos_underline]
        # Check if fixed effect
        if len([k for k in data_pop_values["parameter"]
                if key_dict in k]) > 1:
            dict_param[key_dict] = j
            j += 1
        i += 1
    
    # Fill the covariance matrix
    for i in range(n_parameters, data_pop_values.shape[0]):
        # Get the parameter
        parameter = data_pop_values.iloc[i, 0]
        # The case for a diagonal entry
        if parameter[0:5] == "omega":
            which_param = parameter[parameter.find("_") + 1:]
            index = dict_param[which_param]
            cov_mat[index, index] = data_pop_values.iloc[i, 1] 
            
            # The case for correlation
        elif parameter[0:4] == "corr":
            # Get all underlines 
            index_underline = find_ch(parameter, "_")
            # Get the numbers for the matrix
            i1 = dict_param[parameter[index_underline[0]+1:index_underline[-1]]]
            i2 = dict_param[parameter[index_underline[-1]+1:]]
            # Calculate the covariance
            cov_mat[i1, i2] = data_pop_values.iloc[i, 1]*cov_mat[i1, i1]*cov_mat[i2, i2]
            cov_mat[i2, i1] = data_pop_values.iloc[i, 1]*cov_mat[i1, i1]*cov_mat[i2, i2]
            
    # Ensure variance along diagonal
    for i in range(n_parameters):
        cov_mat[i, i] *= cov_mat[i, i]
        
    parameter_names = dict_param.keys()
    data_to_save = pd.DataFrame(cov_mat, columns = parameter_names)
    path_save = dir_pop_param + "Cov_mat.csv"
    data_to_save.to_csv(path_save)
    
    return cov_mat, parameter_names


# Important that dir_pop_param ends with /
dir_pop_param = "../Monolix_code/Simple_feedback/Simple_feedback/"
cov_mat = create_cor_mat(dir_pop_param, n_var=1)

