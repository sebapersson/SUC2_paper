using CSV
using DifferentialEquations
using Distributions
using Printf
using DataFrames
using ProgressMeter


# This file contains the functions used for simulating cells given estimated
# population parameters in Monolix. The cells are simulated from a
# multivariate normal distribution, where the covariance matrix
# is created in the Create_cov_mat.py file.
# The file is divided into the following sections:
#   structs
#   functions


# ---------------------------------------------------------------------------
# Shared structs
# ---------------------------------------------------------------------------

# Struct that holds the parameter values for a model.
# Args:
#   rate_constants, the rate rate constants
#   initial_values, the initial values for a model
#   delays, the time-delays used for a model
mutable struct ParameterValues
    rate_constants
    initial_values
    delays
end


# Struct that holds the model information required for simulation cells
# Args:
#   states, a list of strings with the state names
#   initial_info, information of the initial value, float -> fixed, u#
#       means unknown and is matched from estimated parameters
#   n_states, the number of states in the model
#   model_name, the function name of the model to be estimated
struct ModelInfo
    states
    initial_info
    n_states
    model_name
end


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

# Function that will generate the inital and paramter values required for
# simulating cells. Note each column is a cell
# Args
#    n_cells_simulate, the number of cells to simulate
#    cov_mat, the covariance matrix
#    pop_param , the population parameters
# Returns
#    parameters, a matrix where each column is a new cell
function generate_parameters(n_cells_simulate, cov_mat, pop_param)

    # Generate the inital value and parameter matrix
    dist_mult_normal = MvNormal(cov_mat)
    rand_effects = rand(dist_mult_normal, n_cells_simulate)
    # For sake of speed apply exponential function once
    rand_effects = exp.(rand_effects)

    # Fixing the parameters
    n_param = length(pop_param)
    param_matrix = zeros(n_param, n_cells_simulate)
    for i in 1:n_param
        param_matrix[i, :] = pop_param[i] .* rand_effects[i, :]
    end

    return param_matrix
end


# Function that maps the initial values from the parameter matrix. Overall
# a vector with inital values is returned
# Args:
#   param_ind, a parameter vector for an individual (column param_matrix)
#   i_init, the indices for the initial values in param_matrix
#   model_info, the model info object
# Returns
#   init_vector, an init vector that can be added to ParameterInfo object
function map_initial_values(param_ind, i_init, model_info)

    init_unknown = param_ind[i_init]
    init_val = zeros(model_info.n_states)
    for i in 1:model_info.n_states
        # The case where the initial value is a scalar
        if isa(model_info.initial_info[i], Float64)
            init_val[i] = model_info.initial_info[i]
        # When the initial value is known
        else
            index_map = parse(Int64, model_info.initial_info[i][2:end])
            init_val[i] = init_unknown[index_map]
        end
    end

    return init_val
end


# Function that will solve the DDE-system for a given model.
# If it fails to solve the dde-system an error messege is
# produced, which will prompt the optimisation to exit.
# Args:
#   model, the ode-model
#   t_span, the time span for the model
#   model_param, the model parameter struct
# Returns
#   dde_solution object if success
#   "exit" if failed to solve the ode-system
function solve_dde_system(model, t_span, param_model, tau)

    # Defining the delay funciton
    h(p, t) = param_model.initial_values
    # Check whetever or not a delay should be estimated
    if param_model.delays == empty
        problem_param = vcat(param_model.rate_constants, tau,
            param_model.initial_values)
    else
        problem_param = vcat(param_model.rate_constants, tau,
            param_model.delays, param_model.initial_values)
    end
    success_symbol = Symbol("Success")
    dde_problem = DDEProblem(model, param_model.initial_values, h,
        t_span, problem_param)

    # Solve the problem
    alg = MethodOfSteps(Tsit5())
    dde_solution = solve(dde_problem, alg, verbose = false)

    if dde_solution.retcode != success_symbol
        return "exit"
    else
        return dde_solution
    end
end


# Function that, given population parameters for a certain model,
# simulates a user provided number of cells. The result is then
# written to file and saved in the Monolix result folder, from where it
# retreives the population parameters.
# Args:
#   result_dir, the result directory for monolix
#   n_cells_simulate, the number of cells to simulate
#   tau, the time-delay for SUC2
#   model_info, model info object about current model
#   cov_mat_name, name of cov-mat if testing different structures
# Returns:
#   0 upon succes
function simulate_cells_nlme(result_dir, n_cells_simulate,
    tau, model_info; tag="", fixed=[], time_span=(0.0, 500.0))


    path_cov = result_dir * "Cov_mat.csv"
    path_pop_effects = result_dir * "populationParameters.txt"
    cov_mat_raw = CSV.read(path_cov)
    cov_mat = convert(Matrix, cov_mat_raw[1:end, 2:end])
    df = CSV.read(path_pop_effects)
    df = df[occursin.(r"s*_pop", df[:, 1]), :]
    not_fixed = []
    for i in 1:length(df[:, 1])
        if df[i, 1] in fixed
            continue
        end
        not_fixed = append!(not_fixed, i)
    end
    df = df[not_fixed, :]
    pop_effects = convert(Array, df[:, 2])

    # Index for matching population parameters to a ParameterInfo objects
    i_delays = findall(x-> x != 0, occursin.(r"^tau", df[:, 1]))
    i_rates = findall(x-> x != 0, occursin.(r"^k", df[:, 1]))
    i_init = .!occursin.(r"^tau", df[:, 1]) + .!occursin.(r"^k", df[:, 1]) .- 1
    i_init = findall(x->x != 0, i_init)

    # Draw the random effects for a cells
    param_matrix = generate_parameters(n_cells_simulate, cov_mat, pop_effects)

    # Details for solving ODE:s
    time_span_int = range(time_span[1], time_span[2], length = 200)

    # Matrix to hold all the simulation results, last column is an id
    mat_sim_result = zeros(n_cells_simulate * 200, model_info.n_states+2)
    mat_sim_result[:, (end-1)] = repeat(time_span_int, inner=1, outer = n_cells_simulate)
    mat_sim_result[:, end] = repeat(1:n_cells_simulate, inner = 200)

    # Simulate cells
    @showprogress 1 "Simulating cells ... " for i in 1:n_cells_simulate
        # Operations for individual cells
        param_test = param_matrix[:, i]
        init_val = map_initial_values(param_test, i_init, model_info)
        param_ind = ParameterValues(param_test[i_rates], init_val, param_test[i_delays])
        dde_sol = solve_dde_system(model_info.model_name, time_span, param_ind, tau)
        if dde_sol == "exit"
            continue
        end
        for j in 1:200
            i_mat = (i-1) * 200 + j
            mat_sim_result[i_mat, 1:(end-2)] = dde_sol(time_span_int[j])
        end
    end

    # Write the result to disk
    result_to_save = convert(DataFrame, mat_sim_result)
    names_save = vcat(model_info.states, "t", "id")
    rename!(result_to_save, Symbol.(names_save))
    if time_span[2] != 500
        path_save = result_dir * "Simulated_cells_e" *
        string(time_span[2]) * ".csv"
    else
        path_save = result_dir * "Simulated_cells" * tag * ".csv"
    end
    CSV.write(path_save, result_to_save)

    return 0
end
