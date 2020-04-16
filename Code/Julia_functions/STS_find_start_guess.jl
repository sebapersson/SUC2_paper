using NLopt
using Printf
using DifferentialEquations
using DataFrames
using CSV
using ProgressMeter
using Statistics
using Random
using Distributions
using LinearAlgebra

# This file contains the functions used for the STS optmisation and the
# generate starting value procedure (where optimisation is run to the
# mean value). As both these procedures basically optimize a likelihood,
# they share many different functions, hence the functions are lumped into
# one file.
# The file is divided into the following sections:
#   shared structs
#   shared function
#   functions for STS
#   functions for starting values

# ---------------------------------------------------------------------------
# Shared structs
# ---------------------------------------------------------------------------
# Structs to hold parameter values for a model.
# Args:
#   rate_constants, the model rate constants
#   initial_values, the inital values for the model
#   delays, the values for potential time-delays (can be empty)
#   noise, the noise value (an additivate noise is assumed )
mutable struct ParameterValues
    rate_constants
    initial_values
    delays
    noise
end


# Struct that contains information of the parameters to estimate
# Args:
#   n_... number of rate constants etc
#   map_init_opt_vector, a vector that maps the unknown initial values
#       to the initial values in the parameter-values object from the
#       opt-vector used in the optmisation.
#   tau, the time delay for maturation (default 32 minutes in all functions)
#   time_span, the time_span for solving the DDE-system (default 0-500 min)
struct ParameterInfo
    n_rate_constants
    n_initial_values
    n_delays
    n_unknown_init
    map_init_opt_vector
    tau
    time_span
end


# Struct that contains information about the states in the model
#   states, a vector with state-names
#   init_status, vector that conveys information about the initial value.
#       A value can be u# (unknown will be estimated, here # = matches the
#       the unknown inital to the start-guess object). If a value is a float,
#       e.g 4.0, the initial value is fixed to that value
#   i_SUC2, (mostly for speed) the index of SUC2 (output) in the state-vector
struct StateInfo
    states
    init_status
    i_SUC2
end

# Struct that contains the start-guess for a model.
# Args:
#   rate_const--noise, start guesses for the different model parameters. Note
#       that delays can be empty delays = []
#   names_opt_vec, names (in order of struct) for the unknown parameters.
struct StartGuess
    rate_const
    delays
    init_values
    noise
    names_opt_vec
end


# ---------------------------------------------------------------------------
# Shared functions
# ---------------------------------------------------------------------------


# Function that will map the opt_vector (vector to optimise) to the
# model parameter object.
# Args:
#   param_model, the model parameters
#   opt_vector, the current value of the opt-vector (is updated)
#   param_info, the parameter info
# Returns
#   param_model, an updated model parameter vector
function map_opt_vec_to_model(param_model, opt_vector, param_info)

    n_delays = param_info.n_delays
    n_rates = param_info.n_rate_constants
    n_unknown_init = param_info.n_unknown_init
    map_init_opt_vector = param_info.map_init_opt_vector

    # Map the rate constants
    param_model.rate_constants = opt_vector[1:n_rates]
    param_model.noise = opt_vector[end]

    # Mapping unknown initial initial_values and potential delays
    if param_model.delays == empty
        param_model.initial_values[map_init_opt_vector] =
            opt_vector[n_rates+1:end-1]
    else
        param_model.initial_values[map_init_opt_vector] =
            opt_vector[(n_rates+n_delays+1):end-1]
        param_model.delays = opt_vector[n_rates+1:n_rates+n_delays]
    end

    return param_model
end


# Function that will solve the dde-system for a given model.
# If it fails to solve the ode-system an error messege is
# produced, which will prompt the optimisation to exit.
# Args:
#   model, the ode-model
#   t_span, the time span for the model
#   model_param, the model parameter struct
# Returns
#   ode_solution object if success
#   "exit" if failed to solve the ode-system
function solve_dde_system(model, param_model, param_info)

    # Defining the delay funciton
    h(p, t) = param_model.initial_values
    # Check whetever or not a delay should be estimated
    tau = param_info.tau
    if param_model.delays == empty
        problem_param = vcat(param_model.rate_constants, tau,
            param_model.initial_values)
    else
        problem_param = vcat(param_model.rate_constants, tau,
            param_model.delays, param_model.initial_values)
    end
    success_symbol = Symbol("Success")
    dde_problem = DDEProblem(model, param_model.initial_values, h,
        param_info.time_span, problem_param)

    # Solve the problem
    alg = MethodOfSteps(Tsit5())
    dde_solution = solve(dde_problem, alg, verbose = false)

    if dde_solution.retcode != success_symbol
        @printf("Failed solving dde system\n")
        return "exit"
    else
        return dde_solution
    end
end


# Function that will calculate the sum of square cost value
# for a model.
# Args:
#   ode_sol, the ode solution
#   SUC2_data, (n x 2) array, first column is time, second is data values
#   state_info, struct containing information of all states
# Returns:
#   cost, the cost function value for the current iteration
function calc_cost(dde_sol, SUC2_data, state_info, param_model)

    l_SUC2 = length(SUC2_data[:, 1])
    simulated_SUC2 = zeros(l_SUC2)
    for i in 1:l_SUC2
        t = SUC2_data[i, 1]
        simulated_SUC2[i] = dde_sol(t)[state_info.i_SUC2]
    end
    cost = sum((SUC2_data[:, 2] .- simulated_SUC2).^2)
    cost /= param_model.noise^2
    cost += l_SUC2 * log(param_model.noise^2)

    return cost
end


function target_function(opt_vector, grad, param_info, state_info,
    SUC2_data, param_model, model)

    if length(grad) > 0
        @printf("Cannot calculate gradient\n")
    end

    param_model = map_opt_vec_to_model(param_model, opt_vector, param_info)
    dde_sol = solve_dde_system(model, param_model, param_info)
    if dde_sol == "exit"
        return 1e8
    end
    cost = calc_cost(dde_sol, SUC2_data, state_info, param_model)

    return cost
end


# Function that given a string of the states will create a StateInfo struct
# Args:
#   state_info, list with all the states
#   observed, which state is observed
#   init_status, a string vector about the status of the initial values
# Returns:
#   states, object with information of the states
function produce_state_info(states, observed, init_status)

    # Define the indices for SUC2
    if "SUC2" in observed
        i_SUC2 = findall(x->x=="SUC2", states)[1]
    else
        i_SUC2 = empty
    end

    state_info = StateInfo(states, init_status,i_SUC2)

    return state_info
end


# Function that will create the parameter info object, param_values objects
# and the opt_vector used in the optmisation. This function sets upp
# all the parameter structs and array used for a model when
# running the optimsation. For sake of speed this function
# collects many operations that other functions could perform
# (i.e n_unkown_init)
# Args:
#   start_guess, a start-guess object with the start-guess for each cell
#   state_info, the state-info object for the model
#   tau (opt), the time-delay for SUC2-expression
#   time_span (opy), time span for solving the dde-system
# Returns:
#   opt_vec, the vector used in the optmisation
#   param_info, the parameter info object for the model
#   param_values, initial parameter values for a model
#   error, upon bad input
function create_param_obj(start_guess, state_info; tau=32.0, time_span=(0.0, 500.0))

    init_status = state_info.init_status

    # Get parameter information from the the start-guess and state ifno
    n_init_values = length(state_info.states)
    n_delays = length(start_guess.delays)
    n_rates = length(start_guess.rate_const)
    n_unkown_init = length(start_guess.init_values)

    # Sanity check
    n_unknown_init_provided = length(findall(x->x[1]=='u', init_status))
    if n_unknown_init_provided != n_unkown_init
        @printf("Error: The init status array does not match the \n")
        @printf("the number of unknowns provided by the start-guess\n")
        return "Error"
    end

    # Create the parameter values object, and set-up mapping from
    # initial values to opt-vector
    map_init_opt_vector = zeros(Int64, n_unkown_init)
    init_val = zeros(n_init_values)
    for (init_info, i) = zip(state_info.init_status, 1:n_init_values)
        if isa(init_info, Float64) || isa(init_info, Int64)
            init_val[i] = convert(Float64, init_info)
        elseif init_info[1] == 'u'
            if length(init_info) == 1
                @printf("Error: Wrong format on provided initial values\n")
                @printf("Should be u#, where # = digit that matches the\n")
                @printf("unknown initial to the start-guess object \n")
                return "Error"
            end
            digit = parse(Int64, init_info[2:end])
            if digit > n_unkown_init
                @printf("Error: Bad initial unknown argument for u#\n")
                @printf("# most lie between 1 - n_unknown init\n")
                @printf("n_unknown_init = %d, # = %d\n", n_unkown_init, digit)
                return "Error"
            end
            init_val[i] = start_guess.init_values[digit]
            map_init_opt_vector[digit] = convert(Int64, i)
        else
            @printf("Error: Unknown format on initial status in state-info\n")
            @printf("Should be a float/int, or u#, where # = integer\n")
            return "Error"
        end
    end

    # Create the objects, note that the last entry is the noise
    param_info = ParameterInfo(n_rates, n_init_values, n_delays, n_unkown_init,
        map_init_opt_vector, tau, time_span)
    param_values = ParameterValues(start_guess.rate_const, init_val,
        start_guess.delays, start_guess.noise)
    if n_delays == 0
        opt_vec = vcat(start_guess.rate_const, start_guess.init_values,
            start_guess.noise)
    else
        opt_vec = vcat(start_guess.rate_const, start_guess.delays,
            start_guess.init_values, start_guess.noise)
    end

    return param_values, param_info, opt_vec
end


# Function that will deepcopy a struct of type
# parameter-values. Note that deepcopy means that an entirely
# new object is created.
# Args:
#   param_val, a parameter value object
# Returns:
#   A deep copy of the parameter value object
function deepcopy_param_val(param_val)
    param_val_copy = ParameterValues(
        deepcopy(param_val.rate_constants),
        deepcopy(param_val.initial_values),
        deepcopy(param_val.delays),
        deepcopy(param_val.noise))

    return param_val_copy
end


# ---------------------------------------------------------------------------
# STS functions
# ---------------------------------------------------------------------------


# Function that will solve the dde and return the individual
# predictions in certain time-pints.
# Args:
#   dde_sol, solution of the dde-system
#   time-points, time-points where to interpolate solution
#   state_info, state info object for the current model
# Returns:
#   Array with simulated data
#   Error, if the input time-points are bigger than the simulated
function calc_ind_pred(dde_sol, time_points, state_info)

    if dde_sol.t[end] < time_points[end]
        @printf("Error: The provided time-points are bigger than\n")
        @printf("the simulated time-points.\n")
        return "Error"
    end

    l_data = length(time_points)
    simulated_obs = zeros(l_data)
    for i in 1:l_data
        t = time_points[i]
        simulated_obs[i] = dde_sol(t)[state_info.i_SUC2]
    end

    return simulated_obs
end


# Function that performs the STS optmisation for a user provided model, with
# a user provided algorithm (BOBYQA). The funciton will save the individual
# fits to disk, and the optmised parameters to a directory in the
# intermediate directory.
# Args:
#   dir_save, the directory to save the result in
#   state_info, the state-info object for the model
#   start_guess, the start-guess object for the model
#   model, the model to optmise
#   alg_opt, the opmisation algorithm
#       NELDERMEAD
#       SBPLX
#       BOBYQA (default)
# Returns:
#   opt_mat, a matrix where each row is the parameters for a cell
#   param_start, a parameter-val object that has all the constant
#       values correct.
function estimate_ind_param(dir_save, state_info, start_guess, model, alg_opt)

    # Extract time and observed values for an individual into array
    path_data = "../../Intermediate/Data_files/Data_monolix_SUC2.csv"
    data = CSV.read(path_data)
    id_list = unique(data[:, :id])
    list_check = 1:length(id_list)

    # Set upp the objects used for optmisation
    param_start, param_info, opt_vec = create_param_obj(start_guess, state_info)
    n_param = length(opt_vec)

    # Data structs for storing the result
    diagnostics = zeros(nrow(data), 3)
    diagnostics[:, 1:2] .= data[:, [:time, :id]]
    opt_mat = zeros(length(list_check), length(opt_vec) + 1)
    opt_mat[:, 1] = id_list[list_check]

    # Order out the individual data
    @showprogress 1 "Doing STS ..." for i in list_check
        id_ind = id_list[i]
        data_ind = filter(row->row[:id] == id_ind, data)[:, [:time, :observation]]
        indices = findall(x -> x == id_ind, data[:, :id])
        i_min, i_max = indices[1], indices[end]
        data_ind = convert(Array{Float64}, data_ind)

        # Solve the optmisation problem
        opt = Opt(alg_opt, n_param)
        opt.lower_bounds = 0.0
        opt.upper_bounds = Inf
        opt.maxeval = 10000

        min_objective!(opt, (opt_vector, grad) -> target_function(opt_vector, grad,
            param_info, state_info, data_ind, param_start, model))
        (minf, minx, ret) = optimize(opt, opt_vec)
        param_final = map_opt_vec_to_model(param_start, minx, param_info)

        # Map the parameter values
        opt_mat[i, 2:end] = minx

        # Calculate the individual predictions (later for diagnostics)
        dde_sol = solve_dde_system(model, param_final, param_info)
        ind_pred = calc_ind_pred(dde_sol, data_ind[:, 1], state_info)
        diagnostics[i_min:i_max, 3] .= ind_pred
    end

    # Save the parameters and diagnostics
    data_diag = convert(DataFrame, diagnostics)
    rename!(data_diag, [:id, :time, :ind_fit])
    data_param = convert(DataFrame, opt_mat)
    rename!(data_param, Symbol.(vcat("id", start_guess.names_opt_vec)))
    CSV.write(dir_save * "Individual_fits.csv", data_diag)
    CSV.write(dir_save * "Estimated_param.csv", data_param)

    return opt_mat, param_start, param_info
end


# Function that given an opt-matrix (from STS optimisation) calculates
# a ML-estimate of the covariance matrix and population parameters
# under the assumption that the parameters follow a log-normal distribution.
# The result is saved as csv-files in the provided dir-save
# Args:
#   dir_save, the directory to save the result in
#   start_guess, the start-guess object for the model
#   opt_mat, the opt-matrix obtained from the optmisation
# Returns:
#   pop_param, the estimated population parameters (excluding noise)
#   cov_mat, the estimated covariance matrix
function ml_est_cov_and_mean(dir_save, start_guess, opt_mat)

    # Separate noise parameter from the other parameters
    opt_mat = opt_mat[:, 2:end]
    col_noise = findall(x -> x[1] == 'a', start_guess.names_opt_vec)
    n_cells = size(opt_mat)[1]
    param_mat = log.(opt_mat[:, Not(col_noise)] .+ 1e-7)
    noise_vec = opt_mat[:, col_noise]

    # Calculate the population values, use likelihood for cov-mat
    pop_param = mean(param_mat, dims = 1)
    noise_mean = mean(noise_vec)
    param_mean = hcat(pop_param, noise_mean)
    cov_mat = Statistics.cov(param_mat, corrected=false)

    # Save the result as csv-files
    names_not_noise = filter(x -> x[1] != 'a', start_guess.names_opt_vec)
    data_mean = convert(DataFrame, param_mean)
    rename!(data_mean, Symbol.(start_guess.names_opt_vec))
    data_cov = convert(DataFrame, cov_mat)
    rename!(data_cov, Symbol.(names_not_noise))
    CSV.write(dir_save * "Mean_param.csv", data_mean)
    CSV.write(dir_save * "Cov_mat.csv", data_cov)

    return pop_param, cov_mat
end


# Function that will simulate cells given estimated covariance matrix
# and population parameters. The function is a part of the perform
# STS pipeline as it relies on structs used for the estimation. If the
# covariance matrix is not positive definite an error is returned
# Args:
#   pop_param, the estimated (logged) population parameters
#   cov_mat, the estimated covariance matrix
#   param_start, the parameter start object
#   param_info, the parameter info object for the model
#   model, the model used
#   state_info, the state info object of the model
#   dir_save, the directory to save the result in
#   n_cells_simulate, optional with default 10 000
# Returns
#   Error if covaraince matrix not positive definite
function simulate_cells(pop_param, cov_mat, param_start, param_info, model,
    state_info, dir_save; n_cells_simulate=10000)

    if !isposdef(cov_mat)
        @printf("Error: Estimated covariance matrix is not positive definite\n")
        return "Error"
    end

    # Simulate cells
    dist_mult_normal = MvNormal(vec(pop_param), cov_mat)
    param_cells = transpose(exp.(rand(dist_mult_normal, n_cells_simulate)))

    # Define matrix where to store the result
    time_span = param_info.time_span
    time_span_int = range(time_span[1], time_span[2], length = 200)
    n_states = length(state_info.states)
    mat_sim_result = zeros(n_cells_simulate*200, n_states+2)
    mat_sim_result[:, (end-1)] = repeat(time_span_int, inner=1, outer = n_cells_simulate)
    mat_sim_result[:, end] = repeat(1:n_cells_simulate, inner = 200)

    @showprogress 1 "Simulating cells ... " for i in 1:n_cells_simulate
        # Note, add noise to use the same struct as for parmeter estimation
        param_cell = vcat(param_cells[i, :], 1.0)
        param_vec = map_opt_vec_to_model(param_start, param_cell, param_info)
        dde_sol = solve_dde_system(model, param_vec, param_info)

        # Add interpolated soluton
        if dde_sol == "exit"
            continue
        end
        # Add interpolated solution
        for j in 1:200
            i_mat = (i-1) * 200 + j
            mat_sim_result[i_mat, 1:(end-2)] = dde_sol(time_span_int[j])
        end
    end

    # Write the result to disk
    result_to_save = convert(DataFrame, mat_sim_result)
    names_save = vcat(state_info.states, "time", "id")
    rename!(result_to_save, Symbol.(names_save))
    path_save = dir_save * "Simulated_cells.csv"
    CSV.write(path_save, result_to_save)
end


# Function that perform the STS approach for a model. It estimates the
# individual parameters (writes to file), the individual predictions
# (writes to file) and the mean population parameters and covaraince matrix
# (write to file). The user can choose from one of three optimisation
# algorithms that do not require a derivative.
# Args:
#   state_info, a state_info object for the model
#   start_guess, a start guess object for the model
#   model, the model top optimise
#   model_name, the name of the model to optmise
#   alg_use, the optmisation algorithm to use
#       1 = LN_NELDERMEAD (default)
#       2 = LN_SBPLX,
#       3 = LN_BOBYQA
# Returns:
#   void
function perform_STS(state_info, start_guess, model; alg_use=1)

    # Sanity check the start guess input
    len_input = length(start_guess.init_values) + length(start_guess.delays) +
        length(start_guess.noise) + length(start_guess.rate_const)
    if length(start_guess.names_opt_vec) != len_input
        @printf("Error: The names for the start guesses have a length \n")
        @printf("that does not match the number of start-guesses")
        return 1
    end

    # The optmisation details
    alg_list = [:LN_NELDERMEAD, :LN_SBPLX, :LN_BOBYQA]

    # Check that the directory where to save the result exists
    model_name = string(model)
    dir_sts = "../../Intermediate/STS/"
    dir_model = dir_sts * model_name * "/"
    dir_save = dir_model * string(alg_list[alg_use]) * "/"
    if !isdir(dir_sts) mkdir(dir_sts) end
    if !isdir(dir_model) mkdir(dir_model) end
    if !isdir(dir_save) mkdir(dir_save) end

    opt_mat, param_start, param_info = estimate_ind_param(dir_save, state_info,
        start_guess, model, alg_list[alg_use])

    pop_param, cov_mat = ml_est_cov_and_mean(dir_save, start_guess, opt_mat)

    simulate_cells(pop_param, cov_mat, param_start, param_info, model,
        state_info, dir_save)

    return pop_param, cov_mat, param_start, param_info
end


# ---------------------------------------------------------------------------
# Generate start guess functions
# ---------------------------------------------------------------------------


# Function that will perturb the start-guess vector and thus generate random
# start guesses using a multiple shooting method for optimising the mean data
# Args:
#   start_guess, the start guess object
#   state_info, the state info objects
#   perturb_vec, vector saying how much a vector is going to be perturbed
# Returns
#   param_start, opt_vec, param_info
function generate_rand_start_guess(start_guess, state_info, perturb_vec)

    param_start, param_info, opt_vec = create_param_obj(start_guess, state_info)
    # Sanity check input
    if length(perturb_vec) != length(opt_vec)-1
        @printf("Error: Perturb vec has a length not matching opt-vec\n")
        return "Error"
    end

    # Generate the random perturbance
    rand_vec = rand(Float64, length(perturb_vec))
    for i in 1:length(rand_vec)
        opt_vec[i] += (rand_vec[i]-0.5) * perturb_vec[i]
    end
    param_start = map_opt_vec_to_model(param_start, opt_vec, param_info)

    return param_start, param_info, opt_vec
end


# Function that will generate a start-guess for a given model by a
# multiple shooting method. These starting values are then used
# for comparing the STS and nlme-approaches
# Args:
#   state_info object for the current model
#   start_guess, a start-guess object for the model
#   perturb_vec, how the start-guesses should be perturbed
#   model, the relevant model
#   alg_choose, which algoruthm to use
#   times_run, the number of runs
# Returns:
#   best_start_guess, the best start guess
function generate_start_guess(state_info, start_guess, perturb_vec, model;
    alg_choose=3, times_run=100)

    Random.seed!(123)
    # Read data
    path_data = "../../Intermediate/Data_files/Data_mean_SUC2.csv"
    data = CSV.read(path_data)
    best_cost = 2000

    # The optmisation details
    alg_list = [:LN_NELDERMEAD, :LN_SBPLX, :LN_BOBYQA]
    alg_opt = alg_list[alg_choose]

    param_best = create_param_obj(start_guess, state_info)

    @showprogress 1 "Finding start guess ..." for i in 1:times_run
        # Generate the random start guess
        param_start, param_info, opt_vec = generate_rand_start_guess(start_guess,
            state_info, perturb_vec)
        n_param = length(opt_vec)

        # Solve the optmisation problem
        opt = Opt(alg_opt, n_param)
        opt.lower_bounds = 0.0
        opt.upper_bounds = Inf
        opt.maxeval = 1000

        min_objective!(opt, (opt_vector, grad) -> target_function(opt_vector, grad,
            param_info, state_info, data, param_start, model))
        (minf, minx, ret) = optimize(opt, opt_vec)
        param_final = map_opt_vec_to_model(param_start, minx, param_info)

        if minf < best_cost
            param_best = deepcopy_param_val(param_final)
            best_cost = minf
        end
    end

    # Print and return the best result
    @printf("\nBest cost = %.3f\n", best_cost)
    @printf("Rate constants\n")
    print(param_best.rate_constants)
    @printf("\nInital values\n")
    print(param_best.initial_values)
    @printf("\nDelays\n")
    print(param_best.delays)
    @printf("\n")

    return param_best
end
