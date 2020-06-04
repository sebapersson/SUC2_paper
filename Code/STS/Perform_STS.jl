# This file will perform the STS optimisation for a given model
# using a derivative free optimisation algorithm. For the simple
# feedback model there are two options, either were parameters are fixed
# or where parameters are not fixed.
# Note that running the STS will output several files (see below)
# Args:
#   model, either simple_feedback or simple_feedback_fixed
#   alg_use, the optmisation algorithm to use
#       1 = LN_NELDERMEAD
#       2 = LN_SBPLX
#       3 = LN_BOBYQA (default, and used in the paper)
# Returns:
#   An error if the input is incorrect
#   Outputing to intermediate:
#       The estimated parameters
#       The individal fits
#       Estimated covariance and population parameters
#       10 000 simulated cells

# Example: ./path_julia simple_feedback 1 will run the STS optimisation
# for the simple feedback model using the nedler-mead algorithm

# To run a model, fill in intervall for start-guesses as below.

include("../Julia_functions/Models.jl")
include("../Julia_functions/STS_find_start_guess.jl")

@printf("\n")
if length(ARGS) != 2
    @printf("Error, not the correct number of input arguments\n")
    exit(1)
end

model_use = ARGS[1]
alg_use = parse(Int64, ARGS[2])


if model_use == "simple_feedback"
    @printf("Simple feedback model STS\n")
    state_info = produce_state_info(["Mig1", "SUC2", "X"], ["SUC2"], ["m", "m", "m"])
    start_guess = StartGuess([0.087, 0.055, 4.83, 12.73, 2.09, 0.87, 0.62, 9.07, 5.17],
        [218.8], [], [1.0],
        ["k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "tau2", "a1"])

    perform_STS(state_info, start_guess, simple_feedback_model,
        alg_use=alg_use, map_init=map_init_simple_feedback)

elseif model_use == "simple_feedback_fixed"
    @printf("Simple feedback fixed parameter STS\n")
    state_info = produce_state_info(["Mig1", "SUC2", "X"], ["SUC2"], ["m", "m", "m"])
    start_guess = StartGuess([0.087, 0.055, 4.83, 12.73, 0.87, 0.62, 5.17],
        [218.8], [], [1.0],
        ["k1", "k2", "k3", "k4", "k6", "k7", "k9", "tau2", "a1"])

    perform_STS(state_info, start_guess, simple_feedback_model_sts_fixed,
        alg_use=alg_use, map_init=map_init_simple_feedback_sts_fixed)

elseif model_use == "simple_feedback_log"
    @printf("Simple feedback model STS log space\n")
    state_info = produce_state_info(["Mig1", "SUC2", "X"], ["SUC2"], ["m", "m", "m"])
    start_guess = StartGuess(log.([0.087, 0.055, 4.83, 12.73, 2.09, 0.87, 0.62, 9.07, 5.17]),
        [log(218.8)], [], [2.0],
        ["k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "tau2", "a1"])

    perform_STS(state_info, start_guess, simple_feedback_model_log,
        alg_use=alg_use, map_init=map_init_simple_feedback_log,
        log_space=true)

elseif model_use == "simple_feedback_log_fixed"
    @printf("Simple feedback fixed parameter log STS\n")
    state_info = produce_state_info(["Mig1", "SUC2", "X"], ["SUC2"], ["m", "m", "m"])
    start_guess = StartGuess(log.([0.087, 0.055, 4.83, 12.73, 0.87, 0.62, 5.17]),
        [log(218.8)], [], [1.0],
        ["k1", "k2", "k3", "k4", "k6", "k7", "k9", "tau2", "a1"])

    perform_STS(state_info, start_guess, simple_feedback_model_log_fixed,
        alg_use=alg_use, map_init=map_init_simple_feedback_log_fixed,
        log_space=true)
else
    @printf("Error: Model provided does not exist\n")
    exit(1)
end
