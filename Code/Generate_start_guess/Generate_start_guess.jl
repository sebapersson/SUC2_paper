# This file will generate a start guess for certain model by using
# a multiple shooting method on the mean-value of the data.
# The optimisation is performed using derivative free algorithms.
# The user can choose from three algorithms (see) below.
# Args:
#   model, either simple_feedback or snf1_feedback
#   alg_use, the optmisation algorithm to use
#       1 = LN_NELDERMEAD
#       2 = LN_SBPLX
#       3 = LN_BOBYQA (default, and used in the paper)
#   times_run, the number of start-guesses
# Returns:
#   An error if the input is incorrect
#   printing the best cost and corresponding parameters

# Example: ./path_julia simple_feedback 1 100 will find start guesses
# for the simple feedback model using the nedler-mead algorithm and
# 100 start guesses.

# To run a model, fill in intervall for start-guesses as below.

include("../Julia_functions/Models.jl")
include("../Julia_functions/STS_find_start_guess.jl")

@printf("\n")
if length(ARGS) != 3
    @printf("Error, not the correct number of input arguments\n")
    exit(1)
end

model_use = ARGS[1]
alg_use = parse(Int64, ARGS[2])
times_run = parse(Int64, ARGS[3])

if model_use == "simple_feedback"
    @printf("\nSimple feedback model\n")
    # Simple feedback model
    # For creating the state info object
    state_info = produce_state_info(["Mig1", "SUC2", "X"], ["SUC2"],
        [1.0, "u1", 0.0])
    # A parameter is perturbed +/- half what is in perturb_vec vector, so here
    # k1 is perturbed 5 +/- 9.8 / 2, tau 112 +/- 40 and SUC2 3.82 +/- 0.5.
    perturb_vec = [9.8, 9.8, 9.8, 9.8, 9.8, 9.8, 9.8, 80, 0.5]
    start_guess = StartGuess([5, 5, 5, 5, 5, 5, 5], [112], [3.82], [0.5],
        ["k1", "k3", "k4", "k5", "k7", "k8", "k9", "tau2", "SUC20", "a1"])
    generate_start_guess(state_info, start_guess, perturb_vec,
        simple_feedback_model, times_run=times_run, alg_choose=alg_use)

elseif model_use == "snf1_feedback"
    @printf("SNF1 feedback model\n")
    state_info = produce_state_info(["SNF1p", "Mig1", "Mig1p", "SUC2", "X"], ["SUC2"], [0.0, 1.0, 0.0, "u1", 0.0])
    perturb_vec = [9.8, 9.8, 9.8, 9.8, 9.8, 9.8, 9.8, 9.8, 0.5]
    start_guess = StartGuess([5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], [], [3.82], [1.0],
        ["k1", "k3", "k4", "k5", "k6", "k8", "k9", "SUC20", "a1"])
    generate_start_guess(state_info, start_guess, perturb_vec,
        snf1_feedback_model, times_run=times_run, alg_choose=alg_use)

else

    @printf("Error: Model provided does not exist\n")
    exit(1)
end
