# This file will simulate cells using the nlme parameter estimates. Thus
# this file is used to simulate the population behaviour, but it is also
# used to simulate the deletion experiments if those are performed
# for a certain model.
# The file outputs a csv-file simulated data for each state. The file
# is tagged if it is a part of a deletion experiment.
# Args:
#   model, simple_feedback, snf1_feedback or snf1_feedback_d (d = deletions)
#   n_cells_simulate, the number of cells to simulate
# Returns:
#   An error if the input is incorrect
#   A csv-file with the simulated data

# Example: ./path_julia simple_feedback 10 will simulate ten cells
# for the simple feedback model

# To run a model, fill in extra model-options as below. If no parameters
# are fixed, argument fixed is not required.

include("../Julia_functions/Models.jl")
include("../Julia_functions/Simulate_NLME.jl")


@printf("\n")
if length(ARGS) != 2
    @printf("Error, not the correct number of input arguments\n")
    exit(1)
end

model_use = ARGS[1]
n_cells_simulate = parse(Int64, ARGS[2])


# u1, see documentation functions, note using relative paths
if model_use == "simple_feedback"
    @printf("Simple feedback model\n")
    # Simulation parameters
    result_dir = "../Monolix_code/Simple_feedback/Simple_feedback/"
    model_info = ModelInfo(["SNF1p", "SUC2", "X"],
        ["m", "m", "m"], 3, simple_feedback_model_nlme)
    tau_m = 32

    simulate_cells_nlme(result_dir, n_cells_simulate, tau_m, model_info,
        map_init_rates=map_init_simple_feedback_nlme, fixed=["k5_pop", "k8_pop"])
elseif model_use == "snf1_feedback"
    @printf("Snf1 feedback model\n")

    # Simulation parameters
    result_dir = "../Monolix_code/Snf1_feedback/Snf1_feedback/"
    model_info = ModelInfo(["SNF1p", "Mig1", "Mig1p", "SUC2", "X"],
        ["m", "m", "m", "m", "m"], 5, snf1_feedback_model)
    tau_m = 32

    simulate_cells_nlme(result_dir, n_cells_simulate, tau_m, model_info,
        map_init_rates=map_init_snf1_feedback)
elseif model_use == "snf1_feedback_d"
    @printf("Deletion experiments Snf1 feedback model\n")

    # Simulation parameters
    result_dir = "../Monolix_code/Snf1_feedback/Snf1_feedback/"
    model_list = [snf1_feedback_model, snf1_feedback_d_snf1_x_model,
        snf1_feedback_d_x_model, snf1_feedback_d_snf1_model,
        snf1_feedback_small_switch_model]
    tag_list = ["_wt", "_dsnf1_x", "_dx", "_dsnf1", "_small_switch"]

    for i in 1:length(model_list)
        model_info = ModelInfo(["SNF1p", "Mig1", "Mig1p", "SUC2", "X"],
            ["m", "m", "m", "m", "m"], 5, model_list[i])
        simulate_cells_nlme(result_dir, n_cells_simulate, 32.0,
            model_info, tag=tag_list[i],
            map_init_rates=map_init_snf1_feedback)
    end
else
    @printf("Error: Model provided does not exist\n")
    exit(1)
end
