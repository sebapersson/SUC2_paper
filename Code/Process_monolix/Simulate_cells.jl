# This file will simulate cells using the nlme parameter estimates. Thus
# this file is used to simulate the population behaviour, but it is also
# used to simulate the deletion experiments if those are performed
# for a certain model.
# The file outputs a csv-file simulated data for each state. The file
# is tagged if it is a part of a deletion experiment.
# Args:
#   model, simple_feedback, complex_model or complex_model_m (d = deletions)
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

if model_use == "simple_feedback"
    @printf("Simple feedback model\n")
    # Relative path to result dir from Monolix
    result_dir = "../Monolix_code/Simple_feedback/Simple_feedback/"
    # u1, see documentation model_info
    model_info = ModelInfo(["Mig1", "SUC2", "X"],
        [1.0, "u1", 0.0], 3, simple_feedback_model_fixed_nlme)
    tau_m = 32
    simulate_cells_nlme(result_dir, n_cells_simulate, tau_m, model_info,
    fixed = ["k7_pop", "k9_pop"])
else
    @printf("Error: Model provided does not exist\n")
    exit(1)
end
