# This file contains the DDE-models in a format where they can be solved using
# the DifferentialEquations package in Julia.
# The functions are split into two parts. First is the standard functions,
# then follows the model functions used when doing special parts with the
# models. This is for example when deleting a state, or running the second
# run in the STS estimation where certain parameters are kept fixed.


# Function that will map the rate-constants to the initial values
# for the simple-feedback model. Note, here none of the parameters
# are assumed to be fixed
# Args:
#  rate_constants, the rate-constants for the model
#  n_states, a vector with the number of states
# Returns:
#  A vector with initial values
function map_init_simple_feedback(rate_constants)
   init_vec = zeros(Float64, 3)
   k1, k2, k3, k4, k5, k6, k7, k8, k9 = rate_constants
   init_vec[1] = k1 / k2
   init_vec[2] = k4 / ((k5 + init_vec[1]^1)*k6)
   init_vec[3] = 0.0

   return init_vec
end


# Function that will map the rate-constants to the initial values
# for the simple-feedback models. Here k5 and k8 are assumed
# to be fixed according to the NLME fit.
# Args:
#  rate_constants, the rate-constants for the model
#  n_states, a vector with the number of states
# Returns:
#  A vector with initial values
function map_init_simple_feedback_nlme(rate_constants)
   init_vec = zeros(Float64, 3)
   k1, k2, k3, k4, k6, k7, k9 = rate_constants
   k5 = 1.65
   init_vec[1] = k1 / k2
   init_vec[2] = k4 / ((k5 + init_vec[1]^1)*k6)
   init_vec[3] = 0.0

   return init_vec
end


# Function that will map the rate-constants to the initial values
# for the simple-feedback models. Here k5 and k8 are assumed
# to be fixed according to the STS fit.
# Args:
#  rate_constants, the rate-constants for the model
#  n_states, a vector with the number of states
# Returns:
#  A vector with initial values
function map_init_simple_feedback_sts_fixed(rate_constants)
   init_vec = zeros(Float64, 3)
   k1, k2, k3, k4, k6, k7, k9 = rate_constants
   k5 = 1.65
   init_vec[1] = k1 / k2
   init_vec[2] = k4 / ((k5 + init_vec[1]^1)*k6)
   init_vec[3] = 0.0

   return init_vec
end


# The simple three state feedback model.
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function simple_feedback_model(du, u, h, p, t)

    k1, k2, k3, k4, k5, k6, k7, k8, k9, tau1, tau2, SNF1p0, SUC20, X0  = p

    # The time-dealys, tau1 = 32 min
    hist_SNF1p = h(p, t - tau1)[1]
    hist_X = h(p, t - tau2)[3]

    # For making the reading easier
    SNF1p, SUC2, X = u[1], u[2], u[3]

    # Glucose downshift
    if t < 0.0483
       rate_in = k1
       HX = 0
    else
       HX = 1
       rate_in = k1 / 40
    end

    # Dynamics, feedback cascade model
    du[1] = rate_in - k2*SNF1p + k3 * hist_X
    du[2] = k4 / (k5 + hist_SNF1p^1) - k6 * SUC2
    du[3] = HX * k7 / (k8 + SNF1p) - k9 * X

end


# The simple three state feedback model. Here k5 and k8 are fixed
# according to the NLME fitted values
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function simple_feedback_model_nlme(du, u, h, p, t)

    k1, k2, k3, k4, k6, k7, k9, tau1, tau2, SNF1p0, SUC20, X0  = p

    # The time-dealys, tau1 = 32 min
    hist_SNF1p = h(p, t - tau1)[1]
    hist_X = h(p, t - tau2)[3]

    # For making the reading easier
    SNF1p, SUC2, X = u[1], u[2], u[3]

    # Glucose downshift
    if t < 0.0483
       rate_in = k1
       HX = 0
    else
       HX = 1
       rate_in = k1 / 40
    end

    # Fixed NLME values
    k5 = 1.65
    k8 = 20.2

    # Dynamics, feedback cascade model
    du[1] = rate_in - k2*SNF1p + k3 * hist_X
    du[2] = k4 / (k5 + hist_SNF1p) - k6 * SUC2
    du[3] = HX * k7 / (k8 + SNF1p) - k9 * X

end


# The simple three state feedback model. Here k5 and k8 are fixed
# according to the NLME fitted values
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function simple_feedback_model_sts_fixed(du, u, h, p, t)

    k1, k2, k3, k4, k6, k7, k9, tau1, tau2, SNF1p0, SUC20, X0  = p

    # The time-dealys, tau1 = 32 min
    hist_SNF1p = h(p, t - tau1)[1]
    hist_X = h(p, t - tau2)[3]

    # For making the reading easier
    SNF1p, SUC2, X = u[1], u[2], u[3]

    # Glucose downshift
    if t < 0.0483
       rate_in = k1
       HX = 0
    else
       HX = 1
       rate_in = k1 / 40
    end

    # Fixed NLME values
    k5 = 1.65
    k8 = 20.2

    # Dynamics, feedback cascade model
    du[1] = rate_in - k2*SNF1p + k3 * hist_X
    du[2] = k4 / (k5 + hist_SNF1p) - k6 * SUC2
    du[3] = HX * k7 / (k8 + SNF1p) - k9 * X
end


# Function that will map the rate-constants to the initial values
# for the simple-feedback models
# Args:
#  rate_constants, the rate-constants for the model
#  n_states, a vector with the number of states
# Returns:
#  A vector with initial values
function map_init_snf1_feedback(rate_constants)
   init_vec = zeros(Float64, 5)
   k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = rate_constants
   init_vec[1] = 0.0
   init_vec[2] = k2 / k5
   init_vec[3] = 0.0
   init_vec[4] = k6 / ((0.1 + init_vec[2]^1)*k7)
   init_vec[5] = 0.0

   return init_vec
end


# The snf1 feedback model (feedback mediated model)
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_model_v2(du, u, h, p, t)

    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

    # The maturation time delay
    hist_Mig1 = h(p, t - tau1)[2]

    SNF1p, Mig1, Mig1p, SUC2, X = u

    # Glucose down shift
    if t < 0.0483
       rate_down = k1
    else
       rate_down = k1 / 40.0 * SNF1p
    end

    # SUC2 dynamics
    du[1] = k1 - rate_down  - k10 * SNF1p * X
    du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
    du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
    du[4] = k6 / (0.1+ hist_Mig1^1) - k7 * SUC2
    du[5] =  k8 * (SUC2-SUC20) - k9*X
end


# The snf1 feedback model without any deletions (wild-type)
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_small_switch_model(du, u, h, p, t)

   k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

   # The maturation time delay
   hist_Mig1 = h(p, t - tau1)[2]

   SNF1p, Mig1, Mig1p, SUC2, X = u

   # Glucose down shift
   if t < 0.0483
      rate_down = k1
   else
      rate_down = k1 / 2.6667 * SNF1p
   end

   # SUC2 dynamics
   du[1] = k1 - rate_down  - k10 * SNF1p * X
   du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
   du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
   du[4] = k6 / (0.1+ hist_Mig1^1) - k7 * SUC2
   du[5] =  k8 * (SUC2-SUC20) - k9*X
end


# The snf1 feedback model snf1 being deleted
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_d_snf1_model(du, u, h, p, t)

   k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

   # The maturation time delay
   hist_Mig1 = h(p, t - tau1)[2]

   SNF1p, Mig1, Mig1p, SUC2, X = u

   # Glucose down shift
   if t < 0.0483
      rate_down = k1
   else
      rate_down = k1 / 2.6667 * SNF1p
   end

   # SUC2 dynamics
   du[1] = 0
   du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
   du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
   du[4] = k6 / (0.1+ hist_Mig1^1) - k7 * SUC2
   du[5] =  k8 * (SUC2-SUC20) - k9*X
end


# The snf1 feedback model with x and snf1 being deleted
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_d_snf1_x_model(du, u, h, p, t)

   k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

   # The maturation time delay
   hist_Mig1 = h(p, t - tau1)[2]

   SNF1p, Mig1, Mig1p, SUC2, X = u

   # Glucose down shift
   if t < 0.0483
      rate_down = k1
   else
      rate_down = k1 / 2.6667 * SNF1p
   end

   # SUC2 dynamics
   du[1] = 0
   du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
   du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
   du[4] = k6 / (0.1+ hist_Mig1^1) - k7 * SUC2
   du[5] = 0
end


# The snf1 feedback model with x being deleted
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_d_x_model(du, u, h, p, t)

   k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

   # The maturation time delay
   hist_Mig1 = h(p, t - tau1)[2]

   SNF1p, Mig1, Mig1p, SUC2, X = u

   # Glucose down shift
   if t < 0.0483
      rate_down = k1
   else
      rate_down = k1 / 2.6667 * SNF1p
   end

   # SUC2 dynamics
   du[1] = k1 - rate_down  - k10 * SNF1p * X
   du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
   du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
   du[4] = k6 / (0.1+ hist_Mig1^1) - k7 * SUC2
   du[5] = 0.0
end
