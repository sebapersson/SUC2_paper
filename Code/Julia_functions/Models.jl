# This file contains the DDE-models in a format where they can be solved using
# the DifferentialEquations package in Julia.
# The functions are split into two parts. First is the standard functions,
# then follows the model functions used when doing special parts with the
# models. This is for example when deleting a state, or running the second
# run in the STS estimation where certain parameters are kept fixed.

# The simple three state feedback model.
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function simple_feedback_model(du, u, h, p, t)

    k1, k3, k5, k6, k7, k8, k9, tau1, tau2, Mig10, SUC20, X0  = p

    # The time-dealys, tau1 = 32 min
    hist_Mig1 = h(p, t - tau1)[1]
    hist_X = h(p, t - tau2)[3]

    # For making the reading easier
    Mig1, SUC2, X = u[1], u[2], u[3]

    # To emulate the glucose cut
    if t < 0.0483
       rate_in = k1
       HX = 0
    else
       HX = 1
       rate_in = k1 / 40
    end

    # Relationships from assuming steady state
    k4 = k3 / ((k7 + Mig10*Mig10) * SUC20)
    k2 = k1 / Mig10

    # The dynamics
    du[1] = rate_in - k2*Mig1 + k5 * hist_X
    du[2] = k3 / (k7 + hist_Mig1^2) - k4 * SUC2
    du[3] = HX * k8 / (k9 + Mig1) - k6 * X
end


# The simple three state feedback model where k7 and k9 are kept fixed
# when simulated. Note, k7 and k9 must be added manually
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function simple_feedback_model_fixed_nlme(du, u, h, p, t)

    k1, k3, k5, k6, k8, tau1, tau2, Mig10, SUC20, X0  = p

    # The time-dealys, tau1 = 32 min
    hist_Mig1 = h(p, t - tau1)[1]
    hist_X = h(p, t - tau2)[3]

    # For making the reading easier
    Mig1, SUC2, X = u[1], u[2], u[3]

    # To emulate the glucose cut
    if t < 0.0483
       rate_in = k1
       HX = 0
    else
       HX = 1
       rate_in = k1 / 40
    end

    # Fixed parameters
    k7 = 0.65911102273653
    k9 = 14.1641123428758

    # Relationships from assuming steady state
    k4 = k3 / ((k7 + Mig10*Mig10) * SUC20)
    k2 = k1 / Mig10

    # The dynamics
    du[1] = rate_in - k2*Mig1 + k5 * hist_X
    du[2] = k3 / (k7 + hist_Mig1^2) - k4 * SUC2
    du[3] = HX * k8 / (k9 + Mig1) - k6 * X
end


# The snf1 feedback model without any deletions (wild-type)
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_model(du, u, h, p, t)

    k1, k3, k4, k5, k6, k8, k9, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

    # The maturation time delay
    hist_Mig1 = h(p, t - tau1)[2]

    SNF1p, Mig1, Mig1p, SUC2, X = u

    # Glucose down shift
    if t < 0.0483
       rate_down = k1
    else
       rate_down = k1 / 40.0 * SNF1p
    end

    # From assuming initial steady state
    k2 = -k4 * Mig1p0 + k5 * Mig10
    k7 = k6 / ((0.1 + Mig10^2) * SUC20)

    # SUC2 dynamics
    du[1] = k1 - rate_down - k9*SNF1p*X
    du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
    du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
    du[4] = k6 / (0.1 + hist_Mig1^2) - k7*SUC2
    du[5] = k8 * (SUC2 - SUC20) - k9*X
end


# The snf1 feedback model without any deletions (wild-type)
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_small_switch_model(du, u, h, p, t)

    k1, k3, k4, k5, k6, k8, k9, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

    # The maturation time delay
    hist_Mig1 = h(p, t - tau1)[2]

    SNF1p, Mig1, Mig1p, SUC2, X = u

    # Glucose down shift
    if t < 0.0483
       rate_down = k1
    else
       rate_down = k1 / 2.677 * SNF1p
    end

    # From assuming initial steady state
    k2 = -k4 * Mig1p0 + k5 * Mig10
    k7 = k6 / ((0.1 + Mig10^2) * SUC20)

    # SUC2 dynamics
    du[1] = k1 - rate_down - k9*SNF1p*X
    du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
    du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
    du[4] = k6 / (0.1 + hist_Mig1^2) - k7*SUC2
    du[5] = k8 * (SUC2 - SUC20) - k9*X
end


# The snf1 feedback model snf1 being deleted
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_d_snf1_model(du, u, h, p, t)

   k1, k3, k4, k5, k6, k8, k9, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

   # The maturation time delay
   hist_Mig1 = h(p, t - tau1)[2]

   SNF1p, Mig1, Mig1p, SUC2, X = u

   # Glucose down shift
   if t < 0.0483
      rate_down = k1
   else
      rate_down = k1 / 40.0 * SNF1p
   end

   # From assuming initial steady state
   k2 = -k4 * Mig1p0 + k5 * Mig10
   k7 = k6 / ((0.1 + Mig10^2) * SUC20)

   # SUC2 dynamics
   du[1] = 0.0
   du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
   du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
   du[4] = k6 / (0.1 + hist_Mig1^2) - k7*SUC2
   du[5] = k8 * (SUC2 - SUC20) - k9*X
end


# The snf1 feedback model with x and snf1 being deleted
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_d_snf1_x_model(du, u, h, p, t)

   k1, k3, k4, k5, k6, k8, k9, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

   # The maturation time delay
   hist_Mig1 = h(p, t - tau1)[2]

   SNF1p, Mig1, Mig1p, SUC2, X = u

   # Glucose down shift
   if t < 0.0483
      rate_down = k1
   else
      rate_down = k1 / 40.0 * SNF1p
   end

   # From assuming initial steady state
   k2 = -k4 * Mig1p0 + k5 * Mig10
   k7 = k6 / ((0.1 + Mig10^2) * SUC20)

   # SUC2 dynamics
   du[1] = 0.0
   du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
   du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
   du[4] = k6 / (0.1 + hist_Mig1^2) - k7*SUC2
   du[5] = 0.0
end


# The snf1 feedback model with x being deleted
# Args:
#  du, the derivates (act as output)
#  u, the model states
#  h, the time-delay function
#  p, the model paraemters
#  t, the time
function snf1_feedback_d_x_model(du, u, h, p, t)

   k1, k3, k4, k5, k6, k8, k9, tau1, SNF1p0, Mig10, Mig1p0, SUC20, X0  = p

   # The maturation time delay
   hist_Mig1 = h(p, t - tau1)[2]

   SNF1p, Mig1, Mig1p, SUC2, X = u

   # Glucose down shift
   if t < 0.0483
      rate_down = k1
   else
      rate_down = k1 / 40.0 * SNF1p
   end

   # From assuming initial steady state
   k2 = -k4 * Mig1p0 + k5 * Mig10
   k7 = k6 / ((0.1 + Mig10^2) * SUC20)

   # SUC2 dynamics
   du[1] = k1 - rate_down - k9*SNF1p*X
   du[2] = k2 - k3*SNF1p*Mig1 + k4*Mig1p - k5*(1 + 1 /(1 + exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1
   du[3] = k3*SNF1p*Mig1 - k4*Mig1p - k5*(1 + 1 /(1 +  exp(-3.0 * (SNF1p-4.5/3) ) ) ) * Mig1p
   du[4] = k6 / (0.1 + hist_Mig1^2) - k7*SUC2
   du[5] = 0.0
end
