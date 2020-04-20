library(tidyverse)
library(RColorBrewer)
library(latex2exp)
library(stringr)
library(xtable)
library(ggthemes)
library(RColorBrewer)
RNGkind("L'Ecuyer-CMRG")

# This file contains the functions that take the exported data from monolix and creates graphs from it. 
# It also plots simulated monolix data, and potential deletion data produced from the monolix 
# models. Example on how to run it using the main function, see the bottom of the file. 
# Output:
#   eta-shrinkage (printed)
#   IPRED-plot
#   IWPRES-plots
#   qq-plots of the estimated parameters 
#   Randomly selected individial plots 
#   qq-plot of IWRES
#   Population parameters (printed and exported to LaTeX table)
#   Simulated cells 


# General plotting parameters 
my_theme <- theme_tufte(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                               plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_colors <- c("#1f78b4", "#33a02c", "#08519c", "#006d2c", "#ff7f00", "#a6cee3")
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0


# Function that will calculate, and print, the eta-shrinkage for a model. 
# Inputs:
#   path_to_result, the path to the result folder
#   out_signals, the out-signals from a model i.e y1 and y2 
# Reaturns:
#   void 
calculate_shrinkage_model <- function(path_to_result, out_signals)
{
  # Get the population values for calculations 
  out_signals <- "observation"
  path_pop_val <- str_c(path_to_result, "/populationParameters.txt")
  pop_val_raw <- read_csv(path_pop_val, col_types = cols())
  var_values <- pop_val_raw %>%
    filter(grepl("^a", parameter)) 
  est_sigma <- var_values$value
  
  # Calculate the epsilon shrinkage
  IWRES_list <- lapply(1:length(out_signals), function(i){
    path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/", out_signals[i], "_obsVsPred.txt")
    obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
      select(-censored, -split, -color, -filter) %>%
      filter(!(time < 0.0)) %>%
      rename("observed" = out_signals[i]) %>%
      rename("IPRED" = "indivPredMode") %>%
      mutate(IWRES = (observed - IPRED) / est_sigma[i])
    return(obs_vs_pred_data$IWRES)})
  
  IWRES_vec <- do.call(c, IWRES_list)
  eps_shrinkage <- 1 - sd(IWRES_vec)
  print(sprintf("Epsilon shrinkage = %.3f", eps_shrinkage))
  
  return(0)
}


# Function that will plot the IPRED results for a given model, note that a specific 
# plot is created for each output. 
# Args:
#   path_to_result, the path to the result folder 
#   outputs, the outputs on a string format 
#   out_signal_name, the real name of the outputs 
#   dir_save, the directory to save the result in 
# Returns:
#   void 
plot_IPRED <- function(path_to_result, out_signals, out_signals_name, dir_save)
{
  
  # Get the population values for calculations 
  path_pop_val <- str_c(path_to_result, "/populationParameters.txt")
  pop_val_raw <- read_csv(path_pop_val, col_types = cols())
  var_values <- pop_val_raw %>%
    filter(grepl("^a", parameter)) 
  est_sigma <- var_values$value
  
  # Calculate the epsilon shrinkage
  IPRED_list <- lapply(1:length(out_signals), function(i){
    path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/", out_signals[i], "_obsVsPred.txt")
    obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
      select(-censored, -split, -color, -filter) %>%
      filter(!(time < 0.0)) %>%
      rename("observed" = out_signals[i]) %>%
      rename("IPRED" = "indivPredMode") 
    return(obs_vs_pred_data)})
  
  
  # Plotting observation plot 
  plot_list <- lapply(IPRED_list, function(data){
    p <- ggplot(data, aes(y = observed, x = IPRED)) + 
      geom_point(alpha = 0.6, size = 0.8) + 
      geom_smooth(aes(color = "Trend line"), method = "lm", size = 1.0, se = F, show.legend = NA) + 
      geom_abline(aes(slope = 1, intercept = 0, color = "y = x"), size = 1.0, show.legend = F, linetype = 2) + 
      labs(x = "Individual predictions", y = "Observed values") +
      scale_color_manual(name = "Line", values = my_colors[c(1, 5)]) + 
      geom_rangeframe() + 
      my_theme + theme(legend.position = "none")
    return(p)})
  
  # Save the figure 
  fig_names <- str_c(dir_save, "IPRED_", out_signals_name, ".pdf")
  for(i in 1:length(fig_names)) ggsave(fig_names[i], plot=plot_list[[i]], height = BASE_HEIGHT, width = BASE_WIDTH)
  
  return(0)
}


# Function that will plot the IWRES (individually weighted residuals)
# for a model. Each residual plot is saved to a pdf-file, and is plotted
# against time and observed value 
# Args:
#   path_to_result, the path to the result folder 
#   out_signals, the outputs on a string format 
#   out_signal_name, the real name of the outputs 
#   dir_save, the directory to save the result in 
# Returns:
#   void 
plot_IWRES <- function(path_to_result, out_signals, out_signals_name, dir_save)
{
  # Get the population values for calculations 
  path_pop_val <- str_c(path_to_result, "/populationParameters.txt")
  pop_val_raw <- read_csv(path_pop_val, col_types = cols())
  var_values <- pop_val_raw %>%
    filter(grepl("^a", parameter)) 
  est_sigma <- var_values$value
  
  # Calculate the epsilon shrinkage
  IWRES_list <- lapply(1:length(out_signals), function(i){
    path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/", out_signals[i], "_obsVsPred.txt")
    obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
      select(-censored, -split, -color, -filter) %>%
      filter(!(time < 0.0)) %>%
      rename("observed" = out_signals[i]) %>%
      rename("IPRED" = "indivPredMode") %>%
      mutate(IWRES = (observed - IPRED) / est_sigma[i])
    return(obs_vs_pred_data)})
  
  names_save <- lapply(out_signals_name, function(name){
    str_c(dir_save, c("IWRES_time_", "IWRES_obs_"), name, ".pdf")})
  
  # Plot IWRES vs time 
  lapply(1:length(IWRES_list), function(i){
    p_time <- ggplot(IWRES_list[[i]], aes(time, IWRES)) + 
      geom_point(size = 0.8, alpha = 0.6) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
      geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
      geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
      geom_smooth(method = "loess", se = F, color = my_colors[1], size = 1.2) + 
      labs(x = "Time", y = "IWRES") +
      geom_rangeframe() + 
      my_theme
    
    # Plot IWRES vs Suc2 expression 
    p_obs <- ggplot(IWRES_list[[i]], aes(observed, IWRES)) + 
      geom_point(size = 0.8, alpha = 0.6) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
      geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
      geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
      geom_smooth(method = "loess", se = F, color = my_colors[1], size = 1.2, linetype = 2) + 
      labs(x = "SUC2-intensity", y = "IWRES") +
      geom_rangeframe() +
      my_theme
    
    ggsave(names_save[[i]][1], plot = p_time, width = 9, height = 5)
    ggsave(names_save[[i]][2], plot = p_obs, width = 9, height = 5)})
  
  return(0)
}


# Function that will produce a qq-plot of the IWRES-residuals. 
# Args:
#   path_to_result, the path to the result folder 
#   out_signals, the outputs on a string format 
#   dir_save, the directory to save the result in 
# Returns:
#   void
plot_qq_IWRES <- function(path_to_result, out_signals, dir_save)
{
  # Get the population values for calculations 
  path_pop_val <- str_c(path_to_result, "/populationParameters.txt")
  pop_val_raw <- read_csv(path_pop_val, col_types = cols())
  var_values <- pop_val_raw %>%
    filter(grepl("^a", parameter)) 
  est_sigma <- var_values$value
  
  # Get all the IWRES-values 
  IWRES_list <- lapply(1:length(out_signals), function(i){
    path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/", out_signals[i], "_obsVsPred.txt")
    obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
      select(-censored, -split, -color, -filter) %>%
      filter(!(time < 0.0)) %>%
      rename("observed" = out_signals[i]) %>%
      rename("IPRED" = "indivPredMode") %>%
      mutate(IWRES = (observed - IPRED) / est_sigma[i])
    return(obs_vs_pred_data$IWRES)})
  IWRES_vec <- tibble(IWRES = do.call(c, IWRES_list))
  
  # Create a qqplot 
  p <- ggplot(IWRES_vec, aes(sample = IWRES)) + 
    geom_qq() + 
    geom_qq_line() + 
    my_theme
  
  file_save <- str_c(dir_save, "IWRES_qq.pdf")  
  ggsave(file_save, plot = p, width = 9, height = 5)
  
  return(0)
}


# Function that plots the distribution for the individual parameters using 
# samples taken from the conditional distribution of the data. 
# Args:
#   path_to_result, the path to the result folder 
#   dir_save, the directory to save the result in 
# Returns:
#   void 
plot_dist_param <- function(path_to_result, dir_save, param_save = F)
{
  # Read sampled parameters from conditional distribution 
  path_sample_param <- str_c(path_to_result, "/IndividualParameters/simulatedIndividualParameters.txt")
  sampled_param_raw <- read_csv(path_sample_param, col_types = cols())
  
  # General parameters 
  n_param <- dim(sampled_param_raw)[2] - 2
  param_names <- names(sampled_param_raw)[3:(n_param+2)]
  
  # Compute the qq-plots 
  plot_list_qq <- lapply(1:n_param, function(i){
    # Get the sampled values
    sampled_val <- sampled_param_raw %>% 
      select(param_names[i])
    
    names(sampled_val) <- "to_plot"
    p1 <- ggplot(sampled_val, aes(sample = log(to_plot))) + 
      geom_qq_line() + 
      geom_qq() + 
      labs(title = param_names[i]) +
      geom_rangeframe() + 
      my_theme
    
    return(p1)})
  
  # The plot always has five columns 
  n_row <- ceiling(n_param / 5)
  p <- ggpubr::ggarrange(plotlist = plot_list_qq, nrow = n_row, ncol = 5)
  file_save <- str_c(dir_save, "qq_param.pdf")
  ggsave(file_save, plot = p, width = 12, height = 7)
  
  if(param_save != F){
    data_par <- sampled_param_raw %>%
      select(param_save) %>%
      mutate_all(log) %>%
      rename("obs" = param_save)
    data_point <- tibble(x = qqnorm(data_par$obs)$x, y = qqnorm(data_par$obs)$y)
    data_line <- robcbi::QQline(data_par$obs)
    p <- ggplot(data_point, aes(x, y)) + 
      geom_abline(slope=data_line[2], intercept = data_line[1], size = 1.0, color = my_colors[1]) + 
      geom_point(size = 1.0, alpha = 0.8) + 
      geom_rangeframe() + 
      labs(x = "Theorethical", y = "Observed") + 
      my_theme
    file_save <- str_c(dir_save, param_save, "_qq.pdf")
    ggsave(file_save, plot = p, width = BASE_WIDTH, height = BASE_HEIGHT)
  }
  
  return(0)
}


# Function that will plot 2 random individual plots (deafult argument) for each 
# output in the model. 
# Args:
#   path_to_result, the path to the result folder 
#   out_signals, the outputs on a string format 
#   out_signal_name, the real name of the outputs 
#   dir_save, the directory to save the result in 
#   n_samples, the number of cell to plot per individual 
# Returns:
#   void 
plot_individual_fit <- function(path_to_result, out_signals, out_signals_name, dir_save, n_samples=2, seed=123)
{
  set.seed(seed)
  # Individual fits 
  path_fits <- str_c(path_to_result, "/ChartsData/IndividualFits/", out_signals, "_fits.txt")
  fit_data <- lapply(path_fits, function(path) read_csv(path, col_types = cols()))
  
  # Observations 
  path_observations <- str_c(path_to_result, "/ChartsData/IndividualFits/", out_signals, "_observations.txt")
  obs_data <- lapply(1:length(path_observations), function(i){ 
    read_csv(path_observations[i], col_types = cols()) %>%
      rename("observed" = out_signals[i])})
  
  # Get the random observations to plot 
  obs_to_plot <- lapply(obs_data, function(data){
    id_s <- unique(data$id)
    to_plot <- sample(id_s, n_samples, replace = F)})
  
  # Plot the cells 
  for(i in 1:length(fit_data)){
    for(which_cell in obs_to_plot[[i]]){
      data_fit_cell <- fit_data[[i]] %>% filter(id == which_cell)
      data_obs_cell <- obs_data[[i]] %>% filter(id == which_cell)
      
      p1 <- ggplot(data_fit_cell, aes(time, indivPredMode)) + 
        geom_line(size = 1.1, color = my_colors[1]) + 
        geom_point(data = data_obs_cell, aes(time, observed), color = "black", size = 1.6) + 
        labs(x = "Time", y = "Suc2 intensity") + 
        geom_rangeframe() + 
        my_theme
      
      # Save the figure 
      file_save <- str_c(dir_save, "Ind_plot_", out_signals_name[i], "_", which_cell, ".pdf")
      ggsave(file_save, plot = p1, width = 9, height = 6)
    }
  }
  
  return(0)
}


# Function that will print the likelihood values for a model, more specifically it will print the 
# likelhood, AIC and BIC values . 
# Args:
#   path_to_result, path to the result folder 
# Returns:
#   void
print_likelihood_values <- function(path_to_result)
{
  # Process the likelihood result 
  path_likelihood <- str_c(path_to_result, "/LogLikelihood/logLikelihood.txt")
  if(!file.exists(path_likelihood)) return(1)
  likelihood_data <- read_csv(path_likelihood, col_types = cols())
  
  # Print the data 
  print(sprintf("Log-likelhood = %.2f", likelihood_data[1, 2]))
  print(sprintf("AIC = %.2f", likelihood_data[2, 2]))
  print(sprintf("BIC = %.2f", likelihood_data[4, 2]))
  return(0)
}


# Function that will create code for a LaTeX table out of the population parameters. 
# It will also save the population parameters as a csv-file in the result folder, 
# togehter with saving the table in a file 
# Args:
#   path_to_result, the path to the result folder 
#   dir_save, the directory where the result is saved 
# Returns:
#   void 
create_latex_table_pop_param <- function(path_to_result, dir_save)
{
  # Read in the population parameters and make latex-table 
  path_pop_param <- str_c(path_to_result, "/populationParameters.txt")
  pop_param <- read_csv(path_pop_param, col_types = cols()) %>%
    mutate(value = as.character(round(value, 2))) %>%
    mutate(rse_sa = as.character(round(rse_sa, 1))) %>%
    select(-se_sa) %>%
    rename("Parameter" = "parameter") %>%
    mutate(value = str_c(value, " (", rse_sa, ")")) %>%
    select(-rse_sa) %>%
    rename("Value" = "value")
  
  path_save_par <- str_c(dir_save, "Pop_param.csv")
  path_save_table <- str_c(dir_save, "LaTeX_table.txt")
  latex_table <- xtable(pop_param)
  write_csv(pop_param, path_save_par)
  print(latex_table, file = path_save_table)
}


# Function that will process the simulated single-cell data for a given model. The function 
# also plots the observed data if the state is one of the observed states. 
# Args:
#   path_to_result, the path to the result folder
#   dir_save, path to the directory where the result is stored 
# Returns:
#   1 if the file does not exist 
plot_simulated_data <- function(path_to_result, dir_save, extrapolated=F)
{
  
  # Process the data into an aggregated format 
  if(extrapolated == F){
    path_data <- str_c(path_to_result, "/Simulated_cells.csv")
    if(!file.exists(path_data)) return(1)
  }else{
    path_data <- path_to_result
  }
  
  data_sim <- read_csv(path_data, col_types = cols(id = col_factor()))
  state_names <- names(data_sim)[1:(dim(data_sim)[2]-2)]
  
  # The special case a model with a sum is used 
  if("Mig1" %in% state_names && "Mig1p" %in% state_names){
    output_sum <- T
    state_names <- c(state_names, "Mig1t")
    data_sim <- data_sim %>%
      mutate(Mig1t = Mig1 + Mig1p)
  }else{
    output_sum <- F
  }
  
  data_tidy <- data_sim %>%
    gather(all_of(state_names), key = "state", value = "value") %>%
    mutate(state = as.factor(state))
  data_agg_sim <- data_tidy %>%
    group_by(state, t) %>%
    summarise(median = median(value, na.rm = T), 
              quant95 = quantile(value, 0.95, na.rm = T), 
              quant05 = quantile(value, 0.05, na.rm = T), 
              quant80 = quantile(value, 0.80, na.rm = T), 
              quant20 = quantile(value, 0.20, na.rm = T))
  
  # Read and process the observed data 
  data_SUC2 <- read_csv("../../Intermediate/Data_files/Data_monolix_SUC2.csv", 
                       col_types = cols(id = col_factor())) %>% 
    filter(time >= 0)
  
  # SUC2 data   
  data_SUC2_sum <- data_SUC2 %>%
    group_by(time) %>%
    summarise(median = median(observation, na.rm = T), 
              quant95 = quantile(observation, 0.95, na.rm = T), 
              quant05 = quantile(observation, 0.05, na.rm = T), 
              quant80 = quantile(observation, 0.80, na.rm = T), 
              quant20 = quantile(observation, 0.20, na.rm = T))
  
  # Make sure that the result can be saved
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  
  # Plot the different states 
  n_states <- length(state_names)
  plot_list <- lapply(state_names, function(state_name){
    
    data_plot <- data_agg_sim %>%
      filter(state == state_name)
    
    if(state_name == "SUC2" && extrapolated == F){
      data_obs <- data_SUC2
      data_sum <- data_SUC2_sum
      data_sum <- data_sum %>% mutate(type = "observed")
      data_sim_suc2 <- data_agg_sim %>% filter(state == "SUC2") %>% filter(t < 480)
      data_min_max <- tibble(x = c(min(data_sim_suc2$t), max(c(data_sim_suc2$t, 500))), 
                             y = c(min(data_sim_suc2$quant05), max(data_sim_suc2$quant95)))

      p <- ggplot(data_sim_suc2, aes(t, median)) + 
        geom_line(size = 1.2, color = my_colors[1]) + 
        geom_ribbon(aes(ymin = quant05, ymax = quant95), color = NA, alpha = 0.25, fill = my_colors[1]) +
        geom_ribbon(aes(ymin = quant20, ymax = quant80), color = NA, alpha = 0.35, fill = my_colors[1]) +
        geom_line(aes(t, quant05), size = 1.2, color = my_colors[1], alpha = 0.8) +
        geom_line(aes(t, quant95), size = 1.2, color = my_colors[1], alpha = 0.8) +
        geom_line(data_sum, mapping=aes(time, median), color = "black", size = 1.2, linetype = 2) +
        geom_line(data_sum, mapping=aes(time, quant95), color = "black", size = 1.2, linetype = 2) +
        geom_line(data_sum, mapping=aes(time, quant05), color = "black", size = 1.2, linetype = 2) +
        geom_rangeframe(data=data_min_max, mapping=aes(x=x, y=y)) +
        labs(x = "Time [min]", y = TeX("SUC2 intesity \\[A.U $\\times 10^{-2}$\\]")) + 
        my_theme
      
    }else{
      data_plot <- data_agg_sim %>% filter(state == state_name)
      data_min_max <- tibble(x = c(min(data_plot$t), max(data_plot$t)), 
                             y = c(min(data_plot$quant05), max(data_plot$quant95)))
      p <- ggplot(data_plot) + 
        geom_line(aes(t, median), size = 1.2, color = my_colors[1]) + 
        geom_ribbon(aes(x = t, ymin = quant05, ymax = quant95), color = NA, alpha = 0.25, fill = my_colors[1]) +
        geom_ribbon(aes(x = t, ymin = quant20, ymax = quant80), color = NA, alpha = 0.35, fill = my_colors[1]) +
        geom_line(aes(t, quant05), size = 1.2, color = my_colors[1], alpha = 0.8) +
        geom_line(aes(t, quant95), size = 1.2, color = my_colors[1], alpha = 0.8) +
        labs(x = "Time [min]", y = str_c("Simulated ", state_name)) +
        geom_rangeframe(data=data_min_max, aes(x=x, y=y)) + 
        my_theme
    }
    path_save <- str_c(dir_save, "Simulated_", state_name, ".pdf")
    ggsave(path_save, plot = p, width = BASE_WIDTH, height = BASE_HEIGHT)
    return(p)
  })
  
  # Save all plots 
  path_save <- str_c(dir_save, "Simulated_all.pdf")
  p <- ggpubr::ggarrange(plotlist = plot_list, ncol = n_states)
  ggsave(path_save, width = 22, height = 6)
  
  return(0)
}


# Functoin that makes a density plot of tau_X and saves it. 
# Args:
#   path_to_result, path to the reuslt vector
#   dir_save, the directory to save the result in 
#   x_int, the x-interval to plot the solution in 
# Returns:
#   void
plot_dist_tau_x <- function(path_to_result, dir_save, x_int = c(0, 500))
{
  my_blues <- c("#1f78b4", "#a6cee3")
  param <- "tau_X"
  path_ind_param <- str_c(path_to_result, "/populationParameters.txt")
  data <- read_csv(path_ind_param, col_types = cols()) %>%
    filter(grepl(str_c("^", param), parameter) | grepl(str_c("^omega_", param), parameter))
  
  pop_val <- as.numeric(data[1, 2])
  var <- as.numeric(data[2, 2])
  data_plot <- tibble(x = seq(from = x_int[1], to = x_int[2], length.out = 300)) %>%
    mutate(y = dlnorm(x, log(pop_val), var))
  
  p <- ggplot(data_plot, aes(x, y)) + 
    geom_line(color = my_blues[1], size = 1.2) + 
    geom_ribbon(aes(ymin = 0, ymax = y), color = NA, alpha = 0.3, fill = my_blues[1]) + 
    geom_rangeframe() + 
    labs(x = TeX("$\\tau$-feedback"), y = "Probability density") +
    theme_tufte()
  
  path_save <- str_c(dir_save, "Dist_tau_X.pdf")
  ggsave(path_save, plot = p, height = BASE_HEIGHT, width = BASE_WIDTH)
  
}


# Function that processes the deletion experiments for the snf1-feedback model. 
# It plots the end result as a bar-chart, where certain time-points are slected. 
# Args:
#   path_to_result, path to the result folder
#   dir_save, the directory in which to save the result 
# Returns:
#   void 
plot_snf1_model_deletions <- function(path_to_result, dir_save)
{
  
  # Process the mutations 
  delete_list <- c("wt", "dsnf1", "dsnf1_x", "dx")
  file_names <- str_c(path_to_result, "/", "Simulated_cells_", delete_list, ".csv")
  data_list <- lapply(1:length(delete_list), function(i){
    data_raw <- read_csv(file_names[i], col_types = cols()) 
    data <- data_raw %>%
      select(t, id, SUC2) %>%
      group_by(t) %>%
      summarise(mean = mean(SUC2, na.rm = T), 
                sd = sd(SUC2, na.rm = T)) %>%
      mutate(type = delete_list[i])
    
    # Select columns to take 
    time_list <- c(0, 75, 150, 225, 300, 375)
    i_list <- sapply(time_list, function(i) which.min((data$t - i)^2))
    data <- data[i_list, ]
    data$t = time_list
    
    return(data)})
  
  data_exp <- do.call(rbind, data_list) %>%
    mutate(t = as.factor(t), type = as.factor(type)) %>%
    mutate(y_frame = mean + sd)
  data_exp$y_frame[1] <- 0
  
  my_blue_scale <- brewer.pal(9, "Blues")[-c(1, 2)]
  p <- ggplot(data_exp, aes(type, mean, fill = t)) + 
    geom_bar(stat='identity', position='dodge', color ="black") +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=0.2, position=position_dodge(.9)) + 
    scale_x_discrete(limit = c("wt", "dsnf1", "dsnf1_x", "dx"), 
                     labels = c("wt", "Δsnf1", "Δsnf1_x", "Δx")) +
    geom_rangeframe(aes(type, y_frame), sides = "l",) + 
    scale_fill_manual(values = my_blue_scale) +
    labs(y = TeX("SUC2 intensity \\[A.U $\\times 10^{-2}$\\]"), x = "") + 
    my_theme + theme(legend.position = "none", axis.title.x = element_blank())
  
  path_save <- str_c(dir_save, "Deletion_bar.pdf")  
  ggsave(path_save, plot = p, width = BASE_WIDTH, height = BASE_HEIGHT)
  
  return(0)
}


# The main function for processing the monolix result (it will call other funcitons)
# in this file. If certain data is not available, i.e the likelihood values, this 
# function will point out that certain data is lacking. 
# Args:
#   path_to_result, path to where the monolix result is stored 
#   out_signals, out-signal in str-format (i.e c(y1, y2))
#   out_signals_name, the name of the outsignals 
#   model_save, the name of the model to save, i.e Reg1_model3 
#   param_save, if a specific parameter is to be saved when plotting the qq-plots 
#   plot_dist_tau_x, control parameter that decides wheter or not to plot distribution of 
#   plot_Deletions, wheter or not deletions for the SNF1-feedback model should be plotted 
# Returns:
#   void 
process_monolix_result <- function(path_to_result, out_signals, out_signals_name, model_save, 
                                   param_save = F, plot_dist_tau_x = F, plot_deletions = F)
{
  # Where to store the result 
  dir_save <- str_c("../../Result/Monolix_processed/", model_save, "/")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)

  calculate_shrinkage_model(path_to_result, out_signals)
  
  plot_IPRED(path_to_result, out_signals, out_signals_name, dir_save)

  plot_IWRES(path_to_result, out_signals, out_signals_name, dir_save)
  
  plot_qq_IWRES(path_to_result, out_signals, dir_save)
  
  plot_dist_param(path_to_result, dir_save, param_save = param_save)
  
  create_latex_table_pop_param(path_to_result, dir_save)
  
  if(print_likelihood_values(path_to_result) != 0) print("Likelihood not calculated")
  
  plot_individual_fit(path_to_result, out_signals, out_signals_name, dir_save)
  
  if(plot_simulated_data(path_to_result, dir_save) != 0) print("Simulated data not available")
  
  # Should only occur for models with tau_x in them 
  if(plot_dist_tau_x) plot_dist_tau_x(path_to_result, dir_save)
  
  # Should only occur for snf1 model 
  if(plot_deletions) plot_snf1_model_deletions(path_to_result, dir_save)
}


# -----------------------------------------------------------------------
# Process the end monolix result 
# -----------------------------------------------------------------------

path_to_result <- "../Monolix_code/Simple_feedback/Simple_feedback"
process_monolix_result(path_to_result, "observation", "SUC2", "Simple_feedback", 
                       param_save = "k3", plot_dist_tau_x = T)

path_to_result <- "../Monolix_code/Snf1_feedback/Snf1_feedback"
process_monolix_result(path_to_result, "observation", "SUC2", "Snf1_feedback", 
                       plot_deletions = T)
