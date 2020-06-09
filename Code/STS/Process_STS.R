library(tidyverse) 
library(stringr)
library(latex2exp)
library(ggthemes)

# This file will process the STS output from the parameter estimation and produce 
# several diagnostic plots, which can be used to judge the quality of the fits. 
# Outputs:
#   plot of the simulated data for each state 
#   IPRED plot 
#   qq-plot of the individial estimates, user can also choose a specific parameter 
#     to plot. 


# General plotting parameters for gg-plot 
my_theme <- theme_tufte(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                                plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_colors <- c("#1f78b4", "#33a02c", "#08519c", "#006d2c", "#ff7f00", "#a6cee3")
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0


# Function that will plot the simulated data for a STS-model. The graphs will 
# follow the standard color-pallette, and use the tufte-theme. SUC2 
# is assumed to be the only observed data 
plot_simulated_data <- function(path_to_result, dir_save, extrapolated = F)
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
    group_by(state, time) %>%
    summarise(median = median(value, na.rm = T), 
              quant95 = quantile(value, 0.95, na.rm = T), 
              quant05 = quantile(value, 0.05, na.rm = T), 
              quant80 = quantile(value, 0.80, na.rm = T), 
              quant20 = quantile(value, 0.20, na.rm = T))
  
  # Read and process the observed data 
  data_obs <- read_csv("../../Intermediate/Data_files/Data_monolix_SUC2.csv", 
                       col_types = cols(id = col_factor())) %>% 
    filter(time >= 0)
  
  # SUC2 data   
  data_SUC2_sum <- data_obs %>%
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
      data_sum <- data_SUC2_sum
      data_sim_suc2 <- data_agg_sim %>% filter(state == "SUC2") %>% filter(time < 480)
      data_min_max <- tibble(x = c(min(data_sim_suc2$time), max(c(data_sim_suc2$time, 500))), 
                             y = c(2, 15))
                             #y = c(min(data_sim_suc2$quant05), max(data_sim_suc2$quant95)))
      
      p <- ggplot(data_sim_suc2, aes(time, median)) + 
        geom_line(size = 1.2, color = my_colors[1]) + 
        geom_line(aes(time, quant05), size = 1.2, color = my_colors[1], alpha = 0.8) +
        geom_line(aes(time, quant95), size = 1.2, color = my_colors[1], alpha = 0.8) +
        geom_ribbon(aes(ymin = quant05, ymax = quant95), color = NA, alpha = 0.25, fill = my_colors[1]) +
        geom_ribbon(aes(ymin = quant20, ymax = quant80), color = NA, alpha = 0.35, fill = my_colors[1]) +
        geom_line(data_sum, mapping=aes(time, median), color = "black", size = 1.2, linetype = 2) +
        geom_line(data_sum, mapping=aes(time, quant95), color = "black", size = 1.2, linetype = 2) +
        geom_line(data_sum, mapping=aes(time, quant05), color = "black", size = 1.2, linetype = 2) +
        geom_rangeframe(data=data_min_max, mapping=aes(x=x, y=y)) +
        labs(x = "Time [min]", y = TeX("YFP intensity \\[A.U $\\times 10^{-2}$\\]")) + 
        ylim(2.0, 15) + 
        my_theme
      
    }else{
      data_plot <- data_agg_sim %>% filter(state == state_name)
      data_min_max <- tibble(x = c(min(data_plot$time), max(data_plot$time)), 
                             y = c(min(data_plot$quant05), max(data_plot$quant95)))
      p <- ggplot(data_plot) + 
        geom_line(aes(time, median), size = 1.2, color = my_colors[1], alpha = 1.0) + 
        geom_line(aes(time, quant05), size = 0.8, color = my_colors[1], alpha = 0.8) +
        geom_line(aes(time, quant95), size = 0.8, color = my_colors[1], alpha = 0.8) +
        geom_ribbon(aes(x = time, ymin = quant05, ymax = quant95), color = NA, alpha = 0.25, fill = my_colors[1]) +
        geom_ribbon(aes(x = time, ymin = quant20, ymax = quant80), color = NA, alpha = 0.35, fill = my_colors[1]) +
        labs(x = "Time [min]", y = str_c("Simulated ", state_name)) +
        geom_rangeframe(data=data_min_max, aes(x=x, y=y)) + 
        my_theme
    }
    path_save <- str_c(dir_save, "Simulated_STS_", state_name, ".pdf")
    ggsave(path_save, plot = p, width = BASE_WIDTH, height = BASE_HEIGHT)
    return(p)
  })
  
  # Save all plots 
  path_save <- str_c(dir_save, "Simulated_all.pdf")
  p <- ggpubr::ggarrange(plotlist = plot_list, ncol = n_states)
  ggsave(path_save, width = 22, height = 6)
  
  return(0)
}


# Function that will plot an IPRED-plot for a STS model 
plot_IPRED <- function(path_to_result, dir_save)
{
  # Make the IPRED plot 
  path_ind <- str_c(path_to_result, "Individual_fits.csv")
  data_ind <- read_csv(path_ind, col_types = cols()) 
  names(data_ind) <- c("time", "id", "ind_fit")
  
  data_obs <- read_csv("../../Intermediate/Data_files/Data_monolix_SUC2.csv",
                       col_types = cols(id = col_factor()))
  data_plot <- tibble(pred = data_ind$ind_fit, obs = data_obs$observation)
  
  p <- ggplot(data_plot, aes(pred, obs)) + 
    geom_point(alpha = 0.6, size = 0.8) +
    geom_smooth(aes(color = "Trend line"), method = "lm", size = 1.0, se = F, show.legend = NA) + 
    geom_abline(aes(slope = 1, intercept = 0, color = "y = x"), size = 1.0, show.legend = F, linetype = 2) + 
    labs(x = "Individual predictions", y = "Observed values") +
    scale_color_manual(name = "Line", values = my_colors[c(1, 5)]) + 
    geom_rangeframe() + 
    my_theme + theme(legend.position = "none")
  
  path_save <- str_c(dir_save, "IPRED_STS.pdf")
  ggsave(path_save, plot = p, height = BASE_HEIGHT, width = BASE_WIDTH)
  
  return(0)
}


# Function that will plot a qq-plot of each estimated parameter by looking at the parameter 
# in the log-space. 
# Args:
#   path_to_result, path to the result directory 
#   dir_save, path to the directory to save the result in 
#   param_save, an argument for a specific parameter to save 
# Returns:
#   void 
plot_qq_param <- function(path_to_result, dir_save, param_save = F, log_space = F)
{
  # Plot distribution tau-feedback 
  path_param <- str_c(path_to_result, "Estimated_param.csv")
  data_param <- read_csv(path_param, col_types = cols())
  
  data_qq_param <- data_param %>%
    select(-id, -a1) 
  
  if(log_space == F){
    data_qq_param <- data_qq_param %>%
    mutate_all(log)
  }
  names_list <- names(data_qq_param)
  plot_list <- lapply(names_list, function(i){
    data_plot <- data_qq_param %>%
      select(i) %>%
      rename("sample" = i)
    ggplot(data_plot, aes(sample=sample)) + 
      geom_qq() +
      geom_qq_line() + 
      labs(title = i) + 
      my_theme})
  

  n_col <- 5
  n_row <- ceiling(length(plot_list) / n_col)
  p <- ggpubr::ggarrange(plotlist = plot_list, ncol = n_col, nrow = n_row)
  file_save <- str_c(dir_save, "qq_param.pdf")
  ggsave(file_save, plot = p, width = 12, height = 7)

    if(param_save != F){
    data_par <- data_param %>%
      select(param_save) %>%
      rename("obs" = param_save)
    
    if(log_space == F){
      data_par <- data_par %>%
        mutate_all(log) 
    }
    
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


# The version with none fixed parameters
path_to_result <- "../../Intermediate/STS/simple_feedback_model_log/LN_BOBYQA/"
dir_save <- "../../Result/STS/Simple_feedback_log/"
# Sanity check that results exist 
if(!dir.exists(path_to_result)) return(1)
r1 <- plot_simulated_data(path_to_result, dir_save)
r2 <- plot_IPRED(path_to_result, dir_save)
r3 <- plot_qq_param(path_to_result, dir_save, param_save = "k3", log_space = T)
if(r1 + r2 + r3 != 0) return(1)


# The version with none fixed parameters
path_to_result <- "../../Intermediate/STS/simple_feedback_model_log_fixed/LN_BOBYQA/"
dir_save <- "../../Result/STS/Simple_feedback_log_fixed/"
# Sanity check that results exist 
if(!dir.exists(path_to_result)) return(1)
r1 <- plot_simulated_data(path_to_result, dir_save)
r2 <- plot_IPRED(path_to_result, dir_save)
r3 <- plot_qq_param(path_to_result, dir_save, param_save = "k6", log_space = T)
if(r1 + r2 + r3 != 0) return(1)
