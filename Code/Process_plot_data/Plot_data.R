library(tidyverse) 
library(stringr)
library(latex2exp)
library(ggthemes)

# The aim with this file is to plot the SUC2-data. The data most have been processed into 
# Monolix-format. Generally, any csv-file with an id-column (named id), time-column (named time)
# and observation (named observation) can be plotted by this file. The file will also plot 
# a dummay example, which is used to produce Fig. 2 in the paper. 
# (1) A plot with:
#   The 0.05, 0.5 and 0.95 quantiles of the observed data 
#   The observed traces (small grey font)
#   Four cell traces choosen randomly 
# (2) A dummy example plot illustrating the data collection process 

# General plotting parameters for gg-plot, theme_tufte is used!
my_theme <- theme_tufte(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                       plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_colors <- c("#1f78b4", "#33a02c", "#08519c", "#006d2c")
BASE_HEIGHT <- 5; BASE_WIDTH <- 7.0


# Function that plots the SUC2-data, where the data should follow the specifications 
# above.
# Args:
#   path_data, relative path to the data set 
#   seed, seed for randomly choosing cells 
# Returns:
#   void 
plot_monolix_data <- function(path_data, seed)
{
  set.seed(seed)
  data <- read_csv(path_data, col_types = cols(id = col_factor()))
  
  # Create directory where to save the data-set 
  dir_save <- "../../Result/Data_set/"
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  path_save <- str_c(dir_save, "Data_set.pdf")
  path_save_png <- str_c(dir_save, "Data_set.png")
  
  # The summarised data for the plot
  data_sum <- data %>%
    group_by(time) %>%
    summarise(median = median(observation, na.rm = T), 
              quant05 = quantile(observation, 0.05, na.rm = T), 
              quant95 = quantile(observation, 0.95, na.rm = T))
  cells_to_plot <- sample(unique(data$id), 4, replace = F)
  data_plot <- data %>% filter(id %in% cells_to_plot)
  
  # Colors for the plot, and range-frame param 
  my_grey <- "#969696"
  my_col <- rep(my_colors[1], 4)  
  data_frame <- tibble(x = c(min(data$time), max(data$time)), 
                       y = c(min(data$observation), max(data$observation)))
  
  p <- ggplot(data_plot, aes(time, observation)) + 
    geom_line(data=data, mapping=aes(time, observation, group=id), size = 0.15, color = my_grey) + 
    geom_line(aes(color = id), size = 1.2) +
    #geom_point(aes(color = id), size = 1.2) +
    scale_color_manual(values = my_col) +
    geom_rangeframe(data=data_frame, mapping=aes(x, y)) +
    geom_line(data_sum, mapping=aes(time, quant05), 
              size = 1.5, linetype = 2) +
    geom_line(data_sum, mapping=aes(time, quant95), 
              linetype = 2, size = 1.5) +
    geom_line(data_sum, mapping=aes(time, median), size = 1.5, linetype = 2) +
    labs(x = "Time [min]", y = TeX("SUC2 intensity \\[A.U.$\\times 10^{-2}$\\]")) +
    my_theme +
    theme(legend.position = "none")
  
  
  ggsave(path_save, plot = p, height = BASE_HEIGHT-0.5, width = BASE_WIDTH)  
  ggsave(path_save_png, plot = p, height = BASE_HEIGHT-0.5, width = BASE_WIDTH)  
  
  return(0)
}


# Function that will plot a dummy-example just to highlight how the data collection
# actually works. 
plot_example_data <- function()
{
  dir_save <- "../../Result/Data_set/"
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  path_save <- str_c(dir_save, "Data_collection_ex.pdf")
  
  x_data <- seq(from = 0, to = 450, by = 50)
  y_data <- c(390, 400, 405, 420, 820, 1100, 700, 480, 450, 455)
  data_plot <- tibble(x = x_data, y = y_data / 100)
  n <- 10
  p <- ggplot(data_plot, aes(x, y)) + 
    geom_point(size = 3.0, color = my_colors[1]) +
    geom_line(data=data.frame(spline(data_plot, n=n*10)), size = 1.2, color = my_colors[1]) +
    geom_rangeframe() + 
    labs(x = "Time [min]", y = "") +
    my_theme 

  
  ggsave(path_save, plot = p, width = 8, height = BASE_HEIGHT-3)
  
  return(0)
}

path_data <- "../../Intermediate/Data_files/Data_monolix_SUC2.csv"
if(!file.exists(path_data)) return(1)
r1 <- plot_monolix_data(path_data, 1)
r2 <- plot_example_data()

if(r1 + r2 != 0) return(1)
