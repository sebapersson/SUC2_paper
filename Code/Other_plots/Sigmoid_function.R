library(tidyverse)
library(ggthemes)
library(stringr)
RNGkind("L'Ecuyer-CMRG")

# This file creates a plot of the sigmoid function used for modelling the sigmoid switch.  
# Ouputs:
#   A sigmoid function graph to describe how the switch in the complex model is callibrated. 

# General plotting parameters 
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0
my_theme <- theme_tufte(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                                plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_colors <- c("#1f78b4", "#33a02c", "#08519c", "#006d2c", "#ff7f00", "#a6cee3")


# Plot the sigmoid function 
data <- tibble(x = seq(from = 0, to = 10, length.out = 500)) %>%
  mutate(y = 1 / (1 + exp(-3*(x - 4.5 /3.5))))
p <- ggplot(data, aes(x, y)) + 
  geom_line(size = 1.2, color = my_colors[1]) + 
  geom_vline(xintercept = 2.667, linetype = 2) +
  geom_rangeframe() + 
  labs(x = "SNF1p", y = "Sigmoid value") + 
  annotate("text", label = "SNF1p = 2.667", x = 3.75, y = 0.5) + 
  my_theme

dir_save <- "../../Result/Paper_figures/Parts/"
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
path_save <- str_c(dir_save, "Sigmoid_func.pdf")
ggsave(path_save, plot = p, width = BASE_WIDTH, height = BASE_HEIGHT)
