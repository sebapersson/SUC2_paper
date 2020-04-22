library(tidyverse)
library(ggthemes)
library(stringr)


# This file produces the plots for the invertase experimental data. To keep a consistent 
# graphical profile through the paper, all results are plotted (including non model results 
# are plotted in R). 
# Outputs:
#   A bar-chart showing invertase activity in high and low glucose 


# General plotting parameters 
my_theme <- theme_tufte(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                               plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_colors <- c("#1f78b4", "#33a02c", "#08519c", "#006d2c", "#ff7f00", "#a6cee3")
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0


# Function that processes the invertase data and outputs the result in a bar-chart 
# Args:
#   void 
# Returns:
#   void 
process_invertase_data <- function()
{
  data_values <- tibble(value = c(158.1394, 714.0953, 33.55298, 31.47122), 
                        sd = c(24.76474, 8.004206629, 19.09195, 8.117907953), 
                        type = as.factor(c("wt", "wt", "dsnf", "dsnf")), 
                        glc = as.factor(c("high", "low", "high", "low"))) %>%
  mutate(frame = c(0, 0, 0, max(value)))
  
  p <- ggplot(data_values, aes(type, value, fill = glc)) + 
    geom_bar(stat='identity', position='dodge', color ="black") +
    geom_errorbar(aes(ymin = value-sd, ymax = value+sd), width=0.2, position=position_dodge(.9)) + 
    scale_x_discrete(limit = c("wt", "dsnf"), 
                     labels = c("wt", "Δsnf1")) +
    geom_rangeframe(aes(type, frame), sides = "l") + 
    scale_fill_manual(values = my_colors[c(1, 6)]) +
    labs(y = "Invertase activity [mU/mg]", x = "") + 
    my_theme + theme(legend.position = "none", 
                     axis.title.x = element_blank(), 
                     axis.text.x = element_blank())
  
  dir_save <- "../../Result/Experimental_data/"
  if(!dir.exists(dir_save)) dir.create(dir_save)
  path_save <- str_c(dir_save, "Invertase_activity.pdf")
  ggsave(path_save, plot = p, width = BASE_WIDTH, height = BASE_HEIGHT)
  
  return(0)
}


process_invertase_data()
