library(tidyverse)

# This file will process the monolix-data into a mean-data set, where the 
# single cell data is aggreated into a population mean value. This mean-value 
# is later used to find starting values for the algorithms. 
# Output:
# (1) CSV-file with mean-value and columns:
#   time, the time-value
#   mean, the mean SUC2 intensity 

path_data <- "../../Intermediate/Data_files/Data_monolix_SUC2.csv"
if(!file.exists(path_data)) return(1)
data_monolix <- read_csv(path_data, col_types = cols())
data_mean <- data_monolix %>%
  group_by(time) %>%
  summarise(mean = mean(observation, na.rm = T))
dir_save <- "../../Intermediate/Data_files/"
if(!dir.exists(dir_save)) dir.create(dir_save)
path_save <- str_c(dir_save, "Data_mean_SUC2.csv")
write_csv(data_mean, path_save)
