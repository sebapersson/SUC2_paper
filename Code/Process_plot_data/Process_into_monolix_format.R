library(tidyverse)
library(readxl)

# This file will process the data raw-data into Monolix friendly format. This means that 
# the two data-sets will be processed into a csv-file.
# Output: 
# (1) CSV-file with individual cells and columns: 
#   id, the unique cell id
#   time, the time-point for the measurment 
#   observation, observed value, here SUC2 intensity.

# Note, cells were removed as they showed unrealstic behaviour. This behaviour can 
# be attributed to cells dying, lumping togehter into clusters or due to the 
# microscope going out of focus. 

# -------------------------------------------------------------------------------------------
# Processing experiments into monolix-format and mean data 
# -------------------------------------------------------------------------------------------
suc2_data1 <- read_csv("../../Data/Suc2_0d1Glc_1.csv", col_types = cols()) %>% 
  select(-X1) %>% 
  rename("X1" = "X1_1")
suc2_data2 <- read_csv("../../Data/Suc2_0d1Glc_2.csv", col_types = cols()) %>% 
  select(-X1) %>% 
  rename("X1" = "X1_1")

# Filter the most obvious crazy observations 
data1_filtered <- suc2_data1 %>%
  select(-matches("^[a-zA-Z]+9$")) %>%
  select(-matches("^[a-zA-Z]+6$")) %>%
  select(-matches("^[a-zA-Z]+30$")) %>%
  select(-matches("^[a-zA-Z]+3$")) %>%
  select(-matches("^[a-zA-Z]+70$")) %>%
  select(-matches("^[a-zA-Z]+31$")) %>%
  select(-matches("^[a-zA-Z]+4$")) %>%
  select(-matches("^[a-zA-Z]+26$")) %>%
  select(-matches("^[a-zA-Z]+29$")) 
data2_filtered <- suc2_data2 %>%
  select(-matches("^[a-zA-Z]+46$")) %>%
  select(-matches("^[a-zA-Z]+45$")) %>%
  select(-matches("^[a-zA-Z]+17$")) %>%
  select(-matches("^[a-zA-Z]+20$")) %>%
  select(-matches("^[a-zA-Z]+88$")) %>%
  select(-matches("^[a-zA-Z]+130$")) %>%
  select(-matches("^[a-zA-Z]+16$")) %>%
  select(-matches("^[a-zA-Z]+15$")) %>%
  select(-matches("^[a-zA-Z]+25$")) %>%
  select(-matches("^[a-zA-Z]+47$")) %>%
  select(-matches("^[a-zA-Z]+83$")) %>%
  select(-matches("^[a-zA-Z]+89$")) 

# Make data tidy and remove cells that were noticed to be 
# outliers (crazy) when the data had been imported into 
# monolix. 
data1_tidy_tmp <- data1_filtered %>%
  select(starts_with("Mean")) %>% 
  mutate(t = seq(from = 0, by = 5, length.out = 97)) %>%
  select(t, everything()) %>%
  gather(starts_with("Mean"), value = "Mean", key = "Cell_index") 
data1_tidy_area_tmp <- data1_filtered %>%
  select(starts_with("Area")) %>% 
  gather(starts_with("Area"), value = "Area", key = "Cell_index_area") 
n_cells <- as.integer(dim(data1_tidy_tmp)[1] / 97)
data1_tidy <- data1_tidy_tmp %>%
  bind_cols(data1_tidy_area_tmp) %>%
  select(-"Cell_index_area") %>%
  mutate(Cell = as.factor(rep(1:n_cells, each = 97))) %>%
  mutate(Intensity = Mean / 100) 

data2_tidy_tmp <- data2_filtered %>%
  select(starts_with("Mean")) %>% 
  mutate(t = seq(from = 0, by = 5, length.out = 97)) %>%
  select(t, everything()) %>%
  gather(starts_with("Mean"), value = "Mean", key = "Cell_index") 
data2_tidy_area_tmp <- data2_filtered %>%
  select(starts_with("Area")) %>% 
  gather(starts_with("Area"), value = "Area", key = "Cell_index_area") 
n_cells <- as.integer(dim(data2_tidy_tmp)[1] / 97)
data2_tidy <- data2_tidy_tmp %>%
  bind_cols(data2_tidy_area_tmp) %>%
  select(-"Cell_index_area") %>%
  mutate(Cell = as.factor(rep(1:n_cells, each = 97))) %>%
  mutate(Intensity = Mean / 100) 

data1_tidy_filtered <- data1_tidy %>%
  filter(Cell != 7) %>%
  filter(Cell != 16) %>%
  filter(Cell != 21) %>%
  filter(Cell != 38) %>%
  filter(Cell != 44) %>%
  filter(Cell != 49) %>%
  filter(Cell != 54) %>%
  filter(Cell != 59) 
data2_tidy_filtered <- data2_tidy %>%
  filter(Cell != 3) %>%
  filter(Cell != 4)%>%
  filter(Cell != 5) %>%
  filter(Cell != 7) %>%
  filter(Cell != 11) %>%
  filter(Cell != 12) %>%
  filter(Cell != 15) %>%
  filter(Cell != 16) %>%
  filter(Cell != 17) %>%
  filter(Cell != 18) %>%
  filter(Cell != 24) %>%
  filter(Cell != 28) %>%
  filter(Cell != 29) %>%
  filter(Cell != 30) %>%
  filter(Cell != 31) %>%
  filter(Cell != 32) %>%
  filter(Cell != 33) %>%
  filter(Cell != 39) %>%
  filter(Cell != 40) %>%
  filter(Cell != 41) %>%
  filter(Cell != 48) %>%
  filter(Cell != 53) %>%
  filter(Cell != 56) %>%
  filter(Cell != 60) %>%
  filter(Cell != 65) %>%
  filter(Cell != 68) %>%
  filter(Cell != 73) %>%
  filter(Cell != 74) %>%
  filter(Cell != 75) %>%
  filter(Cell != 76) %>%
  filter(Cell != 77) %>%
  filter(Cell != 78) %>%
  filter(Cell != 79) %>%
  filter(Cell != 81) %>%
  filter(Cell != 82) %>%
  filter(Cell != 84) %>%
  filter(Cell != 86) %>%
  filter(Cell != 90) %>%
  filter(Cell != 95) %>%
  filter(Cell != 98) %>%
  filter(Cell != 101) %>%
  filter(Cell != 108) %>%
  filter(Cell != 110) %>%
  filter(Cell != 113) %>%
  filter(Cell != 114) %>%
  filter(Cell != 116) %>%
  filter(Cell != 117) %>%
  filter(Cell != 64)

# Join the two data sets 
n_cells1 <- dim(data1_tidy_filtered)[1] / 97
n_cells2 <- dim(data2_tidy_filtered)[1] / 97
n_obs1 <- dim(data1_tidy_filtered)[1] 
n_obs2 <- dim(data2_tidy_filtered)[1]
whole_data_filt_tidy <- data1_tidy_filtered %>%
  bind_rows(data2_tidy_filtered) %>%
  mutate(Experiment = as.factor(c(rep(1, n_obs1), rep(2, n_obs2)))) %>%
  mutate(ID = as.factor(rep(1:(n_cells1 + n_cells2), each = 97))) %>%
  rename("Suc2" = "Intensity")

# Reanme columns for monolix 
data_monolix <- whole_data_filt_tidy %>%
  select(ID, t, Suc2) %>%
  rename("id" = "ID", "observation" = "Suc2", "time" = "t")

# Ensure that the directory to save the result in exists 
dir_save <- "../../Intermediate/Data_files/"
if(!dir.exists(dir_save)) dir.create(dir_save)
path_save <- str_c(dir_save, "Data_monolix_SUC2.csv")
write_csv(data_monolix, path_save)
