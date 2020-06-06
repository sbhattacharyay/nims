library(readxl)
library(tidyverse)

sensor_time_info<- read_xlsx('../accel_sensor_data/accel_time_info.xlsx') %>% dplyr::select(4:24)
hours_durations <- sensor_time_info %>% dplyr::select(ends_with("hours"))
sensor_start_times <- sensor_time_info %>% dplyr::select(ends_with("start"))
sensor_end_times <- sensor_time_info %>% dplyr::select(ends_with("end"))
max_duration_idx<- max.col(hours_durations,ties.method="first")

start_times <- c()
end_times <- c()

for (i in 1:nrow(sensor_time_info)) {
  ref_start_time <- sensor_time_info[[i,(max_duration_idx[i]-1)*3+1]]
  ref_end_time <- sensor_time_info[[i,(max_duration_idx[i]-1)*3+2]]
  
  
  
  start_times <- c(start_times,ref_start_time)
  end_times <- c(end_times,ref_end_time)
}