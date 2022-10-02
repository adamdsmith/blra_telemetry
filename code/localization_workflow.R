##################################################
# LOAD REQUIRED PACKAGES 
##################################################
if (!requireNamespace("pacman", quietly = TRUE)) 
  install.packages("pacman", quiet = TRUE)
pacman::p_load(pbapply, DBI, RPostgres, dplyr, lubridate, 
               ggplot2, mgcv, lme4, geosphere, parallel, googledrive)

##################################################
# LOAD PROJECT-SPECIFIC FUNCTIONALITY
##################################################
## Modified CTT functions
source("code/functions/CTT/api_postgres.R")
source("code/functions/CTT/data_manager.R")

## Our functions
source("code/functions/get_calibrations.R")
source("code/functions/get_tag_deployments.R")
source("code/functions/get_node_deployments.R")

##################################################
# DOWNLOAD MOST RECENT DATA FROM CTT SERVERS
##################################################
new_beep_data <- FALSE
if (new_beep_data) source("code/functions/get_data_from_CTT.R")

##################################################
# READ IN TELEMETRY DATA ASSOCIATED WITH CALIBRATIONS
##################################################
# Multiple SensorStations were used in final calibrations
cal_sensorstations <- c("68E8BC02DD51", "4F627F78EC34")

if (new_beep_data) {
  beep_data <- lapply(cal_sensorstations, function(ss) {
    infile <- file.path("data_raw/GA Black Rail", ss)
    # Read in only the telemetry ("beep") data
    tmp_beep <- load_beep_data(infile) #start_time, end_time, tags
  })
  beep_data <- bind_rows(beep_data)
  saveRDS(beep_data, "data_derived/beep_data.rds")
} else beep_data <- readRDS("data_derived/beep_data.rds")

##################################################
# Before we can estimate the parameters of the RF signal loss as a function
# of distance from nodes *in situ*, we have to put all telemetry detections
# on equal footing by correcting for possible differences in (1) received
# signal strength by individual nodes and (2) emitted signal strength by
# individual tags
##################################################

##################################################
# CALCULATE NODE RSSI ADJUSTMENTS FROM CALIBRATIONS
##################################################
update_node_cal <- FALSE
# Create or load `node_corrections` object
source("code/calc_node_calibrations.R")

##################################################
# CALCULATE TAG RSSI ADJUSTMENTS FROM CALIBRATIONS
##################################################
# Create or load `tag_corrections` object
source("code/calc_tag_calibrations.R")

##################################################
# CALCULATE TAG RSSI VS DISTANCE FUNCTION
##################################################
# Create or load `rssi_dist_m_mn` object
# This is a statistical model (class `merMod`)
# from which we will estimate distance from nodes
# prior to multilateration
source("code/calc_rssi_distance_calibration.R")

##################################################
# Process beep data by deployments
##################################################
# Create or load `processed_beep_dat` object
# This filters beep data based on tag and node
# deployment information and prepares the data
# for multilateration (localization) on 3-minute
# intervals, including adding predicted detection
# distance and its associated uncertainty from the 
# RSSI vs distance calibration
source("code/process_beep_data.R")

##################################################
# Multilaterate individual positions
##################################################
# Create or load `bird_loc_dat` object
# This uses multilateration to estimate an 
# individual's location (in 3-minute intervals) 
# based on detections by at least 3 nodes.
# Distance estimates from a node are weighted
# by the inverse of the estimated uncertainty of
# the distance estimate (at least currently)
source("code/functions/multilaterate_detections.R")
# tag_est_locs <- multilaterate_beep_dat(processed_beep_dat,
#                                        Tags = c("781E664C", "78334B2D"), 
#                                        min_node_beeps = 4, min_n_nodes = 4)

# Estimate locations using full node grids
if (new_beep_data) {
  blra <- multilaterate_beep_dat(processed_beep_dat,
                                 Tags = c("78076161"), 
                                 min_node_beeps = 4, min_n_nodes = 3)
  saveRDS(blra, "data_derived/localized_beep_data.rds")
  gd_file <- googledrive::drive_get("localized_beep_data.rds")
  googledrive::drive_update(gd_file, "data_derived/localized_beep_data.rds")
  
  # Here we estimate locations for our single example BLRA using the reduced node grid
  blra_red_nodes <- get_localization_validation_nodes() %>%
    filter(reduced_grid) %>%
    pull(NodeId)
  blra_ex <- filter(processed_beep_dat, 
                    TagId == "78076161",
                    NodeId %in% blra_red_nodes)
  blra_red <- multilaterate_beep_dat(blra_ex,
                                     Tags = c("78076161"), 
                                     min_node_beeps = 4, min_n_nodes = 3)
  saveRDS(blra_red, "data_derived/blra_locs_reduced.rds")
} else {
  blra <- readRDS("data_derived/localized_beep_data.rds")
  blra_red <- readRDS("data_derived/blra_locs_reduced.rds")
}