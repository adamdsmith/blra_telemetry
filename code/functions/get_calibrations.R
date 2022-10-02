get_tag_calibration <- function() {
  
  if (!requireNamespace("googlesheets4", quietly = TRUE))
    pacman::p_load("googlesheets4")
  
  out_fn <- file.path("data_raw", "tag_calibration.rds")
  
  # Googlesheets access code for archive only, restricted access
  gs_title <- "CTT PowerTag Calibration"
  
  # Set up column types for call # (or call type for BLRA) columns
  col_types <- "cccccciiccl"
  
  out <- googledrive::drive_get(gs_title) %>%
    googlesheets4::read_sheet(col_types = col_types) %>%
    mutate(CalStartDT = ymd_hm(paste(Date, stringr::str_pad(StartTime, 4, pad = "0")), tz = "America/New_York"),
           CalEndDT = ymd_hm(paste(Date, stringr::str_pad(EndTime, 4, pad = "0")), tz = "America/New_York"),
           NodeId1 = toupper(NodeId_Pos1),
           NodeId2 = toupper(NodeId_Pos2),
           NodeId3 = toupper(NodeId_Pos3),
           NodeId4 = toupper(NodeId_Pos4),
           TagId = toupper(TagId)) %>%
    select(TagId, NodeId1:NodeId4, CalStartDT:CalEndDT, Use_in_calibration) %>%
    pivot_longer(cols = NodeId1:NodeId4,
                 names_to = NULL,
                 values_to = "NodeId")
  saveRDS(out, file = out_fn)
  out
}

get_node_calibration <- function() {
  
  if (!requireNamespace("googlesheets4", quietly = TRUE))
    pacman::p_load("googlesheets4")
  
  out_fn <- file.path("data_raw", "node_calibration.rds")
  
  # Googlesheets access code for archive only, restricted access
  gs_title <- "CTT Node Calibration"
  
  # Set up column types for call # (or call type for BLRA) columns
  col_types <- "ccicciiccl"
  
  out <- googledrive::drive_get(gs_title) %>%
    googlesheets4::read_sheet(col_types = col_types) %>%
    mutate(CalStartDT = ymd_hm(paste(Date, stringr::str_pad(StartTime, 4, pad = "0")), 
                               tz = "America/New_York"),
           CalEndDT = ymd_hm(paste(Date, stringr::str_pad(EndTime, 4, pad = "0")), 
                             tz = "America/New_York"),
           NodeId = toupper(NodeId),
           TagId = toupper(TagId)) %>%
    select(NodeId:Date, CalStartDT:CalEndDT, Use_in_calibration)
  saveRDS(out, file = out_fn)
  out
}

get_rssi_calibration <- function() {
  
  if (!requireNamespace("googlesheets4", quietly = TRUE))
    pacman::p_load("googlesheets4")
  
  out_fn <- file.path("data_raw", "rssi_distance_calibration.rds")
  
  # Googlesheets access code for archive only, restricted access
  gs_title <- "CTT PowerTag RSSI vs Distance Calibration"
  
  # Set up column types for call # (or call type for BLRA) columns
  col_types <- "ccccciiiiccl"
  
  out <- googledrive::drive_get(gs_title) %>%
    googlesheets4::read_sheet(col_types = col_types)  %>%
    mutate(CalStartDT = ymd_hm(paste(Date, stringr::str_pad(StartTime, 4, pad = "0")), 
                               tz = "America/New_York"),
           CalEndDT = ymd_hm(paste(Date, stringr::str_pad(EndTime, 4, pad = "0")), 
                             tz = "America/New_York"),
           NodeId1 = toupper(NodeId_Pos1),
           NodeId2 = toupper(NodeId_Pos2),
           Dist_Node2 = Internode_dist - Dist_Node1,
           TagId = toupper(TagId)) %>%
    dplyr::select(Trial, TagId, NodeId1:NodeId2, Dist_Node1, Dist_Node2, 
                  CalStartDT:CalEndDT, Use_in_calibration)
  
  # Partition transect setup by beginning and end nodes and combine
  cal_node1 <- dplyr::select(out, -NodeId2, -Dist_Node2) %>%
    rename(NodeId = NodeId1, Dist_node = Dist_Node1)
  cal_node2 <- dplyr::select(out, -NodeId1, -Dist_Node1) %>%
    rename(NodeId = NodeId2, Dist_node = Dist_Node2)
  out <- bind_rows(cal_node1, cal_node2)
  attr(out$CalStartDT, "tzone") <- "UTC"
  attr(out$CalEndDT, "tzone") <- "UTC"
  
  saveRDS(out, file = out_fn)
  out
}

get_localization_validation <- function() {
  
  if (!requireNamespace("googlesheets4", quietly = TRUE))
    pacman::p_load("googlesheets4")
  
  out_fn <- file.path("data_raw", "localization_validation.rds")

  # Googlesheets access code for archive only, restricted access
  gs_title <- "Localization Validation"
  
  # Set up column types for call # (or call type for BLRA) columns
  col_types <- "cciiciiclc"
  
  out <- googledrive::drive_get(gs_title) %>%
    googlesheets4::read_sheet(col_types = col_types) %>%
    filter(Use_in_validation) %>%
    mutate(ValStartDT = ymd_hm(paste(Date, stringr::str_pad(StartTime, 4, pad = "0")), 
                               tz = "America/New_York"),
           ValEndDT = ymd_hm(paste(Date, stringr::str_pad(EndTime, 4, pad = "0")), 
                             tz = "America/New_York"),
           TagId = toupper(TagId)) %>%
    dplyr::select(TagId, Val_pt, Val_X = UTME_17N, Val_Y = UTMN_17N, ValStartDT:ValEndDT)
    
  attr(out$ValStartDT, "tzone") <- attr(out$ValEndDT, "tzone") <- "UTC"
  
  saveRDS(out, file = out_fn)
  out
}

get_localization_validation_nodes <- function() {
  
  if (!requireNamespace("googlesheets4", quietly = TRUE))
    pacman::p_load("googlesheets4")
  
  # Googlesheets access code for archive only, restricted access
  gs_title <- "Localization Validation Nodes"
  
  # Set up column types for call # (or call type for BLRA) columns
  col_types <- "cciilc"
  
  out <- googledrive::drive_get(gs_title) %>%
    googlesheets4::read_sheet(col_types = col_types) %>%
    dplyr::select(NodeId, X = UTME_17N, Y = UTMN_17N, reduced_grid)
  out
}