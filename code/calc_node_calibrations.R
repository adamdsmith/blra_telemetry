# Read in calibration setup metadata
node_correction_path <- "data_derived/Calibration/node_corrections.rds"

if (!update_node_cal) {
  node_corrections <- readRDS(node_correction_path)
} else {
  node_cal_data <- get_node_calibration() 
  node_cal_data <- filter(node_cal_data, Use_in_calibration)
  
  node_cal_nodes <- pull(node_cal_data, NodeId) %>% unique()
  
  attr(node_cal_data$CalStartDT, "tzone") <- "UTC"
  attr(node_cal_data$CalEndDT, "tzone") <- "UTC"
  
  node_cal_beep_data <- beep_data %>%
    filter(!is.na(NodeId),
           NodeId %in% node_cal_nodes) %>%
    select(DetTime = Time, RadioId:NodeId)
  
  # Join beep data with calibration metadata
  node_cal_beeps <- left_join(node_cal_data, node_cal_beep_data, by = c("NodeId", "TagId")) %>%
    filter(between(DetTime, CalStartDT, CalEndDT))
  
  # Node calibration summary 
  ## Expected calibration sets
  (exp_node_cal <- group_by(node_cal_data, NodeId, NodePos, CalStartDT) %>% distinct())
  
  ## Actual calibration sets
  (act_node_cal <- group_by(node_cal_beeps, NodeId, CalStartDT) %>% tally(sort = TRUE))
  
  missing_node_cal <- exp_node_cal %>% anti_join(act_node_cal, by = c("NodeId", "CalStartDT"))
  if (nrow(missing_node_cal) > 0) {message("Missing node calibrations"); print(missing_node_cal)}
  
  # Calculate the node correction summary (mean and sd)
  message("Calculating naive node corrections (not accounting for position during calibration)")
  naive_node_corrections <- node_cal_beeps %>%
    # Use correction from largest sample of first trial (if replicated for position effect)
    group_by(NodeId, NodePos, RadioId, TagId, Date, CalStartDT, CalEndDT) %>%
    summarize(n = n(),
              TagRSSI_mean = mean(TagRSSI),
              TagRSSI_sd = sd(TagRSSI),
              .groups = "drop") %>%
    group_by(NodeId, CalStartDT) %>%
    arrange(NodeId, Date, desc(n), RadioId) %>% 
    slice(1) %>% ungroup() %>%
    mutate(mean_correction = TagRSSI_mean - mean(TagRSSI_mean)) %>% 
    arrange(desc(TagRSSI_sd))
  
  # with(naive_node_corrections, hist(TagRSSI_sd, breaks = 50))
  with(naive_node_corrections, hist(TagRSSI_mean, breaks = 20))

  # Control for potential position effects on RSSI
  node_cal_beeps <- right_join(node_cal_beeps, select(naive_node_corrections, NodeId, RadioId, n, TagRSSI_mean),
                               by = c("NodeId", "RadioId")) %>%
    mutate(NodeId = factor(NodeId),
           TagRSSI_cent = TagRSSI - mean(naive_node_corrections$TagRSSI_mean))
  message("Fitting GAM to control for node position around test tag during calibration.")
  node_cal_gam <- gam(TagRSSI_cent ~ NodeId + s(NodePos, bs = "cc") - 1,
                      method = "ML", data = node_cal_beeps)
  
  # Retrieve Node-level deviations (fixed effect estimates)
  node_coefs <- grep("NodeId", names(coef(node_cal_gam)))
  node_cal_gam.fe <- coef(node_cal_gam)[node_coefs]
  node_corrections <- tibble(NodeId = levels(node_cal_beeps$NodeId),
                             node_rssi_gam_adj = unname(node_cal_gam.fe))
  message("Saving node calibration corrections to '", node_correction_path, ".'")
  saveRDS(node_corrections, file = node_correction_path)
  with(left_join(node_corrections, naive_node_corrections), plot(mean_correction, node_rssi_gam_adj,
                              xlab = "Naive node correction (db)",
                              ylab = "Position-correction node correction (db)"))
}