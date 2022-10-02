# Read in calibration setup metadata
processed_beep_path <- "data_derived/processed_beep_data.rds"

if (!new_beep_data) {
  processed_beep_dat <- readRDS(processed_beep_path)
} else {
  tag_deps <- get_tag_deployments()
  node_deps <- get_node_deployments()
  processed_beep_dat <- beep_data %>%
    # # Only keep data from 
    # filter(NodeId %in% viable_nodes) %>%
    # Get rid of duplicate detections on multiple SS radios
    select(DetTime = Time, TagId:NodeId) %>% distinct() %>%
    right_join(tag_deps, by = "TagId") %>%
    right_join(sf::st_drop_geometry(node_deps), by = c("TagId", "NodeId", "BaseStation")) %>%
    filter(between(DetTime, NodeStartDT, NodeEndDT)) %>%
    # Add calibration information
    left_join(node_corrections, by = "NodeId") %>%
    left_join(tag_corrections, by = "TagId") %>%
    # Round time to nearest three minutes
    mutate(DetTimeR = round_date(DetTime, unit = "3 minutes"),
           # Incorporate calibration corrections to RSSI
           adjTagRSSI = TagRSSI + node_rssi_gam_adj - tag_rssi_gam_adj) %>%
    # Group and calculate average TagRSSI and sample size
    group_by(TagId, SppCode, Age, Sex, DeployDT, NodeId, X, Y, DetTimeR) %>%
    summarize(n = n(),
              adjTagRSSI_mn = mean(adjTagRSSI),
              adjTagRSSI_sd = sd(adjTagRSSI),
              adjTagRSSI_cv = adjTagRSSI_sd/(abs(adjTagRSSI_mn)),
              join_rssi = paste0("RSSI_", round(adjTagRSSI_mn, 1)),
              .groups = "drop")
  
  hist(processed_beep_dat$n)
  hist(processed_beep_dat$adjTagRSSI_mn)
  
  # Filter to detections with average (adjusted) tag RSSI >= -90
  processed_beep_dat <- processed_beep_dat %>%
    filter(adjTagRSSI_mn >= -90) %>%
    arrange(TagId, DetTimeR, NodeId)
  print(ggplot(processed_beep_dat, aes(x = adjTagRSSI_mn)) + 
    geom_histogram(binwidth = 1, color = "black", fill = NA))
  
  # Join with estimated distance from RSSI vs distance model
  processed_beep_dat <- left_join(processed_beep_dat, pred_dist_sd, by = "join_rssi") %>%
    select(-join_rssi)
  saveRDS(processed_beep_dat, processed_beep_path)
}
