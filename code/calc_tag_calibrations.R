# Read in calibration setup metadata
tag_correction_path <- "data_derived/Calibration/tag_corrections.rds"

if (!update_node_cal) {
  tag_corrections <- readRDS(tag_correction_path)
} else {
  tag_cal_data <- get_tag_calibration() %>%
    filter(Use_in_calibration) %>%
    left_join(node_corrections, by = "NodeId")
  tag_cal_nodes <- pull(tag_cal_data, NodeId) %>% unique()
  
  attr(tag_cal_data$CalStartDT, "tzone") <- "UTC"
  attr(tag_cal_data$CalEndDT, "tzone") <- "UTC"

  tag_cal_beep_data <- beep_data %>%
    filter(!is.na(NodeId),
           NodeId %in% tag_cal_nodes) %>%
    select(DetTime = Time, RadioId:NodeId)
  
  # Join beep data with calibration metadata
  tag_cal_interval <- 4
  tag_cal_beeps <- left_join(tag_cal_data, tag_cal_beep_data, by = c("NodeId", "TagId")) %>%
    filter(between(DetTime, CalStartDT, CalEndDT)) %>%
    # Here we're identifying tag antenna orientation during the trial
    # because reading are much more consistent when the tag antenna is
    # oriented oblique to the node (rather than pointed toward/away from node)
    mutate(cal_time = as.numeric(difftime(DetTime, CalStartDT, units = "mins")),
           cal_interval = as.integer(cut(cal_time, breaks = seq(0, 4 * tag_cal_interval, by = tag_cal_interval),
                              labels = 1:4,
                              include.lowest = TRUE, right = FALSE)),
           NodeId = factor(NodeId, levels = c("327585", "33D972", "33D7DB", "339C5B")),
           node_cal_order = as.integer(NodeId),
           tag_ant_orient = case_when(
             abs(cal_interval - node_cal_order) == 0 ~ "towards",
             abs(cal_interval - node_cal_order) == 2 ~ "away",
             (cal_interval - node_cal_order) %in% c(-3, 1) ~ "right",
             (cal_interval - node_cal_order) %in% c(-1, 3) ~ "left",
             TRUE ~ "ERROR"))
  
  # Plot of calibration progression for diagnostics...
  ggplot(tag_cal_beeps, aes(cal_time, TagRSSI, color = tag_ant_orient)) +
    geom_point() + 
    scale_color_viridis_d("Antenna\norientation") +
    facet_grid(TagId ~ NodeId)
  
  # A few calibration detections with bad timing by tag orientiation to be removed
  # NOTE: Only applies to left/right orientation as they're only ones used in calibration
  bad_timing <- tibble(TagId = c("662A7852", "662A7852", "52782A19", "334C0755", "334C4C19"),
                       NodeId = factor(c("327585", "33D7DB", "33D7DB", "327585", "327585"),
                                       levels = levels(tag_cal_beeps$NodeId)),
                       DetTime = ymd_hms(c("2022-02-22 16:33:00", "2022-02-22 16:33:00", 
                                           "2022-03-02 13:51:00", "2022-03-02 14:24:00",
                                           "2022-03-03 15:25:59"),
                                         tz = "UTC"))

  tag_cal_beeps <- anti_join(tag_cal_beeps, bad_timing)
  
  # Tag calibration summary 
  ## Expected calibration sets
  (exp_tag_cal <- select(tag_cal_data, TagId, CalStartDT) %>% distinct())
  
  ## Actual calibration sets
  (act_tag_cal <- group_by(tag_cal_beeps, TagId, CalStartDT) %>% tally(sort = TRUE))
  
  missing_tag_cal <- exp_tag_cal %>% anti_join(act_tag_cal, by = c("TagId", "CalStartDT"))
  if (nrow(missing_tag_cal) > 0) {message("Missing tag calibrations"); print(missing_tag_cal)}

  keep_radios <- tag_cal_beeps %>%
    # Drop calibrations with tag antenna point towards/away from node
    # Data more inconsistent
    filter(tag_ant_orient %in% c("left", "right")) %>%
    group_by(TagId, NodeId, CalStartDT, CalEndDT, RadioId, tag_ant_orient) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(desc(n), RadioId) %>%
    group_by(TagId, NodeId, CalStartDT, CalEndDT, tag_ant_orient) %>%
    slice(1) %>%
    ungroup() %>%
    select(TagId, NodeId, RadioId, tag_ant_orient)
  
  # Reduce to dataset to calculate tag corrections
  tag_cal_beeps <- right_join(tag_cal_beeps, keep_radios,
                              by = c("TagId", "NodeId", "RadioId", "tag_ant_orient")) %>%
    # Correction for node occurs here...
    # Also add in fixed approximate intercept so add'l tags calibrated similarly later
    mutate(adjTagRSSI = TagRSSI - node_rssi_gam_adj - (-30),
           TagId = factor(TagId),
           tag_ant_orient = factor(tag_ant_orient)) 
  tag_cal_gam <- gam(adjTagRSSI ~ TagId + tag_ant_orient - 1 + s(NodeId, bs = "re"),
                     data = tag_cal_beeps, methods = "ML")
  
  # Retrieve tag-level deviations (fixed effect estimates)
  tag_coefs <- grep("TagId", names(coef(tag_cal_gam)))
  tag_cal_gam.fe <- coef(tag_cal_gam)[tag_coefs]
  
  tag_corrections <- tibble(TagId = levels(tag_cal_beeps$TagId),
                            tag_rssi_gam_adj = unname(coef(tag_cal_gam)[tag_coefs]))
  mean(tag_corrections$tag_rssi_gam_adj); sd(tag_corrections$tag_rssi_gam_adj)

  ggplot(tag_corrections, aes(tag_rssi_gam_adj)) +
    geom_histogram(binwidth = 1, color = "black", fill = "gray50") +
    labs(x = "RSS correction", y = NULL) + theme_bw()
  ggsave("output/figures/tag-specific_correction_histogram.png", width = 4.5, height = 4.5)
  message("Saving tag calibration corrections to '", tag_correction_path, ".'")
  saveRDS(tag_corrections, file = tag_correction_path)
}