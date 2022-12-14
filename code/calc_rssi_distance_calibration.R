pacman::p_load(qpcR, nlraa)

# Read in calibration setup metadata
rssi_distance_path <- "data_derived/rssi_distance_calibration.rds"
rssi_distance_sds_path <- "data_derived/rssi_distance_estimate_sds.rds"

if (!update_node_cal) {
  rssi_dist_m_mn <- readRDS(rssi_distance_path)
  pred_dist_sd <- readRDS(rssi_distance_sds_path)
} else {
  rssi_distance_data <- get_rssi_calibration() %>%
    filter(Use_in_calibration) %>%
    left_join(node_corrections, by = "NodeId") %>%
    left_join(tag_corrections, by = "TagId")
  rssi_cal_nodes <- pull(rssi_distance_data, NodeId) %>% unique()

  attr(rssi_distance_data$CalStartDT, "tzone") <- "UTC"
  attr(rssi_distance_data$CalEndDT, "tzone") <- "UTC"

  rssi_cal_beep_data <- beep_data %>%
    filter(!is.na(NodeId),
           NodeId %in% rssi_cal_nodes) %>%
    dplyr::select(DetTime = Time, RadioId:NodeId)
  
  # Join beep data with calibration metadata
  rssi_cal_beeps <- left_join(rssi_distance_data, rssi_cal_beep_data, by = c("NodeId", "TagId")) %>%
    filter(between(DetTime, CalStartDT, CalEndDT))  
  
  rssi_keep_radios <- rssi_cal_beeps %>%
    group_by(Trial, TagId, NodeId, Dist_node, CalStartDT, CalEndDT, RadioId) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(Trial, TagId, NodeId, Dist_node, CalStartDT, CalEndDT) %>%
    arrange(desc(n), RadioId) %>% 
    slice(1) %>% ungroup() %>% 
    dplyr::select(Trial:RadioId)
  
  rssi_cal_beeps <- left_join(rssi_keep_radios, rssi_cal_beeps,
                              by = c("Trial", "TagId", "NodeId", "Dist_node", 
                                     "CalStartDT", "CalEndDT", "RadioId"))
  
  rssi_cal_beeps_i <- rssi_cal_beeps %>%
    mutate(Trial = factor(Trial),
           # Correct TagRSSI based on node and tag calibrations
           adjTagRSSI = TagRSSI + node_rssi_gam_adj - tag_rssi_gam_adj,
           logDist = log10(Dist_node + 1)) %>%
    dplyr::select(Trial:CalEndDT, logDist, TagRSSI, adjTagRSSI)

  ggplot(rssi_cal_beeps_i, aes(-logDist, adjTagRSSI, color = Trial)) + 
    geom_point(alpha = 0.5) + geom_smooth(method = "lm") +
    facet_wrap(~ Trial)

  rssi_cal_beeps_mn <- rssi_cal_beeps %>%
    mutate(Trial = factor(Trial),
           adjTagRSSI = TagRSSI + node_rssi_gam_adj - tag_rssi_gam_adj) %>%
    group_by(Trial, TagId, NodeId, Dist_node, CalStartDT, CalEndDT) %>%
    summarize(adjTagRSSI_mn = mean(adjTagRSSI), 
              min_elapsed = as.numeric(difftime(max(DetTime), min(DetTime), units = "mins")),
              n = n(), .groups = "drop") %>%
    mutate(logDist = log10(Dist_node + 1))
  
  ggplot(rssi_cal_beeps_mn, aes(-logDist, adjTagRSSI_mn, color = Trial)) + 
    geom_point(alpha = 0.5) + geom_smooth(method = "lm") +
    facet_wrap(~ Trial)
  hist(rssi_cal_beeps_mn$n)
  
  # RF power drops with inverse square of distance
  # RSSI is proportional to log10 (RF Power)
  # RSSI is thus proportional to log10(1/Distance^2)
  # Thus RSSI is proportional to -log10(Distance)
  # So we can estimate the relationship between RSSI and distance in our BLRA system
  # with a linear regression to estimate:
  # RSSI = -K * log10(Distance) + RSSI0
  # where K = 10 * n and n is the signal propagation constant or exponent describing 
  # the particular relationship between RSSI and distance in our system, and 
  # RSSI0 is essentially the RSSI at distance 0 (intercept)
  # With estimates of K and RSSI0, and their uncertainty, we can estimate BLRA distance 
  # from a node (with uncertainty) based on RSSI using:
  #  Distance = 10^((RSSI0 - RSSI) / K)
  
  # First, based on individual pulses
  rssi_dist_m_i <- lmer(adjTagRSSI ~ 1 + logDist + (1 | Trial), data = rssi_cal_beeps_i)
 
  # Second, based on 2-minute mean signal strength
  rssi_dist_m_mn <- lmer(adjTagRSSI_mn ~ 1 + logDist + (1 | Trial), data = rssi_cal_beeps_mn)
  
  # Add residuals
  rssi_cal_beeps_mn$isl_resid <- resid(rssi_dist_m_mn)
  
  # Compare with exponential model of Paxton et al. 2022
  exp_mod_mn <- nls(adjTagRSSI_mn ~ SSasymp(Dist_node, Asym, R0, lrc), data = rssi_cal_beeps_mn)
  exp_mod_coef <- coef(exp_mod_mn)
  # Extract starting values
  a <- exp_mod_coef[["R0"]] - exp_mod_coef[["Asym"]]
  S <- exp(exp_mod_coef[["lrc"]])
  K <- exp_mod_coef[["Asym"]]
  exp_mod_mn <- nls(adjTagRSSI_mn ~ a * exp(-S * Dist_node) + K, 
                      start = list(a = a, S = S, K = K), 
                      data = rssi_cal_beeps_mn)
  # Add residuals
  rssi_cal_beeps_mn$ed_resid <- resid(exp_mod_mn)
  
  # Compare our physical model (inverse squared distance) vs.
  # Paxton et al. exponential decay model
  pred_dat <- data.frame(Dist_node = 0:60) %>% 
    mutate(logDist = log10(Dist_node + 1))
  pred_dat$InvSqLaw <- predict(rssi_dist_m_mn, newdata = pred_dat, re.form = NA)
  pred_dat$ExpDecay <- predict(exp_mod_mn, newdata = pred_dat)
  pred_dat <- pivot_longer(pred_dat,
                           cols = c("InvSqLaw", "ExpDecay"), 
                           names_to = "Model",
                           values_to = "pred") %>%
    mutate(Model = factor(Model, levels = c("InvSqLaw", "ExpDecay"),
                          labels = c("Inverse Square Law", "Exponential Decay")))

  ggplot(data = rssi_cal_beeps_mn, aes(Dist_node, adjTagRSSI_mn)) + 
    geom_line(data = pred_dat, aes(y = pred, color = Model), lwd = 2) +
    geom_point() +
    labs(x = "Distance (m) from node", y = "Calibration-adjusted RSS (2 min mean)") +
    theme_bw() +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c(1, 1),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))
  ggsave("output/figures/RSSI_distance_model_comparison.png", width = 5, height = 5)
  
  resid_dat <- dplyr::select(rssi_cal_beeps_mn, Dist_node, isl_resid, ed_resid) %>%
    pivot_longer(!Dist_node, names_to = "Model", values_to = "resid") %>%
    mutate(Model = factor(Model, levels = c("isl_resid", "ed_resid"),
                          labels = c("Inverse Square Law", "Exponential Decay")))
  ggplot(data = resid_dat, aes(Dist_node, resid)) + 
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_point(aes(color = Model)) +
    labs(x = "Distance (m) from node", y = "Residual RSS") +
    theme_bw() 
  ggsave("output/figures/RSSI_distance_model_residual_comparison.png", width = 6.5, height = 4)
  
  pred_dist <- function(merMod, rssi = -100:-25) {
    calc_D <- function(RSSI, RSSI0, K) 10^((RSSI0 - RSSI)/K)
    coefs <- unname(fixef(merMod))
    k <- coefs[2] / -1; rssi0 <- coefs[1]
    out <- calc_D(rssi, rssi0, k)
    return(out)
  }
  
  boot_i <- t(bootMer(rssi_dist_m_i, pred_dist, 100)$t) %>%
    as.data.frame() %>%
    mutate(rssi = -100:-25) %>%
    pivot_longer(-rssi, names_to = "sim", values_to = "dist")
  ggplot(boot_i, aes(rssi, dist, group = sim)) + geom_line(alpha = 0.1, lwd=2)
  
  boot_mn <- t(bootMer(rssi_dist_m_mn, pred_dist, 100)$t) %>%
    as.data.frame() %>%
    mutate(rssi = -100:-25) %>%
    pivot_longer(-rssi, names_to = "sim", values_to = "dist")
  
  # Not run here, but rerun boot_mn with 1000 iterations and
  # calculate 95% confidence interval for predicted distance at RSS = -90
  # ci <- filter(boot_mn, rssi == -90)
  # quantile(ci$dist, probs = c(0.025, 0.975))
  
  ggplot(boot_mn, aes(rssi, dist)) +
    geom_line(aes(group = sim), alpha = 0.1, lwd=2) +
    labs(x = "Calibration-adjusted RSS (2 min mean)", y = "Estimated distance (m)") +
    theme_bw()
  ggsave("output/figures/RSSI_distance_model_uncertainty.png", width = 5, height = 4)
  
  # Illustrate how predicted distance varies with RSSI
  boot_mn_sum <- boot_mn %>% group_by(rssi) %>%
    summarize(dist_mn = mean(dist),
              dist_sd = sd(dist),
              dist_cv = dist_sd / dist_mn)
  ggplot(boot_mn_sum, aes(rssi, dist_sd)) +
    geom_line() +
    labs(x = "Calibration-adjusted RSS (2 min mean)", 
         y = "Standard deviation of estimated distance (m)") +
    theme_bw()
  ggsave("output/figures/RSSI_distance_model_Sd_vs_RSS.png", width = 5, height = 4)
  
  # RMSE comparisons out of curiosity (packages not loaded by default)
  phys_mod_rmse <- bootMer(rssi_dist_m_mn, qpcR::RMSE, nsim = 500)$t
  exp_mod_rmse <- nlraa::boot_nls(exp_mod_mn, f = qpcR::RMSE, R = 500)$t
  rmse_comparison <- data.frame(Model = rep( c("Inverse Square Law", "Exponential Decay"), each = 500),
                                RMSE = c(phys_mod_rmse, exp_mod_rmse))
  ggplot(rmse_comparison, aes(x = Model, y = RMSE)) + geom_violin(fill = "gray75") +
    theme_bw()
  ggsave("output/figures/RSSI_distance_models_RMSE_comparison.png", width = 5, height = 4)
  
  saveRDS(rssi_dist_m_mn, rssi_distance_path)
  
  # Save lookup table of distance estimate uncertainty by RSSI value
  # Bootstrapped from final model (n = 1000 iterations)
  # nsim = 1000
  out_rssi <- seq(-90, -25, by = 0.1)
  # pred_dist_fin <- function(merMod) pred_dist(merMod, rssi = out_rssi)
  # pred_dist_sd_mn <- t(bootMer(rssi_dist_m_mn, pred_dist_fin, nsim)$t) %>%
  #   as.data.frame() %>%
  #   mutate(rssi = out_rssi,
  #          join_rssi = paste0("RSSI_", round(rssi, 1))) %>%
  #   pivot_longer(starts_with("V"), names_to = "sim", values_to = "dist")
  # pred_dist_sd <- group_by(pred_dist_sd_mn, join_rssi) %>%
  #   summarize(dist_pred = mean(dist),
  #             dist_sd = sd(dist))
  pred_dist_sd <- data.frame(join_rssi = paste0("RSSI_", out_rssi),
                             dist_pred = round(pred_dist(rssi_dist_m_mn, out_rssi), 2))
  saveRDS(pred_dist_sd, rssi_distance_sds_path)
}
