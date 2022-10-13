# NOTE: Requires `code/localization_workflow` to have been run first 
source("code/localization_workflow.R")
pacman::p_load(MASS, DHARMa)

loc_val_path <- "data_derived/localization_validation.rds"
loc_est_path <- "data_derived/localization_validation_estimates.rds"

# Load some functions for retrieving 
source("code/functions/get_calibrations.R")
loc_val_nodes <- get_localization_validation_nodes()

if (!new_beep_data) {
  loc_val_data <- readRDS(loc_val_path)
} else {
  loc_val_meta <- get_localization_validation()
   # Join beep data with validation metadata
  loc_val_data <- beep_data %>%
    filter(between(Time, min(loc_val_meta$ValStartDT), max(loc_val_meta$ValEndDT))) %>%
    # Get rid of duplicate detections on multiple SS radios
    select(DetTime = Time, TagId:NodeId) %>% distinct() %>%
    right_join(loc_val_meta, by = "TagId") %>%
    filter(between(DetTime, ValStartDT, ValEndDT)) %>%
    # Add calibration information
    left_join(node_corrections, by = "NodeId") %>%
    left_join(tag_corrections, by = "TagId") %>%
    # Incorporate calibration corrections to RSSI
    mutate(adjTagRSSI = TagRSSI + node_rssi_gam_adj - tag_rssi_gam_adj) %>%
    # Group and calculate average TagRSSI and sample size
    group_by(TagId, Val_pt, Val_X, Val_Y, NodeId) %>%
    summarize(n = n(),
              adjTagRSSI_mn = mean(adjTagRSSI),
              adjTagRSSI_sd = sd(adjTagRSSI),
              adjTagRSSI_cv = adjTagRSSI_sd/(abs(adjTagRSSI_mn)),
              join_rssi = paste0("RSSI_", round(adjTagRSSI_mn, 1)),
              .groups = "drop") %>%
    # Filter to detections with average (adjusted) tag RSSI >= -90
    filter(adjTagRSSI_mn >= -90) %>%
    # Join with estimated distance from RSSI vs distance model
    left_join(pred_dist_sd, by = "join_rssi") %>%
    select(-join_rssi)   %>%
    left_join(loc_val_nodes, by = "NodeId")
  
  saveRDS(loc_val_data, loc_val_path)
}

if (!new_beep_data) {
  out_loc_ests <- readRDS(loc_est_path)
} else {
  naive_ests <- bind_rows(
    # Full grid
    loc_val_data %>% group_by(Val_pt) %>% slice(which.max(adjTagRSSI_mn)) %>%
      select(Val_pt, X_naive = X, Y_naive = Y) %>% ungroup() %>%
      mutate(grid_config = "full"),
    loc_val_data %>% filter(reduced_grid) %>% group_by(Val_pt) %>% slice(which.max(adjTagRSSI_mn)) %>%
      select(Val_pt, X_naive = X, Y_naive = Y) %>% ungroup() %>%
      mutate(grid_config = "reduced"))

  # Center X and Y to ease computation
  x_mn <- mean(loc_val_nodes$X)
  y_mn <- mean(loc_val_nodes$Y)
  
  # Base list of all validation points with naive estimated locations
  all_val_pts <- loc_val_data %>%
    mutate(dist_grid_center = sqrt((Val_X - x_mn)^2 + (Val_Y - y_mn)^2)) %>%
    select(TagId, Val_pt, Val_X, Val_Y, Val_dist_to_center = dist_grid_center) %>% distinct() %>%
    left_join(naive_ests, by = "Val_pt")
  
  # Specify some minimum detecting node requirements to use localization estimate
  # Note these needn't apply to naive estimates; by def, just need one node detection to estimate location
  min_node_beeps <- 4
  min_n_nodes <- 3
  
  # For location estimation, reduce to detections from required minimum number of nodes
  tmp_beep <- filter(loc_val_data, n >= min_node_beeps) %>%
    group_by(TagId, Val_pt) %>%
    mutate(Xc = X - x_mn,
           Yc = Y - y_mn) 
  tmp_beep_list <- group_split(tmp_beep)
  
  loc_ests <- pbapply::pblapply(tmp_beep_list, function(dat) {
    ests <- lapply(c("full", "reduced"), function(grid_config) {
      if (grid_config == "reduced") dat <- filter(dat, reduced_grid)
      n_nodes <- nrow(dat)
      if (n_nodes >= min_n_nodes) {
        Xstart <- dat$Xc[which.max(dat$adjTagRSSI_mn)]
        Ystart <- dat$Yc[which.max(dat$adjTagRSSI_mn)]
        loc_mod = nls(dist_pred ~ EucDist(cbind(Xc, Yc), cbind(x_est, y_est)), data = dat,
                      start = list(x_est = Xstart, y_est = Ystart),
                      # Weight longer range detections less, by inverse of estimated distance
                      weights = 1/dist_pred,
                      control=nls.control(warnOnly = T, minFactor=1/30000, maxiter = 100))
        loc_coefs <- coef(loc_mod)
        loc_vcov <- vcov(loc_mod)
        loc_sd <- sqrt(diag(loc_vcov))
        X_est = loc_coefs[1] + x_mn
        Y_est = loc_coefs[2] + y_mn
        X_est_sd = loc_sd[1]
        Y_est_sd = loc_sd[2]  
        df = stats::df.residual(loc_mod)
        converged = summary(loc_mod)$convInfo$isConv
        varcov = list(loc_vcov)
      } else {
        X_est = Y_est = X_est_sd = Y_est_sd = df = NA
        varcov = list()
        converged = FALSE
      }
      mod_summary <- tibble::tibble(Val_pt = unique(dat$Val_pt), 
                                    n_nodes = n_nodes,
                                    X_est, Y_est, X_est_sd, Y_est_sd,
                                    RSSI_sd_wt = stats::weighted.mean(dat$adjTagRSSI_sd, dat$n),
                                    RSSI_cv_wt = stats::weighted.mean(dat$adjTagRSSI_cv, dat$n),
                                    df, converged, varcov,
                                    grid_config = grid_config)
    })
    ests <- bind_rows(ests)
    ests
  })
  loc_ests <- bind_rows(loc_ests)
  
  out_loc_ests <- all_val_pts %>%
    left_join(loc_ests, by = c("Val_pt", "grid_config")) %>%
    replace_na(list(converged = FALSE)) %>%
    mutate(xy_dist = distEuc(cbind(Val_X, Val_Y), cbind(X_est, Y_est)),
           xy_dist_naive = distEuc(cbind(Val_X, Val_Y), cbind(X_naive, Y_naive)),
           within_grid = case_when(between(Val_X, min(loc_val_nodes$X), max(loc_val_nodes$X)) & 
                                     between(Val_Y, min(loc_val_nodes$Y), max(loc_val_nodes$Y)) ~ TRUE, 
                                   TRUE ~ FALSE),
           max_loc_sd = pmax(X_est_sd, Y_est_sd)) %>%
    arrange(Val_pt)
  saveRDS(out_loc_ests, loc_est_path)
}

# Compare convergence of estimation model by full vs. reduced grid configuration
with(out_loc_ests, table(converged, grid_config))

# Use only converged estimates hereafter
all_converged_ests <- filter(out_loc_ests, converged)
(est_medians <- out_loc_ests %>%
  group_by(grid_config) %>%
  summarize(med_xy_dist = median(xy_dist, na.rm = TRUE),
            med_xy_dist_naive = median(xy_dist_naive)))

# Histograms of absolute error for estimated location of validation points (w/median absolute error)
ggplot(all_converged_ests, aes(x = xy_dist)) +
  geom_histogram(binwidth = 2, colour="black", fill="gray") +
  geom_vline(data = est_medians, aes(xintercept = med_xy_dist), color = "red", linetype = "dashed", size = 1) + 
  scale_x_continuous("Absolute error (m)", limits = c(0, 100)) +
  geom_text(data = est_medians, aes(label = paste("Median:", round(med_xy_dist, 1), "m")), 
            x = 100, y = 4, vjust = 1, hjust = 1) +
  labs(y = "# validation points") + facet_wrap(~grid_config, ncol = 1) +
  ggtitle("Multilateration error of validation locations") +
  theme_bw()
ggsave("output/figures/absolute_estimation_error_by_grid_configuration.png", width = 6.5, height = 4.5)

# Histograms of absolute error for naive location (i.e., location of node w/strongest RSSI) of validation points (w/median absolute error)
ggplot(out_loc_ests, aes(x = xy_dist_naive)) +
  geom_histogram(binwidth = 2, colour="black", fill="gray") +
  geom_vline(data = est_medians, aes(xintercept = med_xy_dist_naive), color = "red", linetype = "dashed", size = 1) + 
  scale_x_continuous("Absolute error (m)", limits = c(0, 100)) +
  geom_text(data = est_medians, aes(label = paste("Median:", round(med_xy_dist_naive, 1), "m")), 
            x = 100, y = 5, vjust = 1, hjust = 1) +
  labs(y = "# validation points") + facet_wrap(~grid_config, ncol = 1) +
  ggtitle("Naive error of validation locations") +
  theme_bw()
ggsave("output/figures/absolute_naive_error_by_grid_configuration.png", width = 6.5, height = 4.5)

# Relationship between number of nodes used in location estimation and absolute error of that estimate
## First, fit a linear model evaluating the relationship and comparing among grid arrangments
# We use the negative binomial model for overdispersion; we rounded error estimates to nearest meter
# We exclude localizations outside the grid boundaries as they're less stable
mod_dat <- mutate(all_converged_ests, xy_dist = round(xy_dist)) %>%
  filter(within_grid)
nnode_mod <- glm(xy_dist ~ grid_config + n_nodes + grid_config * n_nodes, family = poisson, data = mod_dat)
sim_out <- simulateResiduals(nnode_mod)
plot(sim_out) # Overdispersed related to Poisson; try NB

nnode_mod <- glm.nb(xy_dist ~ grid_config + n_nodes + grid_config * n_nodes, data = mod_dat)
sim_out <- simulateResiduals(nnode_mod)
plot(sim_out) # Much better
# Check if interaction is supported? Nope
summary(nnode_mod) # Drop interaction, z = -1.44, p = 0.15
# Refit without interaction
nnode_mod <-glm.nb(xy_dist ~ grid_config + n_nodes, data = mod_dat)
summary(nnode_mod) # Little evidence of difference between full and reduced grid (z = -1.16, p = 0.25)
(1 - exp(coef(nnode_mod)[3])) # How much does each extra node reduce estimation error, on average
# Fit # nodes only model for plotting
nnode_mod <-glm.nb(xy_dist ~ n_nodes, data = mod_dat)
new_dat <- data.frame(n_nodes = 3:16)
preds <- predict(nnode_mod, newdata = new_dat, se.fit = TRUE)
new_dat$fit_link <- preds$fit; new_dat$se_link <- preds$se.fit
new_dat <- mutate(new_dat,
                  fit = exp(fit_link),
                  lcl = exp(fit_link - 2 * se_link),
                  ucl = exp(fit_link + 2 * se_link))

ggplot(all_converged_ests, aes(x = n_nodes)) +
  geom_ribbon(data = new_dat, aes(ymin = lcl, ymax = ucl), alpha = 0.5) +
  geom_line(data = new_dat, aes(y = fit)) + 
  geom_point(aes(y = xy_dist, fill = grid_config), shape = 21, size = 3, alpha = 0.5) + 
  labs(x = "# detecting nodes", y = "Absolute error (m)", fill = "Grid") + 
  ggtitle("Multilateration error vs. number of estimating nodes") +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
ggsave("output/figures/absolute_estimation_error_by_number_detecting_nodes.png", width = 5, height = 5)

# Relationship between distance from node grid center and absolute error of that estimate (localization)
ggplot(all_converged_ests, aes(x = Val_dist_to_center, y = log(xy_dist))) +
  geom_smooth(method = "lm") + 
  geom_point(aes(fill = within_grid), shape = 21, size = 3, alpha = 0.5) + 
  labs(x = "Distance to node grid center (m)", y = "log (Absolute error) (m)", color = "Point\nin grid") + 
  facet_wrap(~ grid_config, ncol = 2) + 
  ggtitle("Multilateration error vs. distance from center of node grid") +
  theme_bw()
ggsave("output/figures/multilateration_estimation_error_by_distance_to_grid_center.png", width = 6.5, height = 4.5)

# Relationship between distance from node grid center and absolute error of that estimate (naive)
ggplot(out_loc_ests, aes(x = Val_dist_to_center, y = log(xy_dist_naive))) +
  geom_smooth(method = "lm") + 
  geom_point(aes(fill = within_grid), shape = 21, size = 3, alpha = 0.5) + 
  labs(x = "Distance to node grid center (m)", y = "log (absolute error) (m)", color = "Point\nin grid") + 
  facet_wrap(~ grid_config, ncol = 2) + 
  ggtitle("Naive error vs. distance from center of node grid") +
  theme_bw()
ggsave("output/figures/naive_estimation_error_by_distance_to_grid_center.png", width = 6.5, height = 4.5)

# Estimated absolute errors plotted by location within node grid
## Need to first expand node location data.frame for use in faceted plot
exp_node_locs <- bind_rows(
  mutate(loc_val_nodes, grid_config = "full"),
  mutate(filter(loc_val_nodes, reduced_grid), grid_config = "reduced")
)

## Then plot
ggplot(all_converged_ests, aes(Val_X, Val_Y)) + 
  geom_point(aes(color = xy_dist, size = n_nodes)) + 
  geom_point(data = exp_node_locs, aes(X, Y)) +
  scale_color_gradient(low = "#1a9641", high = "#d7191c") +
  labs(x = NULL, y = NULL, size = "# nodes", color = "Absolute\nerror (m)") +
  facet_wrap(~grid_config) +
  theme_bw() +
  theme(axis.text = element_blank())
ggsave("output/figures/absolute_estimation_error_by_location_in_grid.png", width = 6.5, height = 4.5)

# Absolute error of naive estimates (i.e., node w/strongest RSSI) by location in grid
ggplot(out_loc_ests, aes(Val_X, Val_Y)) + 
  geom_point(aes(color = xy_dist_naive), size = 4) + 
  geom_point(data = exp_node_locs, aes(X, Y)) +
  scale_color_gradient(low = "#1a9641", high = "#d7191c") +
  labs(x = NULL, y = NULL, color = "Absolute\nerror (m)") +
  facet_wrap(~grid_config) +
  theme_bw() +
  theme(axis.text = element_blank())
ggsave("output/figures/absolute_naive_error_by_location_in_grid.png", width = 6.5, height = 4.5)

# Relationship between localization standard deviation and absolute error of estimated location
# NOTE: This was just out of curiosity; doubt we use
ggplot(all_converged_ests, aes(x = max_loc_sd, y = xy_dist)) +
  # geom_smooth(method = "lm", formula = log(y) ~ x, se = FALSE) + 
  geom_point() +
  xlab("Localization SD") +
  ylab("Absolute error (m)") + 
  facet_wrap(~ grid_config) +
  theme_bw()
