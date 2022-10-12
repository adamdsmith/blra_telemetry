# Home range estimation via the amt package, though it uses ctmm for aKDE estimation
pacman::p_load(amt, tidyr, sf)
blra_id <- "78076161"

# Using localizations/naive estimates from full node grid
blra_full <- read_rds("data_derived/localized_beep_data.rds")[[blra_id]] %>% 
  # Use data from end of March into April with complete node grid while bird was alive
  filter(DetTimeR >= ymd_hms("2022-03-31 00:00:00"),
         DetTimeR < ymd_hms("2022-04-14 00:00:00"),
         converged, # Only use estimated locations that produced a solution
         # Use only observations every half hour to reduce data set size for example
         minute(DetTimeR) %in% c(0, 15, 30, 45)) 
blra_full_loc_sf <- select(blra_full, TagId, DetTimeR, X_est, Y_est) %>%
  st_as_sf(coords = c("X_est", "Y_est"), crs = 32617)
blra_full_naive_sf <- select(blra_full, TagId, DetTimeR, X_naive, Y_naive) %>%
  st_as_sf(coords = c("X_naive", "Y_naive"), crs = 32617)
blra_full_loc_tracks <- blra_full %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))
blra_full_naive_tracks <- blra_full %>%
  make_track(X_naive, Y_naive, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))

# Using localizations from reduced node grid
blra_red <- read_rds("data_derived/blra_locs_reduced.rds")[[blra_id]] %>% 
  # Use data from end of March into April with complete node grid while bird was alive
  filter(DetTimeR >= ymd_hms("2022-03-31 00:00:00"),
         DetTimeR < ymd_hms("2022-04-14 00:00:00"),
         converged,
         # Use only observations every quarter hour to reduce data set size for example
         minute(DetTimeR) %in% c(0, 15, 30, 45))
blra_red_loc_sf <-  select(blra_red, TagId, DetTimeR, X_est, Y_est) %>%
  st_as_sf(coords = c("X_est", "Y_est"), crs = 32617)
blra_red_naive_sf <- select(blra_red, TagId, DetTimeR, X_naive, Y_naive) %>%
  st_as_sf(coords = c("X_naive", "Y_naive"), crs = 32617)
blra_red_loc_tracks <- blra_red %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))
blra_red_naive_tracks <- blra_red %>%
  make_track(X_naive, Y_naive, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))

recalc_hrs <- FALSE
if (recalc_hrs) {
  # Takes a bit, even for quarter hourly points (about 5-8 minutes per scenario---full grid, reduced grid)
  system.time({
    blra_hr_loc_full <- blra_full_loc_tracks %>%
      mutate(
        hr_mcp = map(data, hr_mcp, levels = c(0.5, 0.95)),
        hr_kde = map(data, hr_kde, levels = c(0.5, 0.95)),
        hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95)))
      )
  })
  blra_hr_loc_full <- blra_hr_loc_full %>%
    select(-data) %>% 
    pivot_longer(hr_mcp:hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_loc_full, "data_derived/blra_hr_loc_full_estimates.rds")

  system.time({
    blra_hr_naive_full <- blra_full_naive_tracks %>%
      mutate(
        hr_mcp = map(data, hr_mcp, levels = c(0.5, 0.95)),
        hr_kde = map(data, hr_kde, levels = c(0.5, 0.95)),
        hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95)))
      )
  })
  blra_hr_naive_full <- blra_hr_naive_full %>%
    select(-data) %>% 
    pivot_longer(hr_mcp:hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_naive_full, "data_derived/blra_hr_naive_full_estimates.rds")
  
  system.time({
    blra_hr_loc_red <- blra_red_loc_tracks %>%
      mutate(
        hr_mcp = map(data, hr_mcp, levels = c(0.5, 0.95)),
        hr_kde = map(data, hr_kde, levels = c(0.5, 0.95)),
        hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95)))
      )
  })
  blra_hr_loc_red <- blra_hr_loc_red %>%
    select(-data) %>% 
    pivot_longer(hr_mcp:hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_loc_red, "data_derived/blra_hr_loc_reduced_estimates.rds")

  system.time({
    blra_hr_naive_red <- blra_red_naive_tracks %>%
      mutate(
        hr_mcp = map(data, hr_mcp, levels = c(0.5, 0.95)),
        hr_kde = map(data, hr_kde, levels = c(0.5, 0.95)),
        hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95)))
      )
  })
  blra_hr_naive_red <- blra_hr_naive_red %>%
    select(-data) %>% 
    pivot_longer(hr_mcp:hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_naive_red, "data_derived/blra_hr_naive_reduced_estimates.rds")
} else {
  blra_hr_loc_full <- readRDS("data_derived/blra_hr_loc_full_estimates.rds")
  blra_hr_naive_full <- readRDS("data_derived/blra_hr_naive_full_estimates.rds")
  blra_hr_loc_red <- readRDS("data_derived/blra_hr_loc_reduced_estimates.rds")
  blra_hr_naive_red <- readRDS("data_derived/blra_hr_naive_reduced_estimates.rds")
}

# With home ranges estimated, here are some example functionality that may be of use
# Both `blra_hr_full` and `blra_hr_red` have a column "hr" that contains the home range estimations
# They are in the following order: MCP, KDE, and aKDE
# To access them, for now, use the following format: object$hr[[index of estimator of interest]]
# For example, to access the KDE using the full node grid, you would use: blra_hr_full$hr[[2]]

# Plotting
plot(blra_hr_loc_full$hr[[1]]) # MCP, full grid
plot(blra_hr_loc_full$hr[[2]]) # KDE, full grid
plot(blra_hr_loc_full$hr[[3]]) # aKDE, full grid

plot(blra_hr_naive_full$hr[[1]]) # MCP, full grid, naive estimates
plot(blra_hr_naive_full$hr[[2]]) # KDE, full grid, naive estimates
plot(blra_hr_naive_full$hr[[3]]) # aKDE, full grid, naive estimates

plot(blra_hr_loc_red$hr[[1]]) # MCP, reduced grid
plot(blra_hr_loc_red$hr[[2]]) # KDE, reduced grid
plot(blra_hr_loc_red$hr[[3]]) # aKDE, reduced grid

# Area of estimated home range
akde_loc_full_area <- hr_area(blra_hr_loc_full$hr[[3]]) %>% mutate(model = "Localization\n(full grid)")
akde_naive_full_area <- hr_area(blra_hr_naive_full$hr[[3]]) %>% mutate(model = "Naive\n(full grid)")
akde_loc_red_area <- hr_area(blra_hr_loc_red$hr[[3]]) %>% mutate(model = "Localization\n(reduced grid)")
akde_naive_red_area <- hr_area(blra_hr_naive_red$hr[[3]]) %>% mutate(model = "Naive\n(reduced grid)")
akde_areas <- bind_rows(akde_loc_full_area,
                        akde_naive_full_area,
                        akde_loc_red_area,
                        akde_naive_red_area) %>%
  filter(what == "estimate") %>%
  mutate(area_ha = round(area / 10000, 2))
# Home range overlap
hr_overlap(blra_hr_loc_full$hr[[3]], blra_hr_loc_red$hr[[3]]) #aKDE, full vs. reduced grid overlap
hr_overlap(blra_hr_loc_full$hr[[3]], blra_hr_naive_full$hr[[3]]) #aKDE, full vs. reduced grid overlap

# Consolidate HR of interest for plotting
akde_loc_full <- hr_isopleths(blra_hr_loc_full$hr[[3]]) %>% mutate(model = "Localization\n(full grid)")
akde_naive_full <- hr_isopleths(blra_hr_naive_full$hr[[3]]) %>% mutate(model = "Naive\n(full grid)")
akde_loc_red <- hr_isopleths(blra_hr_loc_red$hr[[3]]) %>% mutate(model = "Localization\n(reduced grid)")
akde_naive_red <- hr_isopleths(blra_hr_naive_red$hr[[3]]) %>% mutate(model = "Naive\n(reduced grid)")
hr_pts <- rbind(mutate(blra_full_loc_sf, model = "Localization\n(full grid)"),
                mutate(blra_full_naive_sf, model = "Naive\n(full grid)"),
                mutate(blra_red_loc_sf, model = "Localization\n(reduced grid)"),
                mutate(blra_red_naive_sf, model = "Naive\n(reduced grid)"))
akde_hr_ests <- rbind(akde_loc_full, akde_naive_full, akde_loc_red, akde_naive_red) %>%
  filter(what == "estimate") 
p <- ggplot(akde_hr_ests) + geom_sf(aes(fill = factor(level, levels = c(0.5, 0.95), labels = c("50%", "95%")))) + 
  geom_sf(data = hr_pts, size = 0.5, alpha = 0.5) +
  scale_fill_manual("aKDE level", values = c("gray", NA), na.value = NA) +
  facet_wrap(~model, nrow = 2) + theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
# Get X and Y limits for labelling
xr <- layer_scales(p)$x$range$range
yr <- layer_scales(p)$y$range$range
akde_areas <- akde_areas %>%
  mutate(x = min(xr) + 0.025 * diff(xr),
         y = max(yr) - 0.025 * ifelse(level == 0.5, 1, 5) * diff(yr),
         label = paste(paste0(level *100, "%:"), area_ha, "ha"))
p + geom_text(data = akde_areas, aes(x, y, label = label),
              xjust = 0, yjust = 1, size = 2)
ggsave("output/figures/akde_hr_comparisons.png", height = 6.5, width = 6.5)
