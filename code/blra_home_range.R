# Home range estimation via the amt package, though it uses ctmm for aKDE estimation
pacman::p_load(amt, tidyr, sf)
blra_id <- "78076161"
blra_full <- read_rds("data_derived/localized_beep_data.rds")[[blra_id]] %>% 
  # Use data from end of March into April with complete node grid while bird was alive
  filter(DetTimeR >= ymd_hms("2022-03-31 00:00:00"),
         DetTimeR < ymd_hms("2022-04-14 00:00:00"),
         converged, # Only use estimated locations that produced a solution
         # Use only observations every half hour to reduce data set size for example
         minute(DetTimeR) %in% c(0, 15, 30, 45)) %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))
blra_red <- read_rds("data_derived/blra_locs_reduced.rds")[[blra_id]] %>% 
  # Use data from end of March into April with complete node grid while bird was alive
  filter(DetTimeR >= ymd_hms("2022-03-31 00:00:00"),
         DetTimeR < ymd_hms("2022-04-14 00:00:00"),
         converged,
         # Use only observations every quarter hour to reduce data set size for example
         minute(DetTimeR) %in% c(0, 15, 30, 45)) %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))

recalc_hrs <- FALSE
if (recalc_hrs) {
  # Takes a bit, even for half-hourly points (about 8 minutes per scenario---full grid, reduced grid)
  system.time({
    blra_hr_full <- blra_full %>%
      mutate(
        hr_mcp = map(data, hr_mcp, levels = c(0.5, 0.95)),
        hr_kde = map(data, hr_kde, levels = c(0.5, 0.95)),
        hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95)))
      )
  })
  blra_hr_full <- blra_hr_full %>%
    select(-data) %>% 
    pivot_longer(hr_mcp:hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_full, "data_derived/blra_hr_full_estimates.rds")
  
  system.time({
    blra_hr_red <- blra_red %>%
      mutate(
        hr_mcp = map(data, hr_mcp, levels = c(0.5, 0.95)),
        hr_kde = map(data, hr_kde, levels = c(0.5, 0.95)),
        hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95)))
      )
  })
  blra_hr_red <- blra_hr_red %>%
    select(-data) %>% 
    pivot_longer(hr_mcp:hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_red, "data_derived/blra_hr_reduced_estimates.rds")
} else {
  blra_hr_full <- readRDS("data_derived/blra_hr_full_estimates.rds")
  blra_hr_red <- readRDS("data_derived/blra_hr_reduced_estimates.rds")
}

# With home ranges estimated, here are some example functionality that may be of use
# Both `blra_hr_full` and `blra_hr_red` have a column "hr" that contains the home range estimations
# They are in the following order: MCP, KDE, and aKDE
# To access them, for now, use the following format: object$hr[[index of estimator of interest]]
# For example, to access the KDE using the full node grid, you would use: blra_hr_full$hr[[2]]

# Plotting
plot(blra_hr_full$hr[[1]]) # MCP, full grid
plot(blra_hr_full$hr[[2]]) # KDE, full grid
plot(blra_hr_full$hr[[3]]) # aKDE, full grid

plot(blra_hr_red$hr[[1]]) # MCP, reduced grid
plot(blra_hr_red$hr[[2]]) # KDE, reduced grid
plot(blra_hr_red$hr[[3]]) # aKDE, reduced grid

# Area of estimated home range
hr_area(blra_hr_full$hr[[3]]) # aKDE, full grid
hr_area(blra_hr_red$hr[[3]]) # aKDE, reduced grid

# Home range overlap
hr_overlap(blra_hr_full$hr[[3]], blra_hr_red$hr[[3]]) #aKDE, full vs. reduced grid overlap

