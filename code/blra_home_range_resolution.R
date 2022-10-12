# Home range estimation via the amt package, though it uses ctmm for aKDE estimation
pacman::p_load(amt, tidyr, sf)
blra_id <- "78076161"

# Using localizations/naive estimates from full node grid
blra_full <- read_rds("data_derived/localized_beep_data.rds")[[blra_id]] %>% 
  # Use data from end of March into April with complete node grid while bird was alive
  filter(DetTimeR >= ymd_hms("2022-03-31 00:00:00"),
         DetTimeR < ymd_hms("2022-04-14 00:00:00"),
         converged)

# 15 minute home range
blra_15m <- filter(blra_full,
                   minute(DetTimeR) %in% c(0, 15, 30, 45))

# 1 hour home range
blra_1h <- filter(blra_full,
                  minute(DetTimeR) == 0)

# 4 hour home range
blra_4h <- filter(blra_full,
                  minute(DetTimeR) == 0,
                  hour(DetTimeR) %in% seq(0, 20, by = 4))

# 8 hour home range
blra_8h <- filter(blra_full,
                  minute(DetTimeR) == 0,
                  hour(DetTimeR) %in% seq(0, 16, by = 8))

# Create complete points for all time resolutions
blra_full_loc_sf <- select(blra_full, TagId, DetTimeR, X_est, Y_est) %>%
  mutate(is15m = minute(DetTimeR) %in% c(0, 15, 30, 45),
         is1h = minute(DetTimeR) == 0,
         is4h = minute(DetTimeR) == 0 & hour(DetTimeR) %in% seq(0, 20, by = 4),
         is8h = minute(DetTimeR) == 0 & hour(DetTimeR) %in% seq(0, 16, by = 8)) %>%
  pivot_longer(cols = starts_with("is"),
               names_to = "resolution",
               values_to = "resolution_eval") %>%
  filter(resolution_eval) %>%
  mutate(resolution = factor(resolution, levels = c("is15m", "is1h", "is4h", "is8h"),
                             labels = c("15 minutes", "1 hour", "4 hours", "8 hours"))) %>%
  st_as_sf(coords = c("X_est", "Y_est"), crs = 32617)

blra_tracks_15m <- blra_15m %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))
blra_tracks_1h <- blra_1h %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))
blra_tracks_4h <- blra_4h %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))
blra_tracks_8h <- blra_8h %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))

recalc_hrs <- FALSE
if (recalc_hrs) {
  # Takes a bit
  system.time({
    blra_hr_15m <- blra_tracks_15m %>%
      mutate(hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95))))
  })
  blra_hr_15m <- blra_hr_15m %>%
    select(-data) %>% 
    pivot_longer(hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_15m, "data_derived/blra_hr_15m_estimates.rds")
  
  system.time({
    blra_hr_1h <- blra_tracks_1h %>%
      mutate(hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95))))
  })
  blra_hr_1h <- blra_hr_1h %>%
    select(-data) %>% 
    pivot_longer(hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_1h, "data_derived/blra_hr_1h_estimates.rds")
  
  system.time({
    blra_hr_4h <- blra_tracks_4h %>%
      mutate(hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95))))
  })
  blra_hr_4h <- blra_hr_4h %>%
    select(-data) %>% 
    pivot_longer(hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_4h, "data_derived/blra_hr_4h_estimates.rds")

  system.time({
    blra_hr_8h <- blra_tracks_8h %>%
      mutate(hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"), levels = c(0.5, 0.95))))
  })
  blra_hr_8h <- blra_hr_8h %>%
    select(-data) %>% 
    pivot_longer(hr_akde, names_to = "estimator", 
                 values_to = "hr")
  saveRDS(blra_hr_8h, "data_derived/blra_hr_8h_estimates.rds")
} else {
  blra_hr_15m <- readRDS("data_derived/blra_hr_15m_estimates.rds")
  blra_hr_1h <- readRDS("data_derived/blra_hr_1h_estimates.rds")
  blra_hr_4h <- readRDS("data_derived/blra_hr_4h_estimates.rds")
  blra_hr_8h <- readRDS("data_derived/blra_hr_8h_estimates.rds")
}

# Area of estimated home range
akde_15m_area <- hr_area(blra_hr_15m$hr[[1]]) %>% mutate(resolution = "15 minutes")
akde_1h_area <- hr_area(blra_hr_1h$hr[[1]]) %>% mutate(resolution = "1 hour")
akde_4h_area <- hr_area(blra_hr_4h$hr[[1]]) %>% mutate(resolution = "4 hours")
akde_8h_area <- hr_area(blra_hr_8h$hr[[1]]) %>% mutate(resolution = "8 hours")
akde_areas <- bind_rows(akde_15m_area,
                        akde_1h_area,
                        akde_4h_area,
                        akde_8h_area) %>%
  filter(what == "estimate") %>%
  mutate(area_ha = round(area / 10000, 2),
         resolution = factor(resolution, levels = c("15 minutes", "1 hour", "4 hours", "8 hours")))

# Consolidate HR of interest for plotting
akde_15m <- hr_isopleths(blra_hr_15m$hr[[1]]) %>% mutate(resolution = "15 minutes")
akde_1h <- hr_isopleths(blra_hr_1h$hr[[1]]) %>% mutate(resolution = "1 hour")
akde_4h <- hr_isopleths(blra_hr_4h$hr[[1]]) %>% mutate(resolution = "4 hours")
akde_8h <- hr_isopleths(blra_hr_8h$hr[[1]]) %>% mutate(resolution = "8 hours")
akde_hr_ests <- rbind(akde_15m, akde_1h, akde_4h, akde_8h) %>%
  filter(what == "estimate") %>%
  mutate(resolution = factor(resolution, levels = c("15 minutes", "1 hour", "4 hours", "8 hours")))
  
# Plot aKDE home range comparison, with 50% and 95% estimated areas
p <- ggplot(akde_hr_ests) + geom_sf(aes(fill = factor(level, levels = c(0.5, 0.95), labels = c("50%", "95%")))) + 
  geom_sf(data = blra_full_loc_sf, size = 0.5, alpha = 0.5) +
  scale_fill_manual("aKDE level", values = c("gray", NA), na.value = NA) +
  facet_wrap(~resolution, nrow = 2) + theme_bw() + labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
# Get X and Y limits for labelling
xr <- layer_scales(p)$x$range$range
yr <- layer_scales(p)$y$range$range
akde_areas <- akde_areas %>%
  mutate(x = min(xr) + 0.01 * diff(xr),
         y = max(yr) - 0.01 * ifelse(level == 0.5, 1, 10) * diff(yr),
         label = paste(paste0(level *100, "%:"), area_ha, "ha"))
p + geom_text(data = akde_areas, aes(x, y, label = label),
              hjust = 0, vjust = 1, size = 3)
ggsave("output/figures/akde_hr_resolution_comparisons.png", height = 6.5, width = 6.5)
