# Home range estimation via the amt package, though it uses ctmm for aKDE estimation
pacman::p_load(amt, tidyverse, sf, lubridate)
blra1 <- "78076161"
blra_full <- read_rds("data_derived/localized_beep_data.rds")[[blra1]] %>% 
  # Use data from APril with complete node grid
  filter(DetTimeR >= ymd_hms("2022-04-01 00:00:00"),
         DetTimeR < ymd_hms("2022-05-01 00:00:00"),
         # Use only observations every quarter hour to reduce data set size for example
         minute(DetTimeR) %in% c(0, 15, 30, 45)) %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))
blra_red <- read_rds("data_derived/blra_locs_reduced.rds")[[blra1]] %>% 
  # Use data from APril with complete node grid
  filter(DetTimeR >= ymd_hms("2022-04-01 00:00:00"),
         DetTimeR < ymd_hms("2022-05-01 00:00:00"),
         # Use only observations every quarter hour to reduce data set size for example
         minute(DetTimeR) %in% c(0, 15, 30, 45)) %>%
  make_track(X_est, Y_est, DetTimeR, id = TagId, crs = 32617) %>%
  # Nest so we can get different HR estimates in tidy format
  nest(data = c(x_, y_, t_))

# Takes quite a while due to 19K+ points ()
system.time({
  blra_hr <- blra_full %>%
    mutate(
      hr_mcp = map(data, hr_mcp),
      hr_kde = map(data, hr_kde),
      hr_locoh = map(data, ~ hr_locoh(., n = ceiling(sqrt(nrow(.))))))
      # hr_akde = map(data, ~ hr_akde(., fit_ctmm(., "auto"))))
})
saveRDS(blra_hr, "data/blra_hr_estimates.rds")

# @knitr hr-one-animal
mcp1 <- hr_mcp(blra, levels = c(0.5, 0.95))
kde1 <- hr_kde(blra, levels = c(0.5, 0.95))


# @knitr plothr
plot(kde1)
plot(mcp1, add.relocations = FALSE, add = TRUE, border = "red")

# @knitr hr-one-animal-area
hr_area(mcp1)
hr_area(kde1)

# @knitr hr-one-animal-iso
hr_isopleths(mcp1)

# @knitr hr-one-animal-overlap1
hr_overlap(mcp1, kde1)


# @knitr hr-one-animal-overlap2
hr_overlap(kde1, mcp1)

# @knitr nest-1
dat1 <- dat %>% 
  nest(data = c(x_, y_, t_)) 

# @knitr nest-2
dat1

# @knitr nest-3
dat1$data[[1]]

# @knitr hr-many-1.1
hr1 <- dat1 %>%
  mutate(
    hr_mcp = map(data, hr_mcp), 
    hr_kde = map(data, hr_kde), 
    hr_locoh = map(data, ~ hr_locoh(., n = ceiling(sqrt(nrow(.))))), 
    hr_akde_iid = map(data, ~ hr_akde(., fit_ctmm(., "iid"))),
    hr_akde_ou = map(data, ~ hr_akde(., fit_ctmm(., "ou")))
  )

# @knitr hr-many-1.2
str(hr1, 2)

# @knitr hr-many-1.3
hr2 <- hr1 %>% select(-data) %>% 
  pivot_longer(hr_mcp:hr_akde_ou, names_to = "estimator", 
               values_to = "hr")

# @knitr hr-many-1.4
str(hr2, 2)

# @knitr hr-many-1.5
hr2.area <- hr2 %>% 
  mutate(hr_area = map(hr, hr_area)) %>% 
  unnest(cols = hr_area) 

# @knitr hr-many-1.6
head(hr2.area, 2)

# @knitr hrssex
library(patchwork)

hr2.area1 <- hr2.area %>% 
  mutate(sex = str_sub(id, 1, 1), 
         area = area / 1e4) 
hr2.area1$est_lab <- factor(hr2.area1$estimator, 
                            levels = c("hr_mcp", "hr_kde", "hr_locoh", "hr_akde_iid", "hr_akde_ou"), 
                            labels = c("MCP", "KDE", "LoCoH", "aKDE (iid)", "aKDE (ou)"))
ci <- hr2.area1 %>% group_by(est_lab, sex) %>% 
  summarise(m = mean(area), se = sd(area) / sqrt(n()), 
            me = qt(0.975, n() - 1) * se, lci = m - me, uci = m + me)

p1 <- hr2.area1 %>% ggplot(aes(est_lab, area, col = sex)) + 
  geom_jitter(alpha = 0.5, position = position_dodge2(width = 0.5)) +
  geom_pointrange(aes(x = est_lab, y = m, ymin = lci, ymax = uci, col = sex), data = ci, inherit.aes = FALSE, 
                  position = position_dodge2(width = 0.5)) +
  theme_light() + 
  theme(axis.title.x = element_blank()) +
  labs(y = expression(paste("HRS [", ha, "]")))

p2 <- ci %>% nest(data = -est_lab) %>% 
  mutate(f = map_dbl(data, ~ .$m[.$sex == "M"] - .$m[.$sex == "F"] )) %>% 
  ggplot(aes(est_lab, f)) + geom_point() + 
  theme_light() +
  theme(axis.title.x = element_blank()) +
  labs(y = expression(HRS[m] - HRS[f]))

p3 <- ci %>% nest(data = -est_lab) %>% 
  mutate(f = map_dbl(data, ~ .$m[.$sex == "M"] / .$m[.$sex == "F"] )) %>% 
  ggplot(aes(est_lab, f)) + geom_point() + geom_hline(yintercept = 1, lty = 2) +
  theme_light() +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.9, 2.9)) +
  labs(y = expression(HRS[m] / HRS[f]))

(p1 / p2 / p3) + plot_layout(heights = c(0.65, 0.3, 0.3)) + plot_annotation(tag_levels = "A", tag_suffix = ")")

# @knitr hr-many-2.1
env <- raster("data/forest.tif")
hr1.env <- hr1 %>% select(-data) %>% 
  pivot_longer(hr_mcp:hr_akde_ou, names_to = "estimator", 
               values_to = "hr") %>% 
  mutate(forest = map(hr, ~ raster::extract(env, hr_isopleths(.))))


# @knitr hr-many-2.2
hr1.env1 <- hr1.env %>% mutate(
  prop_forest = map_dbl(forest, ~ mean(unlist(.))), 
  area = map(hr, hr_area)) %>% 
  select(id, estimator, prop_forest, area) %>% 
  unnest(cols = area) 


# @knitr hrsforest
hr1.env1 %>% 
  mutate(area = area / 1e4, 
         est_lab = factor(estimator, 
                          levels = c("hr_mcp", "hr_kde", "hr_locoh", "hr_akde_iid", "hr_akde_ou"), 
                          labels = c("MCP", "KDE", "LoCoH", "aKDE (iid)", "aKDE (ou)"))) %>% 
  ggplot(aes(prop_forest, area, col = est_lab)) + geom_point() +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_light() + 
  theme(legend.position = "top") +
  labs(col = "Estimator", y = "Home-range size (ha)", x = "Proportion forest")


# @knitr hr-many-3.1
hr2 <- dat %>%  
  mutate(yday = yday(t_)) %>% 
  nest(data = x_:t_) %>% 
  mutate(n = map_int(data, nrow)) %>% 
  filter(n > 10) %>% 
  mutate(
    hr_mcp = map(data, hr_mcp), 
    hr_kde = map(data, hr_kde),
    hr_locoh = map(data, ~ hr_locoh(., n = sqrt(nrow(.)))),
    hr_od_iid = map(data, ~ hr_od(., model = fit_ctmm(., "iid"))),
    hr_od_ou = map(data, ~ hr_od(., model = fit_ctmm(., "ou")))
  )

# @knitr hrsday
hr2 %>% select(-n, -data) %>% 
  pivot_longer(-c(id, yday, sex)) %>% 
  mutate(area = map(value, hr_area)) %>% select(-value) %>% 
  unnest(area) %>% 
  mutate(area = area / 1e4,
         est_lab = factor(
           name, 
           levels = c("hr_mcp", "hr_kde", "hr_locoh", "hr_od_iid", "hr_od_ou"), 
           labels = c("MCP", "KDE", "LoCoH", "OD (iid)", "OD (ou)"))) %>% 
  ggplot(aes(yday, area, col = est_lab)) + 
  geom_point(alpha = 0.6) + geom_smooth(se = FALSE) + 
  scale_y_continuous(trans = "log10") +
  facet_wrap(~ id, scale = "free") +
  theme_light() +
  theme(legend.position = "top") +
  labs(col = "Estimator", y = "Home-range size (ha)", x = "year day")

