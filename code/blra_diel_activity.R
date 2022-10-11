# Diel activity
blra_id <- "78076161"
blra_act <- read_rds("data_derived/localized_beep_data.rds")[[blra_id]]
attr(blra_act$DetTimeR, "tzone") <- "America/New_York"
blra_act <- blra_act %>% 
  # Use data from end of March into April with complete node grid while bird was alive
  filter(DetTimeR >= ymd_hms("2022-03-31 00:00:00", tz = "America/New_York"),
         DetTimeR < ymd_hms("2022-04-14 00:00:00", tz = "America/New_York")) %>%
  mutate(doy = yday(DetTimeR),
         hour = hour(DetTimeR))
# Put on correct timezone to may sunrise/sunset join easier
blra_act <- mutate(blra_act,
                   date_str = as.character(as_date(DetTimeR)))

sun <- nrsmisc::get_sun(-80.90303, 28.55893, start = as_date("2022-03-30"),
                        end = as_date("2022-04-15"), direction = c("dawn", "dusk"),
                        out_tz = "America/New_York")

blra_act <- left_join(blra_act, sun, by = "date_str") %>%
  mutate(daytime = DetTimeR >= dawn & DetTimeR <= dusk,
         DetTimeQtr = hms::as_hms(floor_date(DetTimeR, "30 mins")))

ggplot(blra_act, aes(DetTimeR, RSSI_sd_wt)) + 
  geom_point(aes(fill = daytime), shape = 21) +
  scale_fill_manual(values = c("black", "white")) +
  scale_x_datetime("Time (America/New York)",
                   date_breaks = "4 hour", date_labels = "%R") +
  labs(y = "Weighted SD of RSSI among detecting nodes") +
  facet_wrap(~ date_str, scales = "free_x", ncol = 2) + 
  theme_bw()
ggsave("output/figures/diel_activity.png", height = 9, width = 6.5)

ggplot(blra_act, aes(DetTimeQtr, RSSI_sd_wt, group=DetTimeQtr)) + geom_boxplot() +
  labs(y = "Weighted SD of RSSI among detecting nodes") +
  theme_bw() + scale_x_time("Time (America/New York)", 
                            breaks = scales::breaks_width("4 hour"))
ggsave("output/figures/diel_activity_30min.png", height = 6.5, width = 9)
