multilaterate_beep_dat <- function(beep_dat, Tags = NULL, min_node_beeps = 4, min_n_nodes = 4) {
  if (is.null(Tags)) stop("You must specify at least one Tag to multilaterate.", call. = FALSE)
  all_tag_loc_ests <- lapply(Tags, function(Tag) {
    tmp_beep <- filter(beep_dat, TagId == Tag, n >= min_node_beeps)
    
    # Center X and Y to ease computation
    x_mn <- mean(tmp_beep$X)
    y_mn <- mean(tmp_beep$Y)
    tmp_beep <- tmp_beep %>%
      # Reduce to detections from at minimum number of nodes
      group_by(TagId, SppCode, Age, Sex, DeployDT, DetTimeR) %>%
      mutate(n_nodes = n(),
             Xc = X - x_mn,
             Yc = Y - y_mn) %>% 
      filter(n_nodes >= min_n_nodes) %>%
      group_by(n_nodes, .add = TRUE) %>%
      mutate(loc_grp = cur_group_id())
    # Split by group to facilitate parallel processing
    tmp_beep_list <- group_split(tmp_beep)
    n_locs <- length(tmp_beep_list)
    n_cores <- min(detectCores() - 1, as.integer(max(1, floor(n_locs/10000))))
    message("Estimating ", n_locs, " locations for TagId ", Tag, " on ", n_cores, " cores/threads.")
    cl <- NULL
    if (n_cores > 1) {
      cl <- makeCluster(n_cores)
      clusterExport(cl, c("tmp_beep_list", "x_mn", "y_mn", "EucDist", "distEuc"),
                    envir=environment())
    }
    loc_ests <- pbapply::pblapply(tmp_beep_list, function(dat) {
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
      mod_summary <- tibble::tibble(loc_grp = unique(dat$loc_grp),
                                    X_est = loc_coefs[1] + x_mn,
                                    Y_est = loc_coefs[2] + y_mn,
                                    X_est_sd = loc_sd[1],
                                    Y_est_sd = loc_sd[2],
                                    RSSI_sd_wt = stats::weighted.mean(dat$adjTagRSSI_sd, dat$n),
                                    RSSI_cv_wt = stats::weighted.mean(dat$adjTagRSSI_cv, dat$n),
                                    df = stats::df.residual(loc_mod),
                                    converged = summary(loc_mod)$convInfo$isConv,
                                    varcov = list(loc_vcov))
      mod_summary
    }, cl = cl)
    if (n_cores > 1) stopCluster(cl)
    loc_ests <- bind_rows(loc_ests)
    out_loc_ests <- tmp_beep %>% group_by(loc_grp, .add = TRUE) %>%
      summarize(.groups = "drop") %>%
      left_join(loc_ests, by = "loc_grp") %>%
      arrange(DetTimeR)
    out_loc_ests
  })
  names(all_tag_loc_ests) <- Tags
  all_tag_loc_ests
}

# Modification of geosphere::distm to not treat X-Y data at lat/lon
# and use Eucldean distance by default
EucDist <- function (x, y, fun = distEuc) {
  x <- geosphere:::.pointsToMatrix(x, checkLonLat = FALSE)
  y <- geosphere:::.pointsToMatrix(y, checkLonLat = FALSE)
  n = nrow(x)
  m = nrow(y)
  dm = matrix(ncol = m, nrow = n)
  for (i in 1:n) {
    dm[i, ] = fun(x[i, ], y)
  }
  return(dm)
}

distEuc <- function (p1, p2) {
  p1 <- geosphere:::.pointsToMatrix(p1, checkLonLat = FALSE)
  p2 <- geosphere:::.pointsToMatrix(p2, checkLonLat = FALSE)
  p = cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2])
  dLat <- p[, 4] - p[, 2]
  dLon <- p[, 3] - p[, 1]
  dist <- sqrt(dLat^2 + dLon^2)
  return(as.vector(dist))
}
