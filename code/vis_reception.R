pacman::p_load(sf, raster, viridis)

vis_reception <- function(x_range, y_range, sigma, spacing = sigma, min_nodes = 3, res = 5) {
  cellsize <- spacing
  
  # We want center coordinate @ nest so we use some trickery
  halfx <- ceiling(x_range / 2)
  halfy <- ceiling(y_range / 2)
  sfc <- st_sfc(st_polygon(list(rbind(c(0, 0), c(0,halfy), c(halfx, halfy), c(halfx, 0), c(0, 0)))))
  quad1 <- st_make_grid(sfc, cellsize = cellsize, square = FALSE)
  pts1 <- round(st_coordinates(st_centroid(quad1)))
  # restrict to area of interest
  pts1 <- pts1[(pts1[, "X"] >= 0 & pts1[, "X"] <= halfx) & (pts1[, "Y"] >= 0 & pts1[, "Y"] <= halfy), ]
  pts2 <- cbind(pts1[, "X"] * -1, pts1[, "Y"])
  pts3 <- cbind(pts1[, "X"] * -1, pts1[, "Y"] * -1)
  pts4 <- cbind(pts1[, "X"], pts1[, "Y"] * -1)
  xy <- rbind(pts1, pts2, pts3, pts4)
  xy <- xy[!duplicated(xy), ]
  rownames(xy) <- NULL
  
  reception <- function(D, sigma) {
    w <- 1/(ceiling(D / sigma))
    w[w < 1] <- 0
    w
  }
  
  # Define study area
  sa <- raster::raster(nrows=320, ncols=320, res = res,
                       xmn = -halfx - sigma, xmx = halfx + sigma, 
                       ymn = -halfy - sigma, ymx = halfy + sigma,
                       crs = "+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
                       vals = 0)
  
  for (i in seq_len(nrow(xy))) {
    dr <- raster::distanceFromPoints(sa, xy[i, ])
    dr[] <- reception(dr[], sigma)
    if (i == 1) st <- stack(dr) else st <- stack(st, dr)
  }
  out <- floor(calc(st, sum))
  out[out < min_nodes] <- NA
  nl <- length(unique(out[]))
  plot(floor(out), col = viridis(nl), main = paste("Detection:", sigma, "m\n",
                                                   "Spacing: X =", spacing, "m; Y =", mean(diff(sort(xy[xy[, 1] == 0, 2]))), "m\n",
                                                   nrow(xy), "nodes"))
  points(xy, col = "gray60", lwd=8)
  
  # Tidy up coordinates
  xy <- as.data.frame(xy) %>%
    dplyr::arrange(Y, X)
  return(xy)
}
  