node_coverage <- function(sa_sf, node_range = 40) {
  
  sa_sp <- as(sa_sf, "Spatial")
  sa_buff <- st_buffer(sa_sf, dist = node_range) %>% as("Spatial")
  
  # Hexagon grid with inradius of `node_range` m
  hex_pts <- sp::spsample(sa_buff, type="hexagonal", cellsize = node_range, 
                      offset = c(0, 0))
  
  # Keep only hexagons largely contained within study area boundary
  hex_pts <- hex_pts[sa_sp, ]
  
  x <- hex_pts@coords[, "x"]; y <- hex_pts@coords[, "y"]
  
  reception <- function(D, sigma) {
    w <- 1/(ceiling(D / sigma))
    w[w < 1] <- 0
    w
  }
  
  # Define study area
  sa_r <- raster::raster(sa_sf, res = node_range / 25, vals = 0)
  sa_r <- fasterize::fasterize(sa_sf, sa_r)
  
  for (i in seq_along(x)) {
    dr <- raster::distanceFromPoints(sa_r, cbind(x[i], y[i]))
    dr[] <- reception(dr[], node_range)
    if (i == 1) st <- raster::stack(dr) else st <- raster::stack(st, dr)
  }
  
  out <- floor(raster::calc(st, sum))
  out[out == 0] <- NA
  out <- raster::mask(out, sa_r)
  out_df <- as.data.frame(as(out, "SpatialPixelsDataFrame"))
  nl <- length(unique(out_df$layer))
  
  sa_ll <- st_transform(sa_sf, '+proj=longlat +datum=WGS84')
  node_pts <- st_transform(st_as_sf(hex_pts), '+proj=longlat +datum=WGS84')
  out_ll <- raster::projectRaster(out, crs = '+proj=longlat +datum=WGS84', method = "ngb")
  
  pal <- colorFactor(viridis::viridis(nl), raster::values(out_ll), 
                     na.color = "transparent")
  tag.map.title <- htmltools::tags$style(
    htmltools::HTML(".leaflet-control.map-title { 
      transform: translate(-50%,20%);
      position: fixed !important;
      left: 50%;
      text-align: center;
      padding-left: 5px; 
      padding-right: 5px; 
      background: rgba(255,255,255,0.75);
      font-weight: bold;
      font-size: 16px;
      }"))
  title <- 
    htmltools::tags$div(tag.map.title, 
                        htmltools::HTML(paste0("Assumed consistent detection distance: ", node_range, " m<br>",
                                               "Number of nodes required: ", length(x))))
  p <- leaflet() %>%
    addProviderTiles("Esri.WorldImagery") %>%
    addRasterImage(out_ll, opacity = 0.75, colors = pal, project = FALSE,
                   group = "Node summary") %>%
    addPolygons(data = sa_ll, fill = FALSE, color = "black", opacity = 1,
                weight = 3) %>%
    addCircleMarkers(data = node_pts, radius = node_range/50, 
                     fillColor = "red", weight = 1, fillOpacity = 0.5,
                     color = "black", opacity = 1) %>%
    addLayersControl(position = "topright", overlayGroups = "Node summary",
                     options = layersControlOptions(collapse = FALSE)) %>%
    addLegend("topright", pal = pal, values = raster::values(out_ll),
              title = "Node detections", opacity = 0.75) %>%
    addScaleBar(position = "bottomright") %>%
    addControl(title, position = "topleft", className="map-title")
  print(p)
  return(node_pts)
}
