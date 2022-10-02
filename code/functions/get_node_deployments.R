get_node_deployments <- function() {
  
  if (!requireNamespace("googlesheets4", quietly = TRUE))
    pacman::p_load("googlesheets4")
  
  out_fn <- file.path("data_raw", "node_deployments.rds")
  
  # Googlesheets access code for archive only, restricted access
  gs_title <- "CTT Node Position Log"
  
  # Set up column types for call # (or call type for BLRA) columns
  col_types <- "cccciciciic"
  
  out <- googledrive::drive_get(gs_title) %>%
    googlesheets4::read_sheet(col_types = col_types)  %>%
    mutate(NodeStartDT = ymd_hm(paste(StartDate, stringr::str_pad(StartTime, 4, pad = "0")), 
                               tz = "America/New_York"),
           NodeEndDT = ymd_hm(paste(EndDate, stringr::str_pad(EndTime, 4, pad = "0")), 
                                 tz = "America/New_York")) %>%
    dplyr::select(TagId:BaseStation, NodeStartDT, NodeEndDT, X = UTME_17N, Y = UTMN_17N) %>%
    # If node still active, set NodeEndDT to current time
    replace_na(list(NodeEndDT = now())) %>% 
    mutate(Xcoord = X, Ycoord = Y) %>%
    ### TEMPORARY
    filter(!is.na(X), !is.na(Y))
  attr(out$NodeStartDT, "tzone") <- attr(out$NodeEndDT, "tzone") <- "UTC"
  out <- sf::st_as_sf(out, coords = c("Xcoord", "Ycoord"), crs = 32617)
  saveRDS(out, file = out_fn)
  out
}