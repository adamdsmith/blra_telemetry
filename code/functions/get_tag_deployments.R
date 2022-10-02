get_tag_deployments <- function() {
  
  if (!requireNamespace("googlesheets4", quietly = TRUE))
    pacman::p_load("googlesheets4")
  
  out_fn <- file.path("data_raw", "tag_deployments.rds")
  
  # Googlesheets access code for archive only, restricted access
  gs_title <- "CTT PowerTag deployment log"
  
  # Set up column types for call # (or call type for BLRA) columns
  col_types <- "ccncccccccc"
  
  out <- googledrive::drive_get(gs_title) %>%
    googlesheets4::read_sheet(col_types = col_types)  %>%
    mutate(DeployDT = ymd_hm(paste(DeployDate, stringr::str_pad(DeployTime, 4, pad = "0")), 
                               tz = "America/New_York")) %>%
    dplyr::select(TagId:Sex, DeployDT, BaseStation)
  
  attr(out$DeployDT, "tzone") <- "UTC"
  
  saveRDS(out, file = out_fn)
  out
}