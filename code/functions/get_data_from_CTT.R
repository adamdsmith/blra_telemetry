message("Pulling telemetry data from CTT servers")

# Path to store the raw data
outpath <- "data_raw/" 

# CTT API access key
my_token <- "1acdee4778800e09e1d50f2ada8471bf16712c0987731969570be32e1f862ed2"

# At some point, this will be much easier/faster if we can simply add to a local
# SQL database, but I don't have the software on my federal machine yet

# db_name <- "calibration"
# conn <- dbConnect(RPostgres::Postgres(), dbname = db_name)

# Download the data from CTT servers
# The `get_my_data` function has been modified to allow us to pass the BLRA project
# so we don't get all of the USFWS data on CTT servers, and add a progress bar
# Still, I think it brings down all the data each time...
options(readr.show_progress = FALSE)
get_my_data(my_token, outpath, proj = "GA Black Rail")#, conn)

# update_db(conn, outpath)
# dbDisconnect(conn)
