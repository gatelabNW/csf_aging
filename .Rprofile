# Activate env ----------------------------------------------------

source("renv/activate.R")
renv::status()

# Load in internal lib files --------------------------------------

# Load in internal lib
lib_files <- list.files(path = "lib", full.names = TRUE)
for (file in lib_files) {source(file)}

# Print execution statement ---------------------------------------
print(".Rprofile of csf_aging project has been sourced.")
