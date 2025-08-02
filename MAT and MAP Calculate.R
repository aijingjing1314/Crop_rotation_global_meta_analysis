# Load required packages
library(readxl)
library(terra)
library(writexl)

# Select Excel file (e.g., "Coordinates.xlsx")
file_path <- file.choose()

# Read the location data
locations <- read_excel(file_path)
colnames(locations)[1:3] <- c("StudyID", "Latitude", "Longitude")

# Ensure that Latitude and Longitude are numeric
locations$Latitude <- as.numeric(locations$Latitude)
locations$Longitude <- as.numeric(locations$Longitude)

# Read precipitation data (12 monthly layers from WorldClim)
prec_files <- sprintf("WorldClim/wc2.1_30s_prec_%02d.tif", 1:12)
prec_stack <- rast(prec_files)
annual_prec <- app(prec_stack, mean, na.rm = TRUE)  # Mean annual precipitation (MAP)

# Read temperature data (12 monthly layers from WorldClim)
temp_files <- sprintf("WorldClim/wc2.1_30s_tavg_%02d.tif", 1:12)
temp_stack <- rast(temp_files)
annual_temp <- app(temp_stack, mean, na.rm = TRUE)  # Mean annual temperature (MAT)

# Construct coordinate matrix and preserve original row order
coords <- cbind(locations$Longitude, locations$Latitude)

# Extract MAP and MAT values for each coordinate
locations$MAP <- extract(annual_prec, coords)[, "mean"]
locations$MAT <- extract(annual_temp, coords)[, "mean"]

# Export the results to Excel
write_xlsx(locations, "Climate_Data.xlsx")

