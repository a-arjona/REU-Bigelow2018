# Build 10-year Average Data
#
#
# 2007 large melting event. Can be used to divide data

# --------- Create Data -----------------

rm(list = ls()) # clear all variables
setwd("C:/files/Work/Bigelow/Data/")
library("ncdf4")

folder <- "./DINEOF_2018_raw_data/" # address of folder where data is

nweeks <- 46
julianday <- seq(0,365,8)
years <- 1998:2017

for (yrs in 1:2) {
  chl_clim <- c()
  if (yrs == 1) {
    set <- 1998:2007
    data_file <- "./DINEOF_2018_processed_data/CHLGLOB_DINEOF19982007.rdata"
  }
  if (yrs == 2) {
    set <- 2008:2017
    data_file <- "./DINEOF_2018_processed_data/CHLGLOB_DINEOF20072017.rdata"
  }
  for (n in 1:nweeks) {
    # date of satelite images 
    days <- format(as.Date(paste(as.character(set),"-01-01",sep = "")) + julianday[n], "*%Y%m%d")
    
    chl_week <- c()
    for (y in 1:length(days)) {
      filename <- list.files(path = folder, pattern = days[y], full.names = T) # filename
      
      if (length(filename) != 0) { # checks whether there is any data for that week
        nc <- nc_open(filename) # if data --> read file
        
        # open latitudes and longitudes data only one time
        # do only once bc lat and lon are the same for all the files, only need to do once
        if (n == 23 & y == 1) { 
          lat <- ncvar_get(nc, "latitude")
          lon <- ncvar_get(nc, "longitude")
          
          LON <- replicate(length(lat),lon)
          LAT <- t(replicate(length(lon),lat))
          
          LAT <- LAT[,dim(LAT)[2]:1]
          
          # longitudes and latitudes from matrix to vector
          LONv <- c(LON)
          LATv <- c(LAT)
          xo <- match(LONv, LON[,1])
          yo <- match(LATv, LAT[1,])
          
        }
        
        # Subset at the location wanted
        CHL <- ncvar_get(nc, varid = "CHL1_intp")
        CHL <- CHL[,dim(CHL)[2]:1] # reverse latitudes
        chl_week <- cbind(chl_week, c(CHL)) # from matrix to vector
        
      } else {
        chl_week <- cbind(chl_week, array(dim = c(20832,1)))
      }
    }
    
    chl_clim <- cbind(chl_clim, rowMeans(chl_week, na.rm = TRUE)) # build climatological dataset
  }

#----- Import Masks from .txt file

mask = read.table(file = "./txt_files/mask.txt")

# ---- Build data and save ------
data <- c()
data$data <- chl_clim
data$lon <- LON    
data$lat <- LAT
data$xo <- xo
data$yo <- yo
data$LON <- LONv
data$LAT <- LATv
data$mask <- mask[,1]
save(data, file = data_file)

}