# Prepare the data matrix for the analysis with interpolated CHL
#
#
# --- Create data --------------------------------------------------
#
# data == List of 5


rm(list=ls()) # clear all variable
setwd("C:/files/Work/Bigelow/Data/")
library("ncdf4")

folder <- "./DINEOF_2018_raw_data/" # address of folder where data is
# creates a subset of data from global set
# lon_lim <- c(-70,20) 
# lat_lim <- c(59,82)

nweeks <- 46
julianday <- seq(0,365,8) # julian day = one 8 day week
years <- 1998:2017
chl_clim <- c()

for (n in 1:nweeks) {
  # date of satellite images
  days <- format(as.Date(paste(as.character(years),"-01-01",sep = "")) + julianday[n], "*%Y%m%d")
  
  # open files and create a weekly climatology
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


#----- Export lon & lat in .txt file to be open in Matlab, to build masks with Matlab 
loc <- cbind(LONv, LATv)
write.table(loc,file = "./txt_files/lonlat.txt", sep = "\t", row.names = F,col.names = F)


#----- Import Masks from .txt file
# 3 masks (1 = ok, 0 = remove): bathy > 300, SubArcticAtlantic, both
# mask gets rid of data from areas shallower than 300 meters (near the coasts) also gets rid of data outside the lat lon of interest

mask = read.table(file = "./txt_files/mask.txt")


#----- Build data and save
data <- c()
data$data <- chl_clim
data$lon <- LON
data$lat <- LAT
data$xo <- xo
data$yo <- yo
data$LON <- LONv
data$LAT <- LATv
data$mask <- mask[,1]
save(data, file = "./DINEOF_2018_processed_data/CHLGLOB_DINEOF19982017.rdata")
