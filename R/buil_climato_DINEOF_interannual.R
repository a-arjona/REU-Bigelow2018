# Prepare the data matrices for the analysis with interpolated CHL
#
#
# -----------Create data---------------------------------
# Result: Data matrices for each year (20 total)
# 46 wks x ~10k stations

rm(list=ls()) # clear all variables
setwd("C:/files/Work/Bigelow/Data/")
#----- Import Masks from .txt file
# 3 masks (1 = ok, 0 = remove): bathy > 300, SubArcticAtlantic, both
# mask gets rid of data from areas shallower than 300 meters (near the coasts) also gets rid of data outside the lat lon of interest

mask = read.table(file="./txt_files/mask.txt")

library("ncdf4")

folder <- "./DINEOF_2018_raw_data/" # address of folder where the data is

nweeks <- 46
julianday <- seq(0,365,8) # julian day is one 8 day week
years <- 1998:2017


# make matrix for every year
for (y in 1:length(years)) {
  # current year
  c_year <- years[y]
  
  # gives filename for all files (not all files exist or have data) from given year (46)
  weeks <- format(as.Date(paste(as.character(c_year), "-01-01",sep = "")) + julianday, "*%Y%m%d")
  
  # open files and have a yearly climatology
  chl_year <- c()
  
  for (w in 1:length(weeks)) {
    filename <- list.files(path = folder, pattern = weeks[w], full.names = T) # filename
    
    if (length(filename) != 0) { # checks whether there is any data for that week
      nc <- nc_open(filename) # if data ---> open file
      
      if (y == 1 & w == 7) {
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
      
      # subset at the location wanted
      CHL <- ncvar_get(nc, varid = "CHL1_intp")
      CHL <- CHL[,dim(CHL)[2]:1] # reverse latitudes so they go low to high
      chl_year <- cbind(chl_year, c(CHL)) # from matrix to vector
    } else {
      chl_year <- cbind(chl_year, array(dim = c(20832,1)))
    }
  }
  
  #--------- Export lon & lat in .txt file to be opened in MATLAB, to build masks with MATLAB
  #loc <- cbind(LONv, LATv)
  # p sure i only need to do this once
  # new_loc <- paste(folder, "lonlat", as.character(c_year), ".txt", sep="")
  # write.table(loc,file = new_loc, sep = "\t", row.names = F, col.names = F)
  
  #----- Build data and save
  data <- c()
  data$data <- chl_year
  data$lon <- LON
  data$lat <- LAT
  data$xo <- xo
  data$yo <- yo
  data$LON <- LONv
  data$LAT <- LATv
  data$mask <- mask[,1]
  
  save_loc <- "./DINEOF_2018_processed_data/"
  new_file <- paste(save_loc, "CHLGLOB_DINEOF", as.character(c_year), ".rdata", sep = "")
  save(data, file = new_file)
  
}

#--------- Export lon & lat in .txt file to be opened in MATLAB, to build masks with MATLAB
loc <- cbind(LONv, LATv)
write.table(loc,file = "./txt_files/", sep = "\t", row.names = F,col.names = F)

clu_array <- c()
clu_array$xo <- xo
clu_array$yo <- yo
save(clu_array, file = "./DINEOF_2018_processed_data/array_index.rdata")
  
