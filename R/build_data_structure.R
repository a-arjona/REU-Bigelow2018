# ---------- Build Data Frames ------------
#
#
#
# Build 20 files that contain the data in a data frame (year, week, pixel, CHL_data, lat, lon, mask)

rm(list=ls()) # clear all variables

#----- Import Masks from .txt file
# 3 masks (1 = ok, 0 = remove): bathy > 300, SubArcticAtlantic, both
# mask gets rid of data from areas shallower than 300 meters (near the coasts) also gets rid of data outside the lat lon of interest

mask = read.table(file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/mask.txt")

library("ncdf4")
library("dplyr")

folder <- "C:/Users/Ade/Documents/Bigelow/DINEOF-2018/DINEOF/" # address of folder where the data is

nweeks <- 46
julianday <- seq(0,365,8) # julian day is one 8 day week
years <- 1998:2017



# make matrix for every year
for (y in 1:length(years)) {
  # current year
  c_year <- years[y]
  pix_col <- c()
  wk_col <- c()
  LONv <- c()
  LATv <- c()
  CHLv <- c()
  mask_col <- c()
  
  # gives filename for all files (not all files exist or have data) from given year (46)
  weeks <- format(as.Date(paste(as.character(c_year), "-01-01",sep="")) + julianday, "*%Y%m%d")
  
  # open files and have a yearly climatology
  chl_year <- c()
  
  for (w in 1:length(weeks)) {
    filename <- list.files(path = folder, pattern = weeks[w], full.names = T) # filename
    
    if (length(filename) != 0) { # checks whether there is any data for that week
      nc <- nc_open(filename) # if data ---> open file
      
      lat <- ncvar_get(nc, "latitude")
      lon <- ncvar_get(nc, "longitude")
      
     
      n_pix <- length(lat) * length(lon)
      wk_col <- c(wk_col, rep(w, n_pix))
      
      
      LON <- replicate(length(lat),lon) 
      LAT <- t(replicate(length(lon),lat))
      
      LAT <- LAT[,dim(LAT)[2]:1] 
      
      # longitudes and latitudes from matrix to vector
      LONv <- c(LONv, LON)
      LATv <- c(LATv, LAT)
      xo <- match(LONv, LON[,1])
      yo <- match(LATv, LAT[1,])
      
      
      # subset at the location wanted
      CHL <- ncvar_get(nc, varid = "CHL1_intp")
      CHL <- CHL[,dim(CHL)[2]:1] # reverse latitudes so they go low to high
      chl_year <- cbind(chl_year, c(CHL)) # from matrix to vector 
      CHLv <- c(CHLv, CHL)
      pix_col <- c(pix_col, 1:n_pix)
      mask_col <- c(mask_col, mask[,1])
      year_col <- rep(c_year, length(wk_col))
    }
    }
  
  
  
  #--------- Export lon & lat in .txt file to be opened in MATLAB, to build masks with MATLAB
  #loc <- cbind(LONv, LATv)
  # p sure i only need to do this once
  # new_loc <- paste(folder, "lonlat", as.character(c_year), ".txt", sep="")
  # write.table(loc,file = new_loc, sep = "\t", row.names = F, col.names = F)
  
  
  # Build data frame
  data <- data.frame("year" = year_col, "week" = wk_col, "pixel" = pix_col, "CHL_data" = CHLv, "Longitude" = LONv, "Latitude" = LATv, "Mask" = mask_col)
  
  save_loc <- "C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/"
  new_file <- paste(save_loc, "CHLGLOB_DINEOF", as.character(c_year), ".rda", sep="")
  save(data, file = new_file)

}

#--------- Export lon & lat in .txt file to be opened in MATLAB, to build masks with MATLAB
loc <- cbind(LONv, LATv)
write.table(loc,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/lonlat.txt", sep = "\t", row.names=F,col.names = F)
