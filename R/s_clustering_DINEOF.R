# Clustering of Nordic Seas
#   Using fuzzy clustering method
#
#
# ------------ Open Data ------------

# years <- 1998:2017
#folder <- "C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/"

#for (y in 1:length(years)) {
  # filename = paste(folder, "CHLGLOB_DINEOF", as.character(years[y]),".rdata", sep="")
  # load(file = filename)
  load(file = "C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/CHLGLOB_DINEOF1998.rdata")
  lat <- data$lat[1,] #latitude
  lon <- data$lon[,1] # longitude
  xdim <- length(lon)
  ydim <- length(lat)
  xo <- data$xo
  yo <- data$yo
  adresults <- data$data[,JDstart:JDend]
  nweeks <- dim(adresults)[2]
  adresults[which(data$mask == 0),] <- NA # apply mask - bathymetry + subarctic atlantic


#---------- Interpolation and Normalization ----------------
#
#
# Linear interpolation to remove NA and divide each time series by its maximal value

  adresults_norm <- array(dim = dim(adresults)) 
  maxo_row <- array(dim = c(dim(adresults)[1],1))
  
  cpt <- 0
  for (i in 1:nrow(adresults)) {
    cpt <- cpt + 1 
    r <- adresults[i,] # time series
    d <- which(is.finite(r) == T) # Locations NA
    diff_d <- d[2:length(d)] - d[1:length(d) - 1] # number of weeks btw NA
    
    if (length(d) >= round(nweeks/2) & length(which(diff_d >= 5)) == 0) { # at least half the weeks have data and a data point every 4 weeks
      row_value <- approx(1:nweeks, r, xout = 1:nweeks, rule = 2)$z # linearly interpolate
      maxo_row[cpt] <- max(row_value) # save max value
      adresults_norm[cpt,] <- row_value/maxo_row[cpt] # normalization (divide by max)
      }
  }

  adresults <- adresults_norm # save data for clustering
  good_points <- which(is.na(adresults[,1]) == F)
  bad_points <- which(is.na(adresults[,1]) == T)
  xo[bad_points] <- NA
  yo[bad_points] <- NA
  
  #--------------- Clustering --------------------
  #
  # the number of clusters that should be used is ~ 6
  
  
#}