# build data matrices from raw satellite data
# returns data matrix (pixels by weeks)
buildData <- function(start, end = start, nweeks = 46, pixels = 20832) { 
  # type true is mean 
  # type false is annual
  
  library("ncdf4")
  folder <- "C:/files/Work/Bigelow/Data/DINEOF_2018_raw_data/"

  julianday <- seq(0,365,8)
  years <- start:end
  
  # Import mask from .txt file
  mask <- read.table(file = "C:/files/Work/Bigelow/Data/txt_files/mask.txt")
  
  
    chl_clim <- c()
    
  for (n in 1:nweeks) { # repeat for every week
    # date of satellite images
    days <- format(as.Date(paste(as.character(years), "-01-01", sep = "")) + julianday[n], "*%Y%m%d")
    
    # open files and create a weekly climatology
    chl_week <- c()
    
    for (y in 1:length(days)) {
      filename <- list.files(path = folder, pattern = days[y], full.names = T) # filename
      
      if (length(filename) != 0) { # checks whether there is any data for that week
        nc <- nc_open(filename) # if data --> read file
        
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
        
        # subset at the location wanted
        CHL <- ncvar_get(nc, varid = "CHL1_intp")
        CHL <- CHL[,dim(CHL)[2]:1] # reverse the latitudes so they go low to high
        chl_week <- cbind(chl_week, c(CHL)) # add data from this week as column (total 46 col)
        
      }
      else {
        chl_week <- cbind(chl_week, array(dim = c(pixels,1)))
      }
    }
    
    chl_clim <- cbind(chl_clim, rowMeans(chl_week, na.rm = TRUE)) 
  }
  
  # Export lon & lat in .txt file to be open in Matlab
  loc <- cbind(LONv, LATv)
  write.table(loc, file = "C:/files/Work/Bigelow/Data/txt_files/lonlat.txt", sep = "\t", row.names = F, col.names = F)
  
  data <- c()
  data$data <- chl_clim
  data$lon <- LON
  data$lat <- LAT
  data$xo <- xo
  data$yo <- yo
  data$LON <- LONv
  data$LAT <- LATv
  data$mask <- mask[,1]
  
  return(data)
}

# Interpolate and Normalize Data
# Returns processed data. Object has data, good_points, bad_points, xo, yo, lat, lon
intNorm <- function(data, wkstart = 9, wkend = 35, pixels = 20832) {
  # define variables
  lat <- data$lat[1,] # Latitude 
  lon <- data$lon[,1] # Longitude 
  xo <- data$xo 
  yo <- data$yo 
  sub_data <- data$data[,wkstart:wkend] # subset of data from week 9 to week 35 (default)
  nweeks <- dim(sub_data)[2]
  sub_data[which(data$mask == 0),] <- NA # apply mask bathy + subarcticatlantic
  
  data_norm <- array(dim = dim(sub_data)) # Colonne = Time (variables), Line = Pixel (stations)
  maxo_row <- array(dim = c(dim(sub_data)[1],1))
  
  cpt <- 0
  for (i in 1:nrow(sub_data)) {
    cpt <- cpt + 1
    r <- sub_data[i,] # Time series
    d <- which(is.finite(r) == T) # Locations NA
    diff_d <- d[2:length(d)] - d[1:length(d) - 1] # number of weeks btw NA
    
    if (length(d) >= round(nweeks/2) & length(which(diff_d >= 5)) == 0) { # at least a half weeks of data with a data point every 4 weeks 
      row_value <- approx(1:nweeks,r,xout = 1:nweeks,rule = 2)$y # linearly interpolate
      maxo_row[cpt] <- max(row_value) # save max value
      data_norm[cpt,] <- row_value/maxo_row[cpt] # Normalization (divide by max)
    }
  }
  
  pData <- c()
  pData$timeseries <- data_norm # save data for clustering
  pData$good_points <- which(is.na(pData$data[,1]) == F) # pixels without NA to be cluterized
  pData$bad_points <- which(is.na(pData$data[,1]) == T) # pixels with NA
  xo[pData$bad_points] <- NA
  yo[pData$bad_points] <- NA
  pData$xo <- xo
  pData$yo <- yo
  pData$lat <- lat
  pData$lon <- lon
  pData$maxo <- maxo_row
  
  return(pData)
}
  
# Clusters data, returns object with kmeans and fcm
cluster <- function(data, center) {
  library("e1071")
  ts <- data$timeseries
  gp <- data$good_points
  goodData <- ts[gp,]
  
  cclData <- cmeans(goodData, center, iter.max = 50, dist = "euclidean", m = 2)
  
  cclData$timeseries <- goodData
  cclData$good_points <- gp
  cclData$lon <- data$lon
  cclData$lat <- data$lat
  cclData$xo <- data$xo
  cclData$yo <- data$yo
  cclData$maxo <- data$maxo
  
  
  return(cclData)
}

# Finds silhouette value. Returns si, clusters (if two arguments the closest clusters), and timeseries
si <- function(data, data2 = NULL) {
  library("fields") 
  ts1 <- data$timeseries
  ts2 <- data2$timeseries
  
  distX <- rdist(ts1, ts2) # distance matrix
  if (!is.null(data2)) { # if there is a second arguement 
    grps <- data2$cluster
  }
  else{# if no second argument
    grps <- data$cluster
  } 
  
  c = apply(distX, 1, function(x) {b = aggregate(x, list(grps), mean, simplify = T) # gives matrix six clusters by pixels with mean distance to every point within one cluster
                                                 return(b[,2])})
  
  new_cl <- apply(c, 2, which.min) # finds the min distance of the column (column rep a pixel and row being a cluster)
  a <- apply(c, 2, min) # gives the actual min dist value
  b <- apply(c, 2, function(x) min(x[c(1:6) != which.min(x)])) # finds the second minimum value
  si <- (b - a)/pmax(a,b)
  
  newData <- c()
  newData$timeseries <- ts1
  newData$cluster <- new_cl
  newData$si <- si
  newData$good_points <- data$good_points
  newData$xo <- data$xo
  newData$yo <- data$yo
  newData$lat <- data$lat
  newData$lon <- data$lon
  newData$maxo <- data$maxo
  
  return(newData)
  
}

# Finds avg si value per pixel. Returns mean si, xo, yo
avgSi <- function(..., pixels = 20832) {
  data <- list(...)
  xo <- data[1]$xo
  yo <- data[1]$yo
  S <- array(NA, dim = c(pixels, length(data)))
  for (i in 1:length(data)) {
    subData <- data[i]
    good_points <- subData$good_points
    S[good_points,i] <- subData$si
  }
  
  meanSi <- array(NA, dim = c(dim(S)[1],1))
  
  for (r in 1:dim(S)[1]) {
    row <- S[r,]
    if (length(which(is.finite(row) == T)) >= dim(S)[2]/2) {
      meanSi[r] <- mean(row)
    }
  }
  
  bad_points <- which(is.na(mean_si) == T) # pixels with NA
  
  xo[bad_points] <- NA
  yo[bad_points] <- NA
  
  meanSi <- na.omit(mean_si)
  xo <- na.omit(xo)
  yo <- na.omit(yo)
  
  mSi <- c()
  mSi$meanSi <- meanSi
  mSi$xo <- xo
  mSi$yo <- yo
  
  return(mSi)
}

# Finds mode cluster and percent associated for given data
modePerc <- function(..., pixels = 20832) {
  data <- list(...)
  xo_c <- data[1]$xo
  yo_c <- data[1]$yo
  xo_p <- xo_c
  yo_p <- yo_c
  C <- array(NA, dim = c(pixels, length(data)))
  for (i in 1:length(data)) {
    subData <- data[i]
    good_points <- subData$good_points
    C[good_points,i] <- subData$cluster
  }
  
  modeCl <- array(NA, dim = c(dim(C)[1],1))
  percAssoc <- array(NA, dim = c(dim(C)[1],1))
  
  for (r in 1:dim(modeCl)[1]) {
    if (length(which(is.finit(C[r,]) == T)) >= dim(C)[2]/2) {
      row <- C[r,]
      uniqv <- unique(row)
      modeCl[r] <- uniqv[which.max(tabulate(match(row,uniqv)))]
      perc_assoc[r] <- length(which(row == modeCl[r]))/dim(C)[2]
    }
  }
  
  bad_points_p <- which(is.na(percAssoc) == T)
  percAssoc <- na.omit(percAssoc)
  xo_p[bad_points_p] <- NA
  yo_p[bad_points_p] <- NA
  xo_p <- na.omit(xo_p)
  yo_p <- na.omit(yo_p)
  
  bad_points_c <- which(is.na(mclusters) == T)
  xo_c[bad_points_c] <- NA
  yo_c[bad_points_c] <- NA
  xo_c <- na.omit(xo_c)
  yo_c <- na.omit(yo_c)
  modeCl <- na.omit(modeCl)
  
  pa <- c()
  pa$percentAssociated <- percAssoc
  pa$xo <- xo_p
  pa$yo <- yo_p
  
  mclu <- c()
  mclu$modeCl <- modeCl
  mclu$xo <- xo_c
  mclu$yo <- yo_c
  
  return(pa, mclu)
}

# Mean Data
processMean <- function(start, end, nclu) {
  climato <- buildData(start,end) # build data matrix from raw data
  climato <- intNorm(climato) # clean data
  
  clData <- cluster(climato, nclu) # cluster the mean data
  sidata <- si(clData) # get silhouette indexes
  
  
}