### Clusterization of the Red Sea ---------------------------------

# require("fpc") # package use for kmeans function
# require("fields") # package use for mapping



# --- Open data --------------------------------------------------
#
# Here mr.rdata => data == List of 5
#   $ data: num [1:138240, 1:46] NA NA NA NA NA NA NA NA NA NA ...
#   $ lat : num [1:288, 1:480] 30 30 30 30 30 ...
#   $ lon : num [1:288, 1:480] 32 32.1 32.1 32.1 32.2 ...
#   $ xo  : int [1:138240] 1 2 3 4 5 6 7 8 9 10 ...
#   $ yo  : int [1:138240] 1 1 1 1 1 1 1 1 1 1 ...
library("ppclust")
load(file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/CHLGLOB_DINEOF19982017.rdata") # load data from build climato

lat <- data$lat[1,] # Latitude 
lon <- data$lon[,1] # Longitude 
xdim <- length(lon)
ydim <- length(lat)
xo <- data$xo 
yo <- data$yo 
adresults <- data$data[,JDstart:JDend]
adresults <- data$data[,JDstart:JDend]
nweeks <- dim(adresults)[2]
adresults[which(data$mask == 0),] <- NA # apply mask bathy + subarcticatlantic


# --- Interpollation and Normalization ----------------------------------------------
#
#
# Linear interpolation to remove NA, and divide each time series by its maximal value 
# Guesstimates data points using a linear regression but only if missing one or two 
# if missing more data points in a row does not interpolate

adresults_norm <- array(dim = dim(adresults)) # Colonne = Time (variables), Line = Pixel (stations)
maxo_row <- array(dim = c(dim(adresults)[1],1))

cpt <- 0
for (i in 1:nrow(adresults)) {
  cpt <- cpt + 1
  r <- adresults[i,] # Time series
  d <- which(is.finite(r) == T) # Locations NA
  diff_d <- d[2:length(d)] - d[1:length(d)-1] # number of weeks btw NA
  
  if (length(d) >= round(nweeks/2) & length(which(diff_d >= 5)) == 0) { # at least a half weeks of data with a data point every 4 weeks 
    row_value <-approx(1:nweeks,r,xout=1:nweeks,rule=2)$y # linearly interpolate
    maxo_row[cpt] <- max(row_value) # save max value
    adresults_norm[cpt,] <- row_value/maxo_row[cpt] # Normalization (divide by max)
  }
}


adresults <- adresults_norm # save data for clustering
good_points <- which(is.na(adresults[,1]) == F) # pixels without NA to be cluterized
bad_points <- which(is.na(adresults[,1]) == T) # pixels with NA
xo[bad_points] <- NA
yo[bad_points] <- NA



# --- Clustering --------------------------------------------------
#
#
# the number of clusters that should be used will be analyzed later 

# nclu <- 7 # Number of clusters
cl<-kmeansCBI(adresults[good_points,],k=nclu,scaling=F,runs=50) # clustering Kmeans, runs: Number of starts of the k-means. scaling: If scaling is TRUE then scaling is done by dividing the variables by their root-mean-square

# save results
cl$good_points <- good_points
cl$lon <- lon
cl$lat <- lat
cl$xo <- xo
cl$yo <- yo
cl$adresults <- adresults
cl$maxo <- maxo_row
