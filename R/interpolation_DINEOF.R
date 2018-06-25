# Clustering of Nordic Seas
# 
#
#
# ------------ Open Data ------------

load(file = cyear_file) # edit to cycle through each year
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
mino_row <- array(dim = c(dim(adresults)[1],1))

cpt <- 0
for (i in 1:nrow(adresults)) {
  cpt <- cpt + 1 
  r <- adresults[i,] # time series
  d <- which(is.finite(r) == T) # Locations NA
  diff_d <- d[2:length(d)] - d[1:length(d) - 1] # number of weeks btw NA
  
  if (length(d) >= round(nweeks/2) & length(which(diff_d >= 5)) == 0) { # at least half the weeks have data and a data point every 4 weeks
    row_value <- approx(1:nweeks, r, xout = 1:nweeks, rule = 2)$y # linearly interpolate
    
    #---- running average
    lag <- 3 # corresponds to 3 weeks 
    filter_coeff <- rep((1/lag),lag) # Coefficient pour pondérer les points
    row_value <- filter(row_value,filter_coeff, method = c("convolution"),sides = 2, circular = T) # Application du lissage (smooth)
    
    maxo_row[cpt] <- max(row_value) # save max value
    adresults_norm[cpt,] <- row_value/maxo_row[cpt] # normalization (divide by max)
  }
}

adresults <- adresults_norm # save data for clustering 
good_points <- which(is.na(adresults[,1]) == F) # pixels without NA to be clusterized
bad_points <- which(is.na(adresults[,1]) == T) # pixels with NA
xo[bad_points] <- NA
yo[bad_points] <- NA
maxo <- maxo_row

timeseries_annual <- adresults[good_points,]

# ----------------- distance matrix ------------------

dist_timeseries <- rdist(timeseries_annual, timeseries) # distance matrix
grps <- cluster
c = apply(dist_timeseries, 1, function(x) {b = aggregate(x, list(cluster), mean, simplify = T) # gives matrix six clusters by pixels with mean distance to every point within one cluster
  return(b[,2])})

# ----------------- silhouette value -----------------
#
# new cl gives a matrix that gives the cluster int which the min dist for every pixel

new_cl <- apply(c, 2, which.min) # finds the min distance of the column (column rep a pixel and row being a cluster)
a <- apply(c, 2, min) # gives the actual min dist value
b <- apply(c, 2, function(x) min(x[c(1:6) != which.min(x)])) # finds the second min value 
si <- (b-a)/pmax(a,b)


