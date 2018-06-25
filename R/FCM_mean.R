# ------------- Fuzzy C Means Analysis --------------
#
#
# Using 20-year mean climatology data
#
#
rm(list=ls())
past_cl <- NA
JDstart <- 9
JDend <- 35
C <- c()
S <- c()

require("fpc") # package use for kmeans function
require("fields") # package use for mapping
setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/scripts/") # the folder with the data

library("e1071")

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


ccl <- cmeans(adresults[good_points,],6,iter.max = 50, dist = "euclidean", m = 2)

# ---------- export text files -----------------

xo_n <- na.omit(xo)
yo_n <- na.omit(yo)
hard_cl <- ccl$cluster
membership_cl <- format(ccl$membership,digits=3, scientific=F)
mean_center <- ccl$centers

tbl <- cbind(membership_cl, xo_n, yo_n, hard_cl)
write.table(tbl,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/fuzzy_cl.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)

# ---------- Fuzzy C Means for Annual Data ----------

setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/")
years <- 1998:2017

load("~/Bigelow/DINEOF-2018/array_index.rdata")
M <- array(NA, dim = c(20832, 20))
xo_M <- clu_array$xo
yo_M <- clu_array$yo

for (yr in 1:length(years)){
  c_year <- years[yr]
  cyear_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/CHLGLOB_DINEOF",as.character(c_year), ".rdata",sep = "")
  
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
  
  ccl_a <- cmeans(adresults[good_points,], mean_center, iter.max = 50, dist = "euclidean", m = 2)
  
  M[,yr] <-
  xo_a <- na.omit(xo)
  yo_a <- na.omit(yo)
  hard_cl_a <- ccl_a$cluster
  membership_cl_a <- format(ccl_a$membership,digits=3, scientific=F)
  mean_center <- ccl_a$centers
  
  
  tbl <- cbind(membership_cl_a, xo_a, yo_a, hard_cl_a)
  filename <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/txt_files/", as.character(c_year), "FCM.txt", sep = "")
  write.table(tbl,file= filename, sep = "\t", row.names=F,col.names = F,quote = FALSE)
  
}

