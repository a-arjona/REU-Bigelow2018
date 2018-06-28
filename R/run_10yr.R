rm(list = ls())
past_cl <- NA
JDstart <- 9
JDend <- 35
nclu <- 6
setwd("C://Files/Work/Bigelow/Data/")

require("fpc") # package use for kmeans function
require("fields") # package use for mapping
library("ppclust")
library("e1071")

for (yr_set in 1:2) {
  C <- c()
  S <- c()
  if (yr_set == 1) {
    data_file <- "./DINEOF_2018_processed_data/CHLGLOB_DINEOF19982007.rdata"
  }
  if (yr_set == 2) {
    data_file <- "./DINEOF_2018_processed_data/CHLGLOB_DINEOF20072017.rdata"
  }
  
  # -------- Interpolation and Normalization ----------
  
  # LOAD DATA
  
  load(file = data_file) # load data from build climato
  
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
  
  # INTERPOLATE AND NORMALIZE
  
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
  
  # ---------------- Clustering ------------
  
  # FUZZY C MEANS 
  
  ccl <- cmeans(adresults[good_points,],6, iter.max = 50, dist = "euclidean", m = 2)
  
  # export text files
  xo_n <- na.omit(xo)
  yo_n <- na.omit(yo)
  hard_cl <- ccl$cluster
  membership_cl <- format(ccl$membership,digits = 3, scientific = F)
  mean_center <- ccl$centers
  
  tbl <- cbind(membership_cl, xo_n, yo_n, hard_cl)
  fname <- paste("./txt_files/10yr_FCM", yr_set, ".txt", sep = "")
  write.table(tbl,file = fname, sep = "\t", row.names = F,col.names = F,quote = FALSE)
  
  # K MEANS CLUSTER
  
  cl <- kmeansCBI(adresults[good_points,],k = nclu,scaling = F,runs = 50) # clustering Kmeans, runs: Number of starts of the k-means. scaling: If scaling is TRUE then scaling is done by dividing the variables by their root-mean-square
  
  # save results
  cl$good_points <- good_points
  cl$lon <- lon
  cl$lat <- lat
  cl$xo <- xo
  cl$yo <- yo
  cl$adresults <- adresults
  cl$maxo <- maxo_row
  
  # -------- Organize variables and Graph ----------
  cluster     <- cl$result$cluster
  centers     <- cl$result$centers
  timeseries  <- cl$adresults[cl$good_points,]
  good_points <- cl$good_points
  maxo        <- cl$maxo
  xo          <- cl$xo
  yo          <- cl$yo
  lon         <- cl$lon
  lat         <- cl$lat
  
  new_cluster <- array(NA, dim = length(cluster))
  new_centers <- array(NA, dim = dim(centers))
  counts_cl <- c()
  if (nclu == 6) {
    keep_clu <- cluster
  }
  if (nclu >= 3) {
    for (k in 1:(nclu - 1)) {
      index <- which(past_cl == k)
      h     <- hist(cluster[index],0:(nclu),plot = F)
      counts_cl <- rbind(counts_cl,h$counts)
    }
    
    index <- max.col(counts_cl)
    for (k in 1:(nclu - 1)) {
      new_cluster[cluster == index[k]] <- k 
      new_centers[k,] <- centers[index[k],]
    }
    
    new_centers[nclu,] <- centers[cluster[is.na(new_cluster)][1],]
    new_cluster[is.na(new_cluster)] <- nclu
    
    cluster <- new_cluster
    centers <- new_centers
  }
  
  past_cl <- cluster
  
  rm(list = setdiff(ls(), c("yr_set","keep_clu","years","cluster","timeseries","good_points","yo","xo","lat","lon","past_cl","centers","nclu","JDstart","JDend","maxo","C","S"))) # delete object except the one useful
  
  # ------- Make Plots ------------

  #source("plot_clustering_DINEOF.R")
  
  source("./Scripts/R/silhouette_analysis_DINEOF.R")
  
  C <- cbind(C, cluster)
  S <- cbind(S, siX[subindex])
  
  rm(list = setdiff(ls(), c("yr_set","keep_clu","years","past_cl","nclu","JDstart","JDend","C","S","maxo","yo","xo","good_points","timeseries"))) # delete object except the one useful
  
  # ------ Save Text Files --------
  fname <- paste("./txt_files/10yr_clusters", yr_set, ".txt", sep = "")
  write.table(C,file = fname, sep = "\t", row.names = F,col.names = F,quote = FALSE)
  
  fname <- paste("./txt_files/10yr_si", yr_set, ".txt", sep = "")
  S <- format(S,digits = 1, scientific = F)
  write.table(S,file = fname, sep = "\t", row.names = F,col.names = F,quote = FALSE)
  
  fname <- paste("./txt_files/10yr_metadata", yr_set, ".txt", sep = "")
  tbl <- cbind(format(maxo[good_points],digits = 3), xo[good_points], yo[good_points], good_points)
  write.table(tbl,file = fname, sep = "\t", row.names = F,col.names = F,quote = FALSE)

  fname <- paste("./txt_files/10yr_timeseries", yr_set, ".txt", sep = "")
  tbl <- cbind(timeseries)
  write.table(format(tbl,digits = 3),file = fname, sep = "\t", row.names = F,col.names = F,quote = FALSE)
  
  }

