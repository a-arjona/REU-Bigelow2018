# main script 
# requires data from build_climato
# requires script from clustering_DINOEF
# requires script from plot_clustering
# requires script from silhouette

rm(list=ls())
past_cl <- NA
JDstart <- 9
JDend <- 35
C <- c()
S <- c()

require("fpc") # package use for kmeans function
require("fields") # package use for mapping
setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/scripts/") # the folder with the data

for (nclu in 2:9) { # tests different number of clusters
  
  source("clustering_DINEOF.R")
  
  rm(list=setdiff(ls(), c("cl","past_cl","nclu","JDstart","JDend","C","S","maxo","yo","xo","good_points"))) # delete object except the one useful
  
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
  
  if (nclu >= 3) {
    for (k in 1:(nclu-1)) {
      index <- which(past_cl == k)
      h     <- hist(cluster[index],0:(nclu),plot=F)
      counts_cl <-rbind(counts_cl,h$counts)
    }
    
    index <- max.col(counts_cl)
    for (k in 1:(nclu-1)) {
      new_cluster[cluster == index[k]] <- k 
      new_centers[k,] <- centers[index[k],]
    }
    
    new_centers[nclu,] <- centers[cluster[is.na(new_cluster)][1],]
    new_cluster[is.na(new_cluster)] <- nclu
    
    cluster <- new_cluster
    centers <- new_centers
  }
  
  past_cl <- cluster
  
  rm(list=setdiff(ls(), c("cluster","timeseries","good_points","yo","xo","lat","lon","past_cl","centers","nclu","JDstart","JDend","maxo","C","S"))) # delete object except the one useful
  
  source("plot_clustering_DINEOF.R")
  
  source("silhouette_analysis_DINEOF.R")
  
  C <- cbind(C, cluster)
  S <- cbind(S, siX[subindex])
  
  rm(list=setdiff(ls(), c("past_cl","nclu","JDstart","JDend","C","S","maxo","yo","xo","good_points","timeseries"))) # delete object except the one useful
  
}

write.table(C,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/clusters.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)

S <- format(S,digits=1, scientific=F)
write.table(S,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/si.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)

tbl <- cbind(format(maxo[good_points],digits=3), xo[good_points], yo[good_points], good_points)
write.table(tbl,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/metadata.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)

tbl <- cbind(timeseries)
write.table(format(tbl,digits=3),file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/timeseries.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)
