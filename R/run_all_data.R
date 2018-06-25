# main script
#
# 
# 

rm(list = ls())
past_cl <- NA
JDstart <- 9 # what week to start on
JDend <- 35 # what week to end on


years <- 1998:2017
C <- c()
S <- c()

require("fpc")
require("fields")
#setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/") # folder with the scripts being used

# ------------- Get Data from Mean Climato Once Only ---------------
# load time series data for mean_climato and cluster #
setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/scripts/")
for (nclu in 4:6) { # tests different number of clusters
  
  source("clustering_DINEOF.R")
  
  rm(list=setdiff(ls(), c("years","cl","past_cl","nclu","JDstart","JDend","C","S","maxo","yo","xo","good_points"))) # delete object except the one useful
  
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
  
  rm(list=setdiff(ls(), c("keep_clu","years","cluster","timeseries","good_points","yo","xo","lat","lon","past_cl","centers","nclu","JDstart","JDend","maxo","C","S"))) # delete object except the one useful
  
  source("plot_clustering_DINEOF.R")
  
  source("silhouette_analysis_DINEOF.R")
  
  C <- cbind(C, cluster)
  S <- cbind(S, siX[subindex])

  rm(list=setdiff(ls(), c("keep_clu","years","past_cl","nclu","JDstart","JDend","C","S","maxo","yo","xo","good_points","timeseries"))) # delete object except the one useful
  
}

write.table(C,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/clusters.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)

S <- format(S,digits=1, scientific=F)
write.table(S,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/si.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)

tbl <- cbind(format(maxo[good_points],digits=3), xo[good_points], yo[good_points], good_points)
write.table(tbl,file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/metadata.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)

tbl <- cbind(timeseries)
write.table(format(tbl,digits=3),file="C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/txt_files/timeseries.txt", sep = "\t", row.names=F,col.names = F,quote = FALSE)



# ------------ Get Data Matrix for one year ------------------
#
load("~/Bigelow/DINEOF-2018/array_index.rdata")

S <- array(NA, dim = c(20832, 20))
xo_s <- clu_array$xo
yo_s <- clu_array$yo

C <- array(NA, dim = c(20832,20))
xo_c <- clu_array$xo 
yo_c <- clu_array$yo  

cluster <- keep_clu 
# Gives the time series 
setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/")

for (yr in 1:length(years)){
  c_year <- years[yr]
  cyear_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/CHLGLOB_DINEOF",as.character(c_year), ".rdata",sep = "")
  
  # Plot sihouette value, time series for each cluster, map with clusters
  source("interpolation_DINEOF.R")
  
  C[good_points,yr] <- new_cl
  S[good_points,yr] <- si
  
  #source("annual_sil_analysis.R")
  
  #rm(list=setdiff(ls(), c("xo_s","yo_s","xo_c","yo_c","c_year","cluster","maxo","years","timeseries","new_cl","nclu","JDstart","JDend","C","S", "timeseries_annual", "si", "dist_timeseries", "c","col","good_points","adresults","lon","lat","xo","yo")))
  
  #source("plot_ts_clustering.R")
  
  rm(list=setdiff(ls(), c("cluster", "xo_s","yo_s","xo_c","yo_c","new_cl", "timeseries_annual","c_year","maxo","timeseries","past_cl","nclu","JDstart","JDend","C","S","years","xo","yo","good_points","si")))
  
  cluster_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/txt_files/", as.character(c_year), "clusters.txt", sep = "")
  write.table(new_cl,file=cluster_file, sep = "\t", row.names=F,col.names = F,quote = FALSE)
  
  si_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/txt_files/", as.character(c_year), "si.txt", sep = "")
  si <- format(si,digits=1, scientific=F)
  write.table(si,file=si_file, sep = "\t", row.names=F,col.names = F,quote = FALSE)
  
  metadata_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/txt_files/", as.character(c_year), "metadata.txt", sep = "")
  tbl <- cbind(format(maxo[good_points],digits=3), xo[good_points], yo[good_points], good_points)
  write.table(tbl,file=metadata_file, sep = "\t", row.names=F,col.names = F,quote = FALSE)
  
  ats_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/txt_files/", as.character(c_year), "timeseries.txt", sep = "")
  tbl <- cbind(timeseries_annual)
  write.table(format(tbl,digits=3),file=ats_file, sep = "\t", row.names=F,col.names = F,quote = FALSE)
  
  }

source("mode_cluster.R")



