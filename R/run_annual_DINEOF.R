# main script
#
# 
# 

rm(list = ls())
past_cl <- NA
JDstart <- 9 # what week to start on
JDend <- 35 # what week to end on


years <- 1998:2017


require("fpc")
require("fields")
#setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/") # folder with the scripts being used

# ------------- Get Data from Mean Climato Once Only ---------------
# load time series data for mean_climato and cluster #

nclu <- 6 # number of clusters in mean climato data
setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/scripts/")

source("clustering_DINEOF.R")

rm(list = setdiff(ls(), c("cl","past_cl","nclu","JDstart","JDend","C","S","maxo","yo","xo","good_points","years"))) # clear variables except ones still necessary

cluster     <- cl$result$cluster # cluster number (integer)
timeseries  <- cl$adresults[cl$good_points,] # time series for each pixel in a matrix wk by pixel
xo_mean     <- na.omit(cl$xo)
yo_mean     <- na.omit(cl$yo)

rm(list=setdiff(ls(), c("nclu","JDstart","JDend","C","S","years","timeseries","cluster"))) # clear all variables except necessary

# ------------ Get Data Matrix for one year ------------------
#
load("~/Bigelow/DINEOF-2018/array_index.rdata")

S <- array(NA, dim = c(20832, length(years)))
xo_s <- clu_array$xo
yo_s <- clu_array$yo

C <- array(NA, dim = c(20832,length(years)))
xo_c <- clu_array$xo 
yo_c <- clu_array$yo  

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
  
  #rm(list=setdiff(ls(), c("xo_c","yo_c","c_year","cluster","maxo","years","timeseries","new_cl","nclu","JDstart","JDend","C","S", "timeseries_annual", "si", "dist_timeseries", "c","col","good_points","adresults","lon","lat","xo","yo")))
  
  #source("plot_ts_clustering.R")
  
  rm(list=setdiff(ls(), c("xo_c","yo_c","new_cl", "timeseries_annual","c_year","maxo","cluster","timeseries","past_cl","nclu","JDstart","JDend","C","S","years","xo","yo","good_points","si")))
  
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
mclu_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/m_clusters.txt", sep = "")
tbl <- cbind(mclusters, xo_c, yo_c)
write.table(tbl,file=mclu_file, sep = "\t", row.names=F,col.names = F,quote = FALSE)


