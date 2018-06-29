setwd("C:/Files/Work/Bigelow/Data/")
source("./Scripts/R/functions.R")

# ------- 1998 - 2007 ------
#run kmeans analysis a chosen number of times and save metadata, timeseries, si, cluster as .txt files 
z <- runK(1998,2007,filename = "nclu9808")
