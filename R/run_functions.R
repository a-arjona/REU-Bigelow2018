rm(list = ls())

setwd("C:/Files/Work/Bigelow/Data/")
source("./Scripts/R/functions.R")

mean20Climato <- buildData(1998, 2017)
mean9807 <- buildData(1998, 2007)
mean0817 <- builData(2008, 2017)
