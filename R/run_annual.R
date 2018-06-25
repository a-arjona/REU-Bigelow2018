rm(list = ls()) # clear all variables
JDstart <- 9 # starting week
JDend <- 35 # ending week
C <- c()
S <- c()
years <- 1998:2017
nclu <- 6 # number of clusters

require("fpc")
require("fields")

# --------------- Get Data from mean Climato once only ------------
# load time series data for mean climato and cluster number

setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Mean-Data/scripts/")

source("clustering_DINEOF.R")

mean_cluster <- cl$result$cluster # vector giving cluster number for every pixel in time series
mean_timeseries <- cl$adresults[cl$good_points,] # timeseries for each pixel. matrix pixel by week

rm(list=setdiff(ls(), c("JDstart","JDend", "C","S","years","nclu","mean_cluster","mean_timeseries")))

#-------------- Process and View Annual Data ---------------------
#
# set new working directory within which are the scripts and text files and figures will be saved here

setwd("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/")

for (yr in 1:length(years)) { # im realizing now there's a better way to do this but oh well. repeats below process for each year in the range
	c_year <- years[yr] # gives the year of the data file to be processed

	cyear_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/CHLGLOB_DINEOF",as.character(c_year), ".rdata",sep = "") # creates filename to be opened

	# Process annual data files to make figures
	#
	# variables needed to run : JDstart, JDend, cyear_file
	#
	# variables to keep from this file: 
	#		timeseries_annual,
	#		xo, yo, maxo, 
	#		good_points, lat, lon
	source("process_annual.R")

}