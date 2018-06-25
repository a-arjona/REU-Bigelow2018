# Process annual data files to make figures
#
#
# variables needed to run : JDstart, JDend, cyear_file
#
# variables to keep from this file: 
#		timeseries_annual,
#		xo, yo, maxo, 
#		good_points, lat, lon
#
# ------------- Open Data -----------------

load(file = cyear_file) # load annual data file
# annual data file contains Large List called "data" with 8 elements
#		data double[(# of pixels (20832)) x [# of weeks (46)]] raw data value (CHL either NA or value)
#		lon double [dimensions of area of interest (248 x 84)] lon values ie. -43.9 -43.6
#		lat double [dimensions of area of interest (248 x 84)] lat values ie. 60.1 60.4
# 		xo integer [# of pixels (20832)] gives the index of a pixel for the longitude vector
# 		yo integer [# of pixels (20832)] gives the index of a pixel for the latitude vector
#		LON double [# of pixels (20832)] lon values in a vector
#		LAT double [# of pixels (20832)] lat values in a vector
# 		mask integer [# of pixels (20832)] tells which pixels to keep and which data points to trow out ie. 0 0 0 1 1 0 1

# assign these elements to easy access variables
lat <- data$lat[1,] # latitude in 'vector' [1:248]
lon <- data$lon[,1] # longitude in 'vector' [1:84]
xdim <- length(lon) 
ydim <- length(lat)
xo <- data$xo
yo <- data$yo
new_data <- data$data[,JDstart:JDend] # data in between weeks 9 and 35 (majority of where the data is)
nweeks <- dim(new_data)[2] # number of weeks in new data (~27)
new_data[which(data$mask == 0),] <- NA # apply mask. all pixels removed by mask now read NA

# ------------- Interpolate and Normalize --------------
#
#
# Linear interpolation to lessen NA points and divide each time series by its maximal value

new_data_norm <- array(dim = dim(new_data)) # creates empty matrix with same dimensions as new_data (pixels (20832) by weeks (27))
maxo_row <- array(dim = c(dim(new_data)[1],1)) # creates empty vector to store max value of time series for each pixel
mino_row <- array(dim = c(dim(new_data)[1],1)) # creates empty vector to store min value of time series for each pixel

cpt <- 0 # counter variable?
for (i in 1:nrow(new_data)) {
	cpt <- cpt + 1 # increment counter
	r <- new_data[i,] # row for pixel i from new_data -- time series
	d <- which(is.finite(r)==T) # locations of pixels that are not NA
	diff_d <- d[2:length(d)] - d[1:length(d) - 1] # number of weeks between NA in vector form

	if (length(d) >= round(nweeks/2) & length(which(diff_d >= 5)) == 0) { # if at least half the weeks have data and there are no gaps wider than 4 weeks
		row_value <- approx(1:nweeks, r, xout = 1:nweeks, rule = 2)$y # linearly interpolate then only store time series value in vector

		# ------- running average -------------
		# smooth annual data
		lag <- 3 # corresponds to 3 weeks
		filter_coeff <- rep((1/lag),lag)  # coefficient pour pondarer les points
		row_value <- filter(row_value, filter_coeff, method = c("convolution"), sides = 2, circular = T) # apply running average to smooth out data

		maxo_row[cpt] <- max(row_value) # save max value in index indicated by cpt counter which correlates to pixel #
		new_data_norm[cpt,] <- row_value/maxo_row[cpt] # normalization (divide by max) then store in new_data_norm matrix
	}
}

# ----------- save data in easy access variables ---------
processed_data <- new_data_norm # save data for plotting
good_points <- which(is.na(processed_data[,1]) == F) # gives a vector with the index of all pixels without NA
bad_points <- which(is.na(processed_data[,1]) == T) # gives a vector with the index of pixels with NA
xo[bad_points] <- NA # replaces all bad point pixels with NA
yo[bad_points] <- NA # ^^
maxo <- maxo_row # vector with max value of time series for every pixel

timeseries_annual <- processed_data[good_points,]