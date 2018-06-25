# --- Time series & Map ------------
#
#
# plot the results of the clustering
#name <- paste("TS-",as.character(c_year),".png",sep = "")
#png(name, width = 1600, height = 1200, res = 150)

layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(3,3,3,3))

NB_pixel_year <- length(new_cl) # Nb of cluterized pixels
adresults <- timeseries_annual # cluterized time series
nweeks <- length(JDstart:JDend)

for (ccc in 1:nclu) {
  
  # Title with % of pixels per cluster
  NB_pixel_ccc <- length(which(new_cl == ccc))
  text_lab <- paste(round((NB_pixel_ccc/NB_pixel_year)*100), " % des pixels",sep="")
  
  # Plot time series
  idx <- which(new_cl == ccc) # index of all pixels with same cluster
  mean_ts <- colMeans(timeseries_annual[idx,])
  
  plot(mean_ts,col=col[ccc],ylim=c(0,1),type="o",xlab = text_lab) # centroid
  res_sd <- array(dim=c(nweeks,1))
  
  for (n in 1:nweeks) {
    res_sd[n,1] <- sd(timeseries_annual[idx,]) # compute STD
  }
  
  lines(mean_ts+res_sd,col=col[ccc]) # centroid + STD
  lines(mean_ts-res_sd,col=col[ccc]) # centroid - STD
  
}

#dev.off()


clus <- array(dim=c(length(lon),length(lat)))

for (ttt in 1:dim(timeseries_annual)[1]) {
  idx <- good_points[ttt]
  clus[xo[idx],yo[idx]]<-new_cl[ttt]
}

# name <- paste("MAP-",as.character(c_year),".png",sep = "")
# png(name, width = 1600, height = 1200, res = 150)
par(mfrow = c(1,1),mar=c(4,4,4,4))

image.plot(lon,lat,clus,col=col[1:nclu])

#dev.off()