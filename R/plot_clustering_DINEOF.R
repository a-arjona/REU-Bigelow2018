# --- Time series & Map ------------
#
#
# plot the results of the clustering

col <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf","#a65628","#999999","#000000")

# name <- paste("TS-",as.character(nclu),".png",sep = "")
# png(name, width = 1600, height = 1200, res = 150)

layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(3,3,3,3))

NB_pixel_year <- length(cluster) # Nb of cluterized pixels
adresults <- timeseries # cluterized time series
nclu <- max(cluster)
nweeks <- dim(adresults)[2]

for (ccc in 1:nclu) {
  
  # Title with % of pixels per cluster
  NB_pixel_ccc <- length(which(cluster == ccc))
  text_lab <- paste(round((NB_pixel_ccc/NB_pixel_year)*100), " % des pixels",sep="")
  
  # Plot time series
  plot(centers[ccc,],col=col[ccc],ylim=c(0,1),type="o",xlab = text_lab) # centroid
  index <- which(cluster == ccc)
  res_sd <- array(dim=c(nweeks,1))
  
  for (n in 1:nweeks) {
    res_sd[n,1] <- sd(adresults[index,n]) # compute STD
  }
  
  lines(centers[ccc,]+res_sd,col=col[ccc]) # centroid + STD
  lines(centers[ccc,]-res_sd,col=col[ccc]) # centroid - STD
  
}

# dev.off()


clus <- array(dim=c(length(lon),length(lat)))

for (ttt in 1:length(good_points)) {
  index <- good_points[ttt]
  clus[xo[index],yo[index]]<-cluster[ttt]
}

# name <- paste("MAP-",as.character(nclu),".png",sep = "")
# png(name, width = 1600, height = 1200, res = 150)
par(mfrow = c(1,1),mar=c(4,4,4,4))

image.plot(lon,lat,clus,col=col[1:nclu])

# dev.off()