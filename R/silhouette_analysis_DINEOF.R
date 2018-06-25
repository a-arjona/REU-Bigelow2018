X <- timeseries

subset <- 1 # percentage

# subset the dataset
subindex <- c()
subX <- c()
subcluster <- c()

for (ccc in 1:dim(centers)[1]) {
  
  loc <- which(ccc == cluster)  # index of all within cluster ccc
  subloc <- sample(loc, round(subset*length(loc)), replace=F) # why? randomize loc matrix
  
  subX <- rbind(subX, X[subloc,]) # matrix with all the time series values pixel x 27 weeks for one cluster
  subcluster <- c(subcluster, cluster[subloc]) # cluster number for all pixels in subloc (should be the same ie. 1 1 1 1 1)
  subindex <- c(subindex, subloc) # all the locations (indexes) of the pixels belonging to cluster ccc
  
}

subindex <- match(1:length(subindex),subindex) # reindex matrix

# silhouette value
distX = rdist(subX)
siX = c()

for (l in 1:dim(subX)[1]) {
  
  value <- aggregate(distX[l,],list(subcluster),mean,simplify = T)
  a <- value[subcluster[l],2]
  b <- min(value[which(subcluster[l] != c(1:dim(centers)[1])),2])
  s <- (b - a)/max(a,b)
  
  siX <- c(siX,s)
}

# order silhouette value
siX_new <- c()
subindex_new <- c()
subX_new <- c()

for (ccc in 1:dim(centers)[1]) {
  
  idx <- which(ccc == subcluster)
  idx <- idx[order(siX[idx],decreasing = T)]
  
  subX_new <- rbind(subX_new,subX[idx,])
  subindex_new <- c(subindex_new,subindex[idx])
  siX_new <- c(siX_new, siX[idx])
  
}

col <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf","#a65628","#999999","#000000")
colclus <- c(rep(col[1],time=length(which(subcluster == 1))), 
             rep(col[2],time=length(which(subcluster == 2))), 
             rep(col[3],time=length(which(subcluster == 3))), 
             rep(col[4],time=length(which(subcluster == 4))), 
             rep(col[5],time=length(which(subcluster == 5))), 
             rep(col[6],time=length(which(subcluster == 6))), 
             rep(col[7],time=length(which(subcluster == 7))), 
             rep(col[8],time=length(which(subcluster == 8))), 
             rep(col[9],time=length(which(subcluster == 9)))) 



# name <- paste("SIL-",as.character(nclu),".png",sep = "")
# png(name, width = 1600, height = 1200, res = 150)
par(mfrow = c(1,1),mar=c(4,4,4,4))

barplot(siX_new, col=colclus, border=NA)
lines(c(0,length(siX_new)*1.25),c(mean(siX_new),mean(siX_new)),col = "black")

# dev.off()