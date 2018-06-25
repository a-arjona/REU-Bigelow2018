# Plot silhouette analysis
#
#
# --------- order silhouette values -----------

X <- timeseries_annual

subset <- 1 # percentage

# subset the dataset # subset into what?
subindex <- c()
subX <- c()
subX_len <- c()
subcluster <- c()

for (ccc in 1:dim(c)[1]) { #1:6
  
  loc <- which(ccc == new_cl) # index of all from cluster ccc
  subloc <- c(loc)
  
  subX <- rbind(subX, X[loc,])
  
  subcluster <- c(subcluster, new_cl[subloc]) 
  subX_len <- c(subX_len, min(which(subcluster==ccc))) # creates a way to track the number of data points in each cluster
  subindex <- c(subindex, subloc)
  
}

# subindex_new <- match(1:length(subindex),subindex)


# order silhouette value
siX <- c()
past_l <- 1
sub_siX <- c()

for (l in 2:dim(c)[1]) {
  if (l == 6) {
    sub_siX <- tail(si, length(si)-subX_len[l])
    sub_siX <- sort(sub_siX, decreasing = T)
    siX <- c(siX, sub_siX)
  }
  else {
    sub_siX <- si[subX_len[past_l]:subX_len[l]-1]
    sub_siX <- sort(sub_siX, decreasing = T)
    past_l <- l
    siX <- c(siX, sub_siX)
  }
  
}

col <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf","#a65628","#999999","#000000")
colclus <- c(rep(col[1],time=length(which(subcluster == 1))), 
             rep(col[2],time=length(which(subcluster == 2))), 
             rep(col[3],time=length(which(subcluster == 3))), 
             rep(col[4],time=length(which(subcluster == 4))), 
             rep(col[5],time=length(which(subcluster == 5))), 
             rep(col[6],time=length(which(subcluster == 6))))

# name <- paste("SIL-",as.character(c_year),".png",sep = "")
#png(name, width = 1600, height = 1200, res = 150)
par(mfrow = c(1,1),mar=c(4,4,4,4))

barplot(siX, col=colclus, border=NA)
lines(c(0,length(siX)*1.25),c(mean(siX),mean(siX)),col = "black")
#dev.off()
