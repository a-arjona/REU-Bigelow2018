# ----- Distance Matrix ---------------
#
#
# 
library("dplyr")
library("ggplot2")
# make matrix of time series data for one year (pixel by week)
ts_x <- c()
tempX <- c()
for (p in 1:max(clean_data_frame$pixel)) {
  sub_cleanX <- filter(clean_data_frame, pixel == p & week >= JDstart & week <= JDend)
  if (length(sub_cleanX$pixel) != 0) {
    sub_cleanv <- c(sub_cleanX$clean)
    tempX <- ts_x
    ts_x <- rbind(tempX, sub_cleanv)
  }
}

dist_ts <- rdist(ts_x, timeseries) # distance matrix
grps <- cluster
c = apply(dist_ts, 1, function(x) {b = aggregate(x, list(grps), mean, simplify = T) # gives matrix six clusters by pixels with mean distance to every point within one cluster
return(b[,2])})

# new cl gives a matrix that gives the cluster int which the min dist for every pixel

new_cl <- apply(c, 2, which.min)
new_cl <- as.numeric(as.character(new_cl))

# finds the min distance of the column (column rep a pixel and row being a cluster)
a <- apply(c, 2, min) # gives the actual min dist value
b <- apply(c, 2, function(x) min(x[c(1:5) != which.min(x)])) # finds the second min value 
si <- (b-a)/pmax(a,b)

new_cl_col <- c()
si_col <- c()

for (n in 1:length(new_cl)) {
  new_cl_temp <- rep(new_cl[n], 27)
  si_temp <- rep(si[n], 27)
  
  new_cl_col <- c(new_cl_col,new_cl_temp)
  si_col <- c(si_col, si_temp)
}

data_clu <- clean_data_frame %>% 
  filter(week >= JDstart & week <= JDend) %>%
  select(year, week, pixel, clean, Longitude, Latitude) %>%
  mutate(cluster = new_cl_col, si = si_col)

data_mean_clu <- data_clu %>%
  group_by(cluster,week) %>%
  summarize(mean_ts_clu = mean(clean), std = sd(clean)) %>%
  mutate(low = mean_ts_clu - std, high = mean_ts_clu + std)

col <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf")

# ----- Time Series -------

ggplot(subset(data_mean_clu, cluster %in% c(1,2,3,4,5,6)), aes(x = week, y = mean_ts_clu)) +
  geom_line(aes(color = factor(cluster))) +
  geom_line(linetype = "dashed", aes(x = week, y = low, color = factor(cluster))) +
  geom_line(linetype = "dashed", aes(x = week, y = high, color = factor(cluster))) +
  facet_wrap(~cluster) + 
  xlab("Weeks") + ylab("CHL") + ggtitle("Time Series", "1998") +
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf"))

# ------- Silhouette Analysis ------

col <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf")
colclus <- c(rep(col[1],time=length(which(data_si$cluster == 1))), 
             rep(col[2],time=length(which(data_si$cluster == 2))), 
             rep(col[3],time=length(which(data_si$cluster == 3))), 
             rep(col[4],time=length(which(data_si$cluster == 4))), 
             rep(col[5],time=length(which(data_si$cluster == 5))), 
             rep(col[6],time=length(which(data_si$cluster == 6))))

data_si <- data_clu %>%
  group_by(cluster, pixel) %>%
  summarize(siX = median(si)) %>%
  arrange(cluster, desc(siX))

barplot(data_si$siX,
        col=colclus, border = NA, main = "Silhouette Analysis for 1998")

# ----- Map ------------------------
rm(list = setdiff(ls(), c("data_clu","data_mean_clu","col","ts_x", "clean_data_frame")))
map_data <- data_clu %>% 
  group_by(pixel) %>%
  summarize(clu = mean(cluster),lat = mean(Latitude), lon = mean(Longitude))

lonV <- map_data %>% 
  distinct(lon)  
lonV <- lonV$lon

latV <- map_data %>%
  distinct(lat)
latV <- latV$lat

clus <- array(dim=c(length(lonV),length(latV)))

for (ttt in 1:dim(ts_x)[1]) {
  idx_xo <- which(lonV == map_data$lon[ttt])
  idx_yo <- which(latV == map_data$lat[ttt])
  clus[idx_xo,idx_yo] <- map_data$clu[ttt]
}

# name <- paste("MAP-",as.character(c_year),".png",sep = "")
# png(name, width = 1600, height = 1200, res = 150)
par(mfrow = c(1,1),mar=c(4,4,4,4))

image.plot(clus, col=col[1:6])

# dev.off()


