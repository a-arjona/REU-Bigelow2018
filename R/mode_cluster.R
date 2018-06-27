# ------ Mode Cluster per Pixel ----------
xo_p <- xo_s
yo_p <- yo_s

mclusters <- array(NA, dim = c(dim(C)[1],1))
perc_assoc <- array(NA, dim = c(dim(C)[1],1))

for (r in 1:dim(mclusters)[1]) {
  if (length(which(is.finite(C[r,]) == T)) >= dim(C)[2]/2) {
    row <- C[r,]
    uniqv <- unique(row)
    mclusters[r] <- uniqv[which.max(tabulate(match(row,uniqv)))]
    perc_assoc[r] <- length(which(C[r,] == mclusters[r])) / dim(C)[2]
  }
}

bad_points <- which(is.na(perc_assoc) == T) # pixels with NA
perc_assoc <- na.omit(perc_assoc)
xo_p[bad_points] <- NA
yo_p[bad_points] <- NA
xo_p <- na.omit(xo_p)
yo_p <- na.omit(yo_p)

perc_assoc_file <- c("./txt_files/perc_assoc.txt")
tbl <- cbind(perc_assoc, xo_p, yo_p)
write.table(tbl,file = perc_assoc_file, sep = "\t", row.names = F,col.names = F,quote = FALSE)

bad_points <- which(is.na(mclusters) == T) # pixels with NA
xo_c[bad_points] <- NA
yo_c[bad_points] <- NA

xo_c <- na.omit(xo_c)
yo_c <- na.omit(yo_c)
mclusters <- na.omit(mclusters)

mclu_file <- paste("./txt_files/m_clusters.txt", sep = "")
tbl <- cbind(mclusters, xo_c, yo_c)
write.table(tbl,file = mclu_file, sep = "\t", row.names = F,col.names = F,quote = FALSE)

# ---------- Average si ------------
mean_si <- array(NA, dim = c(dim(S)[1],1))

for (r in 1:dim(S)[1]) {
  row <- S[r,]
  if (length(which(is.finite(row) == T)) >= dim(S)[2]/2) {
    mean_si[r] <- mean(row)
  }
}

bad_points <- which(is.na(mean_si) == T) # pixels with NA

xo_s[bad_points] <- NA
yo_s[bad_points] <- NA

mean_si <- na.omit(mean_si)
xo_s <- na.omit(xo_s)
yo_s <- na.omit(yo_s)

meansi_file <- paste("./txt_files/mean_si.txt", sep = "")
tbl <- cbind(mean_si, xo_s, yo_s)
write.table(tbl,file = meansi_file, sep = "\t", row.names = F,col.names = F,quote = FALSE)
