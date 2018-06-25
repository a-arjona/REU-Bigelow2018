rm(list=ls())

library("dplyr") # for tidyverse commands

years <- 1998:2017 


for (yr in 1:length(years)){
  clean_data_frame <- data.frame()
  
  c_year <- years[yr] # gives current year
  cyear_file <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/CHLGLOB_DINEOF",as.character(c_year), ".rda",sep = "")
  
  load(cyear_file) # load data frame from current year
  
  npixels <- max(data$pixel) # number of pixels per week of data
  m_data <- filter(data, Mask == 1) # apply mask 
  
  for (p in 1:npixels) {
    data_sub_df <- filter(m_data, pixel == p) # focus on only data from pixel p
    data_sub <- data_sub_df$CHL_data # list chlorophyll data in vector form
    
    nweeks <- length(data_sub) # number of weeks in that year
    if (nweeks != 0) {
      d <- which(is.finite(data_sub) == T) # locations of non-NA 
      diff_d <- d[2:length(d)] - d[1:length(d) - 1] # list of # of weeks between NA points
      
      if (length(d) >= round(nweeks/2) & length(which(diff_d >= 5)) == 0) { # at least half the weeks have data and a data point every 4 weeks
        row_value <- approx(1:nweeks, data_sub, xout = 1:nweeks, rule = 2)$y # linearly interpolate
        
        #---- running average
        lag <- 3 # corresponds to 3 weeks 
        filter_coeff <- rep((1/lag),lag) # Coefficient pour pondérer les points
        row_value <- stats::filter(row_value,filter_coeff, method = c("convolution"),sides = 2, circular = T) # Application du lissage (smooth)
        
        temp_max <- max(row_value) # save max value
        norm_ts <-  row_value / temp_max # normalization (divide by max)
        
      }
      else {
        norm_ts <- rep(NA, length(data_sub))
      }
    clean_data <- mutate(data_sub_df, clean = norm_ts)
    clean_data_frame <- bind_rows(clean_data_frame, clean_data)
    }
  }
  new_name <- paste("C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/CLEAN_CHLGLOB_DINEOF",as.character(c_year), ".rda",sep = "")
  save(clean_data_frame, file = new_name)
}
    

