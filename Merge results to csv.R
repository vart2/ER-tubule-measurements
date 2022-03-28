##This script takes the individual tubule measurements (.csv) and RoiSet (.zip) files
##and merges them into one data frame based on slice number. Saves to "results.csv"
##To run this, press 'Source' and follow instructions in the Console to 
##enter working directory path and the number of tubules in the stack when requested. 

##Load packages. Use install.packages("pkg_name") in the console if first time. 
library(RImageJROI)
library(tidyverse)

wd <- readline("Enter axon working directory:")
setwd(wd)

no_slices <- as.numeric(readline("Number of slices in this stack:"))
results <- data.frame(slice=1:no_slices)

filenames_csv <- list.files(path=wd, pattern="*.csv")
filenames_roi <- list.files(path=wd, pattern="*.zip")

for(i in 1:length(filenames_csv)){

  feret <- read.csv(filenames_csv[i])["MinFeret"]
  if(nrow(feret)<10) {next}
  roiset <- read.ijzip(filenames_roi[i])
  roi_names <- names(lapply(roiset, '[[', "name"))
  slice <- substr(roi_names, 1,4)
  slice <- as.numeric(slice)
  feret <- cbind(feret, slice)
  colnames(feret)[1] <- str_extract(filenames_csv[i], "\\d+, \\d+")
  results <- left_join(results, feret, "slice")
  print(i)
}


write.csv(results, "results>9.csv")
