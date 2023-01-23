#-------------------------------------------------------------------------------

# Convert, import and correct bite force curves and extract bite force data, 
# using ForceR.
# Here, we simply aim at extracting maximum bite force.


# Required libraries
require(devtools)
devtools::install_github("https://github.com/Peter-T-Ruehr/forceR")
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(readr)
library(forceR)

#-------------------------------------------------------------------------------
# Define input folder and files

input_folder <- "../bites/"

# No need to use the following two lines once files have been created.

# dat.files <- list.files(input_folder, pattern = ".dat")

# for (i in 1:length(dat.files)) {convert_measurement(file = dat.files[i], 
#                                 path.data = input_folder)}


data.folder <- input_folder # For correspondence with ForceR vignette

#-------------------------------------------------------------------------------
# Open csv files containing time series of electrical current, which represent
# the raw bite force measurements, and concatenate them for each individual.

csv.files <- list.files(data.folder, 
                        pattern = "converted.csv", 
                        full.names = T)

specimens <- paste("Sch", 1:30, "_", sep="")

for (spe in specimens) {
  a <- csv.files[grep(spe, csv.files)]
  b <- list()
  
  for (i in 1:length(a)) {
    b[[i]] <- read.table(a[i],
                         h = T,
                         sep= ",",
                         dec = ".")
  }
  
  mat <- do.call(rbind, b) 
  #This might not work if some measurements are only contained in one file
  
  colnames(mat) <- colnames(b[[1]])
  
  write.table(mat, 
              file = paste(data.folder, 
                           spe, 
                           "concat.csv",
                           sep = ""), 
              row.names = F,
              dec = ".",
              sep = ",")
}

#-------------------------------------------------------------------------------
# Read one concatenated file and check that everything is in order

concat.files <- list.files(data.folder, 
                           pattern = "concat.csv", 
                           full.names = F)

file <- file.path(data.folder, 
                  concat.files[1])

plot_measurement(file, 
                 columns = c(1:2))

#-------------------------------------------------------------------------------
# Concatenated series can and should be manually cropped, to time window in
# which actual bites were measured.

# /!!\ One must make a "cropped" subfolder manually on Hard drive

cropped.folder <- "../bites/cropped/"

for (i in 1:length(concat.files)) {
  file <- file.path(data.folder, 
                    concat.files[i])
  
  file.cropped <- crop_measurement(file, 
                                   path.data = cropped.folder)
} 
# First and last point of cropping must be defined manually by clicking along 
# the bite series curve

#-------------------------------------------------------------------------------
# Cropped files must then be corrected for amplification drift

cropped.files <- list.files(cropped.folder, 
                            pattern = "cropped.csv",
                            full.names = T)

layout(mat = matrix(1:length(cropped.files), 
                    ncol=5))
par(mar=c(1,1,1,1))

for (i in 1:length(cropped.files)) {
  
  file <- cropped.files[i]
  
  plot_measurement(file)
} # With this loop one can see that the baseline of the curves moves around

file.list <- list.files(cropped.folder, 
                        pattern = ".csv",
                        full.names = T)

# NB: one should create the "ampdriftcorr" folder manually.

ampdriftcorr.folder <- "../bites/cropped/ampdriftcorr/"

for(filename in file.list){
  print(filename)
  amp_drift_corr(filename = filename,
                 tau = 9400,
                 res.reduction = 10,
                 plot.to.screen = FALSE,
                 write.data = TRUE,
                 write.PDFs = TRUE,
                 write.logs = TRUE,
                 output.folder = ampdriftcorr.folder,
                 show.progress = FALSE)
  print("***********")
}

baselinecorr.folder <- "./cropped/ampdriftcorr/baselinecorr"

file.list <- list.files(ampdriftcorr.folder, 
                        pattern = ".csv",
                        full.names = TRUE)

for(filename in file.list){
  print(filename)
  file.baselinecorr <- baseline_corr(filename = filename, 
                                   corr.type = "auto",  
                                   window.size.mins = 2000,
                                   window.size.means = NULL,
                                   quantile.size = 0.05,
                                   y.scale = 0.5,
                                   res.reduction = 10,
                                   Hz = 100,
                                   plot.to.screen = TRUE,
                                   write.data = TRUE,
                                   write.PDFs = TRUE,
                                   write.logs = TRUE,
                                   output.folder = baselinecorr.folder,
                                   show.progress = FALSE)
  print("***********")
}

#-------------------------------------------------------------------------------
# Extract bite force data with and without corrections
#-------------------------------------------------------------------------------
# First for non-corrected data

input_folder <- "../bites/cropped/"

list.files <- list.files(input_folder, 
                         pattern = ".csv",
                         full.names = T)
ls.bf <- list() # Empty list

for (i in 1:length(list.files)){ # Fill list
  ls.bf[[i]] <- read.table(list.files[i], 
                           h=T, 
                           dec=".", 
                           sep=",")
}

maxBF <- rep(NA, length(list.files)) # Empty vector

for (i in 1:length(list.files)){ # Fill vector
  maxBF[i] <- max(ls.bf[[i]]$y)
}

maxBF_nocorr <- maxBF/2*0.5 # Raw data (volts) must be corrected for the
                            # amplification of amplifier device (2A), and 
                            # mechanical advantage of measuring device (0.5)

#-------------------------------------------------------------------------------
# Second for data with corrected amplification drift

input_folder <- "../bites/cropped/ampdriftcorr/"

list.files.corr <- list.files(input_folder, 
                              pattern = ".csv",
                              full.names = T)
ls.corr.bf <- list()

for (i in 1:length(list.files.corr)){
  ls.corr.bf[[i]] <- read.table(list.files.corr[i], 
                                h=T, 
                                dec=".", 
                                sep=",")
}

maxBF <- rep(NA, length(list.files.corr))

for (i in 1:length(list.files.corr)){
  maxBF[i] <- max(ls.corr.bf[[i]]$y)
}

maxBF_ampcorr <- maxBF/2*0.5

#-------------------------------------------------------------------------------
# Third for data with corrected amplification drift and baseline drift

input_folder <- "../bites/cropped/ampdriftcorr/baselinecorr/"

list.files.corr <- list.files(input_folder, 
                              pattern = ".csv",
                              full.names = T)
ls.corr.bf <- list()

for (i in 1:length(list.files.corr)){
  ls.corr.bf[[i]] <- read.table(list.files.corr[i], h=T, dec=".", sep=",")
}

maxBF <- rep(NA, length(list.files.corr))
for (i in 1:length(list.files.corr)){
  maxBF[i] <- max(ls.corr.bf[[i]]$y)
}

maxBF_ampbasecorr <- maxBF/2*0.5

#-------------------------------------------------------------------------------
# Make data frame out of extracted maximum bite forces and save that data

specimens <- gsub(x=list.files, 
                  pattern = "_concat_cropped.csv",
                  replacement = "")

specimens <- gsub(x=specimens, 
                  pattern = "../bites/cropped/",
                  replacement = "")

dfBF <- data.frame(specimens = specimens,
           maxBF_nocorr, 
           maxBF_ampcorr, 
           maxBF_ampbasecorr)
  # Data frame as is is not in the increasing order to to lack of 0 padding

dfBF <- dfBF[c(1, 12, 23, 25:30, 2:11, 13:22, 24),]
rownames(dfBF) <- 1:30
  # Re-order table with increasing specimen numbers

write.table(dfBF, file = "../bites/maxBF_clean.csv", row.names = F)
  # The table is then pasted in the main measurements csv files 
  # "measurments_Sch.csv", with Sch13 removed: outlier, no good bite measurement