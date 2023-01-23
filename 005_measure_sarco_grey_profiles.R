#-------------------------------------------------------------------------------

# Extraction of sarcomere length from grey value profiles from Fiji

#-------------------------------------------------------------------------------

# Grey value profiles where extracted from fiji by drawing a line
# perpendicular to the sarcomere bands, then Analyze > Plot profile

#-------------------------------------------------------------------------------
# Preparations

# Define files of raw data

input_folder <- "../raw_data_sarcomere/"

fil <- list.files(path = input_folder,
                  pattern = "GRAY_PROFILE",
                  full.names = T)

# Make matrix and list to hold the results
# For each profile, record mean, sd, N, max and min sarcormere length

mat <- matrix(NA, nrow = length(fil), ncol = 5)

ls_lengths <- list()

#-------------------------------------------------------------------------------
# Interactive loop for obtaining the data, one must click on graphics
dev.off()
layout(1)

for (i in 1:length(fil)) {
  dat_gray <- read.csv(fil[i],
                       h = T,
                       dec = ".",
                       sep = ",")
  
  plot(dat_gray[,1], 
       dat_gray[,2],
       type = "l")
  
  coords <- locator() # One must click on the repeated max gray value peaks
                      # on the graph ! Then "Finish" in RStudio or right click.
  lengths <- diff(coords$x, lag = 1)
  
  ls_lengths[[i]] <- lengths
  
  mat[i,] <- c(mean(lengths), 
               sd(lengths), 
               length(lengths),
               min(lengths), 
               max(lengths))
}

rownames(mat) <- names(ls_lengths) <- gsub(fil, 
                                           pattern = input_folder, 
                                           replacement = "")

colnames(mat) <- c("Mean", "SD", "N", "Minimum", "Maximum")

#-------------------------------------------------------------------------------
# Save data as separate files

save(mat, file = "mat_sarcomere_lengths_summary.RData")

save(ls_lengths, file = "list_sarcomere_length_measurements.RData")

# Optional, save as a CSV file
write.csv(mat, file = "mat_sarcomere_lengths_summary.csv")

#-------------------------------------------------------------------------------
# Plot data to check if any outliers may be present (e.g. wrong click)

lab <- rep(NA, length(fil))

lab[which(grepl(fil, pattern = "LEFT") & 
          grepl(fil, pattern = "CT"))] <- "LEFT CT"

lab[which(grepl(fil, pattern = "RIGHT") & 
          grepl(fil, pattern = "CT"))] <- "RIGHT CT"

lab[which(grepl(fil, pattern = "LEFT") & 
          !grepl(fil, pattern = "CT"))] <- "LEFT Axio"

lab[which(grepl(fil, pattern = "RIGHT") & 
          !grepl(fil, pattern = "CT"))] <- "RIGHT Axio"

boxplot(ls_lengths, las = 2, names = lab)

boxplot(unlist(ls_lengths[which(grepl(fil, pattern = "LEFT"))]),
        unlist(ls_lengths[which(grepl(fil, pattern = "RIGHT"))]),
        names = c("LEFT", "RIGHT"))

boxplot(unlist(ls_lengths[which(grepl(fil, pattern = "CT"))]),
        unlist(ls_lengths[which(!grepl(fil, pattern = "CT"))]),
        names = c("CT scan", "Axio Zoom"))

avg_global <- mean(unlist(ls_lengths))

hist(unlist(ls_lengths))
abline(v = avg_global)
