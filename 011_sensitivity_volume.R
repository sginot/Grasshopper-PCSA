#-------------------------------------------------------------------------------
# Sensitivity analysis of muscle 3D volume estimates for muscles

# In Fiji, each grayscale image stack for specimens Sch1-Sch15 (Sch2 removed)
# were masked using 2D3D convex hull stacks for individual muscles. This results
# in grayscale stacks for left and right muscles for each specimen. These 
# grayscale images are then binarized using various algorithms available in Fiji
# following which each slice in the stack has the mean pixel value measured, and
# when black pixels are present (i.e. mean pixel value =/= 0), all black pixels
# are selected, and the total area of selection is measured.

# This script then aims to read all tables for individual muscles, compute the
# muscle's volume (by adding up all areas mutliplied by pixel depth), and 
# combine all results into a table with muscles as lines and binarization
# algorithms as columns.

# Required packages:
library(viridis)

#-------------------------------------------------------------------------------
# Setup input folders and files and data frame to be filled

input_folder <- "output_sensitivity/"

input_subfolders <- list.files(input_folder)

input_files <- list.files(paste(input_folder, input_subfolders[1], sep=""))

mat_volumes <- matrix(NA, 
                      nrow = length(input_files),
                      ncol = length(input_subfolders),
                      dimnames = list(gsub(pattern = ".mhd.mha.csv", 
                                           replacement = "", 
                                           x = input_files),
                                      input_subfolders))

#-------------------------------------------------------------------------------
# Use double loop to go through all files in all subfolders, read them, 
# calculate volume based on voxel size

N_s <- length(input_subfolders)
N_f <- length(input_files)

for (i in 1:N_s) {
  
  fol <- input_subfolders[i]
  
  for (j in 1:N_f) {
    
    fil <- input_files[j]
    
    d <- read.csv(paste(input_folder,
                        fol, 
                        fil, 
                        sep = "/"),
                  h = T,
                  row.names = 1)
    
    voxel_char <- unlist(strsplit(input_files[j], 
                                  split = "_"))[3]
    
    voxel_char <- substr(x = voxel_char, 
                         start =1,
                         stop = 5)
    
    voxel <- as.numeric(gsub(pattern = ",",
                             replacement = ".",
                             x = voxel_char))
    
    voxel <- voxel / 1000
    
    index <- which(d[, 2] != 0)
    
    areas <- d[index, 1]
    
    volume <- sum(areas * voxel)
    
    mat_volumes[j, i] <- volume
  }
}

mat_volumes <- mat_volumes[c(1:2, 15:28, 3:14), ]
  # Order specimens in increasing order

for (i in seq(from = 1, 
              to = dim(mat_volumes)[1] - 1, 
              by = 2)) {
  
  if (mat_volumes[i, 1] > mat_volumes[i + 1, 1]) {
    
    mat_volumes[c(i, (i + 1)), ] <- mat_volumes[c((i + 1), i), ]
    
  }
}
  # Order side values so that the larger muscle is the right one (problem due
  # to reversed 3D reconstructions)

#-------------------------------------------------------------------------------
# Compute PCSA for the variously binarized muscles

dat <- read.csv("measurments_Sch.csv", 
                h = T, 
                dec = ".", 
                sep = ",") 
# This is the main data table, containing experimentally derived values

dat_3D <- dat[which(dat$Dissected==0),] # 3D reconstructed subset

load("fiber_lengths_and_summary_stats.RData")
load("muscle_stress.RData")
load("mat_fiber_lgt_angle.RData")

# These are useful R objects, saved from previous scripts

avcosL <-  mean(dat$cos_ang_L, 
                na.rm = T) 

avcosR <-  mean(dat$cos_ang_R, 
                na.rm = T) # cosines of pennation angles

avMAL <- mean(dat$MA_L, 
              na.rm = T)

avMAR <- mean(dat$MA_R, 
              na.rm = T) # Mechanical advantages

vol_L <- dat_3D$L_closer_vol_mm3[-2]
vol_R <- dat_3D$R_closer_vol_mm3[-2]
  #Sch2 removed here, so that the vol_L and vol_R vectors match the mat_volumes

vol <- c(rbind(vol_L, vol_R))

flL <- dat_3D$L_closer_fib_l_mm[-2]
flR <- dat_3D$R_closer_fib_l_mm[-2]

fl <- c(rbind(flL, flR))

mat_global <- cbind(mat_volumes, vol)

mat_global_PCSA <- mat_global

for (i in seq(from = 1, 
              to = dim(mat_global)[1] - 1, 
              by = 2)) {
  
  mat_global_PCSA[i, ] <- (mat_global[i, ] * avcosL) / fl[i]
  mat_global_PCSA[i + 1, ] <- (mat_global[i + 1, ] * avcosR) / fl[i + 1]
  
}

#-------------------------------------------------------------------------------
# Compute bite forces

sigma <- pred_stress_N_cm2 / 100 # Convert to N/mm2
  # Define muscle stress value

mat_global_forces <- mat_global_PCSA


for (i in seq(from = 1, 
              to = dim(mat_global)[1] - 1, 
              by = 2)) {
  
  mat_global_forces[i, ] <- mat_global_PCSA[i, ] * sigma * avMAL
  mat_global_forces[i + 1, ] <- mat_global_PCSA[i + 1, ] * sigma * avMAR
  
}

mat_bite_sensit <- matrix(NA, 
                          nrow = nrow(mat_global)/2,
                          ncol = ncol(mat_global),
                          dimnames = list(dat_3D$Specimen[-2],
                                          colnames(mat_global)))



for (i in 1:nrow(mat_bite_sensit)) {
  
  submat <- mat_global_forces[(i * 2 - 1):(i * 2),]
  
  mat_bite_sensit[i, ] <- apply(submat, 2, sum)

}

#-------------------------------------------------------------------------------
# Test and plot against in vivo bite force

# Combined plots of bite force estimates based on various binarization algo vs.
# In vivo measurements

cols <- viridis(ncol(mat_bite_sensit))

in_vivo <- dat_3D$maxBF_ampbasecorr[-2]

output_folder <- "Figures/"

pdf(file = paste(output_folder,
                 "sensitivity_binarization.pdf"),
    height = 7,
    width = 9)

layout(matrix(c(1:10, 0, 0), 
              nrow = 3,
              ncol = 4,
              byrow = T))

par(mar = c(4, 4, 1, 1))

for (i in 1:ncol(mat_bite_sensit)) {
  
  plot(mat_bite_sensit[, i], 
       in_vivo,
       pch = 21,
       bg = cols[i],
       lwd = 2,
       cex = 2.5,
       asp = 1,
       ylim = c(0, 2),
       xlim = c(0, 2),
       xlab = colnames(mat_bite_sensit)[i])
  
  abline(0, 1,
         lty = 2,
         lwd = 2)
  
}

dev.off()

# Test whether the estimates have some predictive power over in vivo data
mat_lm_sensit <- matrix(NA, 
                        nrow = ncol(mat_bite_sensit),
                        ncol = 4)

for (i in 1:ncol(mat_bite_sensit)) {
  
  slm <- summary(lm(in_vivo ~ mat_bite_sensit[, i]))
  
  mat_lm_sensit[i, ] <- slm$coefficients[2, ]
  
}

rownames(mat_lm_sensit) <- colnames(mat_bite_sensit)
colnames(mat_lm_sensit) <- colnames(slm$coefficients)

# Test slope against 1 instead of 0

mat_lm_sensit2 <- matrix(NA, 
                        nrow = ncol(mat_bite_sensit),
                        ncol = 4)

for (i in 1:ncol(mat_bite_sensit)) {
  
  slm <- summary(lm(in_vivo - mat_bite_sensit[, i] ~ mat_bite_sensit[, i]))
  
  mat_lm_sensit2[i, ] <- slm$coefficients[2, ]
  
}

rownames(mat_lm_sensit2) <- colnames(mat_bite_sensit)
colnames(mat_lm_sensit2) <- colnames(slm$coefficients)
  
write.csv(file = "results_sensitivity_binarization_test0.csv",
          x = mat_lm_sensit)


write.csv(file = "results_sensitivity_binarization_test1.csv",
          x = mat_lm_sensit2)
