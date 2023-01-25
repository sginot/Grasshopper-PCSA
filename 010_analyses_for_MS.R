#-------------------------------------------------------------------------------
# This script produces all analyses of the manuscript, as well as some shown
# in the supplementary information.

# Required libraries:
library(car)
library(MASS)
library(scales)

#-------------------------------------------------------------------------------
# Load all required data

dat <- read.csv("../measurments_Sch.csv", 
                h = T, 
                dec = ".", 
                sep = ",") 
# This is the main data table, containing experimentally derived values

dat_dis <- dat[which(dat$Dissected==1),] # Dissected subset

dat_3D <- dat[which(dat$Dissected==0),] # 3D reconstructed subset

load("fiber_lengths_and_summary_stats.RData")
load("muscle_stress.RData")
load("mat_fiber_lgt_angle.RData")
load("all_BF.RData")
load("fiber_lengths_3D.RData")
load("models.RData")
load("PCSA_matrices.RData")
load("muscle_stress.RData")
load("pennation_angles.RData")

#-------------------------------------------------------------------------------
# t tests of muscular traits and bite forces

#-------------------------------------------------------------------------------
# Linear models of estimated vs. in vivo bite force.

#-------------------------------------------------------------------------------
# Sensitivity analysis???
