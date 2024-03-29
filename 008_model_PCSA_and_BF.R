#-------------------------------------------------------------------------------

# From data collected in dissected specimens and 3D reconstructed specimens
# and values computed in previous scripts (e.g. muscle stress value),
# Compute individual values for individual muscle PCSAs and bite forces

# Required libraries
library(car)
library(MASS)
library(scales)

#-------------------------------------------------------------------------------
# Open and load data

dat <- read.csv("measurments_Sch.csv", 
                h = T, 
                dec = ".", 
                sep = ",") 
# This is the main data table, containing experimentally derived values

head(dat) # Check it is correctly loaded
str(dat)

dat_dis <- dat[which(dat$Dissected==1),] # Dissected subset

dat_3D <- dat[which(dat$Dissected==0),] # 3D reconstructed subset

load("fiber_lengths_and_summary_stats.RData")
load("muscle_stress.RData")
load("mat_fiber_lgt_angle.RData")

# These are useful R objects, saved from previous scripts

#-------------------------------------------------------------------------------

# Computation of PCSAs

#-------------------------------------------------------------------------------

# This dataset has two sub-sets: 
# Data acquired from dissected individuals,
# Data acquired from 3D reconstructed individuals.
# Because of this, some variables could only be derived from one subset, e.g.
# Mechanical Advantage, Muscle Stress, Pennation Angle.
# Therefore these values will be averaged for each side (left and right), 
# and the side averages used for the entire dataset

avcosL <-  mean(dat$cos_ang_L, 
                na.rm = T) 

avcosR <-  mean(dat$cos_ang_R, 
                na.rm = T) # cosines of pennation angles

avMAL <- mean(dat$MA_L, 
              na.rm = T)

avMAR <- mean(dat$MA_R, 
              na.rm = T) # Mechanical advantages

#-------------------------------------------------------------------------------
# PCSA formula for DISSECTION data:
# PCSA = (Muscle_mass * cos(pennation angle)) / (muscle density * fiber_length)

rho <- 0.00106 # Muscle density in g/mm3
  # This muscle density value is a universal standard for vertebrates. However,
  # it may be different for arthropods. We therefore ran experiments to check
  # whether muscle density might differ in Schistocerca. Our results, reported
  # in supplementary information, show that muscle density does not greatly 
  # differ from that of vertebrate muscle, meaning it is fairly close to that of
  # water.

PCSA_L <- (dat_dis$L_closer_g *       # Individual muscle mass in grams
          avcosL) /                   # average cos pennation angle
          (rho *                      # muscle density
          dat_dis$L_closer_fib_l_mm)  # Individual-muscle average fiber length

PCSA_R <- (dat_dis$R_closer_g * avcosR) / (rho * dat_dis$R_closer_fib_l_mm)

mat_PCSA <- cbind(PCSA_L, PCSA_R)

rownames(mat_PCSA) <- dat_dis$Specimen
  # This matrix is then pasted in the main data csv table, "measurments_Sch.csv"
  # for the corresponding specimens.

#-------------------------------------------------------------------------------
# PCSA formula for !!!3D!!! data:
# PCSA = (Muscle_volume * cos(pennation angle)) / fiber_length

# In this case, several different volume estimates were used:
# Direct volume estimate from 3D surface object
# Convex hull volume for the 3D object.
# So-called 2D-3D convex hull = stack of 2D convex hulls

vol_L <- dat_3D$L_closer_vol_mm3
vol_R <- dat_3D$R_closer_vol_mm3

CH_L <- dat_3D$L_convexhull_vol_mm3
CH_R <- dat_3D$R_convexhull_vol_mm3

CH2D3D_L <- dat_3D$L_2d3d_CH
CH2D3D_R <- dat_3D$R_2d3d_CH

flL <- dat_3D$L_closer_fib_l_mm
flR <- dat_3D$R_closer_fib_l_mm

#-------------------------------------------------------------------------------
# Compute PCSA for different volume estimates

PCSA_vol_L <- (vol_L * avcosL) / flL
PCSA_vol_R <- (vol_R * avcosR) / flR

PCSA_CH_L <- (CH_L * avcosL) / flL
PCSA_CH_R <- (CH_R * avcosR) / flR

PCSA_2D3D_L <- (CH2D3D_L * avcosL) / flL
PCSA_2D3D_R <- (CH2D3D_R * avcosR) / flR

mat_vols_PCSA <- cbind(PCSA_vol_L,
                       PCSA_vol_R,
                       PCSA_CH_L,
                       PCSA_CH_R,
                       PCSA_2D3D_L,
                       PCSA_2D3D_R)

rownames(mat_vols_PCSA) <- dat_3D$Specimen

#-------------------------------------------------------------------------------
# Save both matrices. (Data was also pasted into the main data table.)

save(list = c("mat_PCSA", "mat_vols_PCSA"), 
     file = "PCSA_matrices.RData")

#-------------------------------------------------------------------------------

# Computation of bite forces

#-------------------------------------------------------------------------------

# Formula to compute bite force from PCSA is:
# F_b = muscle_stress * PCSA * Mechanical_advantage
# Here we compute forces for individual muscles (left and right), then add them
# up later on, Each force is computed for two opening angle values (min and max)
# which influences the effective mechanical advantage.

#-------------------------------------------------------------------------------
# Define function to compute bite force

PCSA_2_force <- function(stress = stop("Muscle stress value must be defined"),
                         MA = 1, # By default, full force transfer
                         PCSA = stop("PCSA must be given")) {
  
  F_b <- stress * MA * PCSA
  
  return(F_b)
}

#-------------------------------------------------------------------------------
# Apply function


sigma <- pred_stress_N_cm2 / 100 # Convert to N/mm2

# First for closed mandibles, i.e. minimum opening angle, -assumed- to give the
# apodeme a 90° angle, and therefore in lever = effective in lever.

mat_F_L <- apply(X = mat_vols_PCSA[,c(1, 3, 5)], 
                 MARGIN = 2, 
                 FUN = PCSA_2_force, 
                 stress = sigma,
                 MA = avMAL)
  # Apply to LEFT muscle of 3D RECONSTRUCTED data

mat_F_R <- apply(X = mat_vols_PCSA[,c(2, 4, 6)], 
                 MARGIN = 2, 
                 FUN = PCSA_2_force, 
                 stress = sigma,
                 MA = avMAR)
  # Apply to RIGHT muscle of 3D RECONSTRUCTED data

F_PCSA_L <- PCSA_2_force(stress = sigma, 
                         MA = avMAL, 
                         PCSA = mat_PCSA[,1])
  # Apply to LEFT muscle of DISSECTED data


F_PCSA_R <- PCSA_2_force(stress = sigma, 
                         MA = avMAR, 
                         PCSA = mat_PCSA[,2])
  # Apply to RIGHT muscle of DISSECTED data


# Second, do the same for an opening angle of 85°, which makes the in-lever by
# approximately 40°, therefore reducing the effective in-lever. No difference
# for out-lever however.

effMAL <- avMAL * sin(40*pi/180)
effMAR <- avMAR * sin(40*pi/180)


mat_F_L_open <- apply(X = mat_vols_PCSA[,c(1, 3, 5)], 
                 MARGIN = 2, 
                 FUN = PCSA_2_force, 
                 stress = sigma,
                 MA = effMAL)
  # Apply to LEFT muscle of 3D RECONSTRUCTED data


mat_F_R_open <- apply(X = mat_vols_PCSA[,c(2, 4, 6)], 
                 MARGIN = 2, 
                 FUN = PCSA_2_force, 
                 stress = sigma,
                 MA = effMAR)
  # Apply to RIGHT muscle of 3D RECONSTRUCTED data

F_PCSA_L_open <- PCSA_2_force(stress = sigma, 
                         MA = effMAL, 
                         PCSA = mat_PCSA[,1])
  # Apply to LEFT muscle of DISSECTED data


F_PCSA_R_open <- PCSA_2_force(stress = sigma, 
                         MA = effMAR, 
                         PCSA = mat_PCSA[,2])
  # Apply to RIGHT muscle of DISSECTED data

muscle_forces_dissec <- data.frame(MF_closed_L = F_PCSA_L, 
                                   MF_closed_R = F_PCSA_R,
                                   MF_open_L = F_PCSA_L_open,
                                   MF_open_R = F_PCSA_R_open)
  # Make data frame for muscle force for dissected specimens

muscle_forces_VOL <- data.frame(MF_closed_L = mat_F_L[,1], 
                                 MF_closed_R = mat_F_R[,1],
                                 MF_open_L = mat_F_L_open[,1],
                                 MF_open_R = mat_F_R_open[,1])
  # Make data frame for muscle force based on muscle volume (3D)

muscle_forces_2D3D <- data.frame(MF_closed_L = mat_F_L[,3], 
                                 MF_closed_R = mat_F_R[,3],
                                 MF_open_L = mat_F_L_open[,3],
                                 MF_open_R = mat_F_R_open[,3])
  # Make data frame for muscle force based on muscle 2d3d convex hull (3D)

muscle_forces_CH <- data.frame(MF_closed_L = mat_F_L[,2], 
                               MF_closed_R = mat_F_R[,2],
                               MF_open_L = mat_F_L_open[,2],
                               MF_open_R = mat_F_R_open[,2])
  # Make data frame for muscle force based on muscle convex hull (3D)

#-------------------------------------------------------------------------------
# Calculate estimated bite forces (sum of muscle forces from both sides)

BF_closed_dissec <- muscle_forces_dissec[, 1] + muscle_forces_dissec[, 2]
BF_open_dissec <- muscle_forces_dissec[, 3] + muscle_forces_dissec[, 4]

BF_closed_VOL <- muscle_forces_VOL[, 1] + muscle_forces_VOL[, 2]
BF_open_VOL <- muscle_forces_VOL[, 3] + muscle_forces_VOL[, 4]

BF_closed_2D3D <- muscle_forces_2D3D[, 1] + muscle_forces_2D3D[, 2]
BF_open_2D3D <- muscle_forces_2D3D[, 3] + muscle_forces_2D3D[, 4]

BF_closed_CH <- muscle_forces_CH[, 1] + muscle_forces_CH[, 2]
BF_open_CH <- muscle_forces_CH[, 3] + muscle_forces_CH[, 4]

BF_closed_insertion <- PCSA_2_force(stress = sigma,
                            MA = mean(avMAR, avMAL),
                            PCSA = dat$insert_area_total)
  # The total insertion area includes both muscles (R and L), therefore total
  # muscle force is computed, rather than side by side. Average of left and 
  # right mechanical advantage is used.

BF_open_insertion <- PCSA_2_force(stress = sigma,
                                   MA = mean(effMAL, effMAR),
                                   PCSA = dat$insert_area_total)
#-------------------------------------------------------------------------------
# Save all forces in one R object

save(file = "all_BF.RData",
     list = c("BF_open_dissec", 
              "BF_closed_dissec",
              "BF_open_VOL",
              "BF_closed_VOL",
              "BF_open_2D3D",
              "BF_closed_2D3D",
              "BF_open_CH",
              "BF_closed_CH",
              "BF_open_insertion",
              "BF_closed_insertion"))
