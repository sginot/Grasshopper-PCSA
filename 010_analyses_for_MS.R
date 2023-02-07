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
# Make synthetic table for all traits measured, must contains, mean, SD, and N
# Globally and for left and right separately. And separately for 3D and dissec

BF <- dat$maxBF_ampbasecorr
  # In vivo bite forces

HW <- dat$HW
HH <- dat$HH
PW <- dat$PW
TL <- dat$TL
HL <- dat$HL
  # Linear measurements

Lfib3D <- unlist(list_lgt_fib[seq(1, 29, by = 2)])
Rfib3D <- unlist(list_lgt_fib[seq(2, 30, by = 2)])

revpenn <- c(1, 4, 6, 9:12) # Reversed side individuals

ord.pennat <- ls.pennat # Prepare list

ord.pennat[c(rbind(revpenn*2-1, 
                  revpenn*2))] <- ls.pennat[c(rbind(revpenn*2, 
                                                    revpenn*2-1))]
  # Reorder elements of the list

Lpennat3D <- Rpennat3D <- list()
  # Prepare empty lists

for (i in seq(1, 29, by = 2)) {
  
  Lpennat3D[[i]] <- ord.pennat[[i]][,1]
  Rpennat3D[[i]] <- ord.pennat[[i+1]][,1]
  
}
  # Fill lists only with angles, not cosines

Lpennat3D <- unlist(Lpennat3D)
Rpennat3D <- unlist(Rpennat3D)
  # UNlist all values

Lfibdis <- L_fib
Rfibdis <- R_fib
  # Fiber lengths measured from dissected muscles

LmuscWT <- dat_dis$L_closer_g
RmuscWT <- dat_dis$R_closer_g
  # Weight of dissected muscles

synth.fun <- function(x) {
  
  x_nona <- na.omit(x)
  avg <- mean(x_nona)
  mx <- max(x_nona)
  mi <- min(x_nona)
  sdev <- sd(x_nona)
  N <- length(x_nona)
  
  return(c(avg, sdev, mx, mi, N))
  
}

df_glob <- data.frame(
  BF = synth.fun(BF),
  HH = synth.fun(HH),
  HL = synth.fun(HL),
  HW = synth.fun(HW),
  PW = synth.fun(PW),
  TL = synth.fun(TL),
  FL_dis = synth.fun(c(Lfibdis, Rfibdis)),
  FL_3D = synth.fun(c(Lfib3D, Rfib3D)),
  PA = synth.fun(c(Lpennat3D, Rpennat3D)),
  MW = synth.fun(c(LmuscWT, RmuscWT)),
  MV = synth.fun(c(vol_L, vol_R)),
  MV2D3D = synth.fun(c(CH2D3D_L, CH2D3D_R)),
  MVCH = synth.fun((c(CH_L, CH_R)))
)
  # Global data frame combining left and right sides

df_LR <- data.frame(
  FL_dis_L = synth.fun(Lfibdis),
  FL_dis_R = synth.fun(Rfibdis),
  FL_3D_L = synth.fun(Lfib3D),
  FL_3D_R = synth.fun(Rfib3D),
  PA_L = synth.fun(Lpennat3D),
  PA_R = synth.fun(Rpennat3D),
  MW_L = synth.fun(LmuscWT),
  MW_R = synth.fun(RmuscWT),
  MV_L = synth.fun(vol_L),
  MV_R = synth.fun(vol_R),
  MV2D3D_L = synth.fun(CH2D3D_L),
  MV2D3D_R = synth.fun(CH2D3D_R),
  MVCH_L = synth.fun(CH_L),
  MVCH_R = synth.fun(CH_R)
)
  # Split sides data frame

rownames(df_LR) <- rownames(df_glob) <- c("Mean", 
                                            "Std. Dev.",
                                            "Maximum",
                                            "Minimum",
                                            "N")
  # Add rownames to the data frames

write.table(df_LR, 
            file = "df_LR.csv", 
            sep = ",", 
            dec = ".", 
            row.names = T)

write.table(df_glob, 
            file = "df_global.csv", 
            sep = ",", 
            dec = ".", 
            row.names = T)
  # Write data frames as CSV files

ratio_RoL <- df_LR[1, seq(from = 2, 
                         to = dim(df_LR)[2], 
                         by = 2)] /
            df_LR[1, seq(from = 1, 
                         to = dim(df_LR)[2]-1, 
                         by = 2)]

#-------------------------------------------------------------------------------
# Comparisons of average in vivo and estimated (closed mandibles) forces as
# reported in abstract

mean(BF_closed_dissec, na.rm = T) / 
  mean(dat_dis$maxBF_ampbasecorr, na.rm = T)

mean(BF_closed_2D3D, na.rm = T) / 
  mean(dat_3D$maxBF_ampbasecorr, na.rm = T)

mean(BF_closed_CH, na.rm = T) / 
  mean(dat_3D$maxBF_ampbasecorr, na.rm = T)

mean(BF_closed_insertion, na.rm = T) / 
  mean(dat_3D$maxBF_ampbasecorr, na.rm = T)

#-------------------------------------------------------------------------------
# Linear models of bite force against linear measurements

summary(lm(BF ~ HW)) # Highest R2
summary(lm(BF ~ HH)) # Lowest R2
summary(lm(BF ~ HL))
summary(lm(BF ~ TL))
summary(lm(BF ~ PW))

#-------------------------------------------------------------------------------
# t tests of muscular traits and bite forces

# Sexual dimorphism
Fe <- which(dat$Sex == "F")
Ma <- which(dat$Sex == "M")

t.test(BF[Fe], 
       BF[Ma])

mean(BF[Fe], na.rm = T) / mean(BF[Ma], na.rm = T)


t.test(HH[Fe], 
       HH[Ma])

mean(HH[Fe], na.rm = T) / mean(HH[Ma], na.rm = T)

t.test(HW[Fe], 
       HW[Ma])

mean(HW[Fe], na.rm = T) / mean(HW[Ma], na.rm = T)

t.test(HL[Fe], 
       HL[Ma])

mean(HL[Fe], na.rm = T) / mean(HL[Ma], na.rm = T)

t.test(TL[Fe], 
       TL[Ma])

mean(TL[Fe], na.rm = T) / mean(TL[Ma], na.rm = T)

t.test(PW[Fe], 
       PW[Ma])

mean(PW[Fe], na.rm = T) / mean(PW[Ma], na.rm = T)

Anova(lm(BF ~ HW * dat$Sex))
  # Anova in supplementary information

# Left / Right asymmetry muscles dissected

t.test(dat_dis$R_closer_g, 
       dat_dis$L_closer_g,
       paired = T)

mean(dat_dis$R_closer_g/dat_dis$L_closer_g)

t.test(R_fib,
       L_fib)

mean(R_fib) / mean(L_fib)

t.test(PCSA_R_mm2, 
       PCSA_L_mm2,
       paired = T)

mean(PCSA_R_mm2/PCSA_L_mm2,
     na.rm = T)

# Left / Right asymmetry 3D reconstructed muscles

t.test(vol_R, 
       vol_L,
       paired = T)

mean(vol_R/vol_L,
     na.rm = T)

t.test(CH2D3D_R, 
       CH2D3D_L,
       paired = T)

mean(CH2D3D_R/CH2D3D_L,
     na.rm = T)

t.test(CH_R, 
       CH_L,
       paired = T)

mean(CH_R/CH_L,
     na.rm = T)

t.test(fibR, 
       fibL,
       paired = F)

mean(fibR, na.rm = T)/mean(fibL, na.rm = T)

t.test(RP, 
       LP,
       paired = F)

mean(RP, na.rm = T)/mean(LP, na.rm = T)

t.test(PCSA_R_vol, 
       PCSA_L_vol,
       paired = T)

mean(PCSA_R_vol/PCSA_L_vol,
     na.rm = T)

t.test(PCSA_R_2d3dCH, 
       PCSA_L_2d3dCH,
       paired = T)

mean(PCSA_R_2d3dCH/PCSA_L_2d3dCH,
     na.rm = T)

t.test(PCSA_R_CH, 
       PCSA_L_CH,
       paired = T)

mean(PCSA_R_CH/PCSA_L_CH,
     na.rm = T)


# 3D vs. dissec tests

t.test(L_fib, fibL)

t.test(R_fib, fibR)

mean(R_fib)/mean(fibR)

t.test(PCSA_L_mm2, PCSA_L_vol)

mean(PCSA_L_mm2, na.rm = T) / mean(PCSA_L_vol, na = T)

t.test(PCSA_R_mm2, PCSA_R_vol)

mean(PCSA_R_mm2, na.rm = T) / mean(PCSA_R_vol, na = T)

t.test(PCSA_L_mm2, PCSA_L_2d3dCH)

mean(PCSA_L_2d3dCH, na = T) / mean(PCSA_L_mm2, na.rm = T)

t.test(PCSA_R_mm2, PCSA_R_2d3dCH)

mean(PCSA_R_2d3dCH, na = T) / mean(PCSA_R_mm2, na.rm = T)

t.test(PCSA_L_mm2, PCSA_L_CH)

mean(PCSA_L_CH, na = T) / mean(PCSA_L_mm2, na.rm = T)

t.test(PCSA_R_mm2, PCSA_R_CH)

mean(PCSA_R_CH, na = T) / mean(PCSA_R_mm2, na.rm = T)

mean(PCSA_insert, na = T)

mean(PCSA_insert, na = T) / mean(PCSA_R_mm2, na.rm = T)

sd(PCSA_L_mm2, na.rm = T)
sd(PCSA_R_mm2, na.rm = T)
sd(PCSA_L_vol, na.rm = T)
sd(PCSA_R_vol, na.rm = T)
sd(PCSA_L_2d3dCH, na.rm = T)
sd(PCSA_R_2d3dCH, na.rm = T)
sd(PCSA_L_CH, na.rm = T)
sd(PCSA_R_CH, na.rm = T)
sd(PCSA_insert, na.rm = T)

#-------------------------------------------------------------------------------
# Comparisons bite forces

t.test(BF_dis,
       BF_3D)

t.test(BF,
       BF_open_dissec)

t.test(BF,
       BF_open_VOL)

t.test(BF,
       BF_closed_dissec)

t.test(BF,
       BF_closed_VOL)

mean(BF_open_dissec, na.rm = T) / mean(BF, na.rm = T)

mean(BF_open_VOL, na.rm = T) / mean(BF, na.rm = T)

mean(BF_open_2D3D, na.rm = T) / mean(BF, na.rm = T)

mean(BF_open_CH, na.rm = T) / mean(BF, na.rm = T)

mean(BF_open_insertion, na.rm = T) / mean(BF, na.rm = T)

mean(BF_closed_dissec, na.rm = T) / mean(BF, na.rm = T)

mean(BF_closed_VOL, na.rm = T) / mean(BF, na.rm = T)

mean(BF_closed_2D3D, na.rm = T) / mean(BF, na.rm = T)

mean(BF_closed_CH, na.rm = T) / mean(BF, na.rm = T)

mean(BF_closed_insertion, na.rm = T) / mean(BF, na.rm = T)


#-------------------------------------------------------------------------------
# Linear models of estimated vs. in vivo bite force.

summary(lm(BF_3D ~ BF_closed_VOL))

summary(lm(BF_dis ~ BF_closed_dissec))

summary(lm(BF_3D ~ BF_closed_2D3D))

summary(lm(BF_3D ~ BF_closed_CH))

summary(lm(BF_3D ~ BF_closed_insertion[1:15]))
