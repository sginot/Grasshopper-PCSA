#-------------------------------------------------------------------------------
# This script produces all plots of the manuscript, as well as some plots shown
# in the supplementary information.

# Required libraries:
library(car)
library(MASS)
library(scales)

#-------------------------------------------------------------------------------
# Read main data table and load data from previous scripts

dat <- read.csv("measurments_Sch.csv", 
                h = T, 
                dec = ".",
                sep = ",")

dat_taylor <- read.csv("sarcomere_stress_data_Taylor_2001.csv",
                       h = T,
                       sep = "\t",
                       dec = ".")

load("all_BF.RData")
load("fiber_lengths_and_summary_stats.RData")
load("mat_fiber_lgt_angle.RData")
load("PCSA_matrices.RData")
load("pennation_angles.RData")
load("fiber_lengths_3D.RData")
load("mat_sarcomere_lengths_summary.RData")
load("list_sarcomere_length_measurements.RData")
load("models.RData")


#-------------------------------------------------------------------------------
# Set output folder for figures

output_folder <- "Figures/"

#-------------------------------------------------------------------------------
# Main figure showing boxplots and distributions for fiber lengths, pennation
# angles, PCSAs and bite forces

fibL <- unlist(list_lgt_fib[seq(from = 1, 
                          to = 29,
                          by = 2)])
  # Make vector containing all fiber length values based on 3D landmarks. Here
  # the original list is indexed to gather only values for LEFT muscles

fibR <- unlist(list_lgt_fib[seq(from = 2, 
                                 to = 30,
                                 by = 2)])
  # Here the original list is indexed to gather only values for RIGHT muscles

arr.pennat <- array(unlist(ls.pennat), dim = c(25, 2, 30))
  # Turn list containing matrices of pennation angles and cosines into an array
  # containing the same values

LP <- arr.pennat[, 1, c(2, 3, 5, 8, 9, 12, 13, 15, 18, 20, 22, 24, 25, 27, 29)]
  # Restrict to pennation angle values of the fibers from LEFT muscles. Indices
  # must be manually specified because some sides were reverted in the 3D 
  # reconstructions (specimens Sch1, 4, 6, 9 to 12)

RP <- arr.pennat[, 1, c(1, 4, 6, 7, 10, 11, 14, 16, 17, 19, 21, 23, 26, 28, 30)]
  # Restrict to pennation angle values of the fibers from RIGHT muscles

DRP <- density(RP)
DLP <- density(LP)
  # Estimates kernel densities for right and left pennation angles

PCSA_L_mm2 <- mat_PCSA[, 1]
PCSA_R_mm2 <- mat_PCSA[, 2]
PCSA_L_vol <- mat_vols_PCSA[, 1]
PCSA_R_vol <- mat_vols_PCSA[, 2]
PCSA_L_2d3dCH <- mat_vols_PCSA[, 5]
PCSA_R_2d3dCH <- mat_vols_PCSA[, 6]
PCSA_L_CH <- mat_vols_PCSA[, 3]
PCSA_R_CH <- mat_vols_PCSA[, 4]
PCSA_insert <- dat$insert_area_total / 2
  # Single vectors of PCSA estimates, to be used in boxplot

N_PCSA <- c(length(na.omit(PCSA_L_mm2)),
            length(na.omit(PCSA_R_mm2)),
            length(na.omit(PCSA_L_vol)),
            length(na.omit(PCSA_R_vol)),
            length(na.omit(PCSA_L_2d3dCH)),
            length(na.omit(PCSA_R_2d3dCH)),
            length(na.omit(PCSA_L_CH)),
            length(na.omit(PCSA_R_CH)),
            length(na.omit(PCSA_insert)))
  # Vector containing sample sizes for each type of PCSA estimates

BF_dis <- dat$maxBF_ampbasecorr[which(dat$Dissected == 1)]
BF_3D <- dat$maxBF_ampbasecorr[which(dat$Dissected == 0)]
  # Single vectors for in vivo bite force

N_BF <- c(length(na.omit(BF_dis)),
          length(na.omit(BF_3D)),
          length(na.omit(BF_closed_dissec)),
          length(na.omit(BF_open_dissec)),
          length(na.omit(BF_closed_VOL)),
          length(na.omit(BF_open_VOL)),
          length(na.omit(BF_closed_2D3D)),
          length(na.omit(BF_open_2D3D)),
          length(na.omit(BF_closed_CH)),
          length(na.omit(BF_open_CH)),
          length(na.omit(BF_closed_insertion)),
          length(na.omit(BF_open_insertion)))
  # Vector containing sample sizes for bite forces estimates


# Make the actual plot
pdf(file = paste(output_folder,
                 "Main_figure.pdf"), 
    height = 12, 
    width = 7)

layout(matrix(c(1:2, 3, 3, 4, 4), 
              ncol = 2, 
              nrow = 3, 
              byrow = T))

par(mar=c(4, 4.5, 1, 1))

boxplot(L_fib, R_fib, fibL, fibR,
        col = c("darkorange", "darkorchid4"),
        pch = 20,
        cex=2,
        lwd=2,
        at = c(0.7, 1.3, 2.4,3),
        boxwex = 0.5,
        xaxt = "n",
        ylab = "Fiber length (mm)", 
        cex.lab=1.5,
        ylim = c(0, 4))
  # Boxplots of fiber lengths, with left/right muscles as different colors, and
  # comparison between microdissected fibers and 3D landmarked fibers (separated
  # by a large space).

text(x = c(0.7, 1.3, 2.4,3), 
     y = rep(0.2, 4),
     labels = c(length(L_fib), 
                length(R_fib), 
                length(fibL),
                length(fibR)),
     cex = 1.2,
     font = 2)
  # Add number of individual fibers at the bottom of each box

axis(side = 1,
     at = c(1,2.7),
     labels = F)

text(x = c(1,2.7), y = -0.5,
     labels = c("Dissection", 
                "3D (landmarks)"), 
     cex = 1.5,
     srt = 20,
     xpd = T)

text(x = 0.5, y = 3.9,
     labels = "A.", 
     cex = 2,
     font = 2)
  # Add custom x axis

plot(DRP, 
     main = "", 
     xlab = "Pennation angle (degrees)", 
     ylim = c(0, max(DLP$y)), 
     cex.lab = 1.5)
  # Plot density curves of pennation angles (obtained from 3D landmarked)

polygon(DLP$x, 
        DLP$y, 
        col = alpha("darkorange", 0.5), 
        lwd=3)

polygon(DRP$x, 
        DRP$y, 
        col = alpha("darkorchid4", 0.5), 
        lwd = 3)
  # Add transparent colors to the density curves

abline(v = mean(DRP$x), 
       col = "darkorchid4", 
       lty = 2, 
       lwd = 2)

abline(v = mean(DLP$x), 
       col = "darkorange", 
       lty = 2, 
       lwd = 2)
  # Add vertical lines at the respective left and right averages

text(c(mean(DRP$x), mean(DLP$x)),
     c(0.002, 0.002),
     labels = c(DRP$n, DLP$n),
     pos = c(4, 2),
     font = 2,
     col = c("darkorchid4", "darkorange"),
     cex = 1.2)
  # Add sample size (number of fibers)

text(x = 0, y = 0.042,
     labels = "B.", 
     cex=2,
     font=2)

par(mar=c(5,4.5,1,1))
  # Change margins for the next panel

boxplot(PCSA_L_mm2,
        PCSA_R_mm2,
        PCSA_L_vol,
        PCSA_R_vol,
        PCSA_L_2d3dCH,
        PCSA_R_2d3dCH,
        PCSA_L_CH,
        PCSA_R_CH,
        PCSA_insert,
        col = c(rep(c("darkorange", 
                      "darkorchid4"), 
                    4),
                "grey"),
        ylim = c(0, 19),
        pch = 20,
        cex = 2,
        lwd = 2,
        boxwex = 0.5,
        at = c(1,1.5, 3,3.5, 5,5.5, 7,7.5, 9.25),
        xaxt = "n",
        ylab = "Area (mm^2)", 
        cex.lab = 1.5)
  # Boxplot of the various estimates of PCSA

text(c(1,1.5, 3,3.5, 5,5.5, 7,7.5, 9.25),
     rep(0, 9),
     labels = N_PCSA,
     font = 2,
     cex = 1.2)
  # Add sample sizes

axis(side = 1,
     at = c(1.25, 3.25, 5.25, 7.25, 9.25),
     labels = F)

text(x = c(0.5, 2.5, 4.5, 6.5, 8.5),
     y = -2.8,
     labels = c("PCSA", 
                "PCSA 3D",
                "PCSA 2D-3D",
                "PCSA 3D",
                "PCSA"),
     srt = 20,
     xpd = T, 
     cex = 1.5)

text(x = c(1, 3, 5, 7, 9),
     y = -3.2,
     labels = c("dissection", 
                "volume",
                "convex hull",
                "convex hull",
                "muscle insertion"),
     srt = 20,
     xpd = T, 
     cex = 1.5)
  # Add labels: some manual adjustments are necessary depending on the 
  # dimensions of the figure

abline(v = 2.5, 
       lwd = 2, 
       lty = 2, 
       xpd = F)

text(x = 0.5,
     y = 18,
     labels = "C.", 
     cex = 2,
     font = 2)


boxplot(BF_dis,
        BF_3D,
        BF_closed_dissec,
        BF_open_dissec,
        BF_closed_VOL,
        BF_open_VOL,
        BF_closed_2D3D,
        BF_open_2D3D,
        BF_closed_CH,
        BF_open_CH,
        BF_closed_insertion,
        BF_open_insertion,
        col = c("red", "firebrick", rep(c("grey", "bisque2"),5)),
        pch = 20,
        cex=2,
        lwd=2,
        boxwex = 0.3,
        at = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2),
        xaxt = "n",
        ylab = "Bite force (N)", 
        cex.lab=1.5)
  # Final panel containing bite force estimates

polygon(x = c(0,10,10,0), 
        y = c(max(BF_dis, na.rm = T), 
              max(BF_dis, na.rm = T),
              min(BF_3D, na.rm = T),
              min(BF_3D, na.rm = T)),
        col = alpha("firebrick", alpha = 0.2),
        xpd = F,
        border = T)
  # Add colored area showing the max and min bite force value measured in vivo

text(c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2),
     rep(0.2, 12),
     labels = N_BF,
     cex = 1.2,
     font = 2)
  # Add sample sizes

axis(side = 1,
     at = c(1:6),
     labels = F)

text(x = c(1:6),
     y = -1,
     labels = c("In vivo",
                "Dissection", 
                "3D volume",
                "2D-3D hull",
                "3D convex hull",
                "Muscle insertion"),
     srt = 20,
     xpd = T, 
     cex = 1.5)
  # Add x labels

abline(v = 1.5, 
       lwd = 2, 
       lty = 2,
       xpd = F)
abline(v = 2.5, 
       lwd = 2,
       lty = 2,
       xpd = F)
  # Add vertical lines to separate in vivo from dissection estimates from 
  # 3d estimates

text(x = 0.3, y = 7,
     labels = "D.", 
     cex=2,
     font=2)

dev.off()

#-------------------------------------------------------------------------------
# Bite forces correlations biplot

pdf(file = paste(output_folder, "Correlations_fig.pdf"),
    width = 7, 
    height = 14)

layout(matrix(1:2, nrow = 2))

par(mar = c(5, 5, 1, 1))

plot(BF_closed_dissec, 
     BF_dis,
     xlim = c(0, 8),
     ylim = c(0.2, 2.2), 
     pch = 21, 
     col = "black", 
     bg = "red", 
     cex = 2,
     lwd = 2,
     xlab = "Estimated bite force (N)",
     ylab = "In vivo bite force (N)",
     cex.lab = 1.5)

clip(0, 9, 
     min(BF_dis, na.rm = T),
     max(BF_dis, na.rm = T))

abline(lm(BF_dis ~ BF_closed_dissec), 
       col = "red",
       lwd = 3)

clip(0, 9, 0.2 ,2.2)

abline(0, 1, 
       lty = 2, 
       lwd = 1.5, 
       xpd = F)

points(BF_closed_VOL,
       BF_3D,
       pch = 22,
       bg = "blue",
       cex = 2,
       lwd=2, 
       xpd = NA)

points(BF_closed_2D3D,
       BF_3D,
       pch = 23,
       bg = "purple",
       cex = 2,
       lwd = 2)

points(BF_closed_CH,
       BF_3D, 
       pch = 24,
       bg = "darkorange",
       cex = 2,
       lwd=2)

points(na.omit(BF_closed_insertion),
       BF_3D,
       pch = 25,
       bg = "forestgreen",
       cex = 2,
       lwd=2)

clip(0, 9, 
     min(BF_3D, na.rm = T),
     max(BF_3D, na.rm = T))

abline(lm(BF_3D ~ BF_closed_2D3D),
       col = "purple", 
       lwd = 3)

text(x = mean(BF_closed_2D3D, na.rm = T),
     y = mean(BF_3D, na.rm = T),
     pos = 4,
     col = "purple",
     font = 2,
     cex = 1.5,
     labels = round(lm(BF_3D ~ BF_closed_2D3D)$coef[2],2))

abline(lm(BF_3D ~ BF_closed_CH), 
       col = "darkorange", 
       lwd = 3)

text(x = mean(BF_closed_CH, na.rm = T),
     y = mean(BF_3D, na.rm = T),
     pos = 4,
     col = "darkorange",
     font = 2,
     cex = 1.5,
     labels = round(lm(BF_3D ~ BF_closed_CH)$coef[2],2))

abline(lm(BF_3D ~ na.omit(BF_closed_insertion)), 
       col="forestgreen", 
       lwd = 3)

text(x = mean(BF_closed_insertion, na.rm = T),
     y = mean(BF_3D, na.rm = T),
     pos = 4,
     col = "forestgreen",
     font = 2,
     cex = 1.5,
     labels = round(lm(BF_3D ~ na.omit(BF_closed_insertion))$coef[2],2))

text(x = mean(BF_closed_dissec, na.rm = T),
     y = mean(BF_dis, na.rm = T),
     pos = 4,
     col = "red",
     font = 2,
     cex = 1.5,
     labels = round(lm(BF_dis ~ BF_closed_dissec)$coef[2],2))


text(min(BF_closed_VOL, na.rm = T),
     max(BF_dis, na.rm = T),
     labels = "A.",
     xpd = NA,
     cex = 2,
     font = 2,
     pos = 2)

# Bite forces correlations biplot with OPEN mandible estimates

par(mar = c(5, 5, 1, 1))

plot(BF_open_dissec, 
     BF_dis,
     xlim = c(0, 4.5),
     ylim = c(0.2, 2.2), 
     pch = 21, 
     col = "black", 
     bg = "red", 
     cex = 2,
     lwd = 2,
     xlab = "Estimated bite force (N)",
     ylab = "In vivo bite force (N)",
     cex.lab = 1.5)

clip(0, 9, 
     min(BF_dis, na.rm = T),
     max(BF_dis, na.rm = T))

abline(lm(BF_dis ~ BF_open_dissec), 
       col = "red",
       lwd = 3)

clip(0, 9, 0.2 ,2.2)

abline(0, 1, 
       lty = 2, 
       lwd = 1.5, 
       xpd = F)

points(BF_open_VOL,
       BF_3D,
       pch = 22,
       bg = "blue",
       cex = 2,
       lwd=2, 
       xpd = NA)

points(BF_open_2D3D,
       BF_3D,
       pch = 23,
       bg = "purple",
       cex = 2,
       lwd = 2)

points(BF_open_CH,
       BF_3D, 
       pch = 24,
       bg = "darkorange",
       cex = 2,
       lwd=2)

points(na.omit(BF_open_insertion),
       BF_3D,
       pch = 25,
       bg = "forestgreen",
       cex = 2,
       lwd=2)

clip(0, 9, 
     min(BF_3D, na.rm = T),
     max(BF_3D, na.rm = T))

abline(lm(BF_3D ~ BF_open_2D3D),
       col = "purple", 
       lwd = 3)

abline(lm(BF_3D ~ BF_open_CH), 
       col = "darkorange", 
       lwd = 3)

abline(lm(BF_3D ~ na.omit(BF_open_insertion)), 
       col="forestgreen", 
       lwd = 3)

text(min(BF_open_VOL, na.rm = T),
     max(BF_dis, na.rm = T),
     labels = "B.",
     xpd = NA,
     cex = 2,
     font = 2,
     pos = 2)

clip(-5, 8, -1, 10)

legend(-1, 2.72, 
       legend = c("Volume estim.",
                  "Dissection estim.",
                  "2D-3D convex hull estim.",
                  "Convex hull estim.",
                  "Insertion area estim."),
       pch = c(22,21,23:35),
       pt.bg = c("blue", 
                 "red", 
                 "purple", 
                 "darkorange", 
                 "forestgreen"),
       bty = "o",
       bg = "white")

dev.off()

#-------------------------------------------------------------------------------
# Allometry of in vivo bite force against head width, with sexual dimorphism

lm_F <- lm(maxBF_ampbasecorr ~ HW, 
           data = dat[which(dat$Sex == "F"),])

lm_M <- lm(maxBF_ampbasecorr ~ HW, 
           data = dat[which(dat$Sex == "M"),])

CIF <- confint(lm_F)
CIM <- confint(lm_M)

newF <- seq(from = min(dat$HW[which(dat$Sex == "F")]),
            to = max(dat$HW[which(dat$Sex == "F")]),
            by = 0.05)
newM <- seq(from = min(dat$HW[which(dat$Sex == "M")], na.rm = T),
            to = max(dat$HW[which(dat$Sex == "M")], na.rm = T),
            by = 0.05)

pCIF <- predict(lm_F, 
                newdata = data.frame(HW = newF),
                interval = "confidence",
                level = 0.95)

pCIM <- predict(lm_M, 
                newdata = data.frame(HW = newM),
                interval = "confidence",
                level = 0.95)

pdf(file = paste(output_folder, "allometry_BF.pdf"),
    width = 7, 
    height = 7)

plot(dat$HW, 
     dat$maxBF_ampbasecorr, 
     bg = c(4, 2)[as.factor(dat$Sex)],
     pch = c(21, 22)[as.factor(dat$Sex)],
     cex = 2,
     lwd = 2,
     xlab = "Head width (mm)",
     ylab = "In vivo bite force (N)",
     xlim =c(5.9, 7.2),
     ylim = c(-0.1, 2.5))

polygon(c(newF, 
          rev(newF)), 
        c(pCIF[,2], 
          rev(pCIF[,3])),
        border = NA,
        col = alpha(4,
                    0.1))

polygon(c(newM, 
          rev(newM)), 
        c(pCIM[,2], 
          rev(pCIM[,3])),
        border = NA,
        col = alpha(2,
                    0.1))

clip(min(dat$HW[which(dat$Sex == "M")], 
         na.rm = T),
     max(dat$HW[which(dat$Sex == "M")], 
         na.rm = T),
     min(dat$maxBF_ampbasecorr[which(dat$Sex == "M")], 
         na.rm = T),
     max(dat$maxBF_ampbasecorr[which(dat$Sex == "M")], 
         na.rm = T))

polygon(x = c(dat$HW[which(dat$Sex == "F")], 
              rev(dat$HW[which(dat$Sex == "F")])),
        y = c(pCIF[,2:3]))

abline(lm_M,
       col = 2, 
       lwd = 3)

clip(min(dat$HW[which(dat$Sex == "F")], 
         na.rm = T),
     max(dat$HW[which(dat$Sex == "F")], 
         na.rm = T),
     min(dat$maxBF_ampbasecorr[which(dat$Sex == "F")], 
         na.rm = T),
     max(dat$maxBF_ampbasecorr[which(dat$Sex == "F")], 
         na.rm = T))

abline(lm_F,
       col = 4,
       lwd = 3)

legend("topleft",
       pch = c(21, 22),
       pt.bg = c(4, 2),
       legend = c(paste("Y =", 
                        round(lm_F[[1]][2], 2),
                        "[",
                        round(CIF[2, 1], 2),
                        ";",
                        round(CIF[2, 2], 2),
                        "]",
                        "* X", 
                        round(lm_F[[1]][1], 2),
                        "[",
                        round(CIF[1, 1], 2),
                        ";",
                        round(CIF[1, 2], 2),
                        "]",
                        "(females)"),
                  paste("Y =", 
                        round(lm_M[[1]][2], 2),
                        "[",
                        round(CIM[2, 1], 2),
                        ";",
                        round(CIM[2, 2], 2),
                        "]",
                        "* X", 
                        round(lm_M[[1]][1], 2),
                        "[",
                        round(CIM[1, 1], 2),
                        ";",
                        round(CIM[1, 2], 2),
                        "]",
                        "(males)")),
       xpd = NA)

dev.off()

#-------------------------------------------------------------------------------
# Sarcomere length model figure

# Define variables of interest

sarco_lgt <- dat_taylor$Mean
  # Sarcomere length data from Taylor 2001

stress <- dat_taylor$Mean.1
  # Muscle stress data from Taylor 2001

sarco_measured <- mat[, 1]
  # Sarcomere length values obtained from Sch1 muscle fibers

pred_stresses <- predict(object = mod, 
                       newdata = data.frame(sarco_lgt = sarco_measured))
  # Predict muscle stress values for measured sarcomere lengths, based on model
  # from script 004 (Taylor's regression re-run).
  # However, length and stress have an allometric, therefore non_linear relation
  # therefore predictions based on linear regression on raw values are expected
  # to be incorrect. Must use logarithmic scale

log_mean_sarco <- log(mean(sarco_measured))
  # Log of the average of measured sarcomere lengths

pred_log_stress <- modlogs$coefficients[2] *
  log_mean_sarco + modlogs$coefficients[1]
  # Predict the log muscle stress values based on log mean measured sarcomere
  # based on log-log regression of Taylor's data. In script 006, this predicted
  # log stress value is converted back to natural scale using exp()

pdf(file = paste(output_folder,
                 "Reanalysis_Taylor_predict_stress.pdf"),
    width = 12,
    height = 7)

layout(matrix(1:2, ncol = 2))

plot(sarco_lgt, 
     stress, 
     pch = 19, 
     cex = 1.5,
     xlab = "Sarcomere lenght (µm)",
     ylab = "Muscle stress (kN.m^-2)")

abline(mod, 
       lwd = 2)

points(sarco_measured, 
       pred_stresses,
       pch = 22,
       bg = "red", 
       cex = 1.5)

points(mean(sarco_measured), 
       mean(pred_stresses),
       pch = 22,
       bg = "red", 
       cex = 2.5,
       lwd = 3)

text(min(sarco_lgt, na.rm = T), 
     max(stress, na.rm = T), 
     labels = paste("Y =", 
                    round(mod$coefficients[2], 2), 
                    "* x +", 
                    round(mod$coefficients[1], 2)),
     pos = 4)

plot(log(sarco_lgt), 
     log(stress), 
     pch = 19, 
     cex = 1.5,
     xlab = "log Sarcomere lenght (µm)",
     ylab = "log Muscle stress (kN.m^-2)")

abline(modlogs, 
       lwd = 2)

points(log_mean_sarco, 
       pred_log_stress,
       pch = 22,
       bg = "red", 
       cex = 2.5,
       lwd = 3)

lines(x = c(0, 
            log_mean_sarco, 
            log_mean_sarco),
      y = c(pred_log_stress, 
            pred_log_stress, 
            0),
      lwd = 2,
      lty = 2,
      col = "red")

text(min(log(sarco_lgt), na.rm = T), 
     max(log(stress), na.rm = T), 
     labels = paste("Y =", 
                    round(modlogs$coefficients[2], 2), 
                    "* x +", 
                    round(modlogs$coefficients[1], 2)),
     pos = 4)

dev.off()
