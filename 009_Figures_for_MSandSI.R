#-------------------------------------------------------------------------------
# This script produces all plots of the manuscript, as well as some plots shown
# in the supplementary information.

# Required libraries:
library(car)
library(MASS)
library(scales)

#-------------------------------------------------------------------------------
# Read main data table and load data from previous scripts

dat <- read.csv("measurments_Sch.csv", h=T, dec=".", sep=",")
load("all_BF.RData")
load("fiber_lengths_and_summary_stats.RData")
load("mat_fiber_lgt_angle.RData")
load("PCSA_matrices.RData")
load("pennation_angles.RData")
load("fiber_legnths_3D.RData")

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
    height = 10, 
    width = 5)

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

text(x = 0.3, y = 3.7,
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

text(x = 0.5, y = 7,
     labels = "D.", 
     cex=2,
     font=2)

dev.off()
