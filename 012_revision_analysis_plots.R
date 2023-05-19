# Script for revisions proposed by D. Labonte
library(scales)
# Load the data from previous scripts

datas <- list.files(pattern = ".RData")

for (i in 1:length(datas)) {load(datas[i])}

#-------------------------------------------------------------------------------
# Calculate stress value from Taylor's data (with confidence interval)

mean_sarco <- mean(mat[,1])

log_sarco <- log(mean_sarco)

pred_logstress <- predict(modlogs, 
                       newdata = data.frame(logsarco=log_sarco),
                       interval = "prediction")

pred_stress_intervals <- exp(pred_stress)

#-------------------------------------------------------------------------------
# Make log-log regressions of RIGHT MUSCLE ONLY PCSA vs in vivo bite force

BF_3D <- dat_3D$maxBF_ampbasecorr
BF_dissec <- dat_dis$maxBF_ampbasecorr

lBF3D <- log(BF_3D)
lBFdi <- log(BF_dissec)

lAdi <- log(PCSA_R)
lA2D3D <- log(PCSA_2D3D_R)
lACH <- log(PCSA_CH_R)
lAvol <- log(PCSA_vol_R)
insert <- dat_3D$insert_area_total/2
lAinsert <- log(insert)

#Regression and log-log regression for dissected specimens

rawdi <- lm(BF_dissec ~ PCSA_R)
CIdi <- confint(rawdi)

loglogdi <- lm(lBFdi ~ lAdi)
CIlogdi <- confint(loglogdi)

#Regression and log-log regression for 3D reconstructed specimens

  #Raw volume estimates
rawvol <- lm(BF_3D ~ PCSA_vol_R)
CIvol <- confint(rawvol)

loglogvol <- lm(lBF3D ~ lAvol)
CIlogvol <- confint(loglogvol)

  #2D3D volume estimates
raw2D3D <- lm(BF_3D ~ PCSA_2D3D_R)
CI2D3D <- confint(raw2D3D)

loglog2D3D <- lm(lBF3D ~ lA2D3D)
CIlog2D3D <- confint(loglog2D3D)

  #ConvexHull volume estimates
rawCH <- lm(BF_3D ~ PCSA_CH_R)
CICH <- confint(rawCH)

loglogCH <- lm(lBF3D ~ lACH)
CIlogCH <- confint(loglogCH)

 #Insertion area estimates
rawinsert <- lm(BF_3D ~ insert)
CIinsert <- confint(rawinsert)

logloginsert <- lm(lBF3D ~ lAinsert)
CIloginsert <- confint(logloginsert)

#Make tables gathering estimates and their CI

intercepts_CIs_reg <- cbind(c(rawdi$coefficients[1], CIdi[1,]),
                        c(rawvol$coefficients[1], CIvol[1,]),
                        c(raw2D3D$coefficients[1], CI2D3D[1,]),
                        c(rawCH$coefficients[1], CICH[1,]),
                        c(rawinsert$coefficients[1], CIinsert[1,]))

slopes_CIs_reg <- cbind(c(rawdi$coefficients[2], CIdi[2,]),
                            c(rawvol$coefficients[2], CIvol[2,]),
                            c(raw2D3D$coefficients[2], CI2D3D[2,]),
                            c(rawCH$coefficients[2], CICH[2,]),
                            c(rawinsert$coefficients[2], CIinsert[2,]))

slopes_CIs_loglog <- cbind(c(loglogdi$coefficients[2], CIlogdi[2,]),
                        c(loglogvol$coefficients[2], CIlogvol[2,]),
                        c(loglog2D3D$coefficients[2], CIlog2D3D[2,]),
                        c(loglogCH$coefficients[2], CIlogCH[2,]),
                        c(logloginsert$coefficients[2], CIloginsert[2,]))

#Plot values and CIs against expectations
#Intercept of raw regression is expected to be 0 (a muscle with PCSA = 0 
# produces a force of 0)
#Slope of raw regression corresponds to muscle stress: should be compared to
#known values from Taylor/literature
#Slope of loglog regression is expected to be 1

plot(intercepts_CIs[1,], 
     ylim = c(min(intercepts_CIs),
              max(intercepts_CIs)),
     type = "n")

for (i in 1:5) {
  
  polygon(x = c(i - 0.1, 
                i + 0.1,
                i + 0.1,
                i - 0.1),
          y = c(intercepts_CIs[c(2,2, 3,3), i]),
          col = alpha("firebrick", 
                                 0.5))
}

points(intercepts_CIs[1,], 
     ylim = c(min(intercepts_CIs),
              max(intercepts_CIs)),
     cex = 5,
     pch = 20)

abline(h = 0,
       lwd = 2,
       lty = 2)

####

plot(slopes_CIs_loglog[1,], 
     ylim = c(min(slopes_CIs_loglog),
              max(slopes_CIs_loglog)),
     type = "n")

for (i in 1:5) {
  
  polygon(x = c(i - 0.1, 
                i + 0.1,
                i + 0.1,
                i - 0.1),
          y = c(slopes_CIs_loglog[c(2,2, 3,3), i]),
          col = alpha("firebrick", 
                                 0.5))
}

points(slopes_CIs_loglog[1,], 
       ylim = c(min(slopes_CIs_loglog),
                max(slopes_CIs_loglog)),
       cex = 5,
       pch = 20)

abline(h = c(0, 1),
       lwd = 2,
       lty = 2)

#-------------------------------------------------------------------------------
# Three panel figure with 
# i) log-log regression slopes and CI against expectation (1)
# ii) raw regression intercept against expectation (0)
# iii) raw regression slope (ie stress expectation) against literature values

output_folder <- "../Reresubmission/Figures/"

pdf(file = paste(output_folder,
                 "Regressions_param_expect.pdf",
                 sep = ""),
    height = 3,
    width = 11)

layout(mat = matrix(1:3, 
                    ncol = 3))

par(mar = c(7, 5, 1, 1))

#Plots of log log slope 

plot(slopes_CIs_loglog[1, ],
     ylim = c(min(slopes_CIs_loglog),
              max(slopes_CIs_loglog)),
     pch = 20,
     cex = 3,
     ylab = "Slope of log-log regression (95% C.I.)",
     xaxt = "n",
     xlab = "")

axis(side = 1,
     at = 1:5,
     labels = c("Dissection", 
                "Volume", 
                "2D3D", 
                "Convex Hull", 
                "Insertion area"),
     las = 2)

for (i in 1:5) {
  
  lines(x = c(i, i),
        y = c(slopes_CIs_loglog[c(2, 3), i]),
        col = "black",
        lwd = 3)
}

abline(h = c(0, 1),
       lty = c(2, 1),
       lwd = 2)

legend("topleft", 
       legend = "A.", 
       bty = "n",
       cex = 1.5)

#Plot intercept raw regression

plot(intercepts_CIs_reg[1, ],
     ylim = c(min(intercepts_CIs_reg),
              max(intercepts_CIs_reg)),
     pch = 20,
     cex = 3,
     ylab = "Intercept of regression (95% C.I.)",
     xaxt = "n",
     xlab = "")

axis(side = 1,
     at = 1:5,
     labels = c("Dissection", 
                "Volume", 
                "2D3D", 
                "Convex Hull", 
                "Insertion area"),
     las = 2)

for (i in 1:5) {
  
  lines(x = c(i, i),
        y = c(intercepts_CIs_reg[c(2, 3), i]),
        col = "black",
        lwd = 3)
}

abline(h = 0,
       lwd = 2)

legend("topleft", 
       legend = "B.", 
       bty = "n",
       cex = 1.5)

# Plot slopes (i.e. stress estimates), and stress estimated from sarcomere 
# length. The latter must be converted to the correct unit (divided by 1000).

stresses <- c(slopes_CIs_reg[1, ], pred_stress_intervals[1]/1000)
CIs <- cbind(slopes_CIs_reg[2:3, ], pred_stress_intervals[2:3]/1000)

plot(stresses, 
     ylim = c(min(CIs),
              max(CIs)),
     type = "n",
     xaxt = "n",
     xlab = "",
     ylab = "Muscle stress (N.mm^-2)")

axis(side = 1,
     at = 1:6,
     labels = c("Dissection", 
                "Volume", 
                "2D3D", 
                "Convex Hull", 
                "Insertion area",
                "Sarcomere"),
     las = 2)

polygon(x = c(0,7,7,0),
        y = c(0.18, 0.18, 1.16, 1.16),
        #Values for closer muscle in leaf cutter ant and stag beetle
        col = "gray80",
        border = NA)

polygon(x = c(0,7,7,0),
        y = c(0.31, 0.31, 0.58, 0.58), 
        #Values for closer muscle in insects from the literature
        col = "gray60",
        border = NA)

for (i in 1:6) {
  
  lines(x = c(i, i),
          y = c(CIs[c(1, 2), i]),
          col = "black",
        lwd = 3)
}

points(stresses, 
       cex = 3,
       pch = 20)

abline(h = 0,
       lwd = 2,
       lty = 2)

legend("bottomright",
       legend = c("Muscle stress estimate, with 95% C.I.", 
                  "0m1 muscle stress values from the literature",
                  "Idem, excluding leaf-cutter ant and stag beetle"), 
       pch = c(20, 22, 22),
       lty = c(1, NA, NA),
       pt.bg = c(NA, "grey80", "grey60"),
       pt.cex = 1.5,
       cex = 0.8)

legend("topleft", 
       legend = "C.", 
       bty = "n",
       cex = 1.5)

dev.off()
