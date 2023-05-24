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

pred_stress_intervals <- exp(pred_logstress)

#-------------------------------------------------------------------------------
# Make log-log regressions of RIGHT MUSCLE ONLY PCSA vs in vivo bite force

meanMA <- mean(dat$MA_R,
               na.rm = T)

meanMAeff <- meanMA * sin(60*pi/180)

MF_closed_3D <- dat_3D$maxBF_ampbasecorr / meanMA
MF_closed_dissec <- dat_dis$maxBF_ampbasecorr /meanMA

MF_open_3D <- dat_3D$maxBF_ampbasecorr / meanMAeff
MF_open_dissec <- dat_dis$maxBF_ampbasecorr / meanMAeff

lMF3D <- log(MF_3D)
lMFdi <- log(MF_dissec)

lMF3D_open <- log(MF_open_3D)
lMFdi_open <- log(MF_open_dissec)

lAdi <- log(PCSA_R)
lA2D3D <- log(PCSA_2D3D_R)
lACH <- log(PCSA_CH_R)
lAvol <- log(PCSA_vol_R)
insert <- dat_3D$insert_area_total/2
lAinsert <- log(insert)

#Regression and log-log regression for dissected specimens

rawdi <- lm(MF_dissec ~ PCSA_R)

CIdi <- confint(rawdi)

newdi <- seq(from = min(PCSA_R, na.rm = T),
            to = max(PCSA_R, na.rm = T),
            by = 0.05)

prawdi <- predict(object = rawdi,
                  newdata = data.frame(PCSA_R = newdi),
                  interval = "confidence")

loglogdi <- lm(lMFdi ~ lAdi)

CIlogdi <- confint(loglogdi)

newlogdi <- seq(from = min(lAdi, na.rm = T),
             to = max(lAdi, na.rm = T),
             by = 0.05)

ploglogdi <- predict(object = loglogdi,
                  newdata = data.frame(lAdi = newlogdi),
                  interval = "confidence")

#Regression and log-log regression for 3D reconstructed specimens

  #Raw volume estimates

rawvol <- lm(MF_3D ~ PCSA_vol_R)

CIvol <- confint(rawvol)

newvol <- seq(from = min(PCSA_vol_R, na.rm = T),
              to = max(PCSA_vol_R, na.rm = T),
              by = 0.05)

prawvol <- predict(object = rawvol,
                  newdata = data.frame(PCSA_vol_R = newvol),
                  interval = "confidence")

loglogvol <- lm(lMF3D ~ lAvol)

CIlogvol <- confint(loglogvol)

newlogvol <- seq(from = min(lAvol, na.rm = T),
              to = max(lAvol, na.rm = T),
              by = 0.05)

plogvol <- predict(object = loglogvol,
                   newdata = data.frame(lAvol = newlogvol),
                   interval = "confidence")

  #2D3D volume estimates

raw2D3D <- lm(MF_3D ~ PCSA_2D3D_R)

CI2D3D <- confint(raw2D3D)

new2D3D <- seq(from = min(PCSA_2D3D_R, na.rm = T),
               to = max(PCSA_2D3D_R, na.rm = T),
               by = 0.05)

praw2D3D <- predict(object = raw2D3D,
                    newdata = data.frame(PCSA_2D3D_R = new2D3D),
                    interval = "confidence")

loglog2D3D <- lm(lMF3D ~ lA2D3D)

CIlog2D3D <- confint(loglog2D3D)

newlog2D3D <- seq(from = min(lA2D3D, na.rm = T),
               to = max(lA2D3D, na.rm = T),
               by = 0.05)

plog2D3D <- predict(object = loglog2D3D,
                    newdata = data.frame(lA2D3D = newlog2D3D),
                    interval = "confidence")

  #ConvexHull volume estimates

rawCH <- lm(MF_3D ~ PCSA_CH_R)

CICH <- confint(rawCH)

newCH <- seq(from = min(PCSA_CH_R, na.rm = T),
               to = max(PCSA_CH_R, na.rm = T),
               by = 0.05)

prawCH <- predict(object = rawCH,
                    newdata = data.frame(PCSA_CH_R = newCH),
                    interval = "confidence")

loglogCH <- lm(lMF3D ~ lACH)

CIlogCH <- confint(loglogCH)

newlogCH <- seq(from = min(lACH, na.rm = T),
                to = max(lACH, na.rm = T),
                by = 0.05)

plogCH <- predict(object = loglogCH,
                  newdata = data.frame(lACH = newlogCH),
                  interval = "confidence")

 #Insertion area estimates

rawinsert <- lm(MF_3D ~ insert)

CIinsert <- confint(rawinsert)

newinsert <- seq(from = min(insert, na.rm = T),
                 to = max(insert, na.rm = T),
                 by = 0.05)

prawinsert <- predict(object = rawinsert,
                  newdata = data.frame(insert = newinsert),
                  interval = "confidence")

logloginsert <- lm(lMF3D ~ lAinsert)

CIloginsert <- confint(logloginsert)

newloginsert <- seq(from = min(lAinsert, na.rm = T),
                 to = max(lAinsert, na.rm = T),
                 by = 0.05)

ploginsert <- predict(object = logloginsert,
                      newdata = data.frame(lAinsert = newloginsert),
                      interval = "confidence")


#-------------------------------------------------------------------------------
#Automatize the whole thing

#Make function
stress.reg.fun <- function(x, 
                           force = BF) {
  
  raw <- lm(force ~ x)
  
  CI <- confint(raw)
  
  newx <- seq(from = min(x, na.rm = T) - 0.1,
              to = max(x, na.rm = T) + 0.1,
              by = 0.05)
  
  predraw <- predict(object = raw,
                     newdata = data.frame(x = newx),
                     interval = "confidence")
  
  logforce <- log(force)
  
  logx <- log(x)
  
  loglog <- lm(logforce ~ logx)
  
  CIlog <- confint(loglog)
  
  newlogx <- seq(from = min(logx, na.rm = T) - 0.1,
                 to = max(logx, na.rm = T) + 0.1,
                 by = 0.05)
  
  predlog <- predict(object = loglog,
                     newdata = data.frame(logx = newlogx),
                     interval = "confidence")
  
as.list(environment())

}

#Define data frame with only PCSA variables
PCSAs <- dat[, c(grep("PCSA", colnames(dat)), grep("insert", colnames(dat)))]

PCSAs <- PCSAs[, -grep("_L_", colnames(PCSAs))]

PCSAs[, 3:4] <- PCSAs[, 4:3]

PCSAs[, 5] <- PCSAs[, 5] / 2


#-------------------------------------------------------------------------------
# Plot the data and regression lines
# 10 panels arranged in two columns
# Add CI and the regression formula and plot CIs along the regression line

pal <- palette()

pal2 <- adjustcolor(palette(),
                    offset = c(0.5,0,0.3,0))

panel <- LETTERS[1:10]

output_folder <- "../Reresubmission/Figures/"

pdf(file = paste(output_folder, 
                 "Force_Area_regressions.pdf", 
                 sep =""),
    height = 15,
    width = 10)

matlayout <- cbind(rep(11, 6),
                   c(1, 3, 5, 7, 9, 12),
                   c(2, 4, 6, 8, 10, 12))

layout(mat = matlayout,
       widths = c(1, 4, 4),
       heights = c(rep(2.5, 5), 1))

par(mar = c(2.5, 2.5, 0.5, 0.5))

for (i in 1:5) {
  
  ls.reg <- stress.reg.fun(x = PCSAs[, i],
                       force = BF / meanMA)
  
  plot(1, 1,
       type = "n",
       xlim = c(min(ls.reg$newx),
                max(ls.reg$newx)),
       ylim = c(min(ls.reg$predraw),
                max(ls.reg$predraw)))
  
  points(ls.reg$x, 
         ls.reg$force,
         pch = 25,
         cex = 2,
         bg = pal[i],
         lwd = 2)
  
  abline(ls.reg$raw,
         col = pal[i])
  
  polygon(c(ls.reg$newx, rev(ls.reg$newx)),
          c(ls.reg$predraw[, 2], rev(ls.reg$predraw[, 3])),
          border = NA,
          col = alpha(pal[i], alpha = 0.2))
  
  legend("topright",
         bty = "n",
         legend = paste("Y = ", 
                              round(ls.reg$raw[[1]][2], 2),
                              "[",
                              round(ls.reg$CI[2, 1], 2),
                              ";",
                              round(ls.reg$CI[2, 2], 2),
                              "]",
                              " * X + ", 
                              round(ls.reg$raw[[1]][1], 2),
                              "[",
                              round(ls.reg$CI[1, 1], 2),
                              ";",
                              round(ls.reg$CI[1, 2], 2),
                              "]",
                        sep = ""),
         pch = 25,
         cex = 1.2)
  
  legend("topleft", 
         bty = "n",
         legend = panel[i*2-1],
         cex = 2)
  
  ls.reg <- stress.reg.fun(x = PCSAs[, i],
                           force = BF / meanMAeff)
  
  points(ls.reg$x, 
         ls.reg$force,
         pch = 24,
         cex = 2,
         bg = pal2[i],
         lwd = 2)
  
  abline(ls.reg$raw,
         col = pal2[i])
  
  polygon(c(ls.reg$newx, rev(ls.reg$newx)),
          c(ls.reg$predraw[, 2], rev(ls.reg$predraw[, 3])),
          border = NA,
          col = alpha(pal2[i], alpha = 0.2))
  
  legend("bottomright",
         bty = "n",
         legend = paste("Y = ", 
                        round(ls.reg$raw[[1]][2], 2),
                        "[",
                        round(ls.reg$CI[2, 1], 2),
                        ";",
                        round(ls.reg$CI[2, 2], 2),
                        "]",
                        " * X + ", 
                        round(ls.reg$raw[[1]][1], 2),
                        "[",
                        round(ls.reg$CI[1, 1], 2),
                        ";",
                        round(ls.reg$CI[1, 2], 2),
                        "]",
                        sep = ""),
         pch = 24,
         cex = 1.2)
  
  ls.reg <- stress.reg.fun(x = PCSAs[, i],
                           force = BF * meanMA)
  
  plot(1,1,
       type = "n",
       xlim = c(min(ls.reg$newlogx),
                max(ls.reg$newlogx)),
       ylim = c(min(ls.reg$predlog),
                max(ls.reg$predlog)))
  
  points(ls.reg$logx, 
         ls.reg$logforce,
         pch = 25,
         cex = 2,
         bg = pal[i], 
         lwd = 2)
  
  abline(ls.reg$loglog, 
         col = pal[i])
  
  polygon(c(ls.reg$newlogx, rev(ls.reg$newlogx)),
          c(ls.reg$predlog[, 2], rev(ls.reg$predlog[, 3])),
          border = NA,
          col = alpha(pal[i], alpha = 0.2))
  
  legend("topright",
         bty = "n",
         legend = paste("Y = ", 
                        round(ls.reg$loglog[[1]][2], 2),
                        "[",
                        round(ls.reg$CIlog[2, 1], 2),
                        ";",
                        round(ls.reg$CIlog[2, 2], 2),
                        "]",
                        " * X + ", 
                        round(ls.reg$loglog[[1]][1], 2),
                        "[",
                        round(ls.reg$CIlog[1, 1], 2),
                        ";",
                        round(ls.reg$CIlog[1, 2], 2),
                        "]",
                        sep = ""),
         pch = 25,
         cex = 1.2)
  
  legend("topleft", 
         bty = "n",
         legend = panel[i*2],
         cex = 2)
  
  ls.reg <- stress.reg.fun(x = PCSAs[, i],
                           force = BF * meanMAeff)
  
  points(ls.reg$logx, 
         ls.reg$logforce,
         pch = 24,
         cex = 2,
         bg = pal2[i],
         lwd = 2)
  
  abline(ls.reg$loglog,
         col = pal2[i])
  
  polygon(c(ls.reg$newlogx, rev(ls.reg$newlogx)),
          c(ls.reg$predlog[, 2], rev(ls.reg$predlog[, 3])),
          border = NA,
          col = alpha(pal2[i], alpha = 0.2))
  
  legend("bottomright",
         bty = "n",
         legend = paste("Y = ", 
                        round(ls.reg$loglog[[1]][2], 2),
                        "[",
                        round(ls.reg$CIlog[2, 1], 2),
                        ";",
                        round(ls.reg$CIlog[2, 2], 2),
                        "]",
                        " * X + ", 
                        round(ls.reg$loglog[[1]][1], 2),
                        "[",
                        round(ls.reg$CIlog[1, 1], 2),
                        ";",
                        round(ls.reg$CIlog[1, 2], 2),
                        "]",
                        sep = ""),
         pch = 24,
         cex = 1.2)
}

plot(1, 1,
     type = "n",
     axes = F)

text(1,1, 
     labels = "Muscle force (N)",
     cex = 3,
     srt = 90)

plot(1, 1,
     type = "n",
     axes = F)

text(1, 1, 
     labels = "PCSA (mm^2)",
     cex = 3,
     las = 1)

dev.off()

#-------------------------------------------------------------------------------
#Make tables gathering parameter estimates and their CI

intercepts_CIs_reg <- 
  slopes_CIs_reg <- 
  slopes_CIs_loglog <- 
  intercepts_CIs_loglog <- matrix(NA,
                                  ncol = 10,
                                  nrow = 3)

for (i in 1:5) {
  
  ls.reg <- stress.reg.fun(x = PCSAs[, i],
                           force = BF / meanMA)
  
  intercepts_CIs_reg[, i] <- c(ls.reg$raw[[1]][1], ls.reg$CI[1, ])
  
  slopes_CIs_reg[, i] <- c(ls.reg$raw[[1]][2], ls.reg$CI[2, ])
  
  intercepts_CIs_loglog[, i] <- c(ls.reg$loglog[[1]][1], ls.reg$CIlog[1, ])
  
  slopes_CIs_loglog[, i] <- c(ls.reg$loglog[[1]][2], ls.reg$CIlog[2, ])
  
}

for (i in 1:5) {
  
  ls.reg <- stress.reg.fun(x = PCSAs[, i],
                           force = BF / meanMAeff)
  
  intercepts_CIs_reg[, i + 5] <- c(ls.reg$raw[[1]][1], ls.reg$CI[1, ])
  
  slopes_CIs_reg[, i + 5] <- c(ls.reg$raw[[1]][2], ls.reg$CI[2, ])
  
  intercepts_CIs_loglog[, i + 5] <- c(ls.reg$loglog[[1]][1], ls.reg$CIlog[1, ])
  
  slopes_CIs_loglog[, i + 5] <- c(ls.reg$loglog[[1]][2], ls.reg$CIlog[2, ])
  
}

#intercepts_CIs_reg <- cbind(c(rawdi$coefficients[1], CIdi[1,]),
#                            c(rawvol$coefficients[1], CIvol[1,]),
#                            c(raw2D3D$coefficients[1], CI2D3D[1,]),
#                            c(rawCH$coefficients[1], CICH[1,]),
#                            c(rawinsert$coefficients[1], CIinsert[1,]))

#slopes_CIs_reg <- cbind(c(rawdi$coefficients[2], CIdi[2,]),
#                        c(rawvol$coefficients[2], CIvol[2,]),
#                        c(raw2D3D$coefficients[2], CI2D3D[2,]),
#                        c(rawCH$coefficients[2], CICH[2,]),
#                        c(rawinsert$coefficients[2], CIinsert[2,]))

#slopes_CIs_loglog <- cbind(c(loglogdi$coefficients[2], CIlogdi[2,]),
 #                          c(loglogvol$coefficients[2], CIlogvol[2,]),
  #                         c(loglog2D3D$coefficients[2], CIlog2D3D[2,]),
   #                        c(loglogCH$coefficients[2], CIlogCH[2,]),
    #                       c(logloginsert$coefficients[2], CIloginsert[2,]))

#intercepts_CIs_loglog <- cbind(c(loglogdi$coefficients[1], CIlogdi[1,]),
 #                              c(loglogvol$coefficients[1], CIlogvol[1,]),
  #                             c(loglog2D3D$coefficients[1], CIlog2D3D[1,]),
   #                            c(loglogCH$coefficients[1], CIlogCH[1,]),
    #                           c(logloginsert$coefficients[1], CIloginsert[1,]))

#-------------------------------------------------------------------------------
#Plot values and CIs against expectations
#Intercept of raw regression is expected to be 0 (a muscle with PCSA = 0 
# produces a force of 0)
#Slope of raw regression corresponds to muscle stress: should be compared to
#known values from Taylor/literature
#Slope of loglog regression is expected to be 1

# Three panel figure with 
# i) log-log regression slopes and CI against expectation (1)
# ii) raw regression intercept against expectation (0)
# iii) raw regression slope (ie stress expectation) against literature values

output_folder <- "../Reresubmission/Figures/"

pdf(file = paste(output_folder,
                 "Regressions_param_expect.pdf",
                 sep = ""),
    height = 5,
    width = 11)

layout(mat = matrix(1:3, 
                    ncol = 3))

par(mar = c(7, 5, 1, 1))

xpos <- c(0.9, 1.9, 2.9, 3.9, 4.9, 1.1, 2.1, 3.1, 4.1, 5.1)
#Plots of log log slope 

plot(x = xpos,
     y = slopes_CIs_loglog[1, ],
     ylim = c(min(slopes_CIs_loglog),
              max(slopes_CIs_loglog)),
     pch = rep(c(25, 24), each = 5),
     cex = 2,
     ylab = "Slope of log-log regression (95% C.I.)",
     xaxt = "n",
     xlab = "",
     lwd = 2,
     bg = 1)

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

plot(x = xpos,
     y = intercepts_CIs_reg[1, ],
     ylim = c(min(intercepts_CIs_reg),
              max(intercepts_CIs_reg)),
     pch = rep(c(25, 24), each = 5),
     cex = 2,
     ylab = "Intercept of regression (95% C.I.)",
     xaxt = "n",
     xlab = "",
     lwd = 2,
     bg = 1)

axis(side = 1,
     at = 1:5,
     labels = c("Dissection", 
                "Volume", 
                "2D3D", 
                "Convex Hull", 
                "Insertion area"),
     las = 2)

for (i in 1:10) {
  
  lines(x = rep(xpos[i], 2),
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

plot(x = c(xpos, 6),
     y = stresses, 
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

for (i in 1:11) {
  
  lines(x = rep(c(xpos, 6)[i], 2),
        y = c(CIs[c(1, 2), i]),
        col = "black",
        lwd = 3)
}

points(x = c(xpos, 6),
       y = stresses, 
       cex = 2,
       pch = c(rep(c(25, 24), each = 5), 20),
       lwd = 2,
       bg = 1)

abline(h = 0,
       lwd = 2,
       lty = 2)

legend("bottomright",
       legend = c("Muscle stress estimate, with 95% C.I.", 
                  "0m1 muscle stress values from the literature"), 
       pch = c(20, 22),
       lty = c(1, NA),
       pt.bg = c(NA, "grey80"),
       pt.cex = 1.5,
       cex = 0.8)

legend("topleft", 
       legend = "C.", 
       bty = "n",
       cex = 1.5)

dev.off()

#-------------------------------------------------------------------------------
#Test for size-dependency of the shrinkage
# Obtain the difference between 2D3D volume and raw volume then correlate the
# difference with size

shrinkR <- dat$PCSA_R_2d3dCH - dat$PCSA_R_vol
shrinkL <- dat$PCSA_L_2d3dCH - dat$PCSA_L_vol


shrinkR2 <- dat$PCSA_R_vol / dat$PCSA_R_2d3dCH
shrinkL2 <- dat$PCSA_L_vol / dat$PCSA_L_2d3dCH

modregR <- lm(shrinkR2 ~ dat$HW)
modregL <- lm(shrinkL2 ~ dat$HW)

summary(modregL)
summary(modregR)


layout(matrix(1:2, ncol = 2))

plot(dat$HW,shrinkR2)
abline(modregR, 
       lwd = 2)

plot(dat$HW,shrinkL2)
abline(modregL,
       lwd = 2)

#-------------------------------------------------------------------------------
# Anova of asymmetry vs measurement method

PCSA_all <- c(dat$PCSA_L_mm2, 
  dat$PCSA_L_vol, 
  dat$PCSA_L_2d3dCH,
  dat$PCSA_L_CH,
  dat$PCSA_R_mm2, 
  dat$PCSA_R_vol, 
  dat$PCSA_R_2d3dCH,
  dat$PCSA_R_CH)

PCSA_relativediff <- c(dat$PCSA_L_mm2/dat$PCSA_R_mm2,
                       dat$PCSA_L_vol/dat$PCSA_R_vol,
                       dat$PCSA_L_2d3dCH/dat$PCSA_R_2d3dCH,
                       dat$PCSA_L_CH/dat$PCSA_R_CH)

sides <- as.factor(c(rep("L", length(PCSA_all)/2),
                     rep("R", length(PCSA_all)/2)))

metho <- as.factor(c(rep("diss", 35),
                     rep("vol", 35),
                     rep("2D3D", 35),
                     rep("CH", 35)))

df.methoasym <- data.frame(PCSA_all, 
                           sides, 
                           metho)

df.metho2 <- data.frame(PCSA_relativediff,
                        metho)

df.methoasym <- na.omit(df.methoasym)

df.metho2 <- na.omit(df.metho2)

AOV.methoasym <- aov(lm(PCSA_all ~ sides * metho, 
                        data = df.methoasym))

TukeyHSD(x = AOV.methoasym)

boxplot(df.methoasym$PCSA_all ~ df.methoasym$metho)

AOV.metho2 <- aov(lm(PCSA_relativediff ~ metho, 
                     data = df.metho2))

TukeyHSD(x = AOV.metho2)

boxplot(df.metho2$PCSA_relativediff ~ df.metho2$metho)
