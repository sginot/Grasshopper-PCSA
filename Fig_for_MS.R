#List of figures that we need in the MS:
#Boxplot with all PCSAs estimates, with left and right split up
#Density curves/boxplot of fiber lengths 3D vs dissect, with sides split up
#Boxplot with all predicted bite forces (sides added up) and in vivo force
#Density curve with Pennation angles, by side
#Make it into one big figure with all boxplots
#Biplot with predicted vs in vivo BF, differentiating 3D from dissected

#In addition, figure not made in R:
#Good dissection picture+fiber measurments picture
#Volume rendering of head, muscle mesh, volume, convex hull with Fiji steps

library(car)
library(MASS)
library(scales)

setwd("C:/Users/sginot/Desktop/Samuel_Projekt_PCSA")
dat <- read.csv("measurments_Sch.csv", h=T, dec=".", sep=",")

#Fb = stress*PCSA*cos(penn)*MA
#To make it comparable between 3d derived and dissect use:
#avMAL + avMAR, avcosL + avcosR, global stress, and individual PCSAs
PCSA <- (dat$PCSA_L_mm2+dat$PCSA_R_mm2)

stress <- dat$maxBF_ampbasecorr/ #in vivo BF
  (mean(dat$PCSA_L_mm2+dat$PCSA_R_mm2, na.rm = T)* #Avg total PCSA
     mean(c(dat$cos_ang_L, dat$cos_ang_R), na.rm = T)* #Avg cos(Phi)
     mean(c(dat$MA_L, dat$MA_R), na.rm = T)) #Avg MA

sigma <- mean(stress*100, na.rm=T)

avcosL <-  mean(dat$cos_ang_L, na.rm = T)
avcosR <-  mean(dat$cos_ang_R, na.rm = T)

avMAL <- mean(dat$MA_L, na.rm = T)
avMAR <- mean(dat$MA_R, na.rm = T)

avcos <- mean(c(dat$cos_ang_L, dat$cos_ang_R), na.rm = T)
avMA <- mean(c(dat$MA_R, dat$MA_L), na.rm = T)

BF <- dat$maxBF_ampbasecorr #Total BF
MF <- BF/(avcos*avMA) #Total Muscle force (both left and right)
sigma2 <- lm(MF~PCSA)$coef[2] #Stress across both muscles


#Pennation angle
matpen <- matrix(NA, ncol=30,nrow=25)
for (i in 1:30) {matpen[,i] <- ls.pennat[[i]][,1]} 
colnames(matpen) <- fil
#for ls.pennat and fil see "script_pennation_angle"
matpenL <- matpen[,grep("left", fil)]
matpenR <- matpen[,grep("right", fil)]
matpenL <- matpenL[, c(1,8:15,2:7)]
matpenR <- matpenR[, c(1,8:15,2:7)]
Rev <- as.logical(dat$Rev_side_3D[1:15])
matpenL.OK <- matpenL
matpenR.OK <- matpenR
matpenL.OK[,Rev] <- matpenR[,Rev] 
matpenR.OK[,Rev] <- matpenL[,Rev]
#These lines to put back the reversed sides to their original true side 
DRP <- density(matpenR.OK)
DLP <- density(matpenL.OK)

setwd("C:/Users/sginot/Desktop/Raws_small/fibre_lenght_landmarks_correct_sides/")
#Run "script_fiber_lgt_3D.R"
fibL <- unlist(list_lgt_fib[grep("L_", tps_files)])
fibR <- unlist(list_lgt_fib[grep("R_", tps_files)])
DL3 <- density(fibL)
DR3 <- density(fibR)

setwd("C:/Users/sginot/Desktop/Samuel_Projekt_PCSA/fibers/")
#Run fibers_script.R

#To have proper PCSAs, all estimates should include pennation angle!
PCSA_L_mm2 <- (dat$L_closer_g*avcosL)/(0.00106*dat$L_closer_fib_l_mm)
PCSA_R_mm2 <- (dat$R_closer_g*avcosR)/(0.00106*dat$R_closer_fib_l_mm)
PCSA_L_vol <- (dat$L_closer_vol_mm3*dat$cos_ang_L)/(dat$L_closer_fib_l_mm)
PCSA_R_vol <- (dat$R_closer_vol_mm3*dat$cos_ang_R)/(dat$R_closer_fib_l_mm)
PCSA_L_2d3dCH <- (dat$L_2d3d_CH*dat$cos_ang_L)/(dat$L_closer_fib_l_mm)
PCSA_R_2d3dCH <- (dat$R_2d3d_CH*dat$cos_ang_R)/(dat$R_closer_fib_l_mm)
PCSA_L_CH <- (dat$L_convexhull_vol_mm3*dat$cos_ang_L)/(dat$L_closer_fib_l_mm)
PCSA_R_CH <- (dat$R_convexhull_vol_mm3*dat$cos_ang_R)/(dat$R_closer_fib_l_mm)


predBF_PCSA <-  (sigma2*avMAL*PCSA_L_mm2)+
                (sigma2*avMAR*PCSA_R_mm2)
predBF_vol <-   (sigma2*avMAL*PCSA_L_vol)+
                (sigma2*avMAR*PCSA_R_vol)
predBF_2d3d <-  (sigma2*avMAL*PCSA_L_2d3dCH)+
                (sigma2*avMAR*PCSA_R_2d3dCH)
predBF_CH <-    (sigma2*avMAL*PCSA_L_CH)+
                (sigma2*avMAR*PCSA_R_CH)
predBF_insert <-(sigma2*avMA*dat$insert_area_total)


#Recalculate the whole thing for when mandibles are open (changing the effective in lever)

ILL <- dat$InLev_L
ILR <- dat$InLev_R

ILLeff <- ILL*sin(40*pi/180)
ILReff <- ILR*sin(40*pi/180)

MAReff <- ILReff/dat$OutLevR
MALeff <- ILLeff/dat$OutLev_L

avMALeff <- mean(MALeff, na.rm = T)
avMAReff <- mean(MAReff, na.rm = T)
avMAeff <- mean(c(MAReff, MALeff), na.rm = T)


MF2 <- BF/(avcos*avMAeff) #Total Muscle force (both left and right)
sigma3 <- lm(MF2~PCSA)$coef[2]

predBF_PCSA2 <- (sigma2*avMALeff*PCSA_L_mm2)+
                (sigma2*avMAReff*PCSA_R_mm2)
predBF_vol2 <-  (sigma2*avMALeff*PCSA_L_vol)+
                (sigma2*avMAReff*PCSA_R_vol)
predBF_2d3d2 <- (sigma2*avMALeff*PCSA_L_2d3dCH)+
                (sigma2*avMAReff*PCSA_R_2d3dCH)
predBF_CH2 <-   (sigma2*avMALeff*PCSA_L_CH)+
                (sigma2*avMAReff*PCSA_R_CH)
predBF_insert2 <-(sigma2*avMAeff*dat$insert_area_total)

#### PCSA boxplots ####
par(mar = c(6,4,1,1), xpd=T)
boxplot(dat$PCSA_L_mm2,
        dat$PCSA_R_mm2,
        dat$PCSA_L_vol,
        dat$PCSA_R_vol,
        dat$PCSA_L_2d3dCH,
        dat$PCSA_R_2d3dCH,
        dat$PCSA_L_CH,
        dat$PCSA_R_CH,
        dat$insert_area_total/2,
        col = c(rep(c("darkorange", "darkorchid4"),4),"grey"),
        ylim = c(0,19),
        pch = 20,
        cex=2,
        lwd=2,
        boxwex = 0.5,
        at = c(1,1.5, 3,3.5, 5,5.5, 7,7.5, 9.25),
        xaxt = "n",
        ylab = "Area (mm^2)"
        )
axis(side = 1,
     at = c(1.25, 3.25, 5.25, 7.25, 9.25),
     labels = F
     )
text(x = c(0.5, 2.5, 4.5, 6.5, 8.5),
     y = -2.5,
      labels = c("PCSA", 
                  "PCSA 3D",
                  "PCSA 2D-3D",
                  "PCSA 3D",
                  "PCSA"),
      srt = 45)
text(x = c(1, 3, 5, 7, 9),
     y = -2.5,
     labels = c("dissection", 
                "volume",
                "convex hull",
                "convex hull",
                "muscle insertion"),
     srt = 45)
abline(v=2.5, lwd=2, lty=2, xpd=F)
text(x=c(1.25,6),
     y=19,
     labels = c("Traditional", "CT scan-derived"),
     font = 2
     )
text(x=c(1,1.5, 3,3.5, 5,5.5, 7,7.5,9.25),
     y=0,
     labels = c(14,15,15,15,14,14,15,15,15),
     col = c(rep(c("darkorange", "darkorchid4"),4),"black"),
     font = 2
)

#### Bite forces boxplots ####

pdf("boxplots_BF.pdf")
par(mar = c(6,4,1,1), xpd=T)
boxplot(BF,
        predBF_PCSA,predBF_PCSA2,
        predBF_vol,predBF_vol2,
        predBF_2d3d,predBF_2d3d2,
        predBF_CH,predBF_CH2,
        predBF_insert,predBF_insert2,
        col = c("firebrick", rep(c("grey", "bisque2"),5)),
        pch = 20,
        cex=2,
        lwd=2,
        boxwex = 0.3,
        at = c(1, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2),
        xaxt = "n",
        ylab = "Bite force (N)"
)
polygon(x = c(0,10,10,0), 
        y = c(max(BF, na.rm = T), 
              max(BF, na.rm = T),
              min(BF, na.rm = T),
              min(BF, na.rm = T)),
        col = alpha("firebrick", alpha = 0.2),
        xpd = F,
        border = T)
axis(side = 1,
     at = c(1:6),
     labels = F
)
text(x = c(1:6),
     y = -0.5,
     labels = c("In vivo",
                "Dissection", 
                "3D volume",
                "2D-3D convex hull",
                "3D convex hull",
                "Muscle insertion"),
     srt = 45)
abline(v=1.5, lwd=2, lty=2, xpd=F)
abline(v=2.5, lwd=2, lty=2, xpd=F)

text(x=c(0.9,2,4.5),
     y=6.5,
     labels = c("Reference", "Traditional", "CT scan-derived"),
     font = 2,
     cex = 0.8
)
dev.off()

####Pennation angle density####

plot(DRP, 
     main="Pennation angle distribution", 
     xlab="Angle (degrees)", 
     ylim=c(0,max(DLP$y)))
polygon(DLP$x, DLP$y, col=alpha("darkorange",0.5), lwd=3)
polygon(DRP$x, DRP$y, col=alpha("darkorchid4",0.5), lwd=3)
abline(v=mean(DRP$x), col = "darkorchid4", lty=2, lwd=2)
abline(v=mean(DLP$x), col = "darkorange", lty=2, lwd=2)

####Fiber length density (from dissections)####

plot(DR, 
     main="Closer muscle fiber lengths (dissection)",
     xlab="Length (mm)", 
     ylim=c(0,max(DL$y)))
polygon(DL$x, DL$y, col=alpha("darkorange",0.5), lwd=3)
polygon(DR$x, DR$y, col=alpha("purple",0.5), lwd=3)
abline(v=mean(L_fib), col="orange", lwd=3, lty=1)
abline(v=mean(R_fib), col="purple", lwd=3, lty=1)

####Fiber length density (from 3D)####

plot(DR3, 
     main="Closer muscle fiber lengths (3D)",
     xlab="Length (mm)", 
     ylim=c(0,max(DL3$y)))
polygon(DL3$x, DL3$y, col=alpha("darkorange",0.5), lwd=3)
polygon(DR3$x, DR3$y, col=alpha("purple",0.5), lwd=3)
abline(v=mean(fibL), col="orange", lwd=3, lty=1)
abline(v=mean(fibR), col="purple", lwd=3, lty=1)

####Fiber length  (from 2D)####

#Boxplot for more easy reading
boxplot(L_fib, R_fib, fibL, fibR,
        col = c("darkorange", "darkorchid4"),
        pch = 20,
        cex=2,
        lwd=2,
        at = c(0.7, 1.3, 2.4,3),
        boxwex = 0.5,
        xaxt = "n",
        ylab = "Fiber length (mm)"
        )
axis(side = 1,
     at = c(1,2.7),
     labels = c("Dissection", "3D (landmarks)"))

####Everything combined####
setwd("C:/Users/sginot/Desktop/Samuel_Projekt_PCSA/")

pdf("Main_figure.pdf", height = 10, width = 5)

layout(matrix(c(1:2,3,3,4,4), ncol=2, nrow=3, byrow = T))
par(mar=c(4,4.5,1,1))

boxplot(L_fib, R_fib, fibL, fibR,
        col = c("darkorange", "darkorchid4"),
        pch = 20,
        cex=2,
        lwd=2,
        at = c(0.7, 1.3, 2.4,3),
        boxwex = 0.5,
        xaxt = "n",
        ylab = "Fiber length (mm)", cex.lab=1.5)

axis(side = 1,
     at = c(1,2.7),
     labels = F)
text(x = c(1,2.7), y = 0.3,
     labels = c("Dissection", "3D (landmarks)"), 
     cex=1.5,
     srt=20,
     xpd=T)
text(x = 0.5, y = 3.7,
     labels = "A.", 
     cex=2,
     font=2)

plot(DRP, 
     main="", 
     xlab="Pennation angle (degrees)", 
     ylim=c(0,max(DLP$y)), 
     cex.lab=1.5)
polygon(DLP$x, DLP$y, col=alpha("darkorange",0.5), lwd=3)
polygon(DRP$x, DRP$y, col=alpha("darkorchid4",0.5), lwd=3)
abline(v=mean(DRP$x), col = "darkorchid4", lty=2, lwd=2)
abline(v=mean(DLP$x), col = "darkorange", lty=2, lwd=2)
text(x = 0, y = 0.042,
     labels = "B.", 
     cex=2,
     font=2)

par(mar=c(5,4.5,1,1))

boxplot(dat$PCSA_L_mm2,
        dat$PCSA_R_mm2,
        dat$PCSA_L_vol,
        dat$PCSA_R_vol,
        dat$PCSA_L_2d3dCH,
        dat$PCSA_R_2d3dCH,
        dat$PCSA_L_CH,
        dat$PCSA_R_CH,
        dat$insert_area_total/2,
        col = c(rep(c("darkorange", "darkorchid4"),4),"grey"),
        ylim = c(0,19),
        pch = 20,
        cex=2,
        lwd=2,
        boxwex = 0.5,
        at = c(1,1.5, 3,3.5, 5,5.5, 7,7.5, 9.25),
        xaxt = "n",
        ylab = "Area (mm^2)", cex.lab=1.5
)
axis(side = 1,
     at = c(1.25, 3.25, 5.25, 7.25, 9.25),
     labels = F
)
text(x = c(0.5, 2.5, 4.5, 6.5, 8.5),
     y = -2.8,
     labels = c("PCSA", 
                "PCSA 3D",
                "PCSA 2D-3D",
                "PCSA 3D",
                "PCSA"),
     srt = 20,
     xpd = T, 
     cex=1.5)
text(x = c(1, 3, 5, 7, 9),
     y = -3.2,
     labels = c("dissection", 
                "volume",
                "convex hull",
                "convex hull",
                "muscle insertion"),
     srt = 20,
     xpd = T, 
     cex=1.5)
abline(v=2.5, lwd=2, lty=2, xpd=F)
text(x = 0.5, y = 18,
     labels = "C.", 
     cex=2,
     font=2)
boxplot(BF,
        predBF_PCSA,predBF_PCSA2,
        predBF_vol,predBF_vol2,
        predBF_2d3d,predBF_2d3d2,
        predBF_CH,predBF_CH2,
        predBF_insert,predBF_insert2,
        col = c("firebrick", rep(c("grey", "bisque2"),5)),
        pch = 20,
        cex=2,
        lwd=2,
        boxwex = 0.3,
        at = c(1, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2),
        xaxt = "n",
        ylab = "Bite force (N)", 
        cex.lab=1.5
)
polygon(x = c(0,10,10,0), 
        y = c(max(BF, na.rm = T), 
              max(BF, na.rm = T),
              min(BF, na.rm = T),
              min(BF, na.rm = T)),
        col = alpha("firebrick", alpha = 0.2),
        xpd = F,
        border = T)
axis(side = 1,
     at = c(1:6),
     labels = F
)
text(x = c(1:6),
     y = -0.5,
     labels = c("In vivo",
                "Dissection", 
                "3D volume",
                "2D-3D hull",
                "3D convex hull",
                "Muscle insertion"),
     srt = 20,
     xpd = T, 
     cex=1.5)
abline(v=1.5, lwd=2, lty=2, xpd=F)
abline(v=2.5, lwd=2, lty=2, xpd=F)
text(x = 0.5, y = 6,
     labels = "D.", 
     cex=2,
     font=2)

dev.off()

###Same as before but with PCSA including pennation angles
setwd("C:/Users/sginot/Desktop/Samuel_Projekt_PCSA/")

pdf("Main_figure2.pdf", height = 10, width = 5)

layout(matrix(c(1:2,3,3,4,4), ncol=2, nrow=3, byrow = T))
par(mar=c(4,4.5,1,1))

boxplot(L_fib, R_fib, fibL, fibR,
        col = c("darkorange", "darkorchid4"),
        pch = 20,
        cex=2,
        lwd=2,
        at = c(0.7, 1.3, 2.4,3),
        boxwex = 0.5,
        xaxt = "n",
        ylab = "Fiber length (mm)", cex.lab=1.5)

axis(side = 1,
     at = c(1,2.7),
     labels = F)
text(x = c(1,2.7), y = 0.3,
     labels = c("Dissection", "3D (landmarks)"), 
     cex=1.5,
     srt=20,
     xpd=T)
text(x = 0.5, y = 3.7,
     labels = "A.", 
     cex=2,
     font=2)

plot(DRP, 
     main="", 
     xlab="Pennation angle (degrees)", 
     ylim=c(0,max(DLP$y)), 
     cex.lab=1.5)
polygon(DLP$x, DLP$y, col=alpha("darkorange",0.5), lwd=3)
polygon(DRP$x, DRP$y, col=alpha("darkorchid4",0.5), lwd=3)
abline(v=mean(DRP$x), col = "darkorchid4", lty=2, lwd=2)
abline(v=mean(DLP$x), col = "darkorange", lty=2, lwd=2)
text(x = 0, y = 0.042,
     labels = "B.", 
     cex=2,
     font=2)

par(mar=c(5,4.5,1,1))

boxplot(PCSA_L_mm2,
        PCSA_R_mm2,
        PCSA_L_vol,
        PCSA_R_vol,
        PCSA_L_2d3dCH,
        PCSA_R_2d3dCH,
        PCSA_L_CH,
        PCSA_R_CH,
        dat$insert_area_total/2,
        col = c(rep(c("darkorange", "darkorchid4"),4),"grey"),
        ylim = c(0,19),
        pch = 20,
        cex=2,
        lwd=2,
        boxwex = 0.5,
        at = c(1,1.5, 3,3.5, 5,5.5, 7,7.5, 9.25),
        xaxt = "n",
        ylab = "Area (mm^2)", cex.lab=1.5
)
axis(side = 1,
     at = c(1.25, 3.25, 5.25, 7.25, 9.25),
     labels = F
)
text(x = c(0.5, 2.5, 4.5, 6.5, 8.5),
     y = -2.8,
     labels = c("PCSA", 
                "PCSA 3D",
                "PCSA 2D-3D",
                "PCSA 3D",
                "PCSA"),
     srt = 20,
     xpd = T, 
     cex=1.5)
text(x = c(1, 3, 5, 7, 9),
     y = -3.2,
     labels = c("dissection", 
                "volume",
                "convex hull",
                "convex hull",
                "muscle insertion"),
     srt = 20,
     xpd = T, 
     cex=1.5)
abline(v=2.5, lwd=2, lty=2, xpd=F)
text(x = 0.5, y = 18,
     labels = "C.", 
     cex=2,
     font=2)
boxplot(BF[which(dat$Dissected==1)], BF[which(dat$Dissected==0)],
        predBF_PCSA,predBF_PCSA2,
        predBF_vol,predBF_vol2,
        predBF_2d3d,predBF_2d3d2,
        predBF_CH,predBF_CH2,
        predBF_insert,predBF_insert2,
        col = c("red", "firebrick", rep(c("grey", "bisque2"),5)),
        pch = 20,
        cex=2,
        lwd=2,
        boxwex = 0.3,
        at = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2),
        xaxt = "n",
        ylab = "Bite force (N)", 
        cex.lab=1.5
)
polygon(x = c(0,10,10,0), 
        y = c(max(BF, na.rm = T), 
              max(BF, na.rm = T),
              min(BF, na.rm = T),
              min(BF, na.rm = T)),
        col = alpha("firebrick", alpha = 0.2),
        xpd = F,
        border = T)
axis(side = 1,
     at = c(1:6),
     labels = F
)
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
     cex=1.5)
abline(v=1.5, lwd=2, lty=2, xpd=F)
abline(v=2.5, lwd=2, lty=2, xpd=F)
text(x = 0.5, y = 7,
     labels = "D.", 
     cex=2,
     font=2)

dev.off()

####Bite forces correlations biplot####

pdf("Correlations_fig.pdf", height = 5, width = 10)
layout(matrix(1:2, ncol=2))
par(mar = c(5,5,1,1))

plot(BF, 
     predBF_PCSA, 
     ylim = c(0,9), 
     pch=21, 
     col = "black", 
     bg = "red", 
     cex = 2,
     lwd=2,
     xlab = "In vivo bite force (N)",
     ylab = "Predicted bite force (N)",
     cex.lab = 1.5)
clip(min(BF[16:30]),max(BF[16:30]), 0, 9)
abline(lm(predBF_PCSA ~ BF), col="red", lwd = 3)
clip(0,4,0,9)

points(BF, 
       predBF_vol,
       pch = 22,
       bg = "blue",
       cex = 2,
       lwd=2)
points(BF, 
       predBF_2d3d,
       pch = 23,
       bg = "purple",
       cex = 2,
       lwd=2)
points(BF, 
       predBF_CH,
       pch = 24,
       bg = "darkorange",
       cex = 2,
       lwd=2)
points(BF, 
       predBF_insert,
       pch = 25,
       bg = "forestgreen",
       cex = 2,
       lwd=2)
text("n.s.", x=2,y=0.5, col="blue")

clip(min(BF[1:15], na.rm = T),max(BF[1:15], na.rm = T), 0, 9)
abline(lm(predBF_vol ~ BF), col="blue", lwd = 3)
abline(lm(predBF_2d3d ~ BF), col="purple", lwd = 3)
abline(lm(predBF_CH ~ BF), col="darkorange", lwd = 3)
abline(lm(predBF_insert ~ BF), col="forestgreen", lwd = 3)

#Do the same with values centered at intercept = 0
par(mar = c(5,1.5,1,1))

inter <- lm(predBF_PCSA ~ BF)$coef[1]
plot(BF, 
     predBF_PCSA-inter, 
     ylim = c(-0.2,2.5), 
     pch=21, 
     #col = alpha("black", 0.5), 
     bg = alpha("red",0.5), 
     cex = 2,
     lwd=2,
     xlab = "In vivo bite force (N)",
     cex.lab = 1.5)
clip(min(BF[16:30]),max(BF[16:30]), 0, 9)
abline(lm(predBF_PCSA-inter ~ BF), col="red", lwd = 3)

clip(0,4,-1,7)
inter <- lm(predBF_vol ~ BF)$coef[1]
points(BF, 
       predBF_vol-inter,
       pch = 22,
       #col = alpha("black", 0.5),
       bg = alpha("blue",0.5),
       cex = 2,
       lwd=2)
text("n.s.", x=2,y=0, col="blue")
clip(min(BF[1:15], na.rm = T),max(BF[1:15], na.rm = T), 0, 9)
abline(lm(predBF_vol-inter ~ BF), col="blue", lwd = 3)

clip(0,4,-1,7)
inter <- lm(predBF_2d3d ~ BF)$coef[1]
points(BF, 
       predBF_2d3d-inter,
       pch = 23,
       #col = alpha("black", 0.5),
       bg = alpha("purple", 0.5),
       cex = 2,
       lwd=2)
clip(min(BF[1:15], na.rm = T),max(BF[1:15], na.rm = T), 0, 9)
abline(lm(predBF_2d3d-inter ~ BF), col="purple", lwd = 3)

clip(0,4,-1,7)
inter <- lm(predBF_CH ~ BF)$coef[1]
points(BF, 
       predBF_CH-inter,
       pch = 24,
       #col = alpha("black", 0.5),
       bg = alpha("darkorange", 0.5),
       cex = 2,
       lwd=2)
clip(min(BF[1:15], na.rm = T),max(BF[1:15], na.rm = T), 0, 9)
abline(lm(predBF_CH-inter ~ BF), col="darkorange", lwd = 3)

clip(0,4,-1,7)
inter <- lm(predBF_insert ~ BF)$coef[1]
points(BF, 
       predBF_insert-inter,
       pch = 25,
       #col = alpha("black", 0.5),
       bg = alpha("forestgreen", 0.5),
       cex = 2,
       lwd=2)
clip(min(BF[1:15], na.rm = T),max(BF[1:15], na.rm = T), 0, 9)
abline(lm(predBF_insert-inter ~ BF), col="forestgreen", lwd = 3)

clip(-5,8,-1 ,10)

legend(0.1,2.7, 
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


####Bite forces correlations biplot X-Y axes inverted####

pdf("Correlations_fig_invert.pdf", width = 7, height = 7)
layout(matrix(1:2, nrow=2))
par(mar = c(2,5,1,1))

plot(predBF_PCSA, 
     BF,
     xlim = c(0,8),
     ylim = c(0.2,2.2), 
     pch=21, 
     col = "black", 
     bg = "red", 
     cex = 2,
     lwd=2,
     xaxt = "n",
     xlab = "",
     ylab = "In vivo bite force (N)",
     cex.lab = 1.5)
clip(0, 9, min(BF[16:30]),max(BF[16:30]))
abline(lm(BF~predBF_PCSA), col="red", lwd = 3)
clip(0,9, 0.2 ,2.2)
abline(0,1, lty=2, lwd = 1.5, xpd=F)
clip(0,9,0,3)

points(predBF_vol,
       BF,
       pch = 22,
       bg = "blue",
       cex = 2,
       lwd=2)
points(predBF_2d3d,
       BF,
       pch = 23,
       bg = "purple",
       cex = 2,
       lwd=2)
points(predBF_CH,
       BF, 
       pch = 24,
       bg = "darkorange",
       cex = 2,
       lwd=2)
points(predBF_insert,
       BF,
       pch = 25,
       bg = "forestgreen",
       cex = 2,
       lwd=2)

clip(0, 9, min(BF[1:15], na.rm = T),max(BF[1:15], na.rm = T))
abline(lm(BF~predBF_2d3d), col="purple", lwd = 3)
abline(lm(BF~predBF_CH), col="darkorange", lwd = 3)
abline(lm(BF~predBF_insert), col="forestgreen", lwd = 3)

#Do the same with values centered at intercept = 0
par(mar = c(5,5,0,1))

inter <- lm(BF~predBF_PCSA)$coef[1]
oo<-BF-inter

plot(predBF_PCSA, 
     oo,
     ylim = c(-0.2,3.5), 
     xlim = c(0,8),
     pch=21, 
     #col = alpha("black", 0.5), 
     bg = alpha("red",0.5), 
     cex = 2,
     lwd=2,
     ylab = "In vivo bite force (N)",
     xlab = "Modeled bite force (N)",
     cex.lab = 1.5)
#clip(-0.2, 8, min(oo[16:30]),max(oo[16:30]))
clip(-0.2, 10, 0, 3.5)
abline(lm(oo~predBF_PCSA), col="red", lwd = 3)

clip(-0.2, 10, 0, 3.5)
inter <- lm(predBF_vol ~ BF)$coef[1]
oo<-BF-inter
points(predBF_vol,
       oo,
       pch = 22,
       #col = alpha("black", 0.5),
       bg = alpha("blue",0.5),
       cex = 2,
       lwd=2)
       
clip(-0.2, 10, 0, 3.5)
inter <- lm(BF~predBF_2d3d)$coef[1]
oo<-BF-inter
points(predBF_2d3d,
       oo,
       pch = 23,
       #col = alpha("black", 0.5),
       bg = alpha("purple", 0.5),
       cex = 2,
       lwd=2)
#clip(min(predBF_2d3d, na.rm=T), max(predBF_2d3d, na.rm=T), 0,2.5)
abline(lm(oo~predBF_2d3d), col="purple", lwd = 3)

clip(-0.2, 10, 0, 3.5)
inter <- lm(BF~predBF_CH)$coef[1]
oo<-BF-inter
points(predBF_CH,
       oo,
       pch = 24,
       #col = alpha("black", 0.5),
       bg = alpha("darkorange", 0.5),
       cex = 2,
       lwd=2)
#clip(min(oo, na.rm=T), max(oo, na.rm=T), 0,2.5)
abline(lm(oo~predBF_CH), col="darkorange", lwd = 3)

clip(-0.2, 10, 0, 3.5)
inter <- lm(BF~predBF_insert)$coef[1]
oo<-BF-inter
points(predBF_insert,
       oo,
       pch = 25,
       #col = alpha("black", 0.5),
       bg = alpha("forestgreen", 0.5),
       cex = 2,
       lwd=2)
#clip(min(oo, na.rm=T), max(oo, na.rm=T), 0,2.5)
abline(lm(oo~predBF_insert), col="forestgreen", lwd = 3)
abline(0,1, lty=2, lwd=2)

legend("bottomright",
       legend = c("Mesh volume estim.",
                  "Dissection estim.",
                  "2D-3D convex hull estim.",
                  "3D convex hull estim.",
                  "Insertion area estim."),
       pch = c(22,21,23:35),
       pt.bg = c("blue", 
                 "red", 
                 "purple", 
                 "darkorange", 
                 "forestgreen"),
       bty = "o",
       bg = "white",
       xpd=T)

dev.off()
#### Allometry biplot BF vs HW

plot(dat$HW, 
     dat$maxBF_ampbasecorr, 
     bg=c(4,2)[as.factor(dat$Sex)],
     pch=c(21,22)[as.factor(dat$Sex)],
     cex=2,
     lwd=2,
     xlab="Head width (mm)",
     ylab = "In vivo bite force (N)")

clip(min(dat$HW[which(dat$Sex=="M")], na.rm = T),
     max(dat$HW[which(dat$Sex=="M")], na.rm = T),
     min(dat$maxBF_ampbasecorr[which(dat$Sex=="M")], na.rm = T),
     max(dat$maxBF_ampbasecorr[which(dat$Sex=="M")], na.rm = T) )
abline(lm(maxBF_ampbasecorr ~ HW, data=dat[which(dat$Sex=="M"),]),
       col=2, 
       lwd = 3)
clip(min(dat$HW[which(dat$Sex=="F")], na.rm = T),
     max(dat$HW[which(dat$Sex=="F")], na.rm = T),
     min(dat$maxBF_ampbasecorr[which(dat$Sex=="F")], na.rm = T),
     max(dat$maxBF_ampbasecorr[which(dat$Sex=="F")], na.rm = T) )
abline(lm(maxBF_ampbasecorr ~ HW, data=dat[which(dat$Sex=="F"),]),
       col=4,
       lwd=3)
