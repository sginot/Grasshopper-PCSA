#-------------------------------------------------------------------------------

# Extracting fiber length data from .CSV files resulting from manual 
# measurements done in Fiji.
# Original images are focus stacked pictures of petri dishes containing fibers
# from S. gregaria mandible closer muscles, fixed in Bouin's solution, 
# dissected, then microdissected to separate fibers from each other.

# Required libraries
library(scales)

#-------------------------------------------------------------------------------
# Setting up input folder and files

input_folder <- "../fibers/"

input_files <- list.files(input_folder,
           pattern="Sch*",
           full.names = T)

#-------------------------------------------------------------------------------
# Importing data from tables into lists

ls <- list()

for (i in 1:length(input_files)) {
  m <- input_files[i]
  
  ls[[i]] <- read.table(m,
                        h=T, 
                        dec=".", 
                        sep=",")
}

ls_lgt <- list()

for (i in 1:length(input_files)){
  ls_lgt[[i]] <- ls[[i]]$Length
}
# The sample is composed of two groups:
# Sch16-Sch25 were fixed in Bouin, and therefore fibers could be properly 
# dissected from each other.
# Sch26-Sch30 were only fixed in ethanol, which made muscles brittle, and fibers
# impossible to individually dissect, therefore only a few measurements are
# available, and none for the Left muscle of Sch26.

ls_lgt[[2]] <- ls_lgt[[2]]*10/3.032 # This one did not have a set scale,
                                    # scale bar was measured a posteriori, and
                                    # this allows to scales the lengths.

#-------------------------------------------------------------------------------
# Summary statistics for the fiber lengths of individual muscles

mat_mean <- matrix(
              unlist(
                lapply(ls_lgt[1:20], mean)), 
              ncol = 2,
              byrow = T)
# 1:20 only, because measurements for Sch26-Sch30 are not trustworthy

mat_var <- matrix(
              unlist(
                lapply(ls_lgt[1:20], var)),
              ncol = 2, 
              byrow = T)

colnames(mat_mean) <- colnames(mat_var) <- c("L","R")

#-------------------------------------------------------------------------------
# Global vector with individual fiber measurements

L_fib <- unlist(ls_lgt[seq(from=1, to=19, by=2)])
R_fib <- unlist(ls_lgt[seq(from=2, to=20, by=2)])

#-------------------------------------------------------------------------------
# Save as R object

save(list = c("mat_mean", "mat_var", "L_fib", "R_fib"),
     file = "fiber_lengths_and_summary_stats.RData")

#-------------------------------------------------------------------------------
# Some plots (redundant with the figures script)

# Histograms
hist(L_fib, 
     xlim = c(0, 4), 
     col = alpha("orange", 0.5))

hist(R_fib,  
     col = alpha("purple", 0.35), 
     add = T)

abline(v = mean(L_fib), 
       col = "orange", 
       lwd = 3,
       lty = 1)

abline(v = mean(R_fib), 
       col = "purple", 
       lwd = 3, 
       lty = 1)

abline(v = median(L_fib), 
       col = "orange", 
       lwd = 3, 
       lty = 2)

abline(v = median(R_fib), 
       col = "purple", 
       lwd = 3, 
       lty = 2)

# Density curves
DR <- density(R_fib)
DL <- density(L_fib)

plot(DR, 
     main = "Closer muscle fiber lengths (dissection)", 
     xlab = "Length (mm)", 
     ylim = c(0,max(DL$y)))

#abline(v=median(L_fib), col="orange", lwd=3, lty=2)
#abline(v=median(R_fib), col="purple", lwd=3, lty=2)

polygon(DL$x, 
        DL$y, 
        col = alpha("darkorange",0.5), 
        lwd = 3)

polygon(DR$x, 
        DR$y, 
        col = alpha("purple",0.5), 
        lwd = 3)

abline(v = mean(L_fib), 
       col = "orange", 
       lwd = 3, 
       lty = 1)
abline(v = mean(R_fib), 
       col = "purple", 
       lwd = 3, 
       lty = 1)

# Barplots

barplot(unlist(lapply(ls_lgt[1:20], mean)), 
        col = alpha(rep(c("darkorange","darkorchid4"), 5), 0.6), 
        space = c(2,0), 
        main = "Average fiber length per muscle", 
        ylab = "Length")

axis(side = 1, 
     at = seq(from=3, to=42, by=4), 
     labels = c(paste("Sch", 16:25, sep="")),
     las = 2)

sex <- as.factor(c("M","F","F","M","F","M","M","M","M","F"))

ordsex_mean_lgt <- c(t(mat_mean[order(sex),]))

barplot(ordsex_mean_lgt, 
        col = alpha(rep(c("darkorange","darkorchid4"),5),0.6), 
        space = c(2,0), 
        main = "Average fiber length per muscle", 
        ylab = "length", 
        ylim = c(0,3.1))

axis(side = 1, 
     at = seq(from=3, to=42, by=4), 
     labels = c(paste("Sch", c(16:25)[order(sex)], sep="")),
     las = 2)

lines(c(3,15), 
      c(2.5,2.5), 
      lwd = 3)

text(9, 2.7, 
     labels = "Females")

lines(c(19,39), 
      c(2.5,2.5), 
      lwd = 3)

text(30,2.7, 
     labels="Males")

legend("topleft",
       border = "black", 
       fill = alpha(c("darkorange","darkorchid4"),0.6), 
       legend = c("Left","Right"), 
       bty = "n")

t.test(mat_mean[order(sex),1],
       mat_mean[order(sex),2], 
       paired=T)

