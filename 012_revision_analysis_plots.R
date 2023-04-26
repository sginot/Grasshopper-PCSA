# Script for revisions proposed by D. Labonte

# Make log-log regressions of PCSA vs in vivo bite force

logreg <- matrix(NA, 
                 ncol = 4, 
                 nrow = 8)

reg <- summary(lm(log(BF_dis) ~ log(mat_PCSA[, 1])))

logreg[1, ] <- c(reg$coefficients[1, 1:2],
                 reg$coefficients[2, 1:2])

reg <- summary(lm(log(BF_dis) ~ log(mat_PCSA[, 2])))

logreg[2, ] <- c(reg$coefficients[1, 1:2],
                 reg$coefficients[2, 1:2])

for (i in 1:6) {
  
  reg <- summary(lm(log(BF_3D) ~ log(mat_vols_PCSA[, i])))
  
  logreg[i + 2, ] <- c(reg$coefficients[1, 1:2],
                       reg$coefficients[2, 1:2])
}

colnames(logreg) <- c("Intercept", 
                      "Intcpt std. err.",
                      "Slope",
                      "Slope std. err.")

rownames(logreg) <- c("Dissec L",
                      "Dissec R",
                      "vol L",
                      "vol R",
                      "CH L",
                      "CH R",
                      "2D3D L",
                      "2D3D R")
