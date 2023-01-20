#-------------------------------------------------------------------------------

# Re-analysis of Taylor 2001 data on muscle stress and sarcomere lenght

#-------------------------------------------------------------------------------

# Load data
input_folder <- "../raw_data_sarcomere/"

input_files <- list.files(input_folder, 
                          pattern = "sarcomere",
                          full.name = T)

dat_taylor <- read.csv(input_files,
                       h = T,
                       sep = "\t",
                       dec = ".")
head(dat_taylor)

#-------------------------------------------------------------------------------
# Define variables of interest

sarco_lgt <- dat_taylor$Mean
stress <- dat_taylor$Mean.1

#-------------------------------------------------------------------------------
# Analyse as done by Taylor 2001

mod <- lm(stress ~ sarco_lgt)
modlogs <- lm(log(stress) ~ log(sarco_lgt))

#-------------------------------------------------------------------------------
# Save models for further analysis

save(file = "models.RData", list = c("mod", "modlogs"))

#-------------------------------------------------------------------------------
# Plots to check

layout(matrix(1:2, ncol = 2))

plot(sarco_lgt, 
     stress, 
     pch = 19, 
     cex = 1.5,
     xlab = "Sarcomere lenght (µm)",
     ylab = "Muscle stress (kN.m^-2)")

abline(mod, 
       lwd = 2)

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

text(min(log(sarco_lgt), na.rm = T), 
     max(log(stress), na.rm = T), 
     labels = paste("Y =", 
                    round(modlogs$coefficients[2], 2), 
                    "* x +", 
                    round(modlogs$coefficients[1], 2)),
     pos = 4)

# Optional, save the figure
dev.copy2pdf(file = "../Manuscript/Figures/sarco_stress_regressions.pdf")
