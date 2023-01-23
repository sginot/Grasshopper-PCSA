#-------------------------------------------------------------------------------

# Combine data from Taylor 2001 with our sarcomere length data to predict
# muscle stress value for our sample.

#-------------------------------------------------------------------------------
# Load important objects

load("mat_sarcomere_lengths_summary.RData")
load("list_sarcomere_length_measurements.RData")
load("models.RData")

#-------------------------------------------------------------------------------
# Use the model with log values, as because only log values will have a truly
# linear relationship.

mean_sarco <- mean(mat[,1])

log_sarco <- log(mean_sarco)

intercept <- modlogs$coefficients[1]

slope <- modlogs$coefficients[2]

pred_log_stress <- slope * log_sarco + intercept  # Predicted value in the log-
                                                  # log linear regression model

pred_stress <- exp(pred_log_stress) # Exponential to convert back to real scale
                                    # value in kN.m^-2

pred_stress_N_cm2 <- pred_stress * 1000 / 10000 # *1000 converts to N
                                                # /10000 converts to cm^-2

save(file = "muscle_stress.RData",
     list = c("pred_stress", "pred_stress_N_cm2"))
# This independently computed average muscle stress value can be used to 
# further compute muscle force estimates. Since only a few individuals of the
# sample could have their sarcomere length measured, the average value will be
# used for all specimens.
