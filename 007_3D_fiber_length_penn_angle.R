#-------------------------------------------------------------------------------

# Extraction of fiber lengths and pennation angles from 3D landmark data
# Landmarks digitized in MorphoDig (Lebrun 2018) and exported as TPS.

# Required libraries:
source("Rfunctions1.txt")

#-------------------------------------------------------------------------------

# Fiber length extraction

#-------------------------------------------------------------------------------
# Define input folder and files

input_folder <- "../Raws_small/fibre_lenght_landmarks_correct_sides/"

input_files <- list.files(input_folder,
                          pattern = ".tps",
                          full.names = T)

#-------------------------------------------------------------------------------
# Define function to extract fiber lengths from TPS files

lgt_extract <- function(x) {
  
  a <-scan(x, 
           what="character") # Scan the text in file
  
  mat <- matrix(as.numeric(a[-c(1, length(a))]), 
                ncol=3, 
                byrow=T) # Make matrix with landmark coordinates
  
  distances <- rep(NA, dim(mat)[1]/2) # Initiate empty vector for lengths
  
  for (i in 1:(dim(mat)[1]/2)) { # dim/2 because each fiber has 2 landmarks
    distances[i] <- dist(mat[(i*2-1):(i*2),]) 
  } # Measure distance between matching pairs of landmarks
  
  return(distances)
}

#-------------------------------------------------------------------------------
# Apply function to all input tps files

list_lgt_fib <- list()

for (i in 1:length(input_files)) {
  list_lgt_fib[[i]] <- lgt_extract(input_files[i])
}

names(list_lgt_fib) <- gsub(pattern = input_folder,
                            input_files,
                            replacement = "")

av.length <- lapply(list_lgt_fib,
                    mean) # average fiber length for each muscle

#-------------------------------------------------------------------------------

# Pennation angle extraction

#-------------------------------------------------------------------------------
# Define input folder and files

input_folder <- "../Raws_small/fibre_lenght_landmarks/"

fil <- list.files(input_folder,
                  pattern = "apodeme",
                  full.names = T) 
# Files contain the same landmarks as before + two landmarks along major axis
# of the apodeme.

ls.pennat <- list()

for (i in 1:length(fil)) {
  
  TPS <- scan(file = fil[i], 
              what = "character") # Scan text in file
  
  coor <- as.numeric(TPS[-c(1,length(TPS))]) # Get coordinates
  
  matcoo <- matrix(coor, 
                   ncol=3, 
                   nrow=length(coor)/3, 
                   byrow = T) # Format coordinates as matrix
  
  Nfib <- nrow(matcoo)/2-1 # Number of fibers
  
  v2 <- matcoo[nrow(matcoo),] - matcoo[(nrow(matcoo)-1),] # apodeme vector
  
  penn <- matrix(NA, ncol = 2, nrow = Nfib) # Empty matrix for pennation angles
  
  for (j in 1:Nfib) {
    
    v1 <- matcoo[(j * 2), ] - matcoo[(j * 2 - 1), ] # Vector for jth fiber
    
    a <- angle(v1, v2) # Angle function from Claude 2008
    
    if (a > pi/2) {a <- -(a-pi)}  # Some fiber landmarks were reversed, e.g.
                                  # top one at bottom, in which case the 
                                  # direction on the vector must be reversed
    
    penn[j,1] <- a*180/pi # Convert radian to degrees
    
    penn[j,2] <- cos(a) # Cosine of pennation angle in rad
  }
  
  colnames(penn) <- c("angle_deg", "cos_angle")
  
  ls.pennat[[i]] <- penn
}


specimens <- gsub(pattern = input_folder, 
                  replacement = "",
                  x = fil)

specimens <- gsub(pattern = "_0md1_withapodeme.tps", 
                  replacement = "",
                  x = specimens) # Specimen names only
  
av.penn <- lapply(ls.pennat, apply, 2, mean) 
# Average across fibers within individual muscles

names(av.penn) <- specimens

#-------------------------------------------------------------------------------
# Make fiber length and pennation matrices, to be pasted in the main
# measurement table

matlength <- matrix(unlist(av.length), 
                    ncol = 2,
                    nrow = 15,
                    byrow = T) # Matrix with individuals as rows

matpenn <- matrix(unlist(av.penn), 
                  ncol = 4,
                  nrow = 15,
                  byrow = T) # Matrix with individuals as rows

rownames(matlength) <- 
  rownames(matpenn) <- 
  gsub(pattern = "_right", 
       replacement = "", 
       specimens)[2*(1:(length(specimens)/2))]

colnames(matlength) <- c("fiber_length_L",
                         "fiber_length_R")

colnames(matpenn) <- c("angle_deg_L",
                        "cos_angle_L",
                        "angle_deg_R",
                        "cos_angle_R")

mat.penn <- matpenn[c(1, 8:15, 2:7), ] # Increasing order of individuals
mat.length <- matlength[c(1, 8:15, 2:7), ]

mat.penn[c(1, 4, 6, 9:12),] <- mat.penn[c(1, 4, 6, 9:12), 
                                        c(3:4,1:2)] # Correct for reversed sides

save(list = c("mat.penn", "mat.length"), 
     file = "mat_fiber_lgt_angle.RData")

# These matrices are also pasted in the main measurement csv table:
# "measurments_Sch.csv"