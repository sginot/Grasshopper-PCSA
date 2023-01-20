#-------------------------------------------------------------------------------

# Calculating lever arms and mechanical advantage for S. gregaria mandibles
# Data acquired by 3D landmarking in MorphoDig (Lebrun 2018)

# Required libraries
source("Rfunctions1.txt") 
# Functions from J. Claude (2008). Morphometrics with R

#-------------------------------------------------------------------------------
# Function for computing lever arms of insect mandible

# Lever arm is computed as the height of triangle formed by anterior condyle, 
# posterior condyle, and tip or insertion point.

# Function uses Heron's formula + classical formula of triangle area.

# Details of MorphoDig Landmarks used here:
# 1 or 5 = anterior condyle (left or right mandible)
# 2 or 6 = posterior condyle (left or right mandible)
# 3 or 7 = tip of mandible main incisivi (left or right mandible)
# 4 or 8 = insertion point of apodeme (left or right mandible)

levers.fun <- function(m, c1=1, c2=2, tip=3, ins=4) {
  # m is the matrix containing 3D coordinates of landmarks, the order in
  # which landmarks were digitized can be redefined with arguments 
  # c1, c2, tip and ins
  a <- dist(m[c(c1,c2),]) #Length of side between condyles (base of triangle)
  b <- dist(m[c(c1,tip),]) #Length of side between condyle 1 and tip
  c <- dist(m[c(c2,tip),]) #Length of side between condyle 2 and tip
  d <- dist(m[c(c1,ins),]) #Length of side between condyle 1 and muscle insertion
  e <- dist(m[c(c2,ins),]) #Length of side between condyle 2 and muscle insertion
  
  s1 <- (a+b+c)/2 #Semi-perimeter of triangle with tip
  s2 <- (a+d+e)/2 #Semi-perimeter of triangle with insertion
  A1 <- sqrt(s1*(s1-a)*(s1-b)*(s1-c)) #Heron's formula of triangle area
  A2 <- sqrt(s2*(s2-a)*(s2-d)*(s2-e))
  h1 <- 2*A1/a #Triangle area formula: A = (1/2)*base*height <-> height = 2*A/base
  h2 <- 2*A2/a
  
return(c(outL = h1,inL = h2))} 

#-------------------------------------------------------------------------------
# Define input folder and input files

input_folder <- "../Raws_small/Lever arms LMs"

fil <- list.files(input_folder, 
                  pattern = "lever_LMs",
                  full.names = T)

# Each file should contain the coordinates for left and right lever LMs
# Left side first, then right side, 4 LMs for each side

# /!!!\ NB: order of file names in 'fil' is not the expected, due to no padding
# with 0's in front of specimen numbers.

#-------------------------------------------------------------------------------
# Compute lever arms from individual .TPS LMs files
# /!\ will certainly not work for formats other than .TPS

# Function to import files, extract data and format it into an array.

clevers <- function(files, c1=1, c2=2, tip=3, ins=4) { 

  arr <- array(NA, 
               dim = c(2,2,length(files)),
               dimnames = list(c("Left", "Right"), 
                               c("OutL", "InL")))
  
  for (i in 1:length(files)) {
    
    o <- scan(files[i], 
              what="character")
    
    m <- matrix(as.numeric(o[2:25]), 
                ncol=3, 
                byrow=T)
    
    arr[,,i] <- matrix(c(levers.fun(m[1:4,]), 
                         levers.fun(m[5:8,])),
                          ncol=2, 
                          byrow=T)
  }
  
return(arr)}  

#-------------------------------------------------------------------------------
# Slightly different function, that will extract the coordinates, without
# transforming them into the lever arms

clevers2 <- function(files, c1=1, c2=2, tip=3, ins=4) { 
  
  arr <- array(NA, 
               dim = c(8,3,length(files)))
  
  for (i in 1:length(files)) {
    
    o <- scan(files[i], 
              what="character")
    
    m <- matrix(as.numeric(o[2:25]), 
                ncol=3, 
                byrow=T)
    
    arr[,,i] <- m
  }
  return(arr)} 
#-------------------------------------------------------------------------------
#Run the function across files, format it as matrix

levers <- clevers(files = fil)

mat.lev <- matrix(NA, 15, 4)

for (i in 1:15) {
  mat.lev[i,] <- c(levers[,,i])
}

colnames(mat.lev) <- c("OutLev_L", 
                       "OutLev_R", 
                       "InLev_L", 
                       "InLev_R")

mat.MA <- cbind(mat.lev[,3]/mat.lev[,1], 
                mat.lev[,4]/mat.lev[,2])
# Mechanical advantage (MA) = in-lever/out-lever

#-------------------------------------------------------------------------------
# As mentioned previously, the order of the files names is not correct, due to
# no padding.
# Reorder in the expected way (i.e. 1 -> 15)

mat.lev <- mat.lev[c(1,8:15,2:7),]
mat.MA <- mat.MA[c(1,8:15,2:7),]

# Some landmarks were also not digitized in the right sequence,
# and must be reordered.

mat.lev[c(4,6,9:12),] <- mat.lev[c(4,6,9:12), c(2,1,4,3)]
mat.MA[c(4,6,9:12),] <- mat.MA[c(4,6,9:12), c(2,1)]

avg_MA_sides <- apply(mat.MA, 2, mean)

#-------------------------------------------------------------------------------
# Measure angle between in- and out-lever (approx) to obtain an average
# /!\ This is NOT the pennation angle!

coords <- clevers2(files = fil)

angles <- matrix(NA, ncol=2, nrow=dim(coords)[3])

for (i in 1:dim(coords)[3]) {
  mid_axis <- apply(coords[1:2,,i], 2, mean)  # Computes the mid-point along
                                              # the axis between condyles
  tip <- coords[3,,i] - mid_axis  # Center coordinates on mid-point of axis
                                  # so resulting vector has origin of (0,0,0)
  ins <- coords[4,,i] - mid_axis  # Same origin for this vector
  
  angles[i,2] <- angle(tip,ins) * 180/pi # Compute angle and convert to degrees
  
  mid_axis <- apply(coords[5:6,,i], 2, mean) # Do the same for the other side
  
  tip <- coords[7,,i] - mid_axis
  
  ins <- coords[8,,i] - mid_axis
  
  angles[i,1] <- angle(tip,ins)*180/pi
}

angles[c(4,6,9:12),] <- angles[c(4,6,9:12), c(2,1)]

avg_mand_angles <- apply(angles, 2, FUN = mean)

avg_global <- mean(angles)
