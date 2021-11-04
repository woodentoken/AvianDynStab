library(alphahull) # need for the convex hull
library(ptinpoly)  # need for determining points outside the convex hull
library(pracma)
source("support_functions.R")

# ------------ Read in data ---------------
### CHECK THAT THE TRIM FORCES ARE ACTUALLY ZERO
# Note that this does not include the wing shapes from the inertial gull
dat_all <- read.csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2020_05_25_OrientedWings.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
dat_all <- subset(dat_all, species == "lar_gla" & sweep == 0 & dihedral == 0)

dat_out <- read.csv('2021_11_04_LongDynStability_Rigid.csv')
S_max = 0.244657  # wings and body reference area from the gull used in inertial study
c_max = 0.2861011  # wing root chord from the gull used in inertial study

# ------------ Clean data ---------------

# remove the unattainable configurations
tmp     <- cut_trueshape(dat_out,unique(dat_all[4:5]),4,5) # cut elbow and wrist to the true shape of the data
dat_cut <- tmp$dat_cut
# restrict to a trustable angle of attack
# ADD NOTE THAT THIS DOESN'T MEAN BIRDS CAN TRIM AT SLOWER SPEEDS JUST THAT MY AERO DATA IN THAT RANGE ISN'T AS TRUSTWORTHY
dat_cut <- subset(dat_cut, alpha < 5)
# Limit so that the wind direction isn't more than 45deg relative to the horizontal
dat_cut <- subset(dat_cut, gamma0 > -45*pi/180)

dat_cut$halft  = 0.69/abs(dat_cut$eig_real) # captures the rate of decay
dat_cut$period = 2*pi/dat_cut$omega_n
dat_cut$sm = -(dat_cut$Cm_alp/dat_cut$CL_alp)*c_max
# Set the pitch rate as the comparable phase 
dat_cut$phase1 = dat_cut$phase1 - dat_cut$phase3
dat_cut$phase2 = dat_cut$phase2 - dat_cut$phase3
dat_cut$phase3 = dat_cut$phase3 - dat_cut$phase3
dat_cut$phase4 = dat_cut$phase4 - dat_cut$phase3
# adjust negative numbers to be in a positive frame
dat_cut$phase1[which(dat_cut$phase1 < 0)] = (2*pi) + dat_cut$phase1[which(dat_cut$phase1 < 0)]
dat_cut$phase2[which(dat_cut$phase2 < 0)] = (2*pi) + dat_cut$phase2[which(dat_cut$phase2 < 0)]
dat_cut$phase3[which(dat_cut$phase3 < 0)] = (2*pi) + dat_cut$phase3[which(dat_cut$phase3 < 0)]
dat_cut$phase4[which(dat_cut$phase4 < 0)] = (2*pi) + dat_cut$phase4[which(dat_cut$phase4 < 0)]

dat_sp <- subset(dat_cut, eignum == 1 & omega_n > 0) # eigen num == 2 just flips phase
dat_ph <- subset(dat_cut, eignum == 3 & omega_n > 0) # eigen num == 4 just flips phase
# need to make sure that the configurations that are unstable in sp mode are removed here as well - why?
dat_ph <- merge(dat_sp[,c("alpha","elbow","manus")], subset(dat_cut, eignum == 3 & omega_n > 0))
