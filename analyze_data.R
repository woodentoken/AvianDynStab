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
dat_cut$alpha_rad = dat_cut$alpha*pi/180
dat_cut$theta_rad = dat_cut$alpha_rad + dat_cut$gamma0 # angle from the horizon to the bird reference line

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
# need to make sure that the configurations that are unstable in sp mode are removed here as well
dat_ph <- merge(dat_sp[,c("alpha","elbow","manus")], subset(dat_cut, eignum == 3 & omega_n > 0))

## --------- Investigate the effect of morphing on stability --------

mod_sp_zeta    = lm(zeta ~ elbow*manus*sweep*dihedral, data = dat_sp)
mod_sp_omega_n = lm(omega_n ~ elbow*manus*sweep*dihedral, data = dat_sp)

mod_ph_zeta    = lm(zeta ~ elbow*manus*sweep*dihedral, data = dat_ph)
mod_ph_omega_n = lm(omega_n ~ elbow*manus*sweep*dihedral, data = dat_ph)

mod_sp_speed_omega_n = lm(omega_n ~ U0 + theta_rad, data = dat_sp)
mod_sp_speed_zeta    = lm(zeta ~ U0 + theta_rad, data = dat_sp)

mod_ph_speed_omega_n = lm(omega_n ~ U0 + theta_rad, data = dat_ph)
mod_ph_speed_zeta    = lm(zeta ~ U0 + theta_rad, data = dat_ph)


## ------------------ Time response to an initial alpha -------------------

dat_time_1 <- read.csv('./outputdata/2021_11_04_elbow90_manus120_sw-10_di20_dalp.csv', header = FALSE)
colnames(dat_time_1) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_2 <- read.csv('./outputdata/2021_11_04_elbow120_manus120_sw-10_di20_dalp.csv', header = FALSE)
colnames(dat_time_2) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_3 <- read.csv('./outputdata/2021_11_04_elbow120_manus150_sw-10_di20_dalp.csv', header = FALSE)
colnames(dat_time_3) <- c("t","del_u","del_alp","del_q","del_theta")

## ------------------ Time response to a ramped speed -------------------

dat_time_4 <- read.csv('./outputdata/2021_11_04_elbow90_manus120_sw-10_di20_uramp.csv', header = FALSE)
colnames(dat_time_4) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_5 <- read.csv('./outputdata/2021_11_04_elbow120_manus120_sw-10_di20_uramp.csv', header = FALSE)
colnames(dat_time_5) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_6 <- read.csv('./outputdata/2021_11_04_elbow120_manus157_sw-10_di20_uramp.csv', header = FALSE)
colnames(dat_time_6) <- c("t","del_u","del_alp","del_q","del_theta")

