library(alphahull) # need for the convex hull
library(ptinpoly)  # need for determining points outside the convex hull
library(pracma)
source("support_functions.R")

# ------------ Read in data ---------------
# Note that this does not include the wing shapes from the inertial gull
load("/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2021_11_21_function_data.RData")

dat_out <- read.csv('2021_11_22_LongDynStability_Rigid.csv')
S_max = 0.244657   # wings and body reference area from the gull used in inertial study (m^2)
c_max = 0.2861011  # wing root chord from the gull used in inertial study (m)
m     = 1.0154     # mass (kg)
W     = m * 9.81   # weight (N)
# ------------ Clean data ---------------

# remove the unattainable configurations
tmp     <- cut_trueshape(dat_out,unique(dat_all[4:5]),4,5) # cut elbow and wrist to the true shape of the data
dat_cut <- tmp$dat_cut
# restrict to angles of attack with trustworthy aerodynamic data
# ADD NOTE THAT THIS DOESN'T MEAN BIRDS CAN TRIM AT SLOWER SPEEDS JUST THAT MY AERO DATA IN THAT RANGE ISN'T AS TRUSTWORTHY
dat_cut <- subset(dat_cut, alpha < 5)
# Limit so that the wind direction isn't more than 45deg relative to the horizontal
dat_cut <- subset(dat_cut, gamma0 > -45*pi/180)
# limit the data to only discuss configurations at a single dihedral angle
dat_cut <- subset(dat_cut, dihedral == 20)

dat_cut$halft  = 0.69/abs(dat_cut$eig_real) # captures the rate of decay
dat_cut$period = 2*pi/dat_cut$omega_n
dat_cut$SM = -(dat_cut$Cm_alp/dat_cut$CL_alp)
dat_cut$NP = -(dat_cut$Cm_alp/dat_cut$CL_alp)*c_max
dat_cut$alpha_rad = dat_cut$alpha*pi/180
dat_cut$theta_rad = dat_cut$alpha_rad + dat_cut$gamma0 # angle from the horizon to the bird reference line
dat_cut$theta_0 = dat_cut$theta_rad*180/pi # angle from the horizon to the bird reference line
dat_cut$cmcl = dat_cut$Cm_alp/dat_cut$CL_alp
dat_cut$n_a = 0.5*1.225*dat_cut$U0^2*dat_cut$CL_alp/W# load factor per angle of attack
dat_cut$w_sp_na = dat_cut$omega_n^2/dat_cut$n_a
# load outputs from the fit function first
dat_cut$xcg <- -predict(mod_xcg_full,dat_cut) # must be negative to give the distance from shoulder to xcg
dat_cut$xac <- dat_cut$xcg + dat_cut$NP

# Set the pitch rate as the comparable phase 
dat_cut$phase1 = dat_cut$phase1 - dat_cut$phase3
dat_cut$phase2 = dat_cut$phase2 - dat_cut$phase3
dat_cut$phase4 = dat_cut$phase4 - dat_cut$phase3
dat_cut$phase3 = dat_cut$phase3 - dat_cut$phase3

# adjust negative numbers to be in a positive frame
dat_cut$phase1[which(dat_cut$phase1 < 0)] = (2*pi) + dat_cut$phase1[which(dat_cut$phase1 < 0)]
dat_cut$phase2[which(dat_cut$phase2 < 0)] = (2*pi) + dat_cut$phase2[which(dat_cut$phase2 < 0)]
dat_cut$phase4[which(dat_cut$phase4 < 0)] = (2*pi) + dat_cut$phase4[which(dat_cut$phase4 < 0)]
dat_cut$phase3[which(dat_cut$phase3 < 0)] = (2*pi) + dat_cut$phase3[which(dat_cut$phase3 < 0)]

#  Set the pitch rate as the comparable magnitude
dat_cut$mag1 = dat_cut$mag1/dat_cut$mag3
dat_cut$mag2 = dat_cut$mag2/dat_cut$mag3
dat_cut$mag4 = dat_cut$mag4/dat_cut$mag3
dat_cut$mag3 = dat_cut$mag3/dat_cut$mag3

dat_sp <- subset(dat_cut, dihedral == 20 & sweep < 1 & eignum == 1 & omega_n > 0) # eigen num == 2 just flips phase
dat_ph <- subset(dat_cut, dihedral == 20 & sweep < 1 & eignum == 3 & omega_n > 0) # eigen num == 4 just flips phase
# need to make sure that the configurations that are unstable in sp mode are removed here as well
dat_ph <- merge(dat_sp[,c("alpha","elbow","manus")], subset(dat_cut, eignum == 3 & omega_n > 0))

## --------- Investigate the effect of morphing on stability --------
# limit our discussion to configurations that had many positions converge
mod_sp_zeta    = lm(zeta ~ elbow*manus*sweep + I(elbow^2) + I(manus^2), data = dat_sp)
mod_sp_omega_n = lm(omega_n ~ elbow*manus + I(elbow^2) + I(manus^2) + sweep, data = dat_sp)

mod_ph_zeta    = lm(zeta ~ elbow*manus + I(elbow^2) + I(manus^2) + sweep, data = dat_ph)
mod_ph_omega_n = lm(omega_n ~ elbow*manus + I(elbow^2) + I(manus^2) + sweep, data = dat_ph)
                    
mod_speed_morph = lm(U0 ~ elbow*manus*sweep, data = dat_sp)
mod_angle_morph = lm(theta_rad ~ elbow*manus*sweep , data = dat_sp)

test <- lm(eig_imag ~ Cm_q + Cm_alp + CL_alp, data = dat_sp) # this is just a theoretical expectation
test <- lm(zeta ~ Cm_q + Cm_alp + CL_alp, data = dat_sp) # this is just a theoretical expectation
## ------------------ Time response to an initial alpha -------------------
dat_time_1 <- read.csv('./outputdata/2021_11_22_elbow126_manus106_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_1) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_2 <- read.csv('./outputdata/2021_11_22_elbow126_manus116_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_2) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_3 <- read.csv('./outputdata/2021_11_22_elbow126_manus126_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_3) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_4 <- read.csv('./outputdata/2021_11_22_elbow126_manus136_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_4) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_5 <- read.csv('./outputdata/2021_11_22_elbow126_manus146_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_5) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_6 <- read.csv('./outputdata/2021_11_22_elbow126_manus156_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_6) <- c("t","del_u","del_alp","del_q","del_theta")

## ------------------ Time response to a ramped speed -------------------

dat_time_7 <- read.csv('./outputdata/2021_11_22_elbow126_manus106_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_7) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_8 <- read.csv('./outputdata/2021_11_22_elbow126_manus116_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_8) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_9 <- read.csv('./outputdata/2021_11_22_elbow126_manus126_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_9) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_10 <- read.csv('./outputdata/2021_11_22_elbow126_manus136_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_10) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_11 <- read.csv('./outputdata/2021_11_22_elbow126_manus146_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_11) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_12 <- read.csv('./outputdata/2021_11_22_elbow126_manus156_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_12) <- c("t","del_u","del_alp","del_q","del_theta")
