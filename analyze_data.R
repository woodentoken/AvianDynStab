library(alphahull) # need for the convex hull
library(ptinpoly)  # need for determining points outside the convex hull
library(pracma)
source("support_functions.R")

# ------------ Read in data ---------------
# Note that this does not include the wing shapes from the inertial gull
load("/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2022_02_07_function_data.RData")

# Note that the experimental data within dat_exp and the previous study numerical data within dat_num do not have a tail - cannot compare directly to predictions as drag is too low
##------------------- 
# allows a comparison against experimental data
dat_pred <- dat_exp
dat_pred$elbow_scale = dat_pred$elbow/1000
dat_pred$manus_scale = dat_pred$manus/1000
dat_pred$alpha_scale = dat_pred$alpha/10
dat_pred$sweep = 0
dat_pred$dihedral = 0
# predict how our models would fit experimental data
dat_pred$CL_adj_exp = dat_pred$L_comp
dat_pred$CL_adj = predict(mod_CL,dat_pred)
dat_pred$Cm_CG_pred = predict(mod_Cm,dat_pred)
dat_pred$CD_adj_pred = predict(mod_CD,dat_pred)

dat_pred$xcg = predict(mod_xcg_full,dat_pred) # the origin must be at the shoulder joint!!
dat_pred$zcg = predict(mod_zcg_full,dat_pred) # the origin must be at the shoulder joint!!
dat_pred$M_CG = dat_pred$m + ((dat_pred$L*cosd(dat_pred$alpha)+dat_pred$D*sind(dat_pred$alpha))*(-dat_pred$xcg) + 
                              (dat_pred$L*sind(dat_pred$alpha)-dat_pred$D*cosd(dat_pred$alpha))*(-dat_pred$zcg))
#needs to be multiplied by half the wing area as the experimental wings were only half and they were scaled 80%
dat_pred$Cm_CG_exp = dat_pred$M_CG/(0.5*dat_pred$rho*dat_pred$U^2*(0.5*0.8^2*max(dat_num$S_max[which(dat_num$WingID == "17_0285")]))*(0.8*max(dat_num$root_c_max[which(dat_num$WingID == "17_0285")])))

##------------------- 
# allow comparison against previous numerical data
dat_num$sweep = 0
dat_num$dihedral = 0
dat_num$xcg = predict(mod_xcg_full,dat_num) # the origin must be at the shoulder joint!!
dat_num$zcg = predict(mod_zcg_full,dat_num) # the origin must be at the shoulder joint!!

dat_num$L_comp = dat_num$CL_adj # make sure that this is the correct lit coefficient to use to predict the drag in the next step
dat_num$CD_adj_exp = predict(mod_CD,dat_num)
dat_num$D_adj_exp  = (0.5*1.225*10^2*dat_num$S_max)*dat_num$CD_adj_exp
dat_num$M_CG = dat_num$Mm + ((dat_num$FL*cosd(dat_num$alpha)+dat_num$D_adj_exp*sind(dat_num$alpha))*(-dat_num$xcg) + 
                               (dat_num$FL*sind(dat_num$alpha)-dat_num$D_adj_exp*cosd(dat_num$alpha))*(-dat_num$zcg))
dat_num$CL_adj_num = dat_num$L_comp
dat_num$CL_adj = predict(mod_CL,dat_num)
dat_num$Cm_CG_pred = predict(mod_Cm,dat_num)

# include supplemental graph to compare the experimentally predicted drag to the MachUpX drag
# non-dimensionalize
dat_num$Cm_CG = dat_num$M_CG/(0.5*1.225*10^2*dat_num$S_max*dat_num$root_c_max)
dat_wtnum <- subset(dat_num, FrameID == "F1380" | FrameID == "F2195" | FrameID == "F3891" | FrameID == "F4352" | FrameID == "F4546" | FrameID == "F4647" | FrameID == "F4849" | FrameID == "F4911" | FrameID == "F6003")
dat_wtnum$xcg <- predict(mod_xcg_full,dat_wtnum)
dat_wtnum$zcg <- predict(mod_zcg_full,dat_wtnum)

##------------------- Read in Python Outputs ------------------
dat_out <- read.csv('2022_02_07_LongDynStability_Rigid.csv')
S_max = 0.267782   # wings and body reference area from the gull used in inertial study (m^2) - from determine_functions.R
c_max = 0.2861011  # maximum wing root chord from the gull used in inertial study (m) - from determine_functions.R
m     = 1.0154     # mass (kg) - from determine_functions.R
W     = m * 9.81   # weight (N)
len   = unique(dat_inertial$full_length[which(dat_inertial$species == "lar_gla")])
# ------------ Clean data ---------------

# remove the unattainable configurations
tmp     <- cut_trueshape(dat_out,unique(dat_all[4:5]),4,5) # cut elbow and wrist to the true shape of the data
dat_cut <- tmp$dat_cut
# restrict to angles of attack with valid aerodynamic data
dat_cut <- subset(dat_cut, alpha < 5)
# Save data for stitched model - only want one row of each state so save first
dat_stitch <- subset(dat_cut, eignum == 1)
dat_stitch$xcg <- predict(mod_xcg_full,dat_stitch)
dat_stitch$zcg <- predict(mod_zcg_full,dat_stitch)

filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_trimstates.csv",sep="")
write.csv(dat_stitch,filename)
# Limit so that the wind direction isn't more than 45deg relative to the horizontal
dat_cut <- subset(dat_cut, gamma0 > -45*pi/180)
# limit the data to only discuss configurations at a single dihedral angle
dat_cut <- subset(dat_cut, dihedral == 20)
# limit the data to only discuss configurations for forward swept wings
dat_cut <- subset(dat_cut, sweep <= 0)

dat_cut$halft      = 0.69/abs(dat_cut$eig_real) # captures the rate of decay
dat_cut$period     = 1/(dat_cut$omega_n/2*pi) # changes rad/s to Hz then calculates period
dat_cut$SM         = -(dat_cut$Cm_alp/dat_cut$CL_alp)
dat_cut$NP_xcg     = -(dat_cut$Cm_alp/dat_cut$CL_alp)*c_max # this is the center of gravity minus the neutral point (if x-axis points forwards)
dat_cut$alpha_rad  = dat_cut$alpha*pi/180
dat_cut$theta_rad  = dat_cut$alpha_rad + dat_cut$gamma0 # angle from the horizon to the bird reference line
dat_cut$gamma0_deg = dat_cut$gamma0*180/pi
dat_cut$theta_0    = dat_cut$theta_rad*180/pi # angle from the horizon to the bird reference line
dat_cut$cmcl       = dat_cut$Cm_alp/dat_cut$CL_alp
dat_cut$n_a        = 0.5*1.225*dat_cut$U0^2*S_max*dat_cut$CL_alp/W # load factor (L/W) per angle of attack
dat_cut$w_sp_na    = dat_cut$omega_n^2/dat_cut$n_a
# load outputs from the fit function first
dat_cut$xcg <- -predict(mod_xcg_full,dat_cut) # must be negative to give the distance from shoulder to xcg as a positive value
dat_cut$xac <- dat_cut$xcg + dat_cut$NP_xcg 
# this works because the input CG is effectively negative (x-axis pointing back see above line) 
# by adding this value we effectively get the negative x_ac which means that the x-axis points backwards 
# This gives a nicer value for plotting

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

# ONLY KEEPS THE STATICALLY STABLE CONFIGURATIONS
dat_sp <- subset(dat_cut, dihedral == 20 & sweep < 1 & eignum == 1 & cmcl < 0) # eigen num == 2 just flips phase
dat_ph <- subset(dat_cut, dihedral == 20 & sweep < 1 & eignum == 3 & cmcl < 0) # eigen num == 4 just flips phase


## --------- Investigate the effect of morphing on stability --------

# Quantify the effect of joint angle on the increasing speed 
# subest to eignum== 1 to avoid having repeat info
mod_trim_speed    = lm(U0 ~ elbow*manus*sweep + gamma0_deg, data = subset(dat_cut, eignum == 1))
mod_trim_angle    = lm(gamma0_deg ~ elbow*manus*sweep + U0, data = subset(dat_cut, eignum == 1))

# limit our discussion to configurations that had many positions converge
mod_sp_zeta    = lm(zeta ~ elbow*manus*sweep, data = dat_sp)
mod_sp_omega_n = lm(omega_n ~ elbow*manus*sweep, data = dat_sp)

# Flying qualities
mod_flying_freq    = lm(w_sp_na ~ elbow + manus + sweep, data = dat_sp)
mod_flying_zeta    = lm(zeta ~ elbow+manus+sweep, data = dat_sp)

#Fosters scaling
scale_fos = sqrt(80)

#Capello's scaling - Using a Cessna 172 from http://jsbsim.sourceforge.net/MassProps.html
v_c = 140*0.44704
c_c = 4.9*0.3048
b_c = 35.8*0.3048
Iy_c = 1346*14.5939*0.3048^2 # converting from slug/ft^2
scale_cap = (10/v_c)*(c_max/c_c)*sqrt((max(dat_inertial$span)/b_c)*(Iy_c/max(dat_inertial$full_Iyy)))

## ------------------ Time response to an initial alpha -------------------
dat_time_a_1 <- read.csv('./outputdata/2022_02_07_elbow130_manus106_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_a_1) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_a_2 <- read.csv('./outputdata/2022_02_07_elbow130_manus116_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_a_2) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_a_3 <- read.csv('./outputdata/2022_02_07_elbow130_manus126_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_a_3) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_a_4 <- read.csv('./outputdata/2022_02_07_elbow130_manus136_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_a_4) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_a_5 <- read.csv('./outputdata/2022_02_07_elbow130_manus146_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_a_5) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_a_6 <- read.csv('./outputdata/2022_02_07_elbow130_manus156_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_a_6) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_a_7 <- read.csv('./outputdata/2022_02_07_elbow130_manus166_sw-15_di20_dalp.csv', header = FALSE)
colnames(dat_time_a_7) <- c("t","del_u","del_alp","del_q","del_theta")

## ------------------ Time response to a ramped speed -------------------

dat_time_r_1 <- read.csv('./outputdata/2022_02_07_elbow130_manus106_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_r_1) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_r_2 <- read.csv('./outputdata/2022_02_07_elbow130_manus116_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_r_2) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_r_3 <- read.csv('./outputdata/2022_02_07_elbow130_manus126_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_r_3) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_r_4 <- read.csv('./outputdata/2022_02_07_elbow130_manus136_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_r_4) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_r_5 <- read.csv('./outputdata/2022_02_07_elbow130_manus146_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_r_5) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_r_6 <- read.csv('./outputdata/2022_02_07_elbow130_manus156_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_r_6) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_r_7 <- read.csv('./outputdata/2022_02_07_elbow130_manus166_sw-15_di20_uramp.csv', header = FALSE)
colnames(dat_time_r_7) <- c("t","del_u","del_alp","del_q","del_theta")


