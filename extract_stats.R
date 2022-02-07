# number of tested configs
40*37*3*7

## ---------- Methods - Trim states ---------
# number of trimmed configurations
nrow(subset(dat_cut, eignum == 1))

# trim glide speed range
min(dat_cut$U0)
max(dat_cut$U0)
# shallowest glide angle
max(dat_cut$gamma0)*180/pi

# speeds within measured range 
length(which(subset(dat_cut, eignum == 1)$U0 <= 19.9))

# iterate through each configuration to see how many configurations tend to have what effect 
# quantify despite interative effects
## -----------------------------------------------------------------------------
## ----------- Effect of the wrist w/ constant sweep and elbow -------------
dat_manus_effect_sp           <- data.frame(matrix(nrow = 106, ncol = 4))
colnames(dat_manus_effect_sp) <- c("sweep","elbow","d_omega_n","d_zeta")
dat_manus_effect_ph           <- data.frame(matrix(nrow = 106, ncol = 4))
colnames(dat_manus_effect_ph) <- c("sweep","elbow","d_omega_n","d_zeta")
count = 1
for(i in 1:5){
  curr_sweep = unique(dat_sp$sweep)[i]
  for (j in 1:40){
    
    curr_elbow = unique(dat_sp$elbow)[j]
    
    if(nrow(subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow)) < 3){next}
    
    test <- lm(omega_n~manus,data = subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow))
    out_o <- summary(test)
    test <- lm(zeta~manus,data = subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow))
    out_z <- summary(test)
    
    test <- lm(U0~manus,data = subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow))
    out_u <- summary(test)
    test <- lm(gamma0~manus,data = subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow))
    out_g <- summary(test)
    
    dat_manus_effect_sp$sweep[count] = curr_sweep
    dat_manus_effect_sp$elbow[count] = curr_elbow
    dat_manus_effect_sp$d_omega_n[count] = out_o$coefficients["manus","Estimate"]
    dat_manus_effect_sp$d_zeta[count] = out_z$coefficients["manus","Estimate"]
    dat_manus_effect_sp$u[count] = out_u$coefficients["manus","Estimate"]
    dat_manus_effect_sp$g[count] = out_g$coefficients["manus","Estimate"]
    
    test <- lm(omega_n~manus,data = subset(dat_ph, sweep == curr_sweep & elbow == curr_elbow))
    out_o <- summary(test)
    test <- lm(zeta~manus,data = subset(dat_ph, sweep == curr_sweep & elbow == curr_elbow))
    out_z <- summary(test)
    
    dat_manus_effect_ph$sweep[count] = curr_sweep
    dat_manus_effect_ph$elbow[count] = curr_elbow
    dat_manus_effect_ph$d_omega_n[count] = out_o$coefficients["manus","Estimate"]
    dat_manus_effect_ph$d_zeta[count] = out_z$coefficients["manus","Estimate"]
    
    count = count + 1
  }
}

## ----------- Effect of the elbow w/ constant sweep and wrist -------------
dat_elbow_effect          <- data.frame(matrix(nrow = 92, ncol = 4))
colnames(dat_elbow_effect) <- c("sweep","manus","u","g")
count = 1
for(i in 1:5){
  curr_sweep = unique(dat_sp$sweep)[i]
  for (j in 1:37){
    
    curr_manus = unique(dat_sp$manus)[j]
    
    if(nrow(subset(dat_sp, sweep == curr_sweep & manus == curr_manus)) < 3){next}
    
    test <- lm(U0~elbow,data = subset(dat_sp, sweep == curr_sweep & manus == curr_manus))
    out_u <- summary(test)
    test <- lm(gamma0~elbow,data = subset(dat_sp, sweep == curr_sweep & manus == curr_manus))
    out_g <- summary(test)
    
    dat_elbow_effect$sweep[count] = curr_sweep
    dat_elbow_effect$manus[count] = curr_manus
    dat_elbow_effect$u[count] = out_u$coefficients["elbow","Estimate"]
    dat_elbow_effect$g[count] = out_g$coefficients["elbow","Estimate"]
    
    count = count + 1
  }
}

## ----------- Effect of the sweep w/ constant wrist and elbow -------------

dat_sweep_effect          <- data.frame(matrix(nrow = 205, ncol = 4))
colnames(dat_sweep_effect) <- c("elbow","manus","u","g")
count = 1
for(i in 1:40){
  curr_elbow = unique(dat_sp$elbow)[i]
  for (j in 1:37){
    
    curr_manus = unique(dat_sp$manus)[j]
    
    if(nrow(subset(dat_sp, elbow == curr_elbow & manus == curr_manus)) < 3){next}
    
    test <- lm(U0~sweep,data = subset(dat_sp, elbow == curr_elbow & manus == curr_manus))
    out_u <- summary(test)
    test <- lm(gamma0~sweep,data = subset(dat_sp, elbow == curr_elbow & manus == curr_manus))
    out_g <- summary(test)
    
    dat_sweep_effect$elbow[count] = curr_elbow
    dat_sweep_effect$manus[count] = curr_manus
    dat_sweep_effect$u[count] = out_u$coefficients["sweep","Estimate"]
    dat_sweep_effect$g[count] = out_g$coefficients["sweep","Estimate"]
    
    count = count + 1
  }
}
## -----------------------------------------------------------------------------
length(which(dat_manus_effect_sp$u < 0)) / nrow(dat_manus_effect_sp)
length(which(dat_elbow_effect$u < 0)) / nrow(dat_elbow_effect)
length(which(dat_sweep_effect$u < 0)) / nrow(dat_sweep_effect)

## ------------------- Static stability about the trim state ----------------------------

# static stability
nrow(subset(dat_cut, eignum == 1 & cmcl < 0))
nrow(subset(dat_cut, eignum == 1 & cmcl > 0))

## ------------------- Dynamic stability about the trim state ----------------------------
# only negative real parts for statically stable configs
nrow(subset(dat_sp, eignum == 1 & eig_real < 0))
# amount of unstable configurations with no imaginary component
nrow(subset(dat_cut, eignum == 1 & cmcl > 0 & eig_imag == 0))

## --- Short period ---
# short period frequency range
min(dat_sp$omega_n) # results are in rad/s
max(dat_sp$omega_n) # results are in rad/s

# significant interactive effects
summary(mod_sp_zeta)
summary(mod_sp_omega_n)

# general trends due to wrist
length(which(dat_manus_effect_sp$d_omega_n > 0)) # wrist tended to decreases natural frequency
plot(dat_manus_effect_sp$elbow,dat_manus_effect_sp$d_omega_n)
length(which(dat_manus_effect_sp$d_zeta < 0)) # wrist tended to increase damping ratio
plot(dat_manus_effect_sp$elbow,dat_manus_effect_sp$d_zeta)

## --- Phugoid ---
min(dat_ph$omega_n) # results are in rad/s
max(dat_ph$omega_n) # results are in rad/s

# general trends due to wrist
plot(dat_manus_effect_ph$elbow,dat_manus_effect_ph$d_zeta)
plot(dat_manus_effect_ph$elbow,dat_manus_effect_ph$d_omega_n)

## ------------------- Flying qualities ----------------------------

# Level 2 flying qualities on the damping ratio
dat_sp$sweep[which(dat_sp$zeta < 0.3)]
length(dat_sp$sweep[which(dat_sp$zeta < 0.3)])
length(dat_sp$sweep[which(dat_sp$zeta > 0.3)])

# MIL-F-8785 
length(which(dat_sp$w_sp_na < 3.6)) # Level 1 upper bound
length(which(dat_sp$w_sp_na < 10)) # Level 2 upper bound
# joints tended to - note this model doesn't capture interactive effects
summary(mod_flying_freq)

#Fosters scaling
length(which(dat_sp$w_sp_na < scale_fos*3.6)) # Level 1 upper bound
length(which(dat_sp$w_sp_na < scale_fos*10)) # Level 2 upper bound

length(which(dat_sp$w_sp_na > scale_fos*0.085)) # Level 1 lower bound
length(which(dat_sp$w_sp_na > scale_fos*0.038)) # Level 2 lower bound

#Capello's scaling - Using a Cessna 172 from http://jsbsim.sourceforge.net/MassProps.html
length(which(dat_sp$w_sp_na < scale_cap*3.6)) # Level 1 upper bound
length(which(dat_sp$w_sp_na < scale_cap*10)) # Level 2 upper bound

length(which(dat_sp$w_sp_na > scale_cap*0.085)) # Level 1 lower bound
length(which(dat_sp$w_sp_na > scale_cap*0.038)) # Level 2 lower bound

mean(dat_ph$zeta)

## ------------------- Simplified gust response ----------------------------

#explain 156deg
unique(subset(dat_out, elbow == 130 & sweep == -15 & manus ==156 & dihedral == 20)$gamma0)*180/(pi)

# time to halve
subset(dat_ph, elbow == 130 & sweep == -15 & manus %in% c(106,116,126,136,146))$halft
subset(dat_sp, elbow == 130 & sweep == -15 & manus %in% c(106,116,126,136,146))$halft

#time to double
dat_out$halft  = 0.69/abs(dat_out$eig_real) # captures the rate of decay
subset(dat_out, eignum == 1 & elbow == 130 & sweep == -15 & dihedral == 20 & manus %in% c(156,166))$halft
subset(dat_out, eignum == 3 & elbow == 130 & sweep == -15 & dihedral == 20 & manus %in% c(156,166))$halft


# Sanity check that the lift and drag parameters are appropriate
plot(dat_exp$L_comp,dat_exp$CD_true)
points(dat_cut$CL,dat_cut$CD, col = "red")

## -----------------------------------------------------------------------------------------------
## -------------------------------- SUPPLEMENTAL METHODS -----------------------------------------
## -----------------------------------------------------------------------------------------------

nrow(subset(dat_aero_all,alpha < 5))
summary(mod_CL)
summary(mod_Cm)
summary(mod_CL_q_ind)
summary(mod_Cm_q_ind)
# Used for investigating the potential of improved model fit through regression training
# library(caret)
# #set up methods for cross-validation
# ctrl<-trainControl(method = "repeatedcv",number = 10,repeats = 10)
# rf.full<-train(CL_adj ~ elbow_scale + manus_scale + alpha_scale + 
#                  I(alpha_scale^2) + I(alpha_scale^3) +
#                  I(elbow_scale^2) + I(manus_scale^2) + I(manus_scale^3) +
#                  sweep + dihedral + 
#                  sweep:alpha_scale + dihedral:alpha_scale,data=subset(dat_aero_all,alpha < 5),method='lm', trControl = ctrl)
# dat_pred$test <- predict(rf.full,dat_pred)
# mean(abs(subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$test-subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$CL_adj_exp))

# Prediction error - Experimental results
mean(abs(subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$Cm_CG_pred-subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$Cm_CG_exp))
mean(abs(subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$CL_adj-subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$CL_adj_exp))
mean(abs(subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$CD_adj_pred-subset(dat_pred, alpha < 5 & alpha >= -10 & U_des == "low")$CD_true))

# Prediction error - Numerical results
mean(abs(subset(dat_num, alpha < 5)$Cm_CG_pred-subset(dat_num, alpha < 5)$Cm_CG))
mean(abs(subset(dat_num, alpha < 5)$CL_adj-subset(dat_num, alpha < 5)$CL_adj_num))

# Prediction error - Inertial results
dat_iner_all$Iyy_pred <- predict(mod_inertia,dat_iner_all)
mean(abs(dat_iner_all$Iyy_pred - dat_iner_all$full_Iyy_adj))

# Effect of varying flight speed
test <- lm(L_comp ~ elbow + manus + alpha + I(alpha^2) + U, data = subset(dat_exp, alpha >= -10 & alpha < 5))
test <- lm(m_comp ~ elbow + manus + alpha + I(alpha^2) + U, data = subset(dat_exp, alpha >= -10 & alpha < 5))
test <- lm(CD_true ~ elbow + manus + alpha + I(alpha^2) + U, data = subset(dat_exp, alpha >= -10 & alpha < 5))
