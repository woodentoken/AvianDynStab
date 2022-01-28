# number of tested configs
40*37*3*7

# number of trimmed configurations
nrow(subset(dat_cut, eignum == 1))

# effect of joints on trim speed
summary(mod_trim_speed)
summary(mod_trim_angle)

min(dat_cut$U0)
max(dat_cut$U0)
length(which(subset(dat_cut, eignum == 1)$U0 < 20))

max(dat_cut$gamma0)*180/pi
# static stability
nrow(subset(dat_cut, eignum == 1 & cmcl < 0))
nrow(subset(dat_cut, eignum == 1 & cmcl > 0))

plot(subset(dat_sp, sweep == -15 & elbow > 90)$manus, subset(dat_sp, sweep == -15 & elbow > 90)$omega_n)
plot(subset(dat_sp, sweep == -15 & elbow > 90)$manus, subset(dat_sp, sweep == -15 & elbow > 90)$zeta)
test <- lm(omega_n~manus,data = subset(dat_sp, sweep == -15 & elbow == 90)) # this is the elbow angle where the effect of the wrist starts to shift
test <- lm(omega_n~manus,data = subset(dat_sp, sweep == -15 & elbow == 92))

min(dat_sp$period)
max(dat_sp$period)

min(dat_sp$omega_n) # results are in rad/s
max(dat_sp$omega_n) # results are in rad/s

min(dat_ph$period)
max(dat_ph$period)

min(dat_ph$omega_n) # results are in rad/s
max(dat_ph$omega_n) # results are in rad/s

dat_manus_effect_sp           <- data.frame(matrix(nrow = 105, ncol = 4))
colnames(dat_manus_effect_sp) <- c("sweep","elbow","d_omega_n","d_zeta")
dat_manus_effect_ph           <- data.frame(matrix(nrow = 105, ncol = 4))
colnames(dat_manus_effect_ph) <- c("sweep","elbow","d_omega_n","d_zeta")
count = 1
for(i in 1:5){
  
  curr_sweep = unique(dat_sp$sweep)[i]
  
  for (j in 1:40){
    
    curr_elbow = unique(dat_sp$elbow)[j]
    
    if(nrow(subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow)) < 4){next}
    
    test <- lm(omega_n~manus,data = subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow))
    out_o <- summary(test)
    test <- lm(zeta~manus,data = subset(dat_sp, sweep == curr_sweep & elbow == curr_elbow))
    out_z <- summary(test)
    
    dat_manus_effect_sp$sweep[count] = curr_sweep
    dat_manus_effect_sp$elbow[count] = curr_elbow
    dat_manus_effect_sp$d_omega_n[count] = out_o$coefficients["manus","Estimate"]
    dat_manus_effect_sp$d_zeta[count] = out_z$coefficients["manus","Estimate"]
    
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

length(which(dat_manus_effect_sp$d_omega_n > 0)) # wrist tended to decreases natural frequency
plot(dat_manus_effect_sp$elbow,dat_manus_effect_sp$d_omega_n)
length(which(dat_manus_effect_sp$d_zeta < 0)) # wrist tended to increase damping ratio
plot(dat_manus_effect_sp$elbow,dat_manus_effect_sp$d_zeta)

plot(dat_manus_effect_ph$elbow,dat_manus_effect_ph$d_zeta)
plot(dat_manus_effect_ph$elbow,dat_manus_effect_ph$d_omega_n)

dat_exp$LD <- dat_exp$L/dat_exp$D

# Flying qualities

# Level 2 flying qualities on the damping ratio
dat_sp$sweep[which(dat_sp$zeta < 0.3)]
length(dat_sp$sweep[which(dat_sp$zeta < 0.3)])
length(dat_sp$sweep[which(dat_sp$zeta > 0.3)])
#Fosters scaling
length(which(dat_sp$w_sp_na < scale_fos*3.6))/nrow(dat_sp) # Level 1 upper bound
length(which(dat_sp$w_sp_na < scale_fos*10))/nrow(dat_sp) # Level 2 upper bound

length(which(dat_sp$w_sp_na > scale_fos*0.085))/nrow(dat_sp) # Level 1 lower bound
length(which(dat_sp$w_sp_na > scale_fos*0.038))/nrow(dat_sp) # Level 2 lower bound

#Capello's scaling - Using a Cessna 172 from http://jsbsim.sourceforge.net/MassProps.html
length(which(dat_sp$w_sp_na < scale_cap*3.6))/nrow(dat_sp) # Level 1 upper bound
length(which(dat_sp$w_sp_na < scale_cap*10))/nrow(dat_sp) # Level 2 upper bound

length(which(dat_sp$w_sp_na > scale_cap*0.085))/nrow(dat_sp) # Level 1 lower bound
length(which(dat_sp$w_sp_na > scale_cap*0.038))/nrow(dat_sp) # Level 2 lower bound

length(which(dat_sp$w_sp_na < 3.6))
length(which(dat_sp$w_sp_na < 10))

mean(dat_ph$zeta)
mean(dat_ph$omega_n)

# Sanity check that the lift and drag parameters are appropriate
plot(dat_exp$L_comp,dat_exp$CD_true)
points(dat_cut$CL,dat_cut$CD, col = "red")

# time to halve
subset(dat_ph, elbow == 130 & sweep == -15 & manus %in% c(106,116,126,136,146))$halft
subset(dat_sp, elbow == 130 & sweep == -15 & manus %in% c(106,116,126,136,146))$halft

#time to double
dat_out$halft  = 0.69/abs(dat_out$eig_real) # captures the rate of decay

subset(dat_out, eignum == 1 & elbow == 130 & sweep == -15 & dihedral == 20 & manus %in% c(156,166))$halft
subset(dat_out, eignum == 3 & elbow == 130 & sweep == -15 & dihedral == 20 & manus %in% c(156,166))$halft
#explain 156deg
unique(subset(dat_out, elbow == 130 & sweep == -15 & manus ==156 & dihedral == 20)$gamma0)*180/(pi)
