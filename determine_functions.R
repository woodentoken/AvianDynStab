library(lme4)
library(pracma)

# saved after running analyse_llt.R and analyze_expresults.R - no changes
load("/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/aerodynamic_data.RData")
#saved after running process_outputdata.R
load("/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/inertial_data.RData")
# load the Q derivative specific runs
dat_q <- read.csv("/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/MachUpX_qderivatives/List_ConvergedWings_derivatives.csv", header = FALSE)
names(dat_q) <- c("species","WingID","TestID","FrameID","elbow","manus","alpha",
                  "U","build_err_max","date","S","ref_c","b_MX","MAC","b",
                  "sweep","dihedral","twist",'relax',
                  "CL","CD","Cm","FL","FD","Mm",'q')

# determine the root chord
dat_all <- read.csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2020_05_25_OrientedWings.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
dat_all <- subset(dat_all, species == "lar_gla" & sweep == 0 & dihedral == 0)
dat_all$root_c = sqrt((dat_all$Pt12X - dat_all$Pt11X)^2 + (dat_all$Pt12Z - dat_all$Pt11Z)^2)
dat_all$FrameID <- paste("F", dat_all$frameID, sep = "")

dat_num = merge(dat_num,dat_all[,c("species","FrameID","TestID","WingID","elbow","manus","root_c")], by = c("species","FrameID","TestID","WingID","elbow","manus"))

dat_num$root_c_max <- 0
dat_num$root_c_max[which(dat_num$WingID == "17_0285")] = max(dat_num$root_c[which(dat_num$WingID == "17_0285")])
dat_num$root_c_max[which(dat_num$WingID == "17_0243")] = max(dat_num$root_c[which(dat_num$WingID == "17_0243")])
dat_num$root_c_max[which(dat_num$WingID == "16_0048")] = max(dat_num$root_c[which(dat_num$WingID == "16_0048")])

coef_all = as.data.frame(matrix(0,nrow = 7, ncol = 24))
colnames(coef_all) <- c("y.model","intercept","elbow","manus","elbow2","manus2",
                        "elbow3","manus3","elbowmanus",
                        "alpha","alpha2","alpha3",
                        "elbowalpha","manusalpha","elbowmanusalpha",
                        "CL","CL2","CL3",
                        "elbowCL","manusCL","elbowmanusCL",
                        "elbowCL2","manusCL2")

## ---------------------------------------------
## ------------- Inertia Data ------------------
## ---------------------------------------------
#  !! CAUTION: INCLUDE NOTE IN PAPER THAT I EFFECTIVELY AM EXTRAPOLATING THE INERTIAL RESULS INTO HIGHER ELBOW AND WRIST ANGLES

# subset to the same range as used in the aerodynamic data
dat_inertial = subset(dat_final,species == "lar_gla" & elbow > 80 & manus > 100)
remove(dat_final)
# need this to be two wings + full body area
2*max(dat_inertial$S_max) + 0.0298 # body area was determined by dividing 0.4/0.41*(0.0305) that number is taken from 3 different wings test .dist outputs from MachUpX
max(dat_inertial$c_max) # root chord approximation - calculated in the same manner as above in processdata.R
max(dat_inertial$b_max) # of full wing
max(dat_inertial$full_m)
mod_inertia <- lm(full_Iyy ~ elbow*manus + I(elbow^2) + I(elbow^3) +
                    I(manus^2)+ I(manus^3), dat_inertial)

coef_all$y.model[1]    = "Iyy"
coef_all$intercept[1]  = coef(mod_inertia)["(Intercept)"]
coef_all$elbow[1]      = coef(mod_inertia)["elbow"]
coef_all$manus[1]      = coef(mod_inertia)["manus"]
coef_all$elbow2[1]     = coef(mod_inertia)["I(elbow^2)"]
coef_all$elbow3[1]     = coef(mod_inertia)["I(elbow^3)"]
coef_all$manus2[1]     = coef(mod_inertia)["I(manus^2)"]
coef_all$manus3[1]     = coef(mod_inertia)["I(manus^3)"]
coef_all$elbowmanus[1] = coef(mod_inertia)["elbow:manus"]

## need to fit model to the full_CGx to allow Cm to be adjusted appropriately.
mod_xcg <- lm(full_CGx ~ elbow*manus + I(elbow^2) + I(elbow^3) +
                I(manus^2)+ I(manus^3), dat_inertial)

mod_zcg <- lm(full_CGz ~ elbow*manus + I(elbow^2) + I(elbow^3) +
                I(manus^2)+ I(manus^3), dat_inertial)
## -------------------------------------------------
## ------------- Aerodynamic Data ------------------
## -------------------------------------------------

# scaling was elbow/1000 manus/1000 alpha/10
# If I want to input the true values into this same model I instead
# need to divide each coefficient by however many inputs are being changed

# ------ CD data ------
# need to fit to experimental data due to the reduced - following same adjustment as the lift in analyse_exp.R
dat_exp$CD_true     <- dat_exp$D_comp/(0.5*max(dat_num$S[which(dat_num$WingID == "17_0285")]))

# minor changes in the coefficients between speeds, will group them
# interactive terms were non-significant leave them out
mod_CD <- lm(CD_true ~ elbow + manus + L_comp + I(L_comp^2) + I(L_comp^3) + 
               elbow:I(L_comp^2) + manus:I(L_comp^2) + 
               elbow:L_comp + manus:L_comp, data = dat_exp)

coef_all$y.model[5]    = "CD"
coef_all$intercept[5]  = coef(mod_CD)["(Intercept)"]
coef_all$elbow[5]      = coef(mod_CD)["elbow"]
coef_all$manus[5]      = coef(mod_CD)["manus"]
coef_all$CL[5]         = coef(mod_CD)["L_comp"]
coef_all$CL2[5]        = coef(mod_CD)["I(L_comp^2)"]
coef_all$CL3[5]        = coef(mod_CD)["I(L_comp^3)"]
coef_all$elbowCL[5]      = coef(mod_CD)["elbow:L_comp"]
coef_all$manusCL[5]      = coef(mod_CD)["manus:L_comp"]
coef_all$elbowCL2[5]      = coef(mod_CD)["elbow:I(L_comp^2)"]
coef_all$manusCL2[5]      = coef(mod_CD)["manus:I(L_comp^2)"]

# ------ CL data ------
mod_CL <- lmer(CL_adj ~ elbow_scale*manus_scale*alpha_scale + I(alpha_scale^2) + I(alpha_scale^3) + 
                 I(elbow_scale^2) + 
                 I(manus_scale^2) + I(manus_scale^3) + (1|WingID), data = subset(dat_num,alpha < 5))

coef_all$y.model[3]    = "CL"
coef_all$intercept[3]  = summary(mod_CL)$coefficients["(Intercept)","Estimate"]
coef_all$elbow[3]      = summary(mod_CL)$coefficients["elbow_scale","Estimate"]/1000
coef_all$manus[3]      = summary(mod_CL)$coefficients["manus_scale","Estimate"]/1000
coef_all$elbow2[3]     = summary(mod_CL)$coefficients["I(elbow_scale^2)","Estimate"]/(1000^2)
coef_all$manus2[3]     = summary(mod_CL)$coefficients["I(manus_scale^2)","Estimate"]/(1000^2)
coef_all$manus3[3]     = summary(mod_CL)$coefficients["I(manus_scale^3)","Estimate"]/(1000^3)
coef_all$elbowmanus[3] = summary(mod_CL)$coefficients["elbow_scale:manus_scale","Estimate"]/(1000^2)
coef_all$alpha[3]      = summary(mod_CL)$coefficients["alpha_scale","Estimate"]/10
coef_all$alpha2[3]     = summary(mod_CL)$coefficients["I(alpha_scale^2)","Estimate"]/(10^2)
coef_all$alpha3[3]     = summary(mod_CL)$coefficients["I(alpha_scale^3)","Estimate"]/(10^2)
coef_all$elbowalpha[3] = summary(mod_CL)$coefficients["elbow_scale:alpha_scale","Estimate"]/(1000*10)
coef_all$manusalpha[3] = summary(mod_CL)$coefficients["manus_scale:alpha_scale","Estimate"]/(1000*10)
coef_all$elbowmanusalpha[3] = summary(mod_CL)$coefficients["elbow_scale:manus_scale:alpha_scale","Estimate"]/((1000^2)*10)

# ------------ Adjust the pitching moment data ---------------
# need to adjust the moment to be calculated about the current xCG
dat_num$xcg = predict(mod_xcg,dat_num)
dat_num$zcg = predict(mod_zcg,dat_num)

dat_num$L_comp = dat_num$CL_adj # make sure that this is the correct lit coefficient
dat_num$CD_adj_exp = predict(mod_CD,dat_num)

# adjust the pitching moment be about the true center of gravity
dat_num$Cm_CG = dat_num$Cm_adj + 
  (1/dat_num$c_max)*((dat_num$CL_adj*cosd(dat_num$alpha)+dat_num$CD_adj_exp*sind(dat_num$alpha))*(-dat_num$xcg) + 
                       (dat_num$CL_adj*sind(dat_num$alpha)-dat_num$CD_adj_exp*cosd(dat_num$alpha))*(-dat_num$zcg))
#adjust to be about the root chord
dat_num$Cm_CG = dat_num$Cm_CG*(dat_num$c_max/dat_num$root_c_max)

# ------ Cm data ------
mod_Cm <- lmer(Cm_CG ~ elbow_scale*manus_scale*CL_adj + I(CL_adj^2) + I(CL_adj^3) + 
                 I(elbow_scale^2) + I(elbow_scale^3) + 
                 I(manus_scale^2) + (1|WingID), data = subset(dat_num,alpha < 5))

coef_all$y.model[4]    = "Cm"
coef_all$intercept[4]  = summary(mod_Cm)$coefficients["(Intercept)","Estimate"]
coef_all$elbow[4]      = summary(mod_Cm)$coefficients["elbow_scale","Estimate"]/1000
coef_all$manus[4]      = summary(mod_Cm)$coefficients["manus_scale","Estimate"]/1000
coef_all$elbow2[4]     = summary(mod_Cm)$coefficients["I(elbow_scale^2)","Estimate"]/(1000^2)
coef_all$elbow3[4]     = summary(mod_Cm)$coefficients["I(elbow_scale^3)","Estimate"]/(1000^3)
coef_all$manus2[4]     = summary(mod_Cm)$coefficients["I(manus_scale^2)","Estimate"]/(1000^2)
coef_all$elbowmanus[4] = summary(mod_Cm)$coefficients["elbow_scale:manus_scale","Estimate"]/(1000^2)
coef_all$CL[4]         = summary(mod_Cm)$coefficients["CL_adj","Estimate"]
coef_all$CL2[4]        = summary(mod_Cm)$coefficients["I(CL_adj^2)","Estimate"]
coef_all$CL3[4]        = summary(mod_Cm)$coefficients["I(CL_adj^3)","Estimate"]
coef_all$elbowCL[4]    = summary(mod_Cm)$coefficients["elbow_scale:CL_adj","Estimate"]/(1000)
coef_all$manusCL[4]    = summary(mod_Cm)$coefficients["manus_scale:CL_adj","Estimate"]/(1000)
coef_all$elbowmanusCL[4] = summary(mod_Cm)$coefficients["elbow_scale:manus_scale:CL_adj","Estimate"]/(1000^2)

# ------ dCm/dCL data ------

dat_wingspec <- unique(dat_num[c("WingID","TestID","FrameID","elbow","manus","species","twist","sweep","dihedral","S_ref","c_max","elbow_scale","manus_scale")])
no_testedconfigs = nrow(dat_wingspec)
dat_stab_adj  <- data.frame(matrix(NA, nrow = no_testedconfigs, ncol = 8))
names(dat_stab_adj) <- c("species","WingID","TestID","FrameID","elbow","manus","cmcl","cm0")

# need to loop through all configurations to re-calculate the static margin
for (m in 1:no_testedconfigs){
  # subset data to be of one wing configuration at a time and subset to only fit angles under 5deg
  dat_curr <- subset(dat_num, 
                     species == dat_wingspec$species[m] & WingID == dat_wingspec$WingID[m] & 
                       TestID == dat_wingspec$TestID[m] & FrameID == dat_wingspec$FrameID[m] & alpha < 5)
  
  # save all wing specific information  
  dat_stab_adj$species[m] <- as.character(dat_wingspec$species[m])
  dat_stab_adj$WingID[m]  <- dat_wingspec$WingID[m]
  dat_stab_adj$TestID[m]  <- dat_wingspec$TestID[m]
  dat_stab_adj$FrameID[m] <- dat_wingspec$FrameID[m]
  dat_stab_adj$elbow[m]   <- dat_wingspec$elbow[m]
  dat_stab_adj$manus[m]   <- dat_wingspec$manus[m]
  
  if(nrow(dat_curr) < 4){next}
  mod.pstab = lm(Cm_CG ~ CL_adj, data = dat_curr)
  
  dat_stab_adj$cm0[m]     <- summary(mod.pstab)$coefficients[1,1]
  dat_stab_adj$cmcl[m]    <- summary(mod.pstab)$coefficients[2,1]
  dat_stab_adj$R2[m]      <- summary(mod.pstab)$r.squared
}
# remove the incomplete cases
dat_stab_adj <- dat_stab_adj[complete.cases(dat_stab_adj[,7]),]
dat_stab_adj$elbow_scale <- dat_stab_adj$elbow/1000
dat_stab_adj$manus_scale <- dat_stab_adj$manus/1000

mod_cmcl <- lmer(cmcl ~ elbow_scale*manus_scale + 
                   I(elbow_scale^2) + I(elbow_scale^3) +
                   I(manus_scale^2) + I(manus_scale^3) + (1|WingID), data = dat_stab_adj)

coef_all$y.model[2]    = "CmCL"
coef_all$intercept[2]  = summary(mod_cmcl)$coefficients["(Intercept)","Estimate"]
coef_all$elbow[2]      = summary(mod_cmcl)$coefficients["elbow_scale","Estimate"]/1000
coef_all$manus[2]      = summary(mod_cmcl)$coefficients["manus_scale","Estimate"]/1000
coef_all$elbow2[2]     = summary(mod_cmcl)$coefficients["I(elbow_scale^2)","Estimate"]/(1000^2)
coef_all$elbow3[2]     = summary(mod_cmcl)$coefficients["I(elbow_scale^3)","Estimate"]/(1000^3)
coef_all$manus2[2]     = summary(mod_cmcl)$coefficients["I(manus_scale^2)","Estimate"]/(1000^2)
coef_all$manus3[2]     = summary(mod_cmcl)$coefficients["I(manus_scale^3)","Estimate"]/(1000^3)
coef_all$elbowmanus[2] = summary(mod_cmcl)$coefficients["elbow_scale:manus_scale","Estimate"]/(1000^2)

# to show that the speed is not a significant effect on the model prediction - can verify from the two speeds tested in tunnel
mod_CD_check <- lm(CD_true ~ elbow*manus + L_comp + I(L_comp^2) + I(elbow^2) + I(manus^2) + U_des, data = subset(dat_exp, alpha > -5 & alpha < 5))
mod_Cm_check <- lm(m_comp ~ elbow*manus + L_comp + I(L_comp^2) + I(elbow^2) + I(manus^2) + U_des, data = subset(dat_exp, alpha > -5 & alpha < 5))
mod_CL_check <- lm(L_comp ~ elbow*manus + alpha + I(alpha^2) + I(elbow^2) + I(manus^2) + U_des, data = subset(dat_exp, alpha > -5 & alpha < 5))

## ----------------- Pitch rate derivatives ------------------
# CLq 
dat_q$xcg = predict(mod_xcg,dat_q)
dat_q$zcg = predict(mod_zcg,dat_q)
dat_q$CL_adj = dat_q$FL/(0.5*1.225*10^2*max(dat_num$S[which(dat_num$WingID == "17_0285")]))
dat_q$L_comp = dat_q$CL_adj
# predict the drag for this configuration
dat_q$CD_adj_exp = predict(mod_CD,dat_q)
dat_q$D_adj_exp = (0.5*1.225*10^2*max(dat_num$S[which(dat_num$WingID == "17_0285")]))*dat_q$CD_adj_exp

# adjust the pitching moment be about the true center of gravity
dat_q$M_CG = dat_q$Mm + ((dat_q$FL*cosd(dat_q$alpha)+dat_q$D_adj_exp*sind(dat_q$alpha))*(-dat_q$xcg) + 
                     (dat_q$FL*sind(dat_q$alpha)-dat_q$D_adj_exp*cosd(dat_q$alpha))*(-dat_q$zcg))

# Note this pitching moment is defined relative to the root chord
dat_q$Cm_CG <- dat_q$M_CG/(0.5*1.225*10^2*max(dat_num$S[which(dat_num$WingID == "17_0285")])*max(dat_num$root_c_max[which(dat_num$WingID == "17_0285")]))

dat_q_ind = as.data.frame(matrix(NA,nrow=61,ncol=8))
names(dat_q_ind) = c("species","WingID","TestID","FrameID","elbow","manus","CL_q","Cm_q")
count = 1

for (i in 1:length(unique(dat_q$FrameID))){
  curr_dat1 = subset(dat_q, FrameID == unique(dat_q$FrameID)[i])
  
  for(j in 1:length(unique(curr_dat1$alpha))){
    curr_dat = subset(dat_q, FrameID == unique(dat_q$FrameID)[i] & alpha == unique(curr_dat1$alpha)[j])
    if(nrow(curr_dat) < 4){next
      count = count + 1}
    mod_CL_q <- lm(CL_adj ~ q, data = curr_dat)
    mod_Cm_q <- lm(Cm_CG ~ q, data = curr_dat)
    
    dat_q_ind$species[count] <- as.character(curr_dat$species[1])
    dat_q_ind$WingID[count]  <- curr_dat$WingID[1]
    dat_q_ind$TestID[count]  <- curr_dat$TestID[1]
    dat_q_ind$FrameID[count] <- curr_dat$FrameID[1]
    dat_q_ind$alpha[count]   <- unique(curr_dat1$alpha)[j]
    dat_q_ind$elbow[count]   <- curr_dat$elbow[1]
    dat_q_ind$manus[count]   <- curr_dat$manus[1]
    
    dat_q_ind$CL_q[count]    <- coefficients(mod_CL_q)["q"]
    dat_q_ind$Cm_q[count]    <- coefficients(mod_Cm_q)["q"]
    dat_q_ind$CL_q_R2[count]    <- summary(mod_CL_q)$r.squared
    dat_q_ind$Cm_q_R2[count]    <- summary(mod_Cm_q)$r.squared
    dat_q_ind$no_sample[count]    <- nrow(curr_dat)
    count = count + 1
  }
}


mod_CL_q_ind <- lm(CL_q ~ elbow*manus + I(elbow^2) + I(manus^2), data = dat_q_ind)

coef_all$y.model[6]    = "dCLdq"
coef_all$intercept[6]  = coef(mod_CL_q_ind)["(Intercept)"]
coef_all$elbow[6]      = coef(mod_CL_q_ind)["elbow"]
coef_all$manus[6]      = coef(mod_CL_q_ind)["manus"]
coef_all$elbowmanus[6] = coef(mod_CL_q_ind)["elbow:manus"]
coef_all$elbow2[6]      = coef(mod_CL_q_ind)["I(elbow^2)"]
coef_all$manus2[6]      = coef(mod_CL_q_ind)["I(manus^2)"]


mod_Cm_q_ind <- lm(Cm_q ~ elbow*manus + I(elbow^2) + I(manus^2), data = dat_q_ind)

coef_all$y.model[7]    = "dCmdq"
coef_all$intercept[7]  = coef(mod_Cm_q_ind)["(Intercept)"]
coef_all$elbow[7]      = coef(mod_Cm_q_ind)["elbow"]
coef_all$manus[7]      = coef(mod_Cm_q_ind)["manus"]
coef_all$elbowmanus[7] = coef(mod_Cm_q_ind)["elbow:manus"]
coef_all$elbow2[7]      = coef(mod_Cm_q_ind)["I(elbow^2)"]
coef_all$manus2[7]      = coef(mod_Cm_q_ind)["I(manus^2)"]

write.csv(coef_all,'/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2021_07_29_coefficients.csv')

