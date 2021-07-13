library(lme4)
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


coef_all = as.data.frame(matrix(0,nrow = 7, ncol = 22))
colnames(coef_all) <- c("y.model","intercept","elbow","manus","elbow2","manus2",
                        "elbow3","manus3","elbowmanus",
                        "alpha","alpha2","alpha3",
                        "elbowalpha","manusalpha","elbowmanusalpha",
                        "CL","CL2","CL3",
                        "elbowCL","manusCL","elbowmanusCL")

## ---------------------------------------------
## ------------- Inertia Data ------------------
## ---------------------------------------------

# subset to the same range as used in the aerodynamic data
dat_inertial = subset(dat_final,species == "lar_gla" & elbow > 80 & manus > 100)
remove(dat_final)

max(dat_inertial$S_max)
max(dat_inertial$c_max)
max(dat_inertial$b_max)
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

## -------------------------------------------------
## ------------- Aerodynamic Data ------------------
## -------------------------------------------------

# scaling was elbow/1000 manus/1000 alpha/10
# If I want to input the true values into this same model I instead
# need to divide each coefficient by however many inputs are being changed

# Cm/CL data
mod_cmcl <- summary(mod_cmcl_num)

coef_all$y.model[2]    = "CmCL"
coef_all$intercept[2]  = mod_cmcl$coefficients["(Intercept)","Estimate"]
coef_all$elbow[2]      = mod_cmcl$coefficients["elbow_scale","Estimate"]/1000
coef_all$manus[2]      = mod_cmcl$coefficients["manus_scale","Estimate"]/1000
coef_all$elbow2[2]     = mod_cmcl$coefficients["I(elbow_scale^2)","Estimate"]/(1000^2)
coef_all$elbow3[2]     = mod_cmcl$coefficients["I(elbow_scale^3)","Estimate"]/(1000^3)
coef_all$manus2[2]     = mod_cmcl$coefficients["I(manus_scale^2)","Estimate"]/(1000^2)
coef_all$manus3[2]     = mod_cmcl$coefficients["I(manus_scale^3)","Estimate"]/(1000^3)
coef_all$elbowmanus[2] = mod_cmcl$coefficients["elbow_scale:manus_scale","Estimate"]/(1000^2)

# CL data
mod_CL <- summary(mod_con_cL_num)

coef_all$y.model[3]    = "CL"
coef_all$intercept[3]  = mod_CL$coefficients["(Intercept)","Estimate"]
coef_all$elbow[3]      = mod_CL$coefficients["elbow_scale","Estimate"]/1000
coef_all$manus[3]      = mod_CL$coefficients["manus_scale","Estimate"]/1000
coef_all$elbow2[3]     = mod_CL$coefficients["I(elbow_scale^2)","Estimate"]/(1000^2)
coef_all$manus2[3]     = mod_CL$coefficients["I(manus_scale^2)","Estimate"]/(1000^2)
coef_all$manus3[3]     = mod_CL$coefficients["I(manus_scale^3)","Estimate"]/(1000^3)
coef_all$elbowmanus[3] = mod_CL$coefficients["elbow_scale:manus_scale","Estimate"]/(1000^2)
coef_all$alpha[3]      = mod_CL$coefficients["alpha_scale","Estimate"]/10
coef_all$alpha2[3]     = mod_CL$coefficients["I(alpha_scale^2)","Estimate"]/(10^2)
coef_all$alpha3[3]     = mod_CL$coefficients["I(alpha_scale^3)","Estimate"]/(10^2)
coef_all$elbowalpha[3] = mod_CL$coefficients["elbow_scale:alpha_scale","Estimate"]/(1000*10)
coef_all$manusalpha[3] = mod_CL$coefficients["manus_scale:alpha_scale","Estimate"]/(1000*10)
coef_all$elbowmanusalpha[3] = mod_CL$coefficients["elbow_scale:manus_scale:alpha_scale","Estimate"]/((1000^2)*10)

# Cm data

mod_Cm <- summary(mod_con_cm_num)

coef_all$y.model[4]    = "Cm"
coef_all$intercept[4]  = mod_Cm$coefficients["(Intercept)","Estimate"]
coef_all$elbow[4]      = mod_Cm$coefficients["elbow_scale","Estimate"]/1000
coef_all$manus[4]      = mod_Cm$coefficients["manus_scale","Estimate"]/1000
coef_all$elbow2[4]     = mod_Cm$coefficients["I(elbow_scale^2)","Estimate"]/(1000^2)
coef_all$manus2[4]     = mod_Cm$coefficients["I(manus_scale^2)","Estimate"]/(1000^2)
coef_all$manus3[4]     = mod_Cm$coefficients["I(manus_scale^3)","Estimate"]/(1000^3)
coef_all$elbowmanus[4] = mod_Cm$coefficients["elbow_scale:manus_scale","Estimate"]/(1000^2)
coef_all$CL[4]         = mod_Cm$coefficients["CL_adj","Estimate"]
coef_all$CL2[4]        = mod_Cm$coefficients["I(CL_adj^2)","Estimate"]
coef_all$CL3[4]        = mod_Cm$coefficients["I(CL_adj^3)","Estimate"]
coef_all$elbowCL[4]    = mod_Cm$coefficients["elbow_scale:CL_adj","Estimate"]/(1000)
coef_all$manusCL[4]    = mod_Cm$coefficients["manus_scale:CL_adj","Estimate"]/(1000)
coef_all$elbowmanusCL[4] = mod_Cm$coefficients["elbow_scale:manus_scale:CL_adj","Estimate"]/(1000^2)

# CD data
# need to fit to experimental data due to the reduced
dat_exp$D_true     <- dat_exp$D_comp/(0.5*max(dat_num$S[which(dat_num$WingID == "17_0285")]))

# minor changes in the coefficients between speeds, will group them
mod_CD <- lm(D_true ~ elbow + manus + L_comp + I(L_comp^2) + I(elbow^2) + I(manus^2), data = subset(dat_exp, alpha > -5 & alpha < 5))

coef_all$y.model[5]    = "CD"
coef_all$intercept[5]  = coef(mod_CD)["(Intercept)"]
coef_all$elbow[5]      = coef(mod_CD)["elbow"]
coef_all$manus[5]      = coef(mod_CD)["manus"]
coef_all$elbow2[5]     = coef(mod_CD)["I(elbow^2)"]
coef_all$manus2[5]     = coef(mod_CD)["I(manus^2)"]
coef_all$CL[5]         = coef(mod_CD)["L_comp"]
coef_all$CL2[5]        = coef(mod_CD)["I(L_comp^2)"]
filename = paste(format(Sys.Date(), "%Y_%m_%d"),"_coefficients.csv",sep="")
write.csv(coef_all,filename)

# to show that the speed is not a significant effect on the model prediction
mod_CD_check <- lm(D_true ~ elbow*manus + L_comp + I(L_comp^2) + I(elbow^2) + I(manus^2) + U_des, data = subset(dat_exp, alpha > -5 & alpha < 5))
mod_Cm_check <- lm(m_comp ~ elbow*manus + L_comp + I(L_comp^2) + I(elbow^2) + I(manus^2) + U_des, data = subset(dat_exp, alpha > -5 & alpha < 5))
mod_CL_check <- lm(L_comp ~ elbow*manus + alpha + I(alpha^2) + I(elbow^2) + I(manus^2) + U_des, data = subset(dat_exp, alpha > -5 & alpha < 5))

# CLq 

dat_q$CL_adj <- dat_q$FL/(0.5*1.225*10^2*max(dat_num$S[which(dat_num$WingID == "17_0285")]))
dat_q$Cm_adj <- dat_q$Mm/(0.5*1.225*10^2*max(dat_num$S[which(dat_num$WingID == "17_0285")])*max(dat_num$ref_c[which(dat_num$WingID == "17_0285")]))

mod_CL_q <- lm(CL_adj ~ I(elbow^2) + I(elbow^3) +
                    I(manus^2)+ I(manus^3) + alpha + I(alpha^2) + I(alpha^3) + q*elbow*manus, dat_q)

coef_all$y.model[6]    = "dCLdq"
coef_all$intercept[6]  = coef(mod_CL_q)["q"]
coef_all$elbow[6]      = coef(mod_CL_q)["q:elbow"]
coef_all$manus[6]      = coef(mod_CL_q)["q:manus"]
coef_all$elbowmanus[6] = coef(mod_CL_q)["q:elbow:manus"]

# Check how close this value is dCm/dq = dCm/dCL*dCL/dq
mod_Cm_q <- lm(Cm_adj ~ I(elbow^2) + I(elbow^3) +
                 I(manus^2)+ I(manus^3) + alpha + I(alpha^2) + I(alpha^3) + q*elbow*manus, dat_q)

coef_all$y.model[7]    = "dCmdq"
coef_all$intercept[7]  = coef(mod_Cm_q)["q"]
coef_all$elbow[7]      = coef(mod_Cm_q)["q:elbow"]
coef_all$manus[7]      = coef(mod_Cm_q)["q:manus"]
coef_all$elbowmanus[7] = coef(mod_Cm_q)["q:elbow:manus"]

write.csv(coef_all,'/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2021_07_13_coefficients.csv')

