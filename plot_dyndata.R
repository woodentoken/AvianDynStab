library(ggplot2)
library(alphahull) # need for the convex hull
library(ptinpoly)  # need for determining points outside the convex hull
library(gridExtra) # for using grid arrange
library(cowplot)   # need for plot_grid()

source("support_functions.R")
# -------- Main theme ------------
th <- theme_classic() +
  theme(
    # Text
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(margin = margin(t = 10, unit = "pt")),
    axis.text.y = element_text(margin = margin(r = 10)),
    # Axis line
    axis.line = element_blank(),
    axis.ticks.length = unit(-5,"pt"),
    # Background transparency
    # Background of panel
    panel.background = element_rect(fill = "transparent"),
    # Background behind actual data points
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# ------------ Read in data ---------------
dat_out <- read.csv('LongDynStability_Rigid.csv')
# remove the unattainable configurations
tmp     <- cut_trueshape(dat_out,unique(subset(dat_all, elbow > 85 & manus > 100)[4:5]),4,5) # cut elbow and wrist to the true shape of the data
dat_cut <- tmp$dat_cut

dat_cut$halft  = 0.69/abs(dat_cut$eig_real)
dat_cut$period = 2*pi/dat_cut$omega_n

dat_all <- read.csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2020_05_25_OrientedWings.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
dat_all <- subset(dat_all, species == "lar_gla" & sweep == 0 & dihedral == 0)
dat_all$root_c = sqrt((dat_all$Pt12X - dat_all$Pt11X)^2 + (dat_all$Pt12Z - dat_all$Pt11Z)^2)
dat_all$FrameID <- paste("F", dat_all$frameID, sep = "")


ggplot() + geom_point(data = subset(dat_cut,alpha==4), aes(x = zeta, y = omega_n, col = elbow))
ggplot() + geom_point(data = subset(dat_cut,alpha==4), aes(x = eig_real, y = eig_imag, col = elbow, shape = as.factor(eignum)))

ggplot() + geom_point(data = subset(dat_cut, alpha==4 & eignum == 1), aes(x = manus, y = elbow, col = theta0*180/pi)) + th

ggplot() + geom_point(data = subset(dat_cut,alpha==4), aes(x = eignum, y = mag1, col = elbow)) + geom_point(data = dat_exp, aes(x = CD_true, y = L_comp, col = manus))

ggplot() + geom_point(data = subset(dat_cut,alpha==4), aes(x = elbow, y = manus, col = 1/tan((-theta0))))
ggplot() + geom_point(data = subset(dat_cut,alpha==4), aes(x = elbow, y = manus, col = CL_q)) + geom_point(data = dat_q_ind, aes(x = elbow, y = manus, col = CL_q), shape = 15)
ggplot() + geom_point(data = subset(dat_cut,alpha==4 & eignum == 1), aes(x = manus, y = eig_imag, col = elbow))

dat_time <- read.csv('./outputdata/2021_07_30_elbow140_manus125_alpha4.csv', header = FALSE)
colnames(dat_time) <- c("t","del_u","del_alp","del_q","del_theta")
plot_del_u     = ggplot() + geom_line(data = dat_time, aes(x = t, y = del_u), col = "red")+ th
plot_del_alp   = ggplot() + geom_line(data = dat_time, aes(x = t, y = del_alp*180/pi), col = "blue") + th
plot_del_q     = ggplot() + geom_line(data = dat_time, aes(x = t, y = del_q*180/pi), col = "green") + th
plot_del_theta = ggplot() + geom_line(data = dat_time, aes(x = t, y = del_theta*180/pi), col = "purple") + th #+ scale_y_continuous(limits = c(-0.05,0.05)) 
  

plot_grid(plot_del_u,plot_del_alp,plot_del_q,plot_del_theta,
                    #arrangement data
                    ncol = 1,
                    #labels
                    labels = c("A","B","C","D"),
                    label_size = 10,
                    label_fontfamily = "sans")

plot(dat_time[,1],dat_time[,2])


# maximum CL/CD
1/tan(min(-dat_cut$theta0))
dat_cut[which.min(-dat_cut$theta0),c("elbow","manus","alpha")]

