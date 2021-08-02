library(ggplot2)
library(alphahull) # need for the convex hull
library(ptinpoly)  # need for determining points outside the convex hull
library(gridExtra) # for using grid arrange
library(cowplot)   # need for plot_grid()
library(ggthemes)  # need for geom_rangeframe

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

lab_elbow      = "Elbow angle (°)"
lab_manus      = "Wrist angle (°)"
lab_phase      = "Phase (°)"
# ------------ Read in data ---------------
dat_all <- read.csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2020_05_25_OrientedWings.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
dat_all <- subset(dat_all, species == "lar_gla" & sweep == 0 & dihedral == 0)

dat_out <- read.csv('LongDynStability_Rigid.csv')

# ------------ Clean data ---------------

# remove the unattainable configurations
tmp     <- cut_trueshape(dat_out,unique(subset(dat_all, elbow > 85 & manus > 100)[4:5]),4,5) # cut elbow and wrist to the true shape of the data
dat_cut <- tmp$dat_cut


# ------ Calculate new parameters -------

dat_cut$halft  = 0.69/abs(dat_cut$eig_real) # captures the rate of decay
dat_cut$period = 2*pi/dat_cut$omega_n


## ------------- Plot the key trim parameters -------------------

# speak to the glide angle in equilibrium and the associated lift to drag ratio
ggplot() + geom_point(data = subset(dat_cut,alpha==2 & eignum == 1), aes(x = elbow, y = manus, col = 1/tan((-theta0)))) + th

plot_trim <- ggplot() + 
  #add data
  geom_point(data = subset(dat_cut,alpha==2 & eignum == 1), 
             aes(x = U0, y = theta0*180/pi, col = manus)) + 
  #theme control
  th +
  # axis control 
  scale_x_continuous(limits = c(10,18), breaks = c(10,12,14,16,18), name = "Trim flight speed (m/s)") +
  scale_y_continuous(limits = c(-20,-5), breaks = c(-20,-15,-10,-5), name = "Trim glide angle (deg)") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 10, xend = 18, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -20, yend = -5)

## ------------ Plot the damping and frequency characteristics ------------------

# Define the axis ticks and associated labels
axis_begin  = -30
axis_end    = 30
total_ticks = 5
tick_frame <- data.frame(ticks = seq(axis_begin, axis_end, length.out = total_ticks), zero=0)
tick_frame <- subset(tick_frame, ticks != 0)
lab_frame <- data.frame(lab = seq(axis_begin, axis_end, length.out = total_ticks),zero = 0)
lab_frame <- subset(lab_frame, lab != 0)
tick_sz <- (tail(lab_frame$lab, 1) -  lab_frame$lab[1]) / 128

locus <- ggplot() + 
  #theme control
  th +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  # x ticks
  geom_segment(data = tick_frame, 
               aes(x = ticks, xend = ticks, 
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  # y ticks
  geom_segment(data = tick_frame, 
               aes(x = zero - tick_sz, xend = zero + tick_sz, 
                   y = ticks, yend = ticks)) + 
  # labels
  geom_text(data=lab_frame, aes(x=lab, y=zero, label=lab), 
            vjust=1.7, size = 8*5/14, family = "sans") +
  geom_text(data=lab_frame, aes(x=zero, y=lab, label=lab),
            hjust=-0.3, size = 8*5/14, family = "sans") +
  geom_rangeframe() +
  annotate(geom = "segment", x = -30, xend = 30, y = 0, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 0, y = -30, yend = 30) +
  scale_x_continuous(name = "Real") +
  scale_y_continuous(name = "Imaginary") + 
  # add data
  geom_point(data = subset(dat_cut,alpha==2), 
             aes(x = eig_real, y = eig_imag, col = manus, shape = as.factor(eignum))) +
  scale_shape_manual(values = c(16,16,2,2))

## ------------------- Plot damping and frequency results ----------------------

# short period mode - Frequency
plot_sp_o <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==2 & eignum == 1), 
             aes(x = manus, y = omega_n, col = elbow)) + th +
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,20), name = expression(paste(omega)["n"])) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 20) + 
  annotate(geom = "segment", x = -log(0), xend = -log(0), y = 0, yend = 2)


# short period mode - Damping
plot_sp_z <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==2 & eignum == 1), 
             aes(x = manus, y = zeta, col = elbow)) + th +
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,1), name = expression(paste(zeta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 1)


# phugoid mode - Frequency
plot_ph_o <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==2 & eignum == 3), 
             aes(x = manus, y = omega_n, col = elbow)) + th  +
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,2), name = expression(paste(omega)["n"])) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 2)

# phugoid mode - Damping
plot_ph_z <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==2 & eignum == 3), 
             aes(x = manus, y = zeta, col = elbow)) + th +
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,1), name = expression(paste(zeta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 1)


plot_oz <- plot_grid(plot_sp_o,plot_ph_o,plot_sp_z,plot_ph_z,
                    #arrangement data
                    ncol = 2,
                    #labels
                    labels = c("A","B","C","D"),
                    label_size = 10,
                    label_fontfamily = "sans")

## ------------- Plot the phase and magnitude -------------------
dat_sp <- subset(dat_cut, alpha==4 & eignum == 1 & manus < 160) # eigen num == 2 just flips phase
# Set the pitch rate as the comparable phase 
dat_sp$phase1 = dat_sp$phase1 - dat_sp$phase3
dat_sp$phase2 = dat_sp$phase2 - dat_sp$phase3
dat_sp$phase3 = dat_sp$phase3 - dat_sp$phase3
dat_sp$phase4 = dat_sp$phase4 - dat_sp$phase3
# adjust negative numbers to be in a positive frame
dat_sp$phase1[which(dat_sp$phase1 < 0)] = (2*pi) + dat_sp$phase1[which(dat_sp$phase1 < 0)]
dat_sp$phase2[which(dat_sp$phase2 < 0)] = (2*pi) + dat_sp$phase2[which(dat_sp$phase2 < 0)]
dat_sp$phase3[which(dat_sp$phase3 < 0)] = (2*pi) + dat_sp$phase3[which(dat_sp$phase3 < 0)]
dat_sp$phase4[which(dat_sp$phase4 < 0)] = (2*pi) + dat_sp$phase4[which(dat_sp$phase4 < 0)]

col_u     = "#2d8183"
col_alpha = "#83d3d4"
col_q     = "#910c07"
col_theta = "#f48153"
  
plot_sp_magphase <- ggplot() + 
  geom_vline(xintercept = seq(0, 1, by = 0.25), colour = "grey70", size = 0.2) +
  geom_hline(yintercept = seq(0, 315, by = 45), colour = "grey70", size = 0.2) +
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag1, y = phase1*180/pi, yend = phase1*180/pi, group = del_z), col = col_u, alpha = 0.3) + 
  geom_point(data = dat_sp, 
             aes(x = mag1, y = phase1*180/pi,  alpha = 0.5), 
             size = 2, col = col_u) + 
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag2, y = phase2*180/pi, yend = phase2*180/pi, group = del_z), col = col_alpha, alpha = 0.3) + 
  geom_point(data = dat_sp, 
             aes(x = mag2, y = phase2*180/pi,  alpha = 0.5), 
             size = 2, col = col_alpha) +
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag3, y = phase3*180/pi, yend = phase3*180/pi, group = del_z), col = col_theta, alpha = 0.3) + 
  geom_point(data = dat_sp, 
             aes(x = mag3, y = phase3*180/pi,  alpha = 0.5), 
             size = 2, col = col_q) + 
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag4, y = phase4*180/pi, yend = phase4*180/pi, group= del_z), col = col_q, alpha = 0.3) + 
  geom_point(data = dat_sp, 
             aes(x = mag4, y = phase4*180/pi,  alpha = 0.5), 
             size = 2, col = col_theta) + 
  th +
  theme_light()+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank()) +
  panel_border(remove = TRUE) + 
  # axis control
  coord_polar(theta = "y", start = -0.5*pi, direction = -1) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), name = "Magnitude") +
  scale_y_continuous(limits = c(0,360), breaks = c(0,90,180,270), name = lab_phase)
  
## ----------- Phugoid mode phase and magnitude -----------------

dat_ph <- subset(dat_cut, alpha==4 & eignum == 3 & manus < 160) # eignum == 4 just flips phase
# Set the pitch rate as the comparable phase 
dat_ph$phase1 = dat_ph$phase1 - dat_ph$phase3
dat_ph$phase2 = dat_ph$phase2 - dat_ph$phase3
dat_ph$phase3 = dat_ph$phase3 - dat_ph$phase3
dat_ph$phase4 = dat_ph$phase4 - dat_ph$phase3
# adjust negative numbers to be in a positive frame
dat_ph$phase1[which(dat_ph$phase1 < 0)] = (2*pi) + dat_ph$phase1[which(dat_ph$phase1 < 0)]
dat_ph$phase2[which(dat_ph$phase2 < 0)] = (2*pi) + dat_ph$phase2[which(dat_ph$phase2 < 0)]
dat_ph$phase3[which(dat_ph$phase3 < 0)] = (2*pi) + dat_ph$phase3[which(dat_ph$phase3 < 0)]
dat_ph$phase4[which(dat_ph$phase4 < 0)] = (2*pi) + dat_ph$phase4[which(dat_ph$phase4 < 0)]

plot_ph_magphase <- ggplot() + 
  geom_vline(xintercept = seq(0, 1, by = 0.25), colour = "grey70", size = 0.2) +
  geom_hline(yintercept = seq(0, 315, by = 45), colour = "grey70", size = 0.2) +
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag1, y = phase1*180/pi, yend = phase1*180/pi, group = del_z), col = col_u, alpha = 0.3) + 
  geom_point(data = dat_ph, 
             aes(x = mag1, y = phase1*180/pi,  alpha = 0.5), 
             size = 2, col = col_u) + 
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag2, y = phase2*180/pi, yend = phase2*180/pi, group = del_z), col = col_alpha, alpha = 0.3) + 
  geom_point(data = dat_ph, 
             aes(x = mag2, y = phase2*180/pi,  alpha = 0.5), 
             size = 2, col = col_alpha) +
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag3, y = phase3*180/pi, yend = phase3*180/pi, group = del_z), col = col_theta, alpha = 0.3) + 
  geom_point(data = dat_ph, 
             aes(x = mag3, y = phase3*180/pi,  alpha = 0.5), 
             size = 2, col = col_q) + 
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag4, y = phase4*180/pi, yend = phase4*180/pi, group= del_z), col = col_q, alpha = 0.3) + 
  geom_point(data = dat_ph, 
             aes(x = mag4, y = phase4*180/pi,  alpha = 0.5), 
             size = 2, col = col_theta) + 
  th +
  theme_light()+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank()) +
  panel_border(remove = TRUE) + 
  # axis control
  coord_polar(theta = "y", start = -0.5*pi, direction = -1) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), name = "Magnitude") +
  scale_y_continuous(limits = c(0,360), breaks = c(0,90,180,270), name = lab_phase)

plot_grid(plot_sp_magphase,plot_ph_magphase,
          #arrangement data
          ncol = 2,
          #labels
          labels = c("A","B"),
          label_size = 10,
          label_fontfamily = "sans")

## ------------------- Phase with respect to elbow and wrist angle ------------
ggplot() + geom_point(data = dat_cut, aes(x = elbow, y = manus, alpha = phase2)) + th



## ------------------ Time response to an initial alpha -------------------
dat_time_1 <- read.csv('./outputdata/2021_08_01_elbow90_manus155_alpha4.csv', header = FALSE)
colnames(dat_time_1) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_2 <- read.csv('./outputdata/2021_08_01_elbow145_manus155_alpha4.csv', header = FALSE)
colnames(dat_time_2) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_3 <- read.csv('./outputdata/2021_08_01_elbow90_manus105_alpha4.csv', header = FALSE)
colnames(dat_time_3) <- c("t","del_u","del_alp","del_q","del_theta")


plot_time1_dalp <- plot_timeseries(dat_time_1,col_u,col_alpha,col_q,col_theta)
plot_time2_dalp <- plot_timeseries(dat_time_2,col_u,col_alpha,col_q,col_theta)
plot_time3_dalp <- plot_timeseries(dat_time_3,col_u,col_alpha,col_q,col_theta)

dat_time_3 <- read.csv('./outputdata/2021_08_02_elbow145_manus155_alpha4_uramp.csv', header = FALSE)
colnames(dat_time_3) <- c("t","del_u","del_alp","del_q","del_theta")
