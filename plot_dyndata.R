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

# blank plot to fill space
blank_plot <- ggplot() + theme_void()

lab_elbow      = "Elbow angle (°)"
lab_manus      = "Wrist angle (°)"
lab_phase      = "Phase (°)"
lab_clcd       = expression(paste("C" [L], "/C" [D]))

col_u     = "#2d8183"
col_alpha = "#83d3d4"
col_q     = "#910c07"
col_theta = "#f48153"

#---- Elbow colour scale
cc4  <- scales::seq_gradient_pal("#6A1B9A", "black", "Lab")(seq(0,1,length.out=16))
cc3  <- scales::seq_gradient_pal("#E85285", "#6A1B9A", "Lab")(seq(0,1,length.out=42))
cc2  <- scales::seq_gradient_pal("#FFECB3", "#E85285", "Lab")(seq(0,1,length.out=60))
cc1  <- scales::seq_gradient_pal("white", "#FFECB3", "Lab")(seq(0,1,length.out=29))
cc_full   <- c(cc1[26:29],cc2[2:60],cc3[2:40],cc4[2:15])
#Done based on the differences between the angles
cc <- c(cc_full[1],cc_full[16],cc_full[19],cc_full[42],cc_full[47],cc_full[49],cc_full[68],cc_full[80],cc_full[82],cc_full[105],cc_full[106],cc_full[116])

cc_elbow <- c("70" = cc_full[36],
             "80" = cc_full[46],
             "90" = cc_full[56],
             "100" = cc_full[66],
             "110" = cc_full[76],
             "120" = cc_full[86],
             "130" = cc_full[96],
             "140" = cc_full[106],
             "150" = cc_full[116],
             "160" = "black",
             "170" = "black")

alp_fixed = 2 
# ------------ Read in data ---------------
dat_all <- read.csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2020_05_25_OrientedWings.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
dat_all <- subset(dat_all, species == "lar_gla" & sweep == 0 & dihedral == 0)

dat_out <- read.csv('LongDynStability_Rigid.csv')

# ------------ Clean data ---------------

# remove the unattainable configurations
tmp     <- cut_trueshape(dat_out,unique(subset(dat_all, elbow > 85 & manus > 100)[4:5]),4,5) # cut elbow and wrist to the true shape of the data
dat_cut <- tmp$dat_cut

dat_cut$halft  = 0.69/abs(dat_cut$eig_real) # captures the rate of decay
dat_cut$period = 2*pi/dat_cut$omega_n

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


dat_sp <- subset(dat_cut, alpha==alp_fixed & eignum == 1 & omega_n > 0) # eigen num == 2 just flips phase
# need to make sure that the configurations that are unstable in sp mode are removed here as well
dat_ph <- merge(dat_sp[,c("alpha","elbow","manus")], subset(dat_cut, alpha==alp_fixed & eignum == 3 & omega_n > 0)) # eignum == 4 just flips phase
## ------------- Plot the key trim parameters -------------------

plot_trim <- ggplot() + 
  #add data
  geom_point(data = subset(dat_cut,alpha==alp_fixed & eignum == 1), 
             aes(x = U0, y = theta0*180/pi, col = elbow)) + 
  #theme control
  th +
  scale_color_gradientn(colours = cc_full, name = lab_elbow) + 
  # axis control 
  scale_x_continuous(limits = c(10,25), breaks = c(10,15,20,25), name = "Trim flight speed (m/s)") +
  scale_y_continuous(limits = c(-50,0), breaks = c(-45,-30,-15,0), name = "Trim glide angle (deg)") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 10, xend = 25, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -45, yend = 0)

plot_clcd <- ggplot() + 
  #add data
  geom_raster(data = subset(dat_cut,alpha==alp_fixed & eignum == 1), 
             aes(x = elbow, y = manus, fill = 1/tan(-theta0))) + 
  #theme control
  th +
  scale_fill_gradient(high = "#D9D9E0", low = "#1D1D21", name = lab_clcd) +
  # axis control 
  coord_fixed() + 
  scale_x_continuous(limits = c(85,170), breaks = c(90,120,150), name = lab_elbow) +
  scale_y_continuous(limits = c(100,180), breaks = c(100,120,140,160,180), name = lab_manus) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 90, xend = 150, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 100, yend = 180) 

plot_fig1 <- plot_grid(blank_plot,plot_trim,plot_clcd,
                     #arrangement data
                     ncol = 3,
                     rel_widths = c(1,1,1),
                     #labels
                     labels = c("A","B","C"),
                     label_size = 10,
                     label_fontfamily = "sans")

## ------------ Plot the damping and frequency characteristics ------------------

# Define the axis ticks and associated labels
axis_begin  = -30
axis_end    = 30
total_ticks = 5
tick_frame <- data.frame(ticks = seq(axis_begin, axis_end, length.out = total_ticks), zero=0)
tick_frame <- subset(tick_frame, ticks != 0)
lab_frame <- data.frame(lab = seq(axis_begin, axis_end, length.out = total_ticks),zero = 0)
lab_frame <- subset(lab_frame, lab != 0)
tick_sz_sp <- (tail(lab_frame$lab, 1) -  lab_frame$lab[1]) / 128

locus_sp <- ggplot() + 
  #theme control
  th +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') +
  scale_color_gradientn(colours = cc_full, name = lab_elbow) + 
  ggtitle("Short period mode") +
  # axis control
  # x ticks
  geom_segment(data = tick_frame, 
               aes(x = ticks, xend = ticks, 
                   y = zero - tick_sz_sp, yend = zero + tick_sz_sp)) +
  # y ticks
  geom_segment(data = tick_frame, 
               aes(x = zero - tick_sz_sp, xend = zero + tick_sz_sp, 
                   y = ticks, yend = ticks)) + 
  # labels
  geom_text(data=lab_frame, aes(x=lab, y=zero, label=lab), 
            vjust=1.7, size = 8*5/14, family = "sans") +
  geom_text(data=lab_frame, aes(x=zero, y=lab, label=lab),
            hjust=-0.5, size = 8*5/14, family = "sans") +
  geom_rangeframe() +
  annotate(geom = "segment", x = -30, xend = 30, y = 0, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 0, y = -30, yend = 30) +
  scale_x_continuous(name = "Real") +
  scale_y_continuous(name = "Imaginary") + 
  # add data
  geom_point(data = subset(dat_cut,alpha==alp_fixed & eignum < 3), 
             aes(x = eig_real, y = eig_imag, col = elbow)) +
  scale_shape_manual(values = c(16,16,2,2))

# Define the axis ticks and associated labels
axis_begin  = -5
axis_end    = 5
total_ticks = 5
tick_frame <- data.frame(ticks = seq(axis_begin, axis_end, length.out = total_ticks), zero=0)
tick_frame <- subset(tick_frame, ticks != 0)
lab_frame <- data.frame(lab = seq(axis_begin, axis_end, length.out = total_ticks),zero = 0)
lab_frame <- subset(lab_frame, lab != 0)
tick_sz <- (tail(lab_frame$lab, 1) -  lab_frame$lab[1]) / 128

locus_ph <- ggplot() + 
  #theme control
  th +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') +
  scale_color_gradientn(colours = cc_full, name = lab_elbow) + 
  ggtitle("Phugoid mode") +
  # axis control
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
            hjust=-0.5, size = 8*5/14, family = "sans") +
  geom_rangeframe() +
  annotate(geom = "segment", x = -5, xend = 5, y = 0, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 0, y = -5, yend = 5) +
  scale_x_continuous(name = "Real") +
  scale_y_continuous(name = "Imaginary") + 
  # add data
  geom_point(data = subset(dat_cut,alpha==alp_fixed & eignum > 2), 
             aes(x = eig_real, y = eig_imag, col = elbow)) +
  scale_shape_manual(values = c(16,16,2,2))


## ------------------- Plot damping and frequency results ----------------------

# short period mode - Frequency
plot_sp_o <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==alp_fixed & eignum == 1), 
             aes(x = manus, y = omega_n, col = elbow)) + 
  th +
  # colour control
  scale_color_gradientn(colours = cc_full, name = lab_elbow) + 
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,20), name = expression(paste(omega)["n"])) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 20) + 
  annotate(geom = "segment", x = -log(0), xend = -log(0), y = 0, yend = 2)

# Extract the legend to plot later
leg_elbow <- g_legend(plot_sp_o)
plot_sp_o <- plot_sp_o + theme(legend.position = 'none')


# short period mode - Damping
plot_sp_z <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==alp_fixed & eignum == 1), 
             aes(x = manus, y = zeta, col = elbow)) + 
  # theme control
  th +
  theme(legend.position = 'none') +
  # colour control
  scale_color_gradientn(colours = cc_full) + 
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,1), name = expression(paste(zeta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 1)


# phugoid mode - Frequency
plot_ph_o <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==alp_fixed & eignum == 3), 
             aes(x = manus, y = omega_n, col = elbow)) + 

  # theme control
  th +
  theme(legend.position = 'none') +
  # colour control
  scale_color_gradientn(colours = cc_full) + 
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,2), name = expression(paste(omega)["n"])) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 2)

# phugoid mode - Damping
plot_ph_z <- ggplot() + 
  geom_point(data = subset(dat_cut,alpha==alp_fixed & eignum == 3), 
             aes(x = manus, y = zeta, col = elbow)) + 
  # theme control
  th +
  theme(legend.position = 'none') +
  # colour control
  scale_color_gradientn(colours = cc_full) + 
  # axis control 
  scale_x_continuous(limits = c(100,180), name = lab_manus) +
  scale_y_continuous(limits = c(0,1), name = expression(paste(zeta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 100, xend = 180, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 1)


plot_oz <- plot_grid(locus_sp, locus_ph, plot_sp_o,plot_ph_o,plot_sp_z,plot_ph_z,
                    #arrangement data
                    ncol = 2,
                    rel_widths = c(1,1),
                    #labels
                    labels = c("A","B","C","D","E","F"),
                    label_size = 10,
                    label_fontfamily = "sans")

## ------------- Plot the phase and magnitude -------------------

  
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
  # theme control
  th +
  theme_light()+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  panel_border(remove = TRUE) + 
  # axis control
  coord_polar(theta = "y", start = -0.5*pi, direction = -1) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), name = "Magnitude") +
  scale_y_continuous(limits = c(0,360), breaks = c(0,90,180,270), name = lab_phase)
  
## ----------- Phugoid mode phase and magnitude -----------------


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
  # theme control
  th +
  theme_light()+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  panel_border(remove = TRUE) + 
  # axis control
  coord_polar(theta = "y", start = -0.5*pi, direction = -1) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), name = "Magnitude") +
  scale_y_continuous(limits = c(0,360), breaks = c(0,90,180,270), name = lab_phase)


## ------------------- Phase with respect to elbow and wrist angle ------------
plot_ph_phase <- ggplot() + geom_raster(data = dat_ph, 
                                        aes(x = elbow, y = manus, alpha = round(phase4,1)*180/pi), fill = col_theta) + 
  # theme control
  th +
  scale_alpha_continuous(name = lab_phase) +
  # axis control 
  coord_fixed() + 
  scale_x_continuous(limits = c(85,156), breaks = c(90,120,150), name = lab_elbow) +
  scale_y_continuous(limits = c(100,161), breaks = c(100,130,160), name = lab_manus) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 90, xend = 150, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 100, yend = 160) 

plot_sp_phase <- ggplot() + 
  geom_raster(data = dat_sp, aes(x = elbow, y = manus, alpha = round(phase4,1)*180/pi), fill = col_theta) + 
  # theme control
  th +  scale_alpha_continuous(name = lab_phase) +
  # axis control 
  coord_fixed() + 
  scale_x_continuous(limits = c(85,156), breaks = c(90,120,150), name = lab_elbow) +
  scale_y_continuous(limits = c(100,161), breaks = c(100,130,160), name = lab_manus) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 90, xend = 150, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 100, yend = 160) 


plot_fig3 <- plot_grid(plot_sp_magphase,plot_ph_magphase,
          plot_sp_phase,plot_ph_phase,
          #arrangement data
          ncol = 2,
          #labels
          labels = c("A","B", "C","D"),
          label_size = 10,
          label_fontfamily = "sans")



## ------------------ Time response to an initial alpha -------------------
dat_time_1 <- read.csv('./outputdata/2021_08_03_elbow90_manus120_alpha2_dalp.csv', header = FALSE)
colnames(dat_time_1) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_2 <- read.csv('./outputdata/2021_08_03_elbow120_manus120_alpha2_dalp.csv', header = FALSE)
colnames(dat_time_2) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_3 <- read.csv('./outputdata/2021_08_03_elbow120_manus157_alpha2_dalp.csv', header = FALSE)
colnames(dat_time_3) <- c("t","del_u","del_alp","del_q","del_theta")

lim_u     = c(-0.05,0.05)
lim_alpha = c(-5,5)
lim_q     = c(-63,63)
break_q = c(-60,-30,0,30,60)
lim_theta = c(-8,8)
break_theta = c(-8,-4,0,4,8)
plot_time1_dalp <- plot_timeseries(dat_time_1,
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, 40)
plot_time2_dalp <- plot_timeseries(dat_time_2,
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, 40)
plot_time3_dalp <- plot_timeseries(dat_time_3,
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, 40)

plot_fig4 <- plot_grid(plot_time1_dalp,plot_time2_dalp,plot_time3_dalp,
                       #arrangement data
                       ncol = 3,
                       #labels
                       labels = c("A","B", "C"),
                       label_size = 10,
                       label_fontfamily = "sans")



lim_u     = c(0,0.25)
lim_alpha = c(-0.01,0.01)
lim_q     = c(-1,1)
break_q = c(-1,-0.5,0,0.5,1)
lim_theta = c(-10,0)
break_theta = c(-10,-5,0)

dat_time_4 <- read.csv('./outputdata/2021_08_03_elbow90_manus120_alpha2_uramp.csv', header = FALSE)
colnames(dat_time_4) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_5 <- read.csv('./outputdata/2021_08_03_elbow120_manus120_alpha2_uramp.csv', header = FALSE)
colnames(dat_time_5) <- c("t","del_u","del_alp","del_q","del_theta")
dat_time_6 <- read.csv('./outputdata/2021_08_03_elbow120_manus157_alpha2_uramp.csv', header = FALSE)
colnames(dat_time_6) <- c("t","del_u","del_alp","del_q","del_theta")

plot_time4_dalp <- plot_timeseries(dat_time_4,
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, 40)
plot_time5_dalp <- plot_timeseries(dat_time_5,
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, 40)
plot_time6_dalp <- plot_timeseries(dat_time_6,
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, 40)

plot_fig5 <- plot_grid(plot_time4_dalp,plot_time5_dalp,plot_time6_dalp,
                       #arrangement data
                       ncol = 3,
                       #labels
                       labels = c("A","B", "C"),
                       label_size = 10,
                       label_fontfamily = "sans")
