library(ggplot2)
library(gridExtra) # for using grid arrange
library(cowplot)   # need for plot_grid()
library(ggthemes)  # need for geom_rangeframe
library(ggalt)     # need for encircle
library(ggstream)  # stream graph
library(dplyr)
library(tidyr)
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
lab_sweep      = "Sweep angle (°)"
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

cc_sweep <- c("0" = "#BFBFFF",
              "-5" = "#9999FF",
              "-10" = "#6066C6",
              "-15" = "#23388F",
              "-20" = "#000F5C")

## ------------- Plot the key trim parameters -------------------
plot_trim <- ggplot() + 
  #add data
  geom_point(data = subset(dat_cut, eignum == 1), 
             aes(x = U0, y = theta_0, col = manus)) + 
  #theme control
  th +
  scale_fill_gradientn(colours = cc_full, name = lab_manus) + 
  scale_color_gradientn(colours = cc_full, name = lab_manus) + 
  scale_shape_manual(values = c(21,22,23,24,25), name = lab_manus) + 
  # axis control 
  scale_x_continuous(limits = c(10,35), breaks = c(10,15,20,25,30), name = "Trim flight speed (m/s)") +
  scale_y_continuous(limits = c(-53,0), breaks = c(-50,-25,0), name = "Trim pitch angle (°)") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 10, xend = 30, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -50, yend = 0)

# Prep data for the stream plots
test <- as_tibble(subset(dat_cut, eignum == 1 & dihedral == 20))
test$U0_r <- round(test$U0, digits = 0)
test$theta_0_r <- round(test$theta_0, digits = 0)

dat_U0_plot <- test %>%
  group_by(U0_r, sweep, .drop=FALSE) %>%
  summarise(count_man=length(unique(manus,elbow)))
dat_U0_plot$sweep <- as.factor(dat_U0_plot$sweep)
dat_U0_plot$count_man <- as.double(dat_U0_plot$count_man)

dat_t0_plot <- test %>%
  group_by(theta_0_r, sweep, .drop=FALSE) %>%
  summarise(count_man=length(unique(manus,elbow)))
dat_t0_plot$sweep <- as.factor(dat_t0_plot$sweep)
dat_t0_plot$count_man <- as.double(dat_t0_plot$count_man)

# ------ Speed stream -----
plot_U0_stream <- ggplot()+
  geom_stream(data = dat_U0_plot, aes(U0_r,count_man,fill = sweep))  + 
  #colour control 
  scale_fill_manual(values = cc_sweep, name = lab_sweep) + 
  #theme
  th +
  # axis control 
  scale_x_continuous(limits = c(10,35), breaks = c(10,15,20,25,30), name = "Trim flight speed (m/s)") +
  scale_y_continuous(limits = c(-100,100), breaks = c(-30,-15,0,15,30), name = "Number of trimmed configurations") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 10, xend = 30, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -30, yend = 30)

# ------ Pitch angle stream -----
plot_t0_stream <- ggplot()+
  # add data
  geom_stream(data = dat_t0_plot, aes(theta_0_r, count_man,fill = sweep)) + 
  #colour control 
  scale_fill_manual(values = cc_sweep, name = lab_sweep) + 
  #theme
  th +
  # axis control 
  coord_flip()+
  scale_x_continuous(limits = c(-53,0), breaks = c(-50,-25,0), name = "Trim pitch angle (°)") +
  scale_y_continuous(limits = c(-100,100), breaks = c(-30,-15,0,15,30), name = "Number of trimmed configurations") +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = -50, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = -30, yend = 30)

plot_x <- ggplot() + 
  #add data
  geom_point(data = subset(dat_cut,eignum == 1 & sweep == -15), 
             aes(x = xcg, y = elbow, col = manus), pch = 10, alpha = 0.9) +
  geom_point(data = subset(dat_cut,eignum == 1 & sweep == -15), 
            aes(x = xac, y = elbow, col = manus), pch = 15, alpha = 0.9) +
  scale_color_gradientn(colours = cc_full, name = lab_manus) +
  #theme
  th +
  # axis control 
  scale_x_continuous(limits = c(0,0.1), breaks = seq(0,0.1,0.02), name = "-x (m)") +
  scale_y_continuous(limits = c(80,150), breaks = seq(80,150,10), name = lab_elbow) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 0.1, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 80, yend = 150)

plot_em <- ggplot() + 
  #add data
  geom_encircle(data = subset(dat_cut,eignum == 1 & cmcl < 0), 
                aes(x = elbow, y = manus, fill = as.factor(sweep), col = as.factor(sweep)), alpha = 0.25, expand = -0.001) + 
  geom_encircle(data = subset(dat_cut,eignum == 1 & cmcl > 0), 
                aes(x = elbow, y = manus, col = as.factor(sweep)), expand = -0.001) + 
  geom_jitter(data = subset(dat_cut,eignum == 1), 
             aes(x = elbow, y = manus, col = as.factor(sweep), pch = as.factor(sign(cmcl))), alpha = 0.9,
             width = 0.5, height = 0.5) +
  #theme
  th +
  scale_fill_manual(values = cc_sweep, name = lab_sweep) + 
  scale_colour_manual(values = cc_sweep, name = lab_sweep) + 
  scale_shape_manual(values = c(15,0)) +
  # axis control 
  scale_x_continuous(limits = c(80,180), breaks = seq(80,170,10), name = lab_elbow) +
  scale_y_continuous(limits = c(100,180), breaks = seq(100,180,10), name = lab_manus) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 80, xend = 170, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 100, yend = 180)


plot_topleft <- plot_grid(plot_trim,plot_t0_stream,
                     #arrangement data
                     ncol = 2,
                     rel_widths = c(1,1),
                     #labels
                     labels = c("A","B"),
                     label_size = 10,
                     label_fontfamily = "sans")
plot_botleft <- plot_grid(plot_U0_stream,blank_plot,
                          #arrangement data
                          ncol = 2,
                          rel_widths = c(1,1),
                          #labels
                          labels = c("C",""),
                          label_size = 10,
                          label_fontfamily = "sans")
plot_left <- plot_grid(plot_topleft,plot_botleft,
                          #arrangement data
                          nrow = 2,
                          rel_widths = c(1,1),
                          #labels
                          labels = c("",""),
                          label_size = 10,
                          label_fontfamily = "sans")
plot_right <- plot_grid(plot_x,plot_em, blank_plot,blank_plot,
                       #arrangement data
                       nrow = 2,
                       rel_widths = c(1,1),
                       #labels
                       labels = c("D",""),
                       label_size = 10,
                       label_fontfamily = "sans")

fig2_full <- plot_grid(plot_left,plot_right,
                       #arrangement data
                       ncol = 2,
                       rel_widths = c(2,2),
                       #labels
                       labels = c("",""),
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
  annotate(geom = "rect", xmin = -1.25, xmax = 1.25, ymin = -1.25, ymax = 1.25, fill = NA, col = "black") +
  scale_x_continuous(name = "Real") +
  scale_y_continuous(name = "Imaginary") + 
  # add data
  geom_point(data = subset(dat_cut, eignum < 3 & (sweep == -15 | sweep == -5)), 
             aes(x = eig_real, y = eig_imag, col = manus, pch = as.factor(sweep)), alpha =0.9)

# Define the axis ticks and associated labels
axis_begin  = -1.25
axis_end    = 1.25
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
  annotate(geom = "segment", x = axis_begin, xend = axis_end, y = 0, yend = 0) +
  annotate(geom = "segment", x = 0, xend = 0, y = axis_begin, yend = axis_end) +
  scale_x_continuous(name = "Real") +
  scale_y_continuous(name = "Imaginary") + 
  # add data
  geom_point(data = subset(dat_cut,eignum > 2 & (sweep == -15 | sweep == -5)), 
             aes(x = eig_real, y = eig_imag, col = manus, pch = as.factor(sweep)), alpha =0.9) 


## ------------------- Plot damping and frequency results ----------------------

# short period mode - Frequency
plot_sp_o <- ggplot() + 
  geom_point(data = subset(dat_cut,eignum == 1 & (sweep == -15 | sweep == -5)), 
             aes(x = elbow, y = omega_n, col = manus, pch = as.factor(sweep))) + 
  th +
  # colour control
  scale_color_gradientn(colours = cc_full, name = lab_manus) + 
  # axis control 
  scale_x_continuous(limits = c(80,170), breaks = seq(80,170,10), name = lab_elbow) +
  scale_y_continuous(limits = c(0,40), name = expression(paste(omega)["n"])) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 80, xend = 170, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 40) + 
  # use this segment for the "zoom in"
  annotate(geom = "segment", x = -log(0), xend = -log(0), y = 0, yend = 2)

# Extract the legend to plot later
leg_elbow <- g_legend(plot_sp_o)
plot_sp_o <- plot_sp_o + theme(legend.position = 'none')


# short period mode - Damping
plot_sp_z <- ggplot() + 
  geom_point(data = subset(dat_cut,eignum == 1 & (sweep == -15 | sweep == -5)), 
             aes(x = elbow, y = zeta, col = manus, pch = as.factor(sweep))) + 
  geom_hline(yintercept = 0.35) + 
  geom_hline(yintercept = 0.25) + 
  # theme control
  th +
  theme(legend.position = 'none') +
  # colour control
  scale_color_gradientn(colours = cc_full) + 
  # axis control 
  scale_x_continuous(limits = c(80,170), breaks = seq(80,170,10), name = lab_elbow) +
  scale_y_continuous(limits = c(0,1), name = expression(paste(zeta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 80, xend = 170, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 1)


# phugoid mode - Frequency
plot_ph_o <- ggplot() + 
  geom_point(data = subset(dat_cut,eignum == 3 & (sweep == -15 | sweep == -5)), 
             aes(x = elbow, y = omega_n, col = manus, pch = as.factor(sweep))) + 
  # theme control
  th +
  theme(legend.position = 'none') +
  # colour control
  scale_color_gradientn(colours = cc_full) + 
  # axis control 
  scale_x_continuous(limits = c(80,170), breaks = seq(80,170,10), name = lab_elbow) +
  scale_y_continuous(limits = c(0,2), name = expression(paste(omega)["n"])) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 80, xend = 170, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 2)

# phugoid mode - Damping
plot_ph_z <- ggplot() + 
  geom_point(data = subset(dat_cut,eignum == 3 & (sweep == -15 | sweep == -5)), 
             aes(x = elbow, y = zeta, col = manus, pch = as.factor(sweep))) + 
  geom_hline(yintercept = 0.04) + 
  # theme control
  th +
  theme(legend.position = 'none') +
  # colour control
  scale_color_gradientn(colours = cc_full) + 
  # axis control 
  scale_x_continuous(limits = c(80,170), breaks = seq(80,170,10), name = lab_elbow) +
  scale_y_continuous(limits = c(0,1), name = expression(paste(zeta))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 80, xend = 170, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 1)


fig3_full <- plot_grid(locus_sp, locus_ph, plot_sp_o,plot_ph_o,plot_sp_z,plot_ph_z,
                    #arrangement data
                    ncol = 2,
                    rel_widths = c(1,1),
                    #labels
                    labels = c("A","B","C","D","E","F"),
                    label_size = 10,
                    label_fontfamily = "sans")


plot_sp <- ggplot() + 
  geom_rect(aes(ymax = 3.6, ymin = 0.28, xmin = 0.35, xmax = 1), alpha = 0.4) + 
  geom_rect(aes(ymax = 10, ymin = 0.16, xmin = 0.25, xmax = 1), alpha = 0.3) + 
  geom_point(data = subset(dat_cut,eignum == 1 & (sweep == -15 | sweep == -5)), 
             aes(x = zeta, y = w_sp_na, col = manus, pch = as.factor(sweep))) + 
  # theme control
  th +
  theme(legend.position = 'none') +
  # colour control
  scale_color_gradientn(colours = cc_full) + 
  # axis control 
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), name = expression(paste(zeta))) +
  scale_y_continuous(limits = c(0,20), name = expression(paste(omega["n"]^2,"/n"[alpha]))) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 0, xend = 1, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 0, yend = 20)

## ------------- Plot the phase and magnitude -------------------

  
plot_sp_magphase <- ggplot() + 
  geom_vline(xintercept = seq(0, 0.25, by = 0.05), colour = "grey70", size = 0.2) +
  geom_hline(yintercept = seq(0, 315, by = 45), colour = "grey70", size = 0.2) +
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag1, y = phase1*180/pi, yend = phase1*180/pi, group = del_z), col = col_u, alpha = 0.3) + 
  geom_point(data = subset(dat_sp, (sweep == -15 | sweep == -5)), 
             aes(x = mag1, y = phase1*180/pi, pch = as.factor(sweep)), 
             size = 2, col = col_u,  alpha = 0.5) + 
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag2, y = phase2*180/pi, yend = phase2*180/pi, group = del_z), col = col_alpha, alpha = 0.3) + 
  geom_point(data = subset(dat_sp, (sweep == -15 | sweep == -5)), 
             aes(x = mag2, y = phase2*180/pi, pch = as.factor(sweep)), 
             size = 2, col = col_alpha,  alpha = 0.5) +
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag3, y = phase3*180/pi, yend = phase3*180/pi, group = del_z), col = col_theta, alpha = 0.3) + 
  #geom_segment(data = dat_sp, aes(x = 0, xend = mag4, y = phase4*180/pi, yend = phase4*180/pi, group= del_z), col = col_q, alpha = 0.3) + 
  geom_point(data = subset(dat_sp, (sweep == -15 | sweep == -5)), 
             aes(x = mag4, y = phase4*180/pi, pch = as.factor(sweep)), 
             size = 2, col = col_theta,  alpha = 0.5) + 
  # theme control
  th +
  theme_light()+
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  panel_border(remove = TRUE) + 
  # axis control
  coord_polar(theta = "y", start = -0.5*pi, direction = -1) +
  scale_x_continuous(limits = c(0,0.25), breaks = seq(0, 0.25, by = 0.05), name = "Magnitude") +
  scale_y_continuous(limits = c(0,360), breaks = c(0,90,180,270), name = lab_phase)
  
## ----------- Phugoid mode phase and magnitude -----------------

plot_ph_magphase <- ggplot() + 
  geom_vline(xintercept = seq(0, 2.5, by = 0.5), colour = "grey70", size = 0.2) +
  geom_hline(yintercept = seq(0, 315, by = 45), colour = "grey70", size = 0.2) +
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag1, y = phase1*180/pi, yend = phase1*180/pi, group = del_z), col = col_u, alpha = 0.3) + 
  geom_point(data = subset(dat_ph, (sweep == -15 | sweep == -5)), 
             aes(x = mag1, y = phase1*180/pi,  alpha = 0.5, pch = as.factor(sweep)), 
             size = 2, col = col_u) + 
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag2, y = phase2*180/pi, yend = phase2*180/pi, group = del_z), col = col_alpha, alpha = 0.3) + 
  geom_point(data = subset(dat_ph, (sweep == -15 | sweep == -5)), 
             aes(x = mag2, y = phase2*180/pi,  alpha = 0.5, pch = as.factor(sweep)), 
             size = 2, col = col_alpha) +
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag3, y = phase3*180/pi, yend = phase3*180/pi, group = del_z), col = col_theta, alpha = 0.3) + 
  #geom_segment(data = dat_ph, aes(x = 0, xend = mag4, y = phase4*180/pi, yend = phase4*180/pi, group= del_z), col = col_q, alpha = 0.3) + 
  geom_point(data = subset(dat_ph, (sweep == -15 | sweep == -5)), 
             aes(x = mag4, y = phase4*180/pi,  alpha = 0.5, pch = as.factor(sweep)), 
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
  scale_x_continuous(limits = c(0,2.5), breaks = seq(0, 2.5, by = 0.5), name = "Magnitude") +
  scale_y_continuous(limits = c(0,360), breaks = c(0,90,180,270), name = lab_phase)


## ------------------- Phase with respect to elbow and wrist angle ------------
plot_ph_phase <- ggplot() + 
  geom_raster(data = subset(dat_ph, sweep == -15), 
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
  geom_raster(data = subset(dat_sp, sweep == -15), 
              aes(x = elbow, y = manus, alpha = round(phase4,1)*180/pi), fill = col_theta) + 
  # theme control
  th +  scale_alpha_continuous(name = lab_phase) +
  # axis control 
  coord_fixed() + 
  scale_x_continuous(limits = c(85,156), breaks = c(90,120,150), name = lab_elbow) +
  scale_y_continuous(limits = c(100,161), breaks = c(100,130,160), name = lab_manus) +
  geom_rangeframe() +
  annotate(geom = "segment", x = 90, xend = 150, y = log(0), yend = log(0)) +
  annotate(geom = "segment", x = log(0), xend = log(0), y = 100, yend = 160) 


fig4_full <- plot_grid(plot_sp_magphase,plot_ph_magphase,
          plot_sp_phase,plot_ph_phase,
          #arrangement data
          ncol = 2,
          #labels
          labels = c("A","B", "C","D"),
          label_size = 10,
          label_fontfamily = "sans")



## ------------------ Time response to an initial alpha -------------------

lim_u     = c(-0.02,0.02)
lim_alpha = c(-2,2)
lim_q     = c(-50,50)
break_q = seq(-50,50, by = 25)
lim_theta = c(-2.5,2.5)
break_theta = seq(-2.5,2.5,by=1.25)
lim_time = 10
plot_dalp <- plot_timeseries(dat_time_1[1:which(dat_time_1$t==lim_time),],
                                   dat_time_2[1:which(dat_time_1$t==lim_time),],
                                   dat_time_3[1:which(dat_time_1$t==lim_time),],
                                   dat_time_4[1:which(dat_time_1$t==lim_time),],
                                   dat_time_5[1:which(dat_time_1$t==lim_time),],
                                   dat_time_6[1:which(dat_time_1$t==lim_time),],
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, lim_time)

## ------------------ Time response to a ramped speed -------------------

lim_u     = c(-0.06,0)
lim_alpha = c(-0.1,0.1)
lim_q     = c(-1,1)
break_q   = seq(-1,1, by= 0.5)
lim_theta = c(-2.5,2.5)
break_theta = seq(-2.5,2.5, by= 1.25)
lim_time = 25
plot_uramp <- plot_timeseries(dat_time_7[1:which(dat_time_1$t==lim_time),],
                                   dat_time_8[1:which(dat_time_1$t==lim_time),],
                                   dat_time_9[1:which(dat_time_1$t==lim_time),],
                                   dat_time_10[1:which(dat_time_1$t==lim_time),],
                                   dat_time_11[1:which(dat_time_1$t==lim_time),],
                                   dat_time_12[1:which(dat_time_1$t==lim_time),],
                                   col_u,col_alpha,col_q,col_theta,
                                   lim_u,lim_alpha,lim_q,lim_theta, 
                                   break_q, break_theta, lim_time)

fig5_full <- plot_grid(plot_dalp,plot_uramp,
                       #arrangement data
                       ncol = 2,
                       #labels
                       labels = c("A","B"),
                       label_size = 10,
                       label_fontfamily = "sans")

