library(ggplot2)
library(alphahull) # need for the convex hull
library(ptinpoly)  # need for determining points outside the convex hull

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

dat_all <- read.csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2020_05_25_OrientedWings.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
dat_all <- subset(dat_all, species == "lar_gla" & sweep == 0 & dihedral == 0)
dat_all$root_c = sqrt((dat_all$Pt12X - dat_all$Pt11X)^2 + (dat_all$Pt12Z - dat_all$Pt11Z)^2)
dat_all$FrameID <- paste("F", dat_all$frameID, sep = "")



ggplot() + geom_point(data = subset(dat_cut,alpha==6), aes(x = eig_real, y = eig_imag, col = elbow))

ggplot() + geom_point(data = subset(dat_cut, alpha==6 & eignum == 1), aes(x = elbow, y = manus, col = zeta)) + th

ggplot() + geom_point(data = dat_out, aes(x = U0, y = theta0*180/pi, col = manus))

ggplot() + geom_point(data = subset(dat_cut,alpha==6), aes(x = elbow, y = manus, col = del_x))
ggplot() + geom_point(data = dat_num, aes(x = CL_adj, y = Cm_adj, col = elbow))
ggplot() + geom_point(data = dat_stab_adj, aes(x = elbow, y = manus, col = cmcl))

dat_time <- read.csv('./outputdata/2021_07_27_elbow160_manus160_alpha6.csv', header = FALSE)
colnames(dat_time) <- c("t","del_u","del_alp","del_q","del_theta")
ggplot() + 
  geom_line(data = dat_time, aes(x = t, y = del_u), col = "red") +
  geom_line(data = dat_time, aes(x = t, y = del_alp), col = "blue") +
  geom_line(data = dat_time, aes(x = t, y = del_q), col = "green") +
  geom_line(data = dat_time, aes(x = t, y = del_theta), col = "purple") + th + scale_y_continuous(limits = c(-0.05,0.05)) 
  

plot(dat_time[,1],dat_time[,2])

## ---------------- Remove points in the predicted that are outside of the convex hull ---------------- 
cut_trueshape <- function(dat,dat_geom,col_elbow,col_manus){
  # fit the convex hull with an alpha factor
  alphashape <- ahull(dat_geom, alpha = 30)
  # save all the given vertices
  vertices <- as.data.frame(dat_geom[alphashape$ashape.obj$alpha.extremes,])
  # Need to order the points appropriately
  # calculate the mean value
  centerpt <- c(mean(vertices[,1]),mean(vertices[,2]))
  # calculate the angle from the mean of each point
  vertices$angle <- atan2(vertices[,2]-centerpt[2],vertices[,1]-centerpt[1])
  # sort by the angle
  vertices <- vertices[order(vertices$angle),]
  
  # cut to be within the polygon
  filtele   <- pip2d(as.matrix(vertices[,c(1:2)]),as.matrix(dat[,c(col_elbow,col_manus)])) 
  dat_cut   <- dat[filtele==1,]  # filtele==1 retains only the points that are inside the hull
  dat_return <- list()
  dat_return$dat_cut  <- dat_cut
  dat_return$vertices <- vertices
  return(dat_return)
}
