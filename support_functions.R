## --------- This script maintains all the supporting functions needed for analysis

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

## --------------------- Return the time series plots
plot_timeseries <- function(dat_time,col_u,col_alpha,col_q,col_theta){
  
  # ------------ u time series ----------
  plot_del_u     = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_u), col = col_u) + 
    th + 
    # axis control
    scale_y_continuous(limits = c(-0.06,0.06), name = "u (m/s)") +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = 60, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = -0.06, yend = 0.06)
  
  # ------------ alpha time series ----------
  plot_del_alp   = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_alp*180/pi), col = col_alpha) + 
    th + 
    # axis control
    scale_y_continuous(limits = c(-5,5), name = expression(paste(Delta,alpha," (°)"))) +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = 60, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = -5, yend = 5)
  
  # ------------ q time series ----------
  plot_del_q     = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_q*180/pi), col = col_q) + 
    th + 
    # axis control
    scale_y_continuous(limits = c(-50,50), name = expression(paste(Delta,"q (°/s)"))) +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = 60, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = -50, yend = 50)
  
  # ------------ theta time series ----------
  plot_del_theta = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_theta*180/pi), col = col_theta) + 
    th + 
    # axis control
    scale_y_continuous(limits = c(-10,10), name = expression(paste(Delta,theta," (°)"))) +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = 60, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = -10, yend = 10)
  
  plot_out <- plot_grid(plot_del_u,plot_del_alp,plot_del_q,plot_del_theta,
                        #arrangement data
                        ncol = 1,
                        #labels
                        labels = c("A","B","C","D"),
                        label_size = 10,
                        label_fontfamily = "sans")
}
