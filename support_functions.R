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
plot_timeseries <- function(dat_time,col_u,col_alpha,col_q,col_theta,lim_u,lim_alpha,lim_q,lim_theta,break_q, break_theta,lim_t){
  
  # ------------ u time series ----------
  plot_del_u     = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_u), col = col_u) + 
    th + 
    # axis control
    scale_x_continuous(limits = c(0,lim_t), name = "Time (s)") + 
    scale_y_continuous(limits = lim_u, name = "u (m/s)") +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = lim_t, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = lim_u[1], yend = lim_u[2])
  
  # ------------ alpha time series ----------
  plot_del_alp   = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_alp*180/pi), col = col_alpha) + 
    th + 
    # axis control
    scale_x_continuous(limits = c(0,lim_t), name = "Time (s)") + 
    scale_y_continuous(limits = lim_alpha, name = expression(paste(Delta,alpha," (°)"))) +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = lim_t, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = lim_alpha[1], yend = lim_alpha[2])
  
  # ------------ q time series ----------
  plot_del_q     = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_q*180/pi), col = col_q) + 
    th + 
    # axis control    
    scale_x_continuous(limits = c(0,lim_t), name = "Time (s)") + 
    scale_y_continuous(limits = lim_q, breaks = break_q, name = expression(paste(Delta,"q (°/s)"))) +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = lim_t, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = break_q[1], yend = break_q[5])
  
  # ------------ theta time series ----------
  plot_del_theta = ggplot() + 
    geom_line(data = dat_time, aes(x = t, y = del_theta*180/pi), col = col_theta) + 
    th + 
    # axis control
    scale_x_continuous(limits = c(0,lim_t), name = "Time (s)") + 
    scale_y_continuous(limits = lim_theta, breaks = break_theta, name = expression(paste(Delta,theta," (°)"))) +
    geom_rangeframe() +
    annotate(geom = "segment", x = 0, xend = lim_t, y = log(0), yend = log(0)) +
    annotate(geom = "segment", x = log(0), xend = log(0), y = lim_theta[1], yend = lim_theta[2])
  
  plot_out <- plot_grid(plot_del_u,plot_del_alp,plot_del_q,plot_del_theta,
                        #arrangement data
                        ncol = 1,
                        #labels
                        label_size = 10,
                        label_fontfamily = "sans")
}


#create function to retrieve the legend as an object
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

shift_Iorigin <- function(input_I,input_origin,input_CG,input_cg_or_a,input_m,new_origin,name,dat_info){
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  colnames(mass_properties) = c("species","BirdID","TestID","FrameID",
                                "component","object","value")
  dat    = list()
  dat$I  = matrix(0, nrow = 3, ncol = 3)
  dat$CG = matrix(0, nrow = 3, ncol = 1)
  
  #I defined about the wing CG
  if(input_cg_or_a == "A"){
    I_CG  = parallelaxis(input_I,(input_origin-input_CG),input_m,"A")
  }else{
    I_CG = input_I
  }
  
  dat$CG = new_origin - input_CG
  
  dat$I  = parallelaxis(I_CG,-dat$CG,input_m,"CG")
  
  dat$m = input_m
  return(dat)
}

### ---------------------------------------------------------------------------------------
### ------------- Compute extremes of the CG position due to shoulder motion --------------
### ---------------------------------------------------------------------------------------

adjust_inertia <- function(sweep,dihedral,dat_in){

  ## ---------------------------- Center of Gravity ----------------------------
  shoulder_motion = dat_inertial
  shoulder_motion$rest_m   = (shoulder_motion$full_m-2*shoulder_motion$wing_m)
  shoulder_motion$rest_CGx = (shoulder_motion$full_m*shoulder_motion$full_CGx - 2*(shoulder_motion$wing_m*shoulder_motion$wing_CGx))/shoulder_motion$rest_m
  shoulder_motion$rest_CGz = (shoulder_motion$full_m*shoulder_motion$full_CGz - 2*(shoulder_motion$wing_m*shoulder_motion$wing_CGz))/shoulder_motion$rest_m
  
  # - Sweep adjustment - this set up ensures that a positive sweep pushes the wing backwards (i.e. more negative x)
  new_wing_CGx1 = (cosd(sweep)*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) - sind(sweep)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_X
  new_wing_CGy1 = (sind(sweep)*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) + cosd(sweep)*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Y
  shoulder_motion$full_CGx = (shoulder_motion$rest_m*shoulder_motion$rest_CGx + 2*shoulder_motion$wing_m*new_wing_CGx1)/shoulder_motion$full_m
  
  shoulder_motion$full_CGx_orgShoulder  = (shoulder_motion$full_CGx-shoulder_motion$pt1_X)
  
  # - Dihedral adjustment - this set up ensures that a positive dihedral raises the wing (i.e. more negative z)
  # caution to also account for the sweep change must read in the previous new_wing_CGy1
  new_wing_CGz2 = (cosd(dihedral)*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z) - sind(dihedral)*(new_wing_CGy1-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Z
  new_wing_CGy2 = (sind(dihedral)*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z) + cosd(dihedral)*(new_wing_CGy1-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Y
  shoulder_motion$full_CGz   = (shoulder_motion$rest_m*shoulder_motion$rest_CGz + 2*shoulder_motion$wing_m*new_wing_CGz2)/shoulder_motion$full_m
  shoulder_motion$full_CGz_orgShoulder = (shoulder_motion$full_CGz-shoulder_motion$pt1_Z)
  
  dat_in$full_CGx_orgShoulder_adj = shoulder_motion$full_CGx_orgShoulder 
  dat_in$full_CGz_orgShoulder_adj = shoulder_motion$full_CGz_orgShoulder 
  
  dat_in$wing_CGx_adj = new_wing_CGx1
  dat_in$wing_CGy_adj = new_wing_CGy2
  dat_in$wing_CGz_adj = new_wing_CGz2
  
  ## ---------------------------- Moment of Inertia ----------------------------
  
  # first need to shift origin from the VRP to the shoulder joint

  for (i in 1:nrow(dat_in)){
    
    wing    = list()
    wing$I  = matrix(0, nrow = 3, ncol = 3)
    wing$CG = matrix(0, nrow = 3, ncol = 1)
    
    #current I is about the VRP
    wing$I[1,1] = dat_in$wing_Ixx[i]
    wing$I[2,2] = dat_in$wing_Iyy[i]
    wing$I[3,3] = dat_in$wing_Izz[i]
    
    wing$I[1,2] = dat_in$wing_Ixy[i]
    wing$I[2,1] = dat_in$wing_Ixy[i]
    wing$I[2,3] = dat_in$wing_Iyz[i]
    wing$I[3,2] = dat_in$wing_Iyz[i]
    wing$I[1,3] = dat_in$wing_Ixz[i]
    wing$I[3,1] = dat_in$wing_Ixz[i]
    
    wing$CG[1] = dat_in$wing_CGx[i]
    wing$CG[2] = dat_in$wing_CGy[i]
    wing$CG[3] = dat_in$wing_CGz[i]
    
    wing$m = dat_in$wing_m[i]
    
    hum_head = c(dat_in$pt1_X[i],
                 dat_in$pt1_Y[i],
                 dat_in$pt1_Z[i])
    
    org_shoulder = shift_Iorigin(wing$I,c(0,0,0),wing$CG,"A",wing$m,hum_head,"wing_hum",dat_final[i,]) # Origin: Humeral head
    
    # Define the rest of the body moment of inertia 
    # Origin: VRP
    
    rest$I[1,1] = dat_in$full_VRP_Ixx[i] - 2*dat_in$wing_Ixx[i]
    rest$I[2,2] = dat_in$full_VRP_Iyy[i] - 2*dat_in$wing_Iyy[i]
    rest$I[3,3] = dat_in$full_VRP_Izz[i] - 2*dat_in$wing_Izz[i]
    
    rest$I[1,2] = dat_in$full_VRP_Ixy[i] - 2*dat_in$wing_Ixy[i]
    rest$I[2,1] = dat_in$full_VRP_Ixy[i] - 2*dat_in$wing_Ixy[i]
    rest$I[2,3] = dat_in$full_VRP_Iyz[i] - 2*dat_in$wing_Iyz[i]
    rest$I[3,2] = dat_in$full_VRP_Iyz[i] - 2*dat_in$wing_Iyz[i]
    rest$I[1,3] = dat_in$full_VRP_Ixz[i] - 2*dat_in$wing_Ixz[i]
    rest$I[3,1] = dat_in$full_VRP_Ixz[i] - 2*dat_in$wing_Ixz[i]
    
    rest$m = dat_in$rest_m[i]
    
    rest$CG[1] = (dat_in$full_m[i]*dat_in$full_CGx[i] - 2*(dat_in$wing_m[i]*dat_in$wing_CGx[i]))/rest$m
    rest$CG[2] = 0
    rest$CG[3] = (dat_in$full_m[i]*dat_in$full_CGz[i] - 2*(dat_in$wing_m[i]*dat_in$wing_CGz[i]))/rest$m
    
    # SET THE SHOULDER JOINT AS THE ORIGIN FOR THE WING ONLY MOMENT OF INERTIA
    
    # ROTATE THE WING MOMENT OF INERTIA FOLLOWING SWEEP AND THEN DIHEDRAL
    
    # ADJUST THE WING MOMENT OF INERTIA TO BE ABOUT IT'S NEW CG
    ##WRONG
    input_I,input_origin,input_CG,input_cg_or_a,input_m,new_origin,name,dat_info
    wing2 = shift_Iorigin(wing_rot$I,hum_head,wing_rot$CG,"A",wing$m,hum_head,"wing_hum",dat_final[i,]) # Origin: Humeral head
    
    # ADJUST THE WING MOMENT OF INERTIA TO BE ABOUT THE VRP
    org_shoulder = shift_Iorigin(wing$I,c(0,0,0),wing$CG,"CG",wing$m,hum_head,"wing_hum",dat_final[i,]) # Origin: Humeral head
    # ADD THE WING MOMENT OF INERTIA TO THE REST MOMENT OF INERTIA
    
    # SHIFT THE FULL MOMENT OF INERTIA TO BE ABOUT THE FINAL CG
    org_shoulder = shift_Iorigin(wing$I,c(0,0,0),wing$CG,"A",full$m,hum_head,"wing_hum",dat_final[i,]) # Origin: Humeral head
  }
  
  
  return(dat_in)
}

