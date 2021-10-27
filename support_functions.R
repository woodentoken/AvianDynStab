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

kronecker_delta <- function(i,j){
  
  if (i == j) {
    return(1)
  }
  else {
    return(0)
  }
}

parallelaxis <- function(I, offset_vec, m, cg_a){
  
  # CAUTION: the parallel axis theorem only works if the I_CG is given about
  #          the object's centroidal axis
  
  I_new = matrix(0, nrow = 3, ncol = 3) # predefine matrix
  
  if(cg_a == "CG"){
    sign = 1
  }
  
  if(cg_a == "A"){
    sign = -1
  }
  
  for (i in 1:3){
    for (j in 1:3){
      I_new[i,j] = I[i,j] +
        sign*m*((kronecker_delta(i,j)*pracma::dot(offset_vec,offset_vec)) -
                  (offset_vec[i]*offset_vec[j]))
    }
  }
  return(I_new)
}

shift_Iorigin <- function(input_I,input_origin,input_CG,input_cg_or_a,input_m,new_origin){
  mass_properties = as.data.frame(matrix(0, nrow = 0, ncol = 7)) # overall data
  colnames(mass_properties) = c("species","BirdID","TestID","FrameID",
                                "component","object","value")
  dat    = list()
  dat$I  = matrix(0, nrow = 3, ncol = 3)
  dat$CG = matrix(0, nrow = 3, ncol = 1)
  
  #Calculate the moment of inertia about the known CG point
  if(input_cg_or_a == "A"){
    I_CG  = parallelaxis(input_I,(input_origin-input_CG),input_m,"A")
  }else{
    I_CG = input_I
  }
  
  dat$CG = new_origin - input_CG
  #Calculate the moment of inertia about the new arbitrary point
  dat$I  = parallelaxis(I_CG,-dat$CG,input_m,"CG")
  
  dat$m = input_m
  return(dat)
}

### ---------------------------------------------------------------------------------------
### ------------- Compute extremes of the CG position due to shoulder motion --------------
### ---------------------------------------------------------------------------------------

adjust_inertia <- function(sweep,dihedral,dat_in){

  ## ---------------------------- Center of Gravity ----------------------------
  shoulder_motion = dat_in
  shoulder_motion$rest_m   = (shoulder_motion$full_m-2*shoulder_motion$wing_m)
  shoulder_motion$rest_CGx = (shoulder_motion$full_m*shoulder_motion$full_CGx - 2*(shoulder_motion$wing_m*shoulder_motion$wing_CGx))/shoulder_motion$rest_m
  shoulder_motion$rest_CGz = (shoulder_motion$full_m*shoulder_motion$full_CGz - 2*(shoulder_motion$wing_m*shoulder_motion$wing_CGz))/shoulder_motion$rest_m
  
  # - Sweep adjustment - this set up ensures that a positive sweep pushes the wing backwards (i.e. more negative x)
  rot_sw = matrix(0, nrow = 3, ncol = 3)
  rot_sw[1,1] = cosd(sweep)
  rot_sw[2,2] = cosd(sweep)
  rot_sw[2,1] = sind(sweep)
  rot_sw[1,2] = -sind(sweep)
  rot_sw[3,3] = 1
  
  new_wing_CGx1 = (rot_sw[1,1]*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) + rot_sw[1,2]*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_X
  new_wing_CGy1 = (rot_sw[1,1]*(shoulder_motion$wing_CGx-shoulder_motion$pt1_X) + rot_sw[2,2]*(shoulder_motion$wing_CGy-shoulder_motion$pt1_Y)) + shoulder_motion$pt1_Y
  shoulder_motion$full_CGx_adj = (shoulder_motion$rest_m*shoulder_motion$rest_CGx + 2*shoulder_motion$wing_m*new_wing_CGx1)/shoulder_motion$full_m
  
  shoulder_motion$full_CGx_orgShoulder  = (shoulder_motion$full_CGx_adj-shoulder_motion$pt1_X)
  
  # - Dihedral adjustment 
  # this set up ensures that a positive dihedral raises the wing (i.e. more negative z) - this is like a negative rotation about x
  # caution to also account for the sweep change must read in the previous new_wing_CGy1
  rot_di = matrix(0, nrow = 3, ncol = 3)
  rot_di[1,1] = 1
  rot_di[2,2] = cosd(-dihedral)
  rot_di[3,3] = cosd(-dihedral)
  rot_di[3,2] = sind(-dihedral)
  rot_di[2,3] = -sind(-dihedral)
  
  new_wing_CGy2 = (rot_di[2,2]*(new_wing_CGy1-shoulder_motion$pt1_Y) + rot_di[2,3]*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z)) + shoulder_motion$pt1_Y
  new_wing_CGz2 = (rot_di[3,2]*(new_wing_CGy1-shoulder_motion$pt1_Y) + rot_di[3,3]*(shoulder_motion$wing_CGz-shoulder_motion$pt1_Z)) + shoulder_motion$pt1_Z
  shoulder_motion$full_CGz_adj   = (shoulder_motion$rest_m*shoulder_motion$rest_CGz + 2*shoulder_motion$wing_m*new_wing_CGz2)/shoulder_motion$full_m
  shoulder_motion$full_CGz_orgShoulder = (shoulder_motion$full_CGz_adj-shoulder_motion$pt1_Z)
  
  # Save the CG data to the main data frame 
  dat_in$full_CGx_orgShoulder_adj = shoulder_motion$full_CGx_orgShoulder 
  dat_in$full_CGz_orgShoulder_adj = shoulder_motion$full_CGz_orgShoulder 
  
  dat_in$wing_CGx_adj = new_wing_CGx1
  dat_in$wing_CGy_adj = new_wing_CGy2
  dat_in$wing_CGz_adj = new_wing_CGz2
  
  dat_in$full_CGx_adj = shoulder_motion$full_CGx_adj
  dat_in$full_CGz_adj = shoulder_motion$full_CGz_adj
  
  ## ---------------------------- Moment of Inertia ----------------------------
  
  # first need to shift origin from the VRP to the shoulder joint
  dat_in$full_Iyy_adj = NA
  
  for (i in 1:nrow(dat_in)){

    wing    = list()
    wing$I  = matrix(0, nrow = 3, ncol = 3)
    wing$CG = matrix(0, nrow = 3, ncol = 1)

    # Pull out the moment of inertia of the wing alone
    # Origin: VRP
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

    # Define the rest of the body moment of inertia
    # Origin: VRP
    
    rest    = list()
    rest$I  = matrix(0, nrow = 3, ncol = 3)
    rest$CG = matrix(0, nrow = 3, ncol = 1)
    
    rest$I[1,1] = dat_in$full_VRP_Ixx[i] - 2*dat_in$wing_Ixx[i]
    rest$I[2,2] = dat_in$full_VRP_Iyy[i] - 2*dat_in$wing_Iyy[i]
    rest$I[3,3] = dat_in$full_VRP_Izz[i] - 2*dat_in$wing_Izz[i]
    
    rest$I[1,2] = dat_in$full_VRP_Ixy[i] - 2*dat_in$wing_Ixy[i]
    rest$I[2,1] = dat_in$full_VRP_Ixy[i] - 2*dat_in$wing_Ixy[i]
    rest$I[2,3] = dat_in$full_VRP_Iyz[i] - 2*dat_in$wing_Iyz[i]
    rest$I[3,2] = dat_in$full_VRP_Iyz[i] - 2*dat_in$wing_Iyz[i]
    rest$I[1,3] = dat_in$full_VRP_Ixz[i] - 2*dat_in$wing_Ixz[i]
    rest$I[3,1] = dat_in$full_VRP_Ixz[i] - 2*dat_in$wing_Ixz[i]
    
    rest$m = dat_in$full_m[i]-2*dat_in$wing_m[i]
    
    rest$CG[1] = (dat_in$full_m[i]*dat_in$full_CGx[i] - 2*(dat_in$wing_m[i]*dat_in$wing_CGx[i]))/rest$m
    rest$CG[2] = 0
    rest$CG[3] = (dat_in$full_m[i]*dat_in$full_CGz[i] - 2*(dat_in$wing_m[i]*dat_in$wing_CGz[i]))/rest$m
    
    # Step 1: Shift the moment of inertia origin to the humeral head
    
    hum_head = c(dat_in$pt1_X[i],dat_in$pt1_Y[i],dat_in$pt1_Z[i])
    
    wing_rot = shift_Iorigin(wing$I,c(0,0,0),wing$CG,"A",wing$m,hum_head) # Output Origin: Humeral head
    
    # Step 2: Rotate the wing moment of inertia about the shoulder joint
    # Origin: Humeral head 
    wing_rot$I_new  =  t(rot_di)%*%(t(rot_sw)%*%wing_rot$I%*%rot_sw)%*%rot_di
    wing_rot$CG_new = c(dat_in$wing_CGx_adj[i], dat_in$wing_CGy_adj[i], dat_in$wing_CGz_adj[i])

    # Step 3: Shift the wing moment of inertia to be about the new CG and then the VRP
    wing_rot_VRP = shift_Iorigin(wing_rot$I_new,hum_head,wing_rot$CG_new,"A",wing_rot$m,c(0,0,0)) # Origin: Output Full Bird VRP
    
    # Step 5: Add the new wing parameters back to the rest
    # Origin: VRP
    full = list()
    full$I  = matrix(0, nrow = 3, ncol = 3)
    full$CG = matrix(0, nrow = 3, ncol = 1)
    
    full$I[1,1] = rest$I[1,1] + 2*wing_rot_VRP$I[1,1]
    full$I[2,2] = rest$I[2,2] + 2*wing_rot_VRP$I[2,2]
    full$I[3,3] = rest$I[3,3] + 2*wing_rot_VRP$I[3,3]
    
    full$I[1,2] = rest$I[1,2] + 2*wing_rot_VRP$I[1,2]
    full$I[2,1] = rest$I[2,1] + 2*wing_rot_VRP$I[2,1]
    full$I[2,3] = rest$I[2,3] + 2*wing_rot_VRP$I[2,3]
    full$I[3,2] = rest$I[3,2] + 2*wing_rot_VRP$I[3,2]
    full$I[1,3] = rest$I[1,3] + 2*wing_rot_VRP$I[1,3]
    full$I[3,1] = rest$I[3,1] + 2*wing_rot_VRP$I[3,1]
    
    # the CG y position doesn't change with this symmetric morphing
    full$CG = c(dat_in$full_CGx_adj[i], dat_in$full_CGy[i], dat_in$full_CGz_adj[i]) 
    
    full$m = dat_in$full_m[i]
    
    # Step 6: Shift the full wing moment of inertia to be about the new CG
    full_CG = shift_Iorigin(full$I, c(0,0,0), full$CG, "A", full$m, full$CG) # Origin: Full bird CG
    
    dat_in$full_Iyy_adj[i] = full_CG$I[2,2]
  }
  return(dat_in)
}

