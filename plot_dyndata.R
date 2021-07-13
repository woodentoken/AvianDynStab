
dat_out <- read.csv('LongDynStability_Rigid.csv')

plot(dat_out$zeta[which(dat_out$alpha == 0)],dat_out$omega_n[which(dat_out$alpha == 0)])
