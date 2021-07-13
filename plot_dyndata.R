
dat_out <- read.csv('LongDynStability_Rigid.csv')

plot(dat_out$eig_real[which(dat_out$alpha == 0)],dat_out$eig_imag[which(dat_out$alpha == 0)])
