# This code subsamples from all the gull wings to limit to one sample for each degree
max_sample_no = 1 # maximum amount of samples to have within one bin
bin_size      = 5 # #deg x #deg bins that will have the max amount of samples
dat_all       = read.csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability/2020_05_25_OrientedWings.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
dat_all       = subset(dat_all, species == "lar_gla" & sweep == 0 & dihedral == 0 & elbow > 80 & manus > 100)
dat_all$elbow_round = bin_size*round(dat_all$elbow/bin_size)
dat_all$manus_round = bin_size*round(dat_all$manus/bin_size)

# only keep the unique values
dat_keep <- dat_all[!allDup(dat_all[,c("elbow_round","manus_round")]),]
# save all duplicated rows so that we can randomly select one value for each bin and add back
dat_dup  <- dat_all[allDup(dat_all[,c("elbow_round","manus_round")]),]
# try to keep mainly 17_0285 if possible
curr_WingID = "17_0285"
wt_wings = c(4849,4911,6003,2195,4647,4352,3891,1380,4546)
dat_subsample = subsample_frames(dat_keep,dat_dup,curr_WingID,max_sample_no)

filename_new <- paste(format(Sys.Date(), "%Y_%m_%d"),"_dyn_subsamplewings.csv",sep = "")
write.csv(dat_subsample,filename_new)

################################## duplication function ###################################

allDup <- function (value) # found online by jholtman at gmail.com
{duplicated(value) | duplicated(value, fromLast = TRUE)}

################################## subsample function ###################################

subsample_frames <- function(dat_keep,dat_dup,curr_WingID,max_sample_no){
  
  dat_dup_uni   <- unique(dat_dup[,c("elbow_round","manus_round")])
  # Loop through each duplicated values and save one random sample
  for (j in 1:nrow(dat_dup_uni)){
    # limit to the current unique bin for the specific individual we are looking at
    tmp = subset(dat_dup, elbow_round == dat_dup_uni$elbow_round[j] &
                   manus_round == dat_dup_uni$manus_round[j] &
                   WingID == curr_WingID)
    
    # if this selection contains a wind tunnel wing then select that one
    if(any(tmp$frameID%in%wt_wings)){
      tmp = subset(tmp, frameID %in% wt_wings)
    }
    
    # if there are more frames than desired samples than choose at random
    if (nrow(tmp) > max_sample_no){
      tmp = subset(tmp, frameID %in% sample(tmp$frameID, max_sample_no))
    } 
    
    # if there are not enough frames than desired samples than keep those but add in randomly sampled frames from other individuals
    if (nrow(tmp) < max_sample_no){
      tmp_oth = subset(dat_dup, elbow_round == dat_dup_uni$elbow_round[j] &
                         manus_round == dat_dup_uni$manus_round[j] &
                         WingID != curr_WingID)
      row_rand = sample(tmp_oth[,c("frameID","WingID")], max_sample_no-nrow(tmp))
      tmp_oth  = subset(tmp_oth, frameID %in% sample(tmp_oth$frameID, max_sample_no-nrow(tmp)))
      tmp      = rbind(tmp,tmp_oth)
    }
    # if there is only one row this will be the final tmp
    
    # Save the selected rows
    if (j == 1){
      dat_subsample = rbind(dat_keep,tmp)
    }else{
      dat_subsample = rbind(dat_subsample,tmp)
    }
    
  } # end of the loop through joint angle bins
  
  return(dat_subsample)
}
