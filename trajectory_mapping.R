#trajectory sampling analysis for paper
setwd("/Users/kasturilele/Documents/SLiM")
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)

fixedmuts <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_comb_16_init_new.csv", header = T)) #use new fixed mutation set

b <- 20 #binsize
#from previous analysis
time_ext_binned <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_extinct_16_new.csv", header = T)) #use new binned dataset

#read the bins established previously
bins <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/fig5_bins_16.csv", header = T))

#code to calculate bins from time_ext_binned
bins <- data_frame(species=integer(),
                   bin=integer(),
                   min=numeric(),
                   max=numeric()) #only for mutation effect
counter <- 1


mk <- 1
for (mk in 1:b) {
  ext_sub <-  subset(time_ext_binned, mk2_bin == mk)
  nt <- nrow(ext_sub)
  if(nt > 0){
    bins[counter,1] <- 2
    bins[counter,2] <- mk
    bins[counter,3] <- min(ext_sub$mut_sp2)
    bins[counter,4] <- max(ext_sub$mut_sp2)
    counter <- counter + 1
  }
}

#proportion of populations extinct
prop_ext <- data_frame(mutation_kernel_1=numeric(),
                       mutation_kernel_2=numeric(),
                       total = integer(),
                       extinct=integer(),
                       end=integer()) #only for mutation effect

mk_1 <- 1
mk_2 <- 1
counter <- 1

for (mk_1 in 1:b) {
  for (mk_2 in 1:b){
    ext_sub <-  subset(time_ext_binned, mk1_bin == mk_1 & mk2_bin == mk_2)
    nt <- nrow(ext_sub)
    if(nt > 0){
      n_ext <- length(which(ext_sub$tick_extinct < 500000))
      n_false <- length(which(ext_sub$end == TRUE)) # only for mutation effect 
      prop_ext[counter,1] <- mk_1
      prop_ext[counter,2] <- mk_2
      prop_ext[counter,3] <- nt
      prop_ext[counter,4] <- n_ext
      prop_ext[counter,5] <- n_false
      counter <- counter + 1
    }
  }
}

frac_extinct <- prop_ext$extinct/prop_ext$total
prop_ext <- cbind(prop_ext, frac_extinct)

#-------functions--------

#write a function to calculate reweights for trajectories
# calc_pxi_qxi <- function(df, mutq_sp1, mutq_sp2, sp1_cols, sp2_cols, mutp_sp1, mutp_sp2) {
#   frac_sp1 <- mutp_sp1 / df[[mutq_sp1]]
#   frac_sp2 <- mutp_sp2 / df[[mutq_sp2]]
#   
#   sp1_terms <- sweep(df[, sp1_cols], 1, frac_sp1, "*")
#   sp2_terms <- sweep(df[, sp2_cols], 1, frac_sp2, "*")
#   
#   reweight_num <- rowSums(sp1_terms) + rowSums(sp2_terms)
#   denom <- rowSums(df[, sp1_cols]) + rowSums(df[, sp2_cols])
#   
#   reweight_num/denom
#   #rowSums(sp1_terms) + rowSums(sp2_terms)
# }
# 
# #function that does this+adjusts for supply of trait-increasing vs trait-decreasing muts?
# #3,4,7,8 increasing, 1,2,5,6 decreasing
# #write a function to calculate reweights for trajectories
# calc_pxi_qxi_new <- function(df, mutq_sp1, mutq_sp2, sp1_cols_inc,sp1_cols_dec, sp2_cols_inc,sp2_cols_dec, mutp_sp1, mutp_sp2) {
#   
#   bias_inc <- 1 #(0.25/0.05)
#   bias_dec <- 1 #(0.25/0.45)
#   
#   frac_sp1_inc <- mutp_sp1 / df[[mutq_sp1]]*bias_inc
#   frac_sp2_inc <- mutp_sp2 / df[[mutq_sp2]]*bias_inc
#   frac_sp1_dec <- mutp_sp1 / df[[mutq_sp1]]*bias_dec
#   frac_sp2_dec <- mutp_sp2 / df[[mutq_sp2]]*bias_dec
#   
#   sp1_terms_inc <- sweep(df[, sp1_cols_inc], 1, frac_sp1_inc, "*")
#   sp2_terms_inc <- sweep(df[, sp2_cols_inc], 1, frac_sp2_inc, "*")
#   sp1_terms_dec <- sweep(df[, sp1_cols_dec], 1, frac_sp1_dec, "*")
#   sp2_terms_dec <- sweep(df[, sp2_cols_dec], 1, frac_sp2_dec, "*")
#   
#   reweight_num <- rowSums(sp1_terms_inc) + rowSums(sp2_terms_inc) + rowSums(sp1_terms_dec) + rowSums(sp2_terms_dec)
#   denom <- rowSums(df[, sp1_cols_inc]) + rowSums(df[, sp2_cols_inc]) + rowSums(df[, sp1_cols_dec]) + rowSums(df[, sp2_cols_dec])
#   
#   reweight_num/denom
#   #rowSums(sp1_terms_inc) + rowSums(sp2_terms_inc) + rowSums(sp1_terms_dec) + rowSums(sp2_terms_dec)
# }
# 
# calc_pxi_qxi_prod <- function(df, mutq_sp1, mutq_sp2, sp1_cols, sp2_cols, mutp_sp1, mutp_sp2) {
#   
#   frac_sp1 <- as.matrix(mutp_sp1 / df[[mutq_sp1]])
#   frac_sp2 <- as.matrix(mutp_sp2 / df[[mutq_sp2]])
#   
#   sp1_inc <- rowSums(df[, sp1_cols])
#   sp2_inc <- rowSums(df[, sp2_cols])
#   
#   sp1_terms_inc <- sweep(frac_sp1, 1, sp1_inc, "^")
#   sp2_terms_inc <- sweep(frac_sp2, 1, sp2_inc, "^")
#   
#   sp1_terms_inc*sp2_terms_inc
#   
# }

#new function, based on updated discussion. now we're considering relative rates of mutations of each type, compared to the total mutation rate(add the rates for each type)
calc_pxi_qxi_prod_new <- function(df, mutq_sp1, mutq_sp2, sp1_cols, sp2_cols, mutp_sp1, mutp_sp2) {
  
  #to simplify, group mutations from each species together
  denom_p <- mutp_sp1 + mutp_sp2 #sum of mutation rates of each species
  denom_q <- df[[mutq_sp1]] + df[[mutq_sp2]] #sum of mutation rates of each species
  
  frac_sp1 <- as.matrix((mutp_sp1/denom_p) / (df[[mutq_sp1]]/denom_q))
  frac_sp2 <- as.matrix((mutp_sp2/denom_p) / (df[[mutq_sp2]]/denom_q))
  
  sp1_inc <- rowSums(df[, sp1_cols])
  sp2_inc <- rowSums(df[, sp2_cols])
  
  sp1_terms_inc <- sweep(frac_sp1, 1, sp1_inc, "^")
  sp2_terms_inc <- sweep(frac_sp2, 1, sp2_inc, "^")
  
  sp1_terms_inc*sp2_terms_inc
  
}

# calc_pxi_qxi_prod_bias <- function(df, mutq_sp1, mutq_sp2, sp1_cols_inc,sp1_cols_dec, sp2_cols_inc,sp2_cols_dec, mutp_sp1, mutp_sp2) {
#   
#   bp_inc <- 0.25
#   bp_dec <- 0.25
#   bq_inc <- 0.05
#   bq_dec <- 0.45
#   
#   denom_p <- mutp_sp1*bp_inc*2 + mutp_sp1*bp_dec*2 + mutp_sp2*bp_inc*2 + mutp_sp2*bp_dec*2
#   denom_q <- df[[mutq_sp1]]*bq_inc*2 + df[[mutq_sp1]]*bq_dec*2 + df[[mutq_sp2]]*bq_inc*2 + df[[mutq_sp2]]*bq_dec*2 
#   
#   frac_sp1_inc <- as.matrix((mutp_sp1*bp_inc/denom_p) / (df[[mutq_sp1]]*bq_inc/denom_q))
#   frac_sp2_inc <- as.matrix((mutp_sp2*bp_inc/denom_p) / (df[[mutq_sp2]]*bq_inc/denom_q))
#   frac_sp1_dec <- as.matrix((mutp_sp1*bp_dec/denom_p) / (df[[mutq_sp1]]*bq_dec/denom_q))
#   frac_sp2_dec <- as.matrix((mutp_sp2*bp_dec/denom_p) / (df[[mutq_sp2]]*bq_dec/denom_q))
#   
#   sp1_inc <- rowSums(df[, sp1_cols_inc])
#   sp2_inc <- rowSums(df[, sp2_cols_inc])
#   sp1_dec <- rowSums(df[, sp1_cols_dec])
#   sp2_dec <- rowSums(df[, sp2_cols_dec])
#   
#   sp1_terms_inc <- sweep(frac_sp1_inc, 1, sp1_inc, "^")
#   sp2_terms_inc <- sweep(frac_sp2_inc, 1, sp2_inc, "^")
#   sp1_terms_dec <- sweep(frac_sp1_dec, 1, sp1_dec, "^")
#   sp2_terms_dec <- sweep(frac_sp2_dec, 1, sp2_dec, "^")
#   
#   sp1_terms_inc*sp2_terms_inc*sp1_terms_dec*sp2_terms_dec
#   
# }

#------- main analysis ------------
#using prop_ext, make a list of bins where ~50% extinction.
#sp1- bins 2,3,4 and sp2- bins 14,15,16
#try new subset: sp1: 18,19,20 and sp2-18,19,20
#use smaller bin: sp1: bin 2, sp2: bin 8
#sp1: 4,5 sp2: 18,19
#make subsets of prop_ext and fixed_muts for these
xh <- 1
xl <- 1
yh <- 10
yl <- 10
mk_sub <- subset(time_ext_binned, mk1_bin >= xl & mk1_bin <= xh & mk2_bin >= yl & mk2_bin <= yh)

#make subset of fixedmuts to work with
fixedmuts_sub <- fixedmuts %>% 
  semi_join(mk_sub, by = c("rep", "mes"))

#check if reps being selected are correct 
#temp <- unique(fixedmuts_sub[,1:2])
#setdiff(mk_sub[,2:3],temp)
#temp2 <- subset(fixedmuts, rep==39 &mes==497)
#checked that the missing values are because no mutations occurred in those trajectories
#start from here for new analysis (bottom-most)

muttypes_grouped <- fixedmuts_sub |> group_by(mes,rep,mutation_type)
muttypes_count <- muttypes_grouped |> tally(sort = F)

mt_count_wide <- muttypes_count |>
  pivot_wider(names_from = mutation_type, values_from = n, values_fill = 0)

#make vector of extinction outcomes (0 or 1) for chosen trajectories
mk_sub$ext_outcome <- as.numeric(mk_sub$tick_extinct < 1000000)

#make subset of mk_sub with only the info needed
mk_subset <- mk_sub[,c(2,3,10,11,14)]
#mk_subset <- mk_sub[,c(1,2,10,11,12)] #new analysis

#join with mt_count_wide to make bigger data frame for all calculation needed
all_mutdata <- merge(mk_subset,mt_count_wide)

#come up with 20 equally spaced mutation rates between 10^-7 and 10^-6, for testing 
mqs_list <- seq(1e-7,1e-6, length.out=20)

#redo this plot with the bin means from the actual sampled data 

bins_mod <- mutate(bins, binmean = rowMeans(select(bins, c("min","max"))))

all_mutdata_ext <- subset(all_mutdata, ext_outcome > 0)
all_mutdata_cox <- subset(all_mutdata, ext_outcome < 1)

new_prob_ext <- data_frame(mutation_kernel_1=numeric(),
                           mutation_kernel_2=numeric(),
                           prob_extinct=numeric()) #only for mutation effect
counter <- 1
for(i in 1:20){
  mutp_sp1 <- bins_mod$binmean[i]
  for(j in 1:20){
    mutp_sp2 <- bins_mod$binmean[j+20]
    #pxi_qxi_ext <- calc_pxi_qxi_prod_bias(all_mutdata_ext,"mut_sp1","mut_sp2",c("7","8"),c("5","6"),c("3","4"),c("1","2"),mutp_sp1,mutp_sp2)
    #pxi_qxi_cox <- calc_pxi_qxi_prod_bias(all_mutdata_cox,"mut_sp1","mut_sp2",c("7","8"),c("5","6"),c("3","4"),c("1","2"),mutp_sp1,mutp_sp2)
    pxi_qxi_ext <- calc_pxi_qxi_prod_new(all_mutdata_ext,"mut_sp1","mut_sp2",c("7","8","5","6"),c("3","4","1","2"),mutp_sp1,mutp_sp2)
    pxi_qxi_cox <- calc_pxi_qxi_prod_new(all_mutdata_cox,"mut_sp1","mut_sp2",c("7","8","5","6"),c("3","4","1","2"),mutp_sp1,mutp_sp2)
    
    ext_prob <- sum(pxi_qxi_ext)/(sum(pxi_qxi_ext) + sum(pxi_qxi_cox))
    new_prob_ext[counter,1] <- i#mutq_sp1
    new_prob_ext[counter,2] <- j#mutq_sp2
    new_prob_ext[counter,3] <- ext_prob
    counter <- counter + 1
    
  }
}

p_h1 <- ggplot(data = new_prob_ext, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = prob_extinct)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1, limits=c(0,1)) + #to make scale comparable to main text figure 6
  geom_rect(aes(xmin = xl - 0.5, xmax = xh + 0.5, ymin = yl - 0.5, ymax = yh + 0.5),
            fill = "transparent", color = "black", size = 1)
p_h1
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_may26/exts_15.pdf", width=4.5, height=3.75, units = "in")

 # p_hold <-  ggplot(data = prop_ext, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = frac_extinct)) +
#   theme_bw()+
#   geom_tile() +
#   #xlim(10,20)+
#   #ylim(10,20)+
#   #scale_fill_viridis_c(option = "viridis", direction = -1, begin = 0.4856479, end = 1) #begin = 0.4856479 for mk15
#   scale_fill_viridis_c(option = "viridis", direction = -1)
# p_hold

#differences between the two plots
probs_old_new <- merge(prop_ext,new_prob_ext)
probs_old_new$diffs <- probs_old_new$frac_extinct - probs_old_new$prob_extinct

p_diffs <-  ggplot(data = probs_old_new, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = diffs)) +
  theme_bw()+
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", limits=c(-1,1))+
  geom_rect(aes(xmin = xl - 0.5, xmax = xh + 0.5, ymin = yl - 0.5, ymax = yh + 0.5),
            fill = "transparent", color = "black", size = 1)

p_diffs
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_may26/exts_15_diffs.pdf", width=4.5, height=3.75, units = "in")

end

#----------including information on length of trajectory:-------
time_data_16 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_exp_16_ext2.csv", header = T))
#calculating extinction for all the data
tmax <- max(time_data_16$tick)

#check for which ones did not end simulation
length(which(time_data_16$tick == tmax))
temp_not_ext <- time_data_16[which(time_data_16$tick == tmax),]

#remove the ones that did not reach end of simulation from our analysis?
temp_rem <- unique(temp_not_ext[,c(1,2)])

exts <- which(time_data_16$num_individuals_species1 < 1 | time_data_16$num_individuals_species2 < 1)
time_data_sub <- time_data_16[exts,]
temp_ext <- unique(time_data_sub[,c(1,2)])
print(nrow(temp_ext))
ext_outcome <- rep(1,nrow(temp_ext))
temp2 <- cbind(temp_ext, time_data_sub[rownames(temp_ext),]$tick, ext_outcome)
colnames(temp2) <- c("rep","mes","tick_extinct","ext_outcome")

temp_all <- unique(time_data_16[,c(1,2)], fromLast = TRUE)
IDs <- as.numeric(rownames(temp_all))-10
tick_extinct <- time_data_16[IDs,]$tick
temp_all <- cbind(temp_all[,c(1,2)],tick_extinct)
temp_cox2 <- anti_join(temp_all, temp_ext)
print(nrow(temp_cox2))
temp_cox <- anti_join(temp_cox2, temp_rem)
ext_outcome <- rep(0,nrow(temp_cox))

temp3 <- cbind(temp_cox, ext_outcome)

time_ext_temp <- rbind(temp2, temp3) #contains length of trajectory info
#write.table(time_ext_temp, file = "~/Documents/SLiM/Rstuff/time_extinct_16_endtick.csv", append = F, sep = ",")
time_ext_temp <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_extinct_16_endtick.csv", header = T))
all_mutdata_large <- inner_join(all_mutdata,time_ext_temp)

#repeating above analysis for this new data frame now
#come up with 20 equally spaced mutation rates between 10^-7 and 10^-6, for testing 
mqs_list <- seq(1e-7,1e-6, length.out=20)

all_mutdata_ext <- subset(all_mutdata_large, ext_outcome > 0)
all_mutdata_cox <- subset(all_mutdata_large, ext_outcome < 1)

#function that does this+adjusts for supply of trait-increasing vs trait-decreasing muts?
#3,4,7,8 increasing, 1,2,5,6 decreasing
#write a function to calculate reweights for trajectories
calc_pxi_qxi_new <- function(df, mutp_sp1, mutp_sp2, sp1_cols_inc, sp1_cols_dec, sp2_cols_inc, sp2_cols_dec, mutq_sp1, mutq_sp2, traj) {
  
  bias_inc <- 1 #(0.25/0.05)
  bias_dec <- 1 #(0.25/0.45)
  
  frac_sp1_inc <- (mutq_sp1 * bias_inc) / (df[[mutp_sp1]]*df[[traj]])
  frac_sp2_inc <- (mutq_sp2 * bias_inc) / (df[[mutp_sp2]]*df[[traj]])
  frac_sp1_dec <- (mutq_sp1 * bias_dec) / (df[[mutp_sp1]]*df[[traj]])
  frac_sp2_dec <- (mutq_sp2 * bias_dec) / (df[[mutp_sp2]]*df[[traj]])
  
  sp1_terms_inc <- sweep(df[, sp1_cols_inc], 1, frac_sp1_inc, "*")
  sp2_terms_inc <- sweep(df[, sp2_cols_inc], 1, frac_sp2_inc, "*")
  sp1_terms_dec <- sweep(df[, sp1_cols_dec], 1, frac_sp1_dec, "*")
  sp2_terms_dec <- sweep(df[, sp2_cols_dec], 1, frac_sp2_dec, "*")
  
  rowSums(sp1_terms_inc) + rowSums(sp2_terms_inc) + rowSums(sp1_terms_dec) + rowSums(sp2_terms_dec)
}

new_prob_ext <- data_frame(mutation_kernel_1=numeric(),
                           mutation_kernel_2=numeric(),
                           prob_extinct=numeric()) #only for mutation effect
counter <- 1
for(i in 1:20){
  mutq_sp1 <- mqs_list[i]
  for(j in 1:20){
    mutq_sp2 <- mqs_list[j]
    pxi_qxi_ext <- calc_pxi_qxi_new(all_mutdata_ext,"mut_sp1","mut_sp2",c("7","8"),c("5","6"),c("3","4"),c("1","2"),mutq_sp1,mutq_sp2,"tick_extinct")
    pxi_qxi_cox <- calc_pxi_qxi_new(all_mutdata_cox,"mut_sp1","mut_sp2",c("7","8"),c("5","6"),c("3","4"),c("1","2"),mutq_sp1,mutq_sp2,"tick_extinct")
    
    ext_prob <- sum(pxi_qxi_ext)/(sum(pxi_qxi_ext) + sum(pxi_qxi_cox))
    new_prob_ext[counter,1] <- mutq_sp1
    new_prob_ext[counter,2] <- mutq_sp2
    new_prob_ext[counter,3] <- ext_prob
    counter <- counter + 1
    
  }
}

p_h1 <- ggplot(data = new_prob_ext, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = prob_extinct)) +
  theme_bw()+
  geom_tile() +
  #xlim(10,20)+
  #ylim(10,20)+
  #scale_fill_viridis_c(option = "viridis", direction = -1, begin = 0.4856479, end = 1) #begin = 0.4856479 for mk15
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h1

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_may26/extra_mut3.pdf", width=4.5, height=3.75, units = "in")

#--------- modify the code from figure 6 in main file, on extra simulations conducted with fixed mutation kernel params------

#load files
time_data_16 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_16_redo.csv", header = T))
mk16 <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutkern_16_redo.csv", header = T)) #loading this even though it's just. one value

#calculating extinction for all the data
tmax <- max(time_data_16$tick)

#check for which ones did not end simulation
length(which(time_data_16$tick == tmax))
temp_not_ext <- time_data_16[which(time_data_16$tick == tmax),]

#remove the ones that did not reach end of simulation from our analysis?
temp_rem <- unique(temp_not_ext[,c(1,2)]) #9.95% trajectories excluded. is this an issue?

exts <- which(time_data_16$num_individuals_species1 < 1 | time_data_16$num_individuals_species2 < 1)
time_data_sub <- time_data_16[exts,]
temp_ext <- unique(time_data_sub[,c(1,2)])
print(nrow(temp_ext))
temp2 <- cbind(temp_ext, time_data_sub[rownames(temp_ext),]$tick)
colnames(temp2) <- c("rep","mes","tick_extinct")

temp_all <- unique(time_data_16[,c(1,2)])
temp_cox2 <- anti_join(temp_all, temp_ext)
print(nrow(temp_cox2))
temp_cox <- anti_join(temp_cox2, temp_rem)
tick_extinct <- rep(500002, nrow(temp_cox))
temp3 <- cbind(temp_cox, tick_extinct)

time_ext_temp <- rbind(temp2, temp3)
time_ext <-  merge(time_ext_temp, mk16, by = c('rep', 'mes'))

#we don't have to bin because it's all one value

write.table(time_ext, file = "~/Documents/SLiM/Rstuff/time_extinct_16_redo_long.csv", append = F, sep = ",")
#read previously saved file rather than having to do all the analysis again
time_ext <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_extinct_16_redo.csv", header = T))

fixedmuts <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_comb_16_redo.csv", header = T))
mk_sub <- time_ext[sample(nrow(time_ext), 150), ] #random sampling
mk_sub <- time_ext
fixedmuts_sub <- fixedmuts %>% 
  semi_join(mk_sub, by = c("rep", "mes"))

#repeat the code from above, except mk_sub <- time_ext (still semi join fixedmuts_sub since a fraction of trajectories are excluded)

#--------- for effect size multiplier -----------
#using the fixed mutation supply parameters data here, even though it is on fixed time trajectories 

# old code based on limited trajectories
# fixed_muts_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_pair2new.csv", header = T))
# colnames_temp <- colnames(fixed_muts_full)
# colnames_temp <- c("rep","mes","mutation_kernel", "species","time_origin", "time_fixed","effect_size","mutation_type","blank") 
# colnames(fixed_muts_full) <- colnames_temp
# #only subset where replicates went extinct: log_combined_extinct.csv
# fixed_muts_extinct <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_extinct.csv", header = T))
# fixed_muts_extinct <- fixed_muts_extinct[,1:9]
# 
# #make combined data frame to facet simulation runs by mutation kernel parameters?
# mks <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutation_kernel_extra.csv", header = T))
# mks_actual <- mks[which(mks$mutr > 0.002),]
# 
# mks_act_list <- unique(mks_actual$mutation_kernel)
# mks_extinct <- unique(time_data_extinct$mutation_kernel)
# 
# mks_rest <- setdiff(mks_act_list, mks_extinct)
# 
# fixed_muts_sub <- subset(fixed_muts_full, mutation_kernel %in% mks_rest)
# fixed_muts_new <- rbind(fixed_muts_sub,fixed_muts_extinct) #data frame combining the data from extinct parameters
# 
# colnames_temp <- colnames(fixed_muts_new)
# colnames_temp <- c("rep","mes","mutation_kernel", "species","time_origin", "time_fixed","effect_size","mutation_type","blank") 
# colnames(fixed_muts_new) <- colnames_temp
# fixed_muts <-  subset(fixed_muts_new, rep == 1)
# fixed_muts$mutation_type <- as.factor(fixed_muts$mutation_type)
# 
# time_ext <- read.csv(file = "~/Documents/SLiM/Rstuff/time_extinct_fig4.csv", header = T) #extinction outcome for trajectories
# time_ext <- merge(time_ext, mks, by = "mutation_kernel")

#repeat the code above, make subset: choose mutation_kernel 47 (see AD_all_analysis for the code to get fraction of trjectories extinct)

mk_sub2 <- subset(time_ext, mutation_kernel == 47)

#make subset of fixedmuts to work with
fixedmuts_sub2 <- fixed_muts %>% 
  semi_join(mk_sub2, by = "mutation_kernel")

muttypes_grouped2 <- fixedmuts_sub2 |> group_by(mes,mutation_type)
muttypes_count2 <- muttypes_grouped2 |> tally(sort = F)

mt_count_wide2 <- muttypes_count2 |>
  pivot_wider(names_from = mutation_type, values_from = n, values_fill = 0)

#make vector of extinction outcomes (0 or 1) for chosen trajectories
mk_sub2$ext_outcome <- as.numeric(mk_sub2$extinct)

#make subset of mk_sub with only the info needed
mk_subset2 <- mk_sub2[,c(1,3,6,10,11,12,14)]

#join with mt_count_wide to make bigger data frame for all calculation needed
all_mutdata2 <- merge(mk_subset2,mt_count_wide2)

all_mutdata_ext2 <- subset(all_mutdata2, ext_outcome > 0)
all_mutdata_cox2 <- subset(all_mutdata2, ext_outcome < 1)

#test the previous functions on this new data? calc_pxi_qxi and calc_pxi_qxi_prod
#compare mutation_kernel 72, 47, 22, mutation probability 0.33, 0.49, 0.24

pxi_qxi_ext2 <- calc_pxi_qxi(all_mutdata_ext2,"mut","mut",c("7","8","5","6"),c("3","4","1","2"),1e-08,1e-08)
pxi_qxi_cox2 <- calc_pxi_qxi(all_mutdata_cox2,"mut","mut",c("7","8","5","6"),c("3","4","1","2"),1e-08,1e-08)

ext_prob2 <- sum(pxi_qxi_ext2)/(sum(pxi_qxi_ext2) + sum(pxi_qxi_cox2))
ext_prob2

#ok, so neither of them are working particularly well.

#new subset of fixed mutations: but not grouped
mk_subfixed2 <- fixedmuts_sub2[,c(2,3,4,7,8)]
mk_subfixed2$effect_size <- abs(mk_subfixed2$effect_size) #make sure all effect sizes are positive

#join with mt_count_wide to make bigger data frame for all calculation needed
all_mutdata2 <- merge(mk_subset2,mk_subfixed2)

all_mutdata_ext2 <- subset(all_mutdata2, ext_outcome > 0)
all_mutdata_cox2 <- subset(all_mutdata2, ext_outcome < 1)

#getting probability densities for all effect sizes at once?
test_effectsizes <- dnorm(all_mutdata2$effect_size, 0.5, 2)
effectsizes_new <- dnorm(all_mutdata2$effect_size, 0.5, 2*(0.1/0.05))
new_weights <- effectsizes_new/test_effectsizes

calc_pxi_qxi_me <- function(df, mutq_sp1, mutq_sp2, sp1_cols, sp2_cols, mutp_sp1, mutp_sp2) {
  
  frac_sp1 <- as.matrix(mutp_sp1 / df[[mutq_sp1]])
  frac_sp2 <- as.matrix(mutp_sp2 / df[[mutq_sp2]])
  
  sp1_inc <- rowSums(df[, sp1_cols])
  sp2_inc <- rowSums(df[, sp2_cols])
  
  sp1_terms_inc <- sweep(frac_sp1, 1, sp1_inc, "^")
  sp2_terms_inc <- sweep(frac_sp2, 1, sp2_inc, "^")
  
  sp1_terms_inc*sp2_terms_inc
  
}

#------------ 27 april 2026 --------
#extra analysis: what was the distribution of excluded simulation runs?

time_data_16 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_exp_16_ext2.csv", header = T))
mk16 <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutkern_16_ext2.csv", header = T))

#calculating extinction for all the data
tmax <- max(time_data_16$tick)

#check for which ones did not end simulation
length(which(time_data_16$tick == tmax))
temp_not_ext <- time_data_16[which(time_data_16$tick == tmax),]

#remove the ones that did not reach end of simulation from our analysis?
temp_rem <- unique(temp_not_ext[,c(1,2)]) # we need these

time_rem <-  merge(temp_rem, mk16, by = c('rep', 'mes'))
#read bins

bins <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/fig5_bins_16.csv", header = T))
bins_sp1 <- bins[1:20,]
bins_sp2 <- bins[21:40,]

#assigning bins to the excluded simulation runs
a <- time_rem$mut_sp1
bins_assigned_sp1 <- findInterval(a, bins_sp1$min)
bins_assigned_sp1_alt <- findInterval(a, bins_sp1$max)
bins_assigned_sp1[which(bins_assigned_sp1 - bins_assigned_sp1_alt < 1)] <- 1

b <- time_rem$mut_sp2
bins_assigned_sp2 <- findInterval(b, bins_sp2$min)
bins_assigned_sp2_alt <- findInterval(b, bins_sp2$max)

time_rem <- cbind(time_rem,bins_assigned_sp1,bins_assigned_sp2)

prop_rem <- data_frame(mutation_kernel_1=numeric(),
                       mutation_kernel_2=numeric(),
                       excluded=integer()) #only for mutation effect

mk_1 <- 1
mk_2 <- 1
counter <- 1
b <- 20 #binsize

for (mk_1 in 1:b) {
  for (mk_2 in 1:b){
    ext_sub <-  subset(time_rem, bins_assigned_sp1 == mk_1 & bins_assigned_sp2 == mk_2)
    nt <- nrow(ext_sub)
    if(nt > 0){
      prop_rem[counter,1] <- mk_1
      prop_rem[counter,2] <- mk_2
      prop_rem[counter,3] <- nt
    }else{
      prop_rem[counter,1] <- mk_1
      prop_rem[counter,2] <- mk_2
      prop_rem[counter,3] <- 0
    }
    counter <- counter + 1
    
  }
}

prop_all <- merge(prop_ext,prop_rem)
frac_excluded <- prop_all$excluded / (prop_all$total+prop_all$excluded)
prop_all <- cbind(prop_all, frac_excluded)

p_excl <- ggplot(data = prop_all, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = frac_excluded)) +
  theme_bw()+
  geom_tile() +
  xlim(4,21)+
  scale_fill_viridis_c(option = "viridis", direction = -1) + 
  geom_rect(aes(xmin = xl - 0.5, xmax = xh + 0.5, ymin = yl - 0.5, ymax = yh + 0.5),
            fill = "transparent", color = "black", size = 1)
p_excl
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_may26/fraction_excluded.pdf", width=4.5, height=3.75, units = "in")

#checking how close the existing trajectories that were removed were to extinction (or coexistence)
remIDs <- rownames(temp_rem)
time_data_rem <- merge(time_data_16,temp_rem)
time_data_rem_sub <- subset(time_data_rem, rep < 30 & tick > 400000)

gc1 <- ggplot()+
  theme_bw()+
  geom_line(data = time_data_rem_sub, mapping = aes(x = tick,y = r1, group = interaction(rep,mes)), alpha = 0.5, colour = "#26828e")+
  geom_hline(yintercept = r1_2sp)

gc2 <- ggplot()+
  theme_bw()+
  geom_line(data = time_data_rem_sub, mapping = aes(x = tick,y = r2, group = interaction(rep,mes)), alpha = 0.5, colour = "#6ece58")+
  geom_hline(yintercept = r2_2sp)
gc1
gc2
