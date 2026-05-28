#plotting trajectories for populations that went extinct
library(ggplot2)
library(dplyr)
library(cowplot)

# #full paired data: log_combined_pair2new.csv
# time_data_paired_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_combined_pair2new.csv", header = T))
# time_data_paired <- subset(time_data_paired_full, rep == 1)
# 
# #only subset where replicates went extinct: log_combined_extinct.csv
# time_data_extinct <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_combined_extinct.csv", header = T))
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
# time_data_paired_sub <- subset(time_data_paired, mutation_kernel %in% mks_rest)
# time_data_paired_new <- rbind(time_data_paired_sub, time_data_extinct) #data frame combining the data from extinct parameters
# write.table(time_data_paired_new, file = "~/Documents/SLiM/Rstuff/time_data_fig4.csv", append = F, sep = ",")
# 
# #replace the fixed muts like we replaced the trajectories
# #read fixed mutations
# fixed_muts_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_pair2new.csv", header = T))
# colnames_temp <- colnames(fixed_muts_full)
# colnames_temp <- c("rep","mes","mutation_kernel", "species","time_origin", "time_fixed","effect_size","mutation_type","blank") 
# colnames(fixed_muts_full) <- colnames_temp
# #only subset where replicates went extinct: log_combined_extinct.csv
# fixed_muts_extinct <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_extinct.csv", header = T))
# fixed_muts_extinct <- fixed_muts_extinct[,1:9]
#   
# fixed_muts_sub <- subset(fixed_muts_full, mutation_kernel %in% mks_rest)
# fixed_muts_new <- rbind(fixed_muts_sub,fixed_muts_extinct) #data frame combining the data from extinct parameters
# 
# colnames_temp <- colnames(fixed_muts_new)
# colnames_temp <- c("rep","mes","mutation_kernel", "species","time_origin", "time_fixed","effect_size","mutation_type","blank") 
# colnames(fixed_muts_new) <- colnames_temp
# fixed_muts <-  subset(fixed_muts_new, rep == 1)
# fixed_muts$mutation_type <- as.factor(fixed_muts$mutation_type)
# write.table(fixed_muts, file = "~/Documents/SLiM/Rstuff/fixed_muts_data_fig4.csv", append = F, sep = ",")

#---------- main analysis ------------
#code above is reading large files and making subsets
#1- #plotting trajectories for populations that went extinct
#use mutation kernel 47 or 48
#time_data_paired_new <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_data_fig4.csv", header = T))
#fixed_muts <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/fixed_muts_data_fig4.csv", header = T))


#fm_subset1 <- subset(fixed_muts, mutation_kernel == 47)
#fm_subset2 <- subset(fixed_muts, mutation_kernel == 48)

#gcsub1 <- subset(time_data_paired_new, mutation_kernel == 47)
#write.table(gcsub1, file = "~/Documents/SLiM/Rstuff/time_data_fig4_subset.csv", append = F, sep = ",")
gcsub1 <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_data_fig4_subset.csv", header = T))
#gcsub2 <- subset(time_data_paired_new, mutation_kernel == 48)

cmut <- c('1' = "#6ece58",
          '2' = "#6ece58",
          '3' = "#6ece58",
          '4' = "#6ece58",
          '5' = "#26828e",
          '6' = "#26828e",
          '7' = "#26828e",
          '8' = "#26828e") #colour names for the different species

#commenting out plots showing distribution of fixed mutation
# pl1 <- ggplot(data = fm_subset1, aes(x = time_origin, y = effect_size, colour = mutation_type)) +
#   theme_bw() +
#   geom_point(alpha = 0.5) +
#   scale_colour_manual(values = cmut) +
#   facet_wrap(~ mutation_type, ncol = 4) +
#   theme(legend.position = "none")
# 
# pl1
# 
# pl2 <- ggplot(data = fm_subset2, aes(x = time_origin, y = effect_size, colour = mutation_type)) +
#   theme_bw() +
#   geom_point(alpha = 0.5) +
#   scale_colour_manual(values = cmut) +
#   facet_wrap(~ mutation_type, nrow = 2) +
#   theme(legend.position = "none")
# 
# pl2
# 
# #pick rep 55 for subset 1 and 56 for subset 2 (extinct)
# #pick rep 16 for subset 1 and 14 for subset 2 (did not go extinct)
# fm_subset1_mes <- subset(fm_subset1, mes == 16)
# fm_subset2_mes <- subset(fm_subset2, mes == 14)
# 
# pl1_spec <- pl1 +
#   geom_point(data = subset(fm_subset1_mes, mutation_type %in% c(3,4,7,8)), aes(x = time_origin, y = effect_size), shape = 2, size = 2, stroke = 0.75, colour = '#000000') +
#   geom_point(data = subset(fm_subset1_mes, mutation_type %in% c(1,2,5,6)), aes(x = time_origin, y = effect_size), shape = 6, size = 2, stroke = 0.75,  colour = '#B40000')
# pl1_spec
# 
# pl2_spec <- pl2 +
#   geom_point(data = subset(fm_subset2_mes, mutation_type %in% c(3,4,7,8)), aes(x = time_origin, y = effect_size), shape = 2, size = 2, stroke = 0.75, colour = '#000000') +
#   geom_point(data = subset(fm_subset2_mes, mutation_type %in% c(1,2,5,6)), aes(x = time_origin, y = effect_size), shape = 6, size = 2, stroke = 0.75,  colour = '#B40000')
# 
# pl2_spec
# 
# p_all <- plot_grid(
#   plotlist = c(pl2_spec,pl1_spec),
#   ncol = 1,
#   byrow = TRUE
# )
# p_all #1048X1048
# #ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig4_supa.pdf", width=10.48, height=10.48, units = "in")
# 
# s <- 1
#growth curve plots

#growth curve subsets
#gcsub1_mes <- subset(gcsub1, mes == 16)
#gcsub2_mes <- subset(gcsub2, mes == 14)
gcsub1_mes1 <- subset(gcsub1, mes == 16)
gcsub1_mes2 <- subset(gcsub1, mes == 36)

#growth curve plots also highlighting the growth curves fron specific measurement reps
# gc1full <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub1, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") +
#   geom_line(data = gcsub1, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") +
#   facet_wrap(~ mes, nrow = 10)
# gc1full
# 
# gc2full <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub2, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") +
#   geom_line(data = gcsub2, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") +
#   facet_wrap(~ mes, nrow = 10)
# gc2full

# gc1_1 <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub1, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
#   geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 7), aes(x = time_origin, y = 0.245), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp1 
#   geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 5), aes(x = time_origin, y = 0.245), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp1 
#   geom_hline(yintercept = r1e)
# #gc1_1
# 
# gc1_2 <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub1, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") + 
#   geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = r2, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 3), aes(x = time_origin, y = 0.20), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp2 
#   geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 1), aes(x = time_origin, y = 0.20), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp2
#   geom_hline(yintercept = r2e)
# #gc1_2
# 
# gc2_1 <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub2, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
#   geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 7), aes(x = time_origin, y = 0.245), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp1 
#   geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 5), aes(x = time_origin, y = 0.245), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp1 
#   geom_hline(yintercept = r1e)
# #gc2_1
# 
# gc2_2 <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub2, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") + 
#   geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = r2, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 3), aes(x = time_origin, y = 0.20), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp2 
#   geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 1), aes(x = time_origin, y = 0.20), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp2 
#   geom_hline(yintercept = r2e)
# #gc2_2
# 
# #plot the same graphs for a11 and a22
# gc1_1a <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub1, mapping = aes(x = tick,y = a11, group = mes), alpha = 0.5, colour = "#26828e") + 
#   geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = a11, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 8), aes(x = time_origin, y = -0.0001), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp1 
#   geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 6), aes(x = time_origin, y = -0.0001), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp1 
#   geom_hline(yintercept = a11e)
# #gc1_1a
# 
# gc1_2a <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub1, mapping = aes(x = tick,y = a22, group = mes), alpha = 0.5, colour = "#6ece58") + 
#   geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = a22, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 4), aes(x = time_origin, y = -0.000035), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp2 
#   geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 2), aes(x = time_origin, y = -0.000035), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp2 
#   geom_hline(yintercept = a22e)
# #gc1_2a
# 
# gc2_1a <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub2, mapping = aes(x = tick,y = a11, group = mes), alpha = 0.5, colour = "#26828e") + 
#   geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = a11, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 8), aes(x = time_origin, y = -0.00008), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp1 
#   geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 6), aes(x = time_origin, y = -0.00008), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp1
#   geom_hline(yintercept = a11e)
# #gc2_1a
# 
# gc2_2a <- ggplot()+
#   theme_bw()+
#   geom_line(data = gcsub2, mapping = aes(x = tick,y = a22, group = mes), alpha = 0.5, colour = "#6ece58") + 
#   geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = a22, group = mes), colour = "#000000") + 
#   geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 4), aes(x = time_origin, y = -0.00003), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp2 
#   geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 2), aes(x = time_origin, y = -0.00003), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp2
#   geom_hline(yintercept = a22e)
# #gc2_2a
# 
# p_gc2 <- plot_grid(
#   plotlist = c(gc2_1,gc2_1a,gc2_2,gc2_2a,gc1_1,gc1_1a,gc1_2,gc1_2a),
#   ncol = 4,
#   byrow = TRUE
# )
# p_gc2#1048X1048
# ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_mar26/fig4b_extinct_alt.pdf", width=10.00, height=4.00, units = "in")

#----------- plotting the actual growth curves -----------
K1 <- 3079.024202
K2 <- 7433.758301
gc1_1 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = num_individuals_species1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub1_mes1, mapping = aes(x = tick,y = num_individuals_species1, group = mes), colour = "#000000")+
  geom_hline(yintercept = K1, colour = "#444444")

gc1_2 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = num_individuals_species2, group = mes), alpha = 0.5, colour = "#6ece58") + 
  geom_line(data = gcsub1_mes1, mapping = aes(x = tick,y = num_individuals_species2, group = mes), colour = "#000000")+
  geom_hline(yintercept = K2, colour = "#666666") 

gc2_1 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = num_individuals_species1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub1_mes2, mapping = aes(x = tick,y = num_individuals_species1, group = mes), colour = "#000000")+
  geom_hline(yintercept = K1, colour = "#444444")

gc2_2 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = num_individuals_species2, group = mes), alpha = 0.5, colour = "#6ece58") + 
  geom_line(data = gcsub1_mes2, mapping = aes(x = tick,y = num_individuals_species2, group = mes), colour = "#000000") +
  geom_hline(yintercept = K2, colour = "#666666") 

p_gc3 <- plot_grid(
  plotlist = c(gc2_1,gc2_2,gc1_1,gc1_2),
  ncol = 2,
  byrow = TRUE
)
p_gc3#1048X1048
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_may26/fig4b_gcs_extinct2.pdf", width=6.00, height=4.00, units = "in")

#part 1: distribution of mutations 

#read file in the future
time_ext <- read.csv(file = "~/Documents/SLiM/Rstuff/time_extinct_fig4.csv", header = T)
mks <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutation_kernel_extra.csv", header = T))
tt_ext_large <- merge(time_ext, mks, by = "mutation_kernel")
tt_ext_large$mut <- as.factor(tt_ext_large$mut)

frac_extinct <- data_frame(mutation_kernel=integer(),
                           extinct = integer())
mks_actual <- mks[which(mks$mutr > 0.002),]
mks_list <- mks_actual$mutation_kernel
counter <- 1
for(mk in mks_list){
  ts <- subset(time_ext, mutation_kernel == mk)
  num_extinct <- length(which(ts$tick_extinct < 50000))
  frac_extinct[counter,1] <- mk
  frac_extinct[counter,2] <- num_extinct
  counter <- counter+1
}
frac_extinct_large <- merge(frac_extinct, mks, by = "mutation_kernel")
frac_extinct_large$mut <- as.factor(frac_extinct_large$mut)
frac_extinct_large$mutr <- as.factor(frac_extinct_large$mutr)
frac_extinct_large$m3 <- as.factor(frac_extinct_large$m3)
frac_extinct_subset <- subset(frac_extinct_large, mut=="1e-07")
p_tt3 <- ggplot(data = frac_extinct_subset, aes(x = m3, y = mutr, fill = extinct)) +
  theme_bw()+
  #geom_boxplot(alpha = 0.6)+
  geom_tile()+
  scale_fill_viridis_c(begin = 0, end = 1, direction = -1, limits=c(0,100)) +
  geom_rect(aes(xmin = 2.55, xmax = 3.45, ymin = 2.55, ymax = 3.45),
            fill = "transparent", color = "#dddddd", size = 1)+
  labs(title = "mutation rate = 1e-07", x = "mutation bias", y = "mutation effect size") +
  theme(plot.title = element_text(hjust = 0.5))
p_tt3
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_may26/fig4a_alt.pdf", width=4.20, height=3.40, units = "in")
