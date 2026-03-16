#putting together code for the figures

setwd("/Users/kasturilele/Documents/SLiM")
library(ggplot2)
library(dplyr)
library(cowplot)

#-------- figure 2 ----------
# Define bt1 (exponential function)
bt1 <- function(b) {-0.0000206 * exp(3.2*b)} #function(b) {-0.0000108 * exp(3.5*b)} function(b) {-0.0000206 * exp(3.2 * b)}
bt2 <- function(b) {-0.0000108 * exp(3.5*b)}

#list of alternative functions bt2 (for supplementary figure)
#bt2 <- function(b) 
  #{-0.00001} #constant
  #{-0.0004*b} #linear 
  #{-0.0001 * b*b} #square
  #{-0.0001 * b*b*b} #cubed
  #{-0.00001 * b /(1 - b)} #to get linear r-K tradeoff
  #{-0.00001/(-1*(b - 1.2))} #inverse, also modified concave boots (10)
  #{-0.00001 * b / ((b - 0.5)*(b-0.5) + 0.01)} #inverse parabola tradeoff between r and K
#{-0.0001 * b / log(10*b + 1)} #to get log relationship between r and K
#{-0.0001/((b + 1))}  #modified convex boots

# Define selection gradient
selection_gradient1 <- function(x, y) {
  y - (x * bt2(y) / bt2(x))
}

# Endpoint for the graph #use (-0.0001, 0.0001) for a and 1 for r (only for supplementary figure)
e <- 1

# Create grid
num_points <- 200
x <- seq(0, e, length.out = num_points)
y <- seq(0, e, length.out = num_points)

grid <- expand.grid(
  x1 = x,
  y1 = y
)

# Evaluate selection gradient on grid
grid_s <- grid %>%
  mutate(Z1 = selection_gradient1(x1, y1))

#tradeoff function plot
p2 <- ggplot() +
  theme_bw() +
  stat_function(fun = bt2, linewidth = 1.5, colour = "#26828e") +
  xlim(0,1) +
  theme(legend.position = 'none')
p2

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/for_figures/fig_supfun_6_fun.pdf", width=2.00, height=2.00, units = "in")

#plot region where selection gradient positive
p1 <- ggplot(grid_s, aes(x = x1, y = y1)) +
  theme_bw()+
  geom_raster(aes(fill = Z1 > 0), alpha = 0.1) +
  scale_fill_manual(values = c("white", "#26828e")) +
  # Boundary where Z1 = 0
  geom_contour(aes(z = Z1),breaks = 0,color = "#26828e",alpha = 0.75) +
  labs(x = "resident", y = "mutant") +
  coord_equal() +
  theme(legend.position = 'none')
p1
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/for_figures/fig_supfun_6.pdf", width=2.00, height=2.00, units = "in")

# Define paired selection gradient
sel_grad_paired <- function(x1, y1, x2, y2) {
  a12 = -2.5e-05
  a21 = -2.5e-05
  den = bt2(x2)*bt1(x1) - a12*a21
  s1m = y1 + bt1(y1)*(a12*x2 - bt2(x2)*x1)/den + (a12*a21*x1 - bt1(x1)*a12*x2)/den
  s2m = y2 + bt2(y2)*(a21*x1 - bt1(x1)*x2)/den + (a12*a21*x2 - bt2(x2)*a21*x1)/den
  
  dftemp <- data.frame(s1m, s2m)
  return (dftemp)
}

grid_p <- expand.grid(
  x2 = x,
  y2 = y
)

grid_p <- cbind(grid, grid_p)

# Evaluate selection gradient on paired grid
sel_paired <- sel_grad_paired(grid_p$x1, grid_p$y1, grid_p$x2, grid_p$y2)

grid_p <- cbind(grid_p, sel_paired)

#plot region where selection gradient positive
p3 <- ggplot(grid_p, aes(x = x1, y = y1)) +
  theme_bw()+
  geom_raster(aes(fill = s1m > 0), alpha = 0.1) +
  scale_fill_manual(values = c("white", "#26828e")) +
  # Boundary where Z1 = 0
  geom_contour(aes(z = s1m),breaks = 0,color = "#26828e",alpha = 0.75) +
  labs(x = "resident", y = "mutant") +
  coord_equal() +
  theme(legend.position = 'none')
p3
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/for_figures/fig2c_1_supr.pdf", width=3.00, height=3.00, units = "in")

p4 <- ggplot(grid_p, aes(x = x2, y = y2)) +
  theme_bw()+
  geom_raster(aes(fill = s2m > 0), alpha = 0.1) +
  scale_fill_manual(values = c("white", "#6ece58")) +
  # Boundary where Z1 = 0
  geom_contour(aes(z = s2m),breaks = 0,color = "#6ece58",alpha = 0.75) +
  labs(x = "resident", y = "mutant") +
  coord_equal() +
  theme(legend.position = 'none')
p4
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/for_figures/fig2c_2_supr.pdf", width=3.00, height=3.00, units = "in")

p5 <- ggplot() +
  theme_bw() +
  stat_function(fun = bt1, linewidth = 1.5, colour = "#26828e") +
  stat_function(fun = bt2, linewidth = 1.5, colour = "#6ece58") +
  xlim(0,1) +
  theme(legend.position = 'none')
p5

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/for_figures/fig2d.pdf", width=3.00, height=3.00, units = "in")

#read file

ESSpred <- as.data.frame(read.csv("~/Documents/SLiM/tradeoff_fun_ESS.csv", header = T))

ESSpred_S1 <- as.data.frame(read.csv("~/Documents/SLiM/tradeoff_fun_ESS1.csv", header = T))
ESSpred_S2 <- as.data.frame(read.csv("~/Documents/SLiM/tradeoff_fun_ESS2.csv", header = T))

p_e1 <- ggplot(data = ESSpred_S1, aes(x = C2, y = C1, fill = outcome_pair)) +
  theme_bw()+
  geom_tile() +
  scale_fill_manual(values = c("#440154", "#fde725")) +
  geom_point(x = 3.2, y = -0.0000206, size = 4, colour = "#26828e" ) +
  theme(legend.position = 'none')
p_e1
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/for_figures/fig2e1.pdf", width=3.20, height=3.20, units = "in")

p_e2 <- ggplot(data = ESSpred_S2, aes(x = C4, y = C3, fill = outcome_pair)) +
  theme_bw()+
  geom_tile() +
  scale_fill_manual(values = c("#440154", "#fde725")) +
  geom_point(x = 3.5, y = -0.0000108, size = 4, colour = "#6ece58") +
  theme(legend.position = 'none')
p_e2
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/for_figures/fig2e2.pdf", width=3.20, height=3.20, units = "in")

#figure 2 supplementary
# Define selection gradient without constraints 
#(sub am = a for muts in r)
selection_gradient1 <- function(x, y) {
  y - x 
}

#(sub rm = r for muts in a)
selection_gradient1 <- function(x, y) {
  1 - y/x 
}

#paired selection gradients without constraints
#(sub am = a for muts in r)
sel_grad_paired <- function(x1, y1, x2, y2) {
  s1m = y1 - x1
  s2m = y2 - x2
  dftemp <- data.frame(s1m, s2m)
  return (dftemp)
}

#(sub rm = r for muts in a)
sel_grad_paired <- function(x1, y1, x2, y2) {
  a12 = -2.5e-05
  a21 = -2.5e-05
  r1 = 0.5
  r2 = 0.5
  den = x2*x1 - a12*a21
  s1m = ((r2*a12 - r1*x2)*(y1 - x1))/den
  s2m = ((r1*a21 - r2*x1)*(y2 - x2))/den
  
  dftemp <- data.frame(s1m, s2m)
  return (dftemp)
}

#-------- figure 3 - single species simulations ----------

time_data_full <- as.data.frame(read.csv("~/Documents/SLiM/output2/log_combined_single3.csv", header = T))
mks_full <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutation_kernel_extra.csv", header = T))
reps <- as.data.frame(read.csv("~/Documents/SLiM/final_params/single_params_3.csv", header = T)) #don't use this for the subset with no parameters

#functions: #function(b) {-0.0000206 * exp(3.2*b)} #function(b) {-0.0000108 * exp(3.5*b)}
fun_exp <- function(b) {-0.0000206 * exp(3.2*b)}
fun_K <- function(r) {-r/fun_exp(r)}
#distance function for thresholding
fun_minim <- function(x, y, px, py){
  sqrt((x - px)*(x - px)/0.01 + (y - py)*(y - py)/1e-10)
}

#equilibrium values of r and alpha (from python)
r1e <- 0.3125 #0.3125 #0.28571429
a11e <- -5.59966057e-05 #-5.59966057e-05 #-2.93574437e-05

time_data <- subset(time_data_full, rep == 2) #do all the reps?

time_threshold <- data_frame(mutation_kernel=character(),
                             rep=integer(),
                             mes=integer(),
                             time=integer(),
                             distance=double())
mk <- 0
r <- 2
m <- 1
counter <- 1

mks_list <- unique(time_data$mutation_kernel)
mks_actual <- mks[which(mks$mutr > 0.002),]

for (mk in mks_list) {
  #for (r in 0:5) {
  for (m in 1:100) {
    time_data_sub <- subset(time_data, rep == r & mes == m & mutation_kernel == mk)
    if(nrow(time_data_sub > 0)){
      #calculate minimum distance from ESS (to minimize oscillations around the ESS?)
      distance_all <- fun_minim(r1e, a11e, time_data_sub$r1, time_data_sub$a11)
      distance <- distance_all[which(distance_all == min(distance_all))][1]
      
      tmax <- max(time_data_sub$tick)
      tt <- 0
      if(length(which(distance_all < 0.5)) > 0){
        tt <- min(time_data_sub$tick[which(distance_all < 0.5)])
      } else {
        tt <- tmax
      }
      
      #save all the values
      time_threshold[counter,1] <- as.character(time_data_sub[1,1])
      time_threshold[counter,2] <- time_data_sub[1,2]
      time_threshold[counter,3] <- time_data_sub[1,3]
      time_threshold[counter,4] <- tt
      time_threshold[counter,5] <- distance 
      counter <- counter + 1
    }
  }
  #}
}

tt_combine <- merge(time_threshold, mks_full, by = "mutation_kernel")
#tt_combine_large <- merge(tt_combine, reps, by.x = "rep", by.y = "IDs")
tt_combine_large <- tt_combine
tt_combine_large$mut <- as.factor(tt_combine_large$mut)

#plot time taken to reach threshold
p_tt <- ggplot(data = tt_combine_large, aes(x = mut, y = time, colour = mut, fill = mut)) +
  theme_bw()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6)+
  scale_colour_manual(values = c("#26828E","#2A5B6F","#203446"))+ #c("#6ece58","#6D8F3D","#576435")
  scale_fill_manual(values = c("#26828E","#2A5B6F","#203446"))+
  #scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  #scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(-mutr~m3)+
  theme(legend.position = "none")
p_tt
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig3a_tt2.pdf", width=7.5, height=5.0, units = "in")

p_dist <-ggplot(data = tt_combine_large, aes(x = mut, y = distance, colour = mut, fill = mut)) +
  theme_bw()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6)+
  #scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  #scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  scale_colour_manual(values = c("#26828E","#2A5B6F","#203446"))+ #c("#6ece58","#6D8F3D","#576435")
  scale_fill_manual(values = c("#26828E","#2A5B6F","#203446"))+
  facet_grid(-mutr~m3)+
  theme(legend.position = "none")
p_dist
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig3a_dist2.pdf", width=7.5, height=5.0, units = "in")

fixed_muts_sing <- as.data.frame(read.csv("~/Documents/SLiM/output2/fixed_combined_single3.csv", header = T))
fixed_muts_sing <-  subset(fixed_muts_sing, rep == 2)
fixed_muts_sing$mutation_type <- as.factor(fixed_muts_sing$mutation_type)

fm_sub1 <- subset(fixed_muts_sing, mutation_kernel == 37)
fm_sub2 <- subset(fixed_muts_sing, mutation_kernel == 70)

gc_sub1 <- subset(time_data, mutation_kernel == 37)
gc_sub2 <- subset(time_data, mutation_kernel == 70)

cmut <- c('1' = "#26828e",'2' = "#26828e",'3' = "#26828e",'4' = "#26828e") #colour names for the different species

#c('1' = "#6ece58",'2' = "#6ece58",'3' = "#6ece58",'4' = "#6ece58",)

pl1 <- ggplot(data = fm_sub1, aes(x = time_origin, y = effect_size, colour = mutation_type)) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = cmut) +
  facet_wrap(~ mutation_type, ncol = 4) +
  theme(legend.position = "none")

pl1

pl2 <- ggplot(data = fm_sub2, aes(x = time_origin, y = effect_size, colour = mutation_type)) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = cmut) +
  facet_wrap(~ mutation_type, ncol = 4) +
  theme(legend.position = "none")

pl2

#growth curve plots also highlighting the growth curves fron specific measurement reps
gc1full <- ggplot()+
  theme_bw()+
  geom_line(data = gc_sub1, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") +
  facet_wrap(~ mes, nrow = 10)
gc1full

gc2full <- ggplot()+
  theme_bw()+
  geom_line(data = gc_sub2, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  facet_wrap(~ mes, nrow = 10)
gc2full

#pick rep 45 for subset 1 and 65 for subset 2 - single4
#pick rep 26 for subset 1 and 35 for subset 2 - single3
fm_sub1_mes <- subset(fm_sub1, mes == 26)
fm_sub2_mes <- subset(fm_sub2, mes == 35)

#growth curve subsets
gc_sub1_mes <- subset(gc_sub1, mes == 26)
gc_sub2_mes <- subset(gc_sub2, mes == 35)

pl1_spec <- pl1 +
  geom_point(data = subset(fm_sub1_mes, mutation_type %in% c(3,4)), aes(x = time_origin, y = effect_size), shape = 2, size = 2, stroke = 0.75, colour = '#000000') +
  geom_point(data = subset(fm_sub1_mes, mutation_type %in% c(1,2)), aes(x = time_origin, y = effect_size), shape = 6, size = 2, stroke = 0.75,  colour = '#B40000')+
  theme(strip.text = element_blank())

pl2_spec <- pl2 +
  geom_point(data = subset(fm_sub2_mes, mutation_type %in% c(3,4)), aes(x = time_origin, y = effect_size), shape = 2, size = 2, stroke = 0.75, colour = '#000000') +
  geom_point(data = subset(fm_sub2_mes, mutation_type %in% c(1,2)), aes(x = time_origin, y = effect_size), shape = 6, size = 2, stroke = 0.75,  colour = '#B40000')+
  theme(strip.text = element_blank())

p_all <- plot_grid(
  plotlist = c(pl2_spec,pl1_spec),
  ncol = 1,
  byrow = TRUE
)
p_all #1048X1048
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig3_supa2.pdf", width=10.48, height=5.12, units = "in")

s <- 1
#growth curve plots
gc1_1 <- ggplot()+
  theme_bw()+
  geom_line(data = gc_sub1, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gc_sub1_mes, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_sub1_mes, mutation_type == 3), aes(x = time_origin, y = r1e), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp1 
  geom_point(data = subset(fm_sub1_mes, mutation_type == 1), aes(x = time_origin, y = r1e), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp1 
  geom_hline(yintercept = r1e)

gc2_1 <- ggplot()+
  theme_bw()+
  geom_line(data = gc_sub2, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gc_sub2_mes, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_sub2_mes, mutation_type == 3), aes(x = time_origin, y = r1e), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp1 
  geom_point(data = subset(fm_sub2_mes, mutation_type == 1), aes(x = time_origin, y = r1e), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp1 
  geom_hline(yintercept = r1e)
#gc2_1

#plot the same graphs for a11 and a22
gc1_1a <- ggplot()+
  theme_bw()+
  geom_line(data = gc_sub1, mapping = aes(x = tick,y = a11, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gc_sub1_mes, mapping = aes(x = tick,y = a11, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_sub1_mes, mutation_type == 4), aes(x = time_origin, y = a11e), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp1 
  geom_point(data = subset(fm_sub1_mes, mutation_type == 2), aes(x = time_origin, y = a11e), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp1 
  geom_hline(yintercept = a11e)
#gc1_1a

gc2_1a <- ggplot()+
  theme_bw()+
  geom_line(data = gc_sub2, mapping = aes(x = tick,y = a11, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gc_sub2_mes, mapping = aes(x = tick,y = a11, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_sub2_mes, mutation_type == 4), aes(x = time_origin, y = a11e), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp1 
  geom_point(data = subset(fm_sub2_mes, mutation_type == 2), aes(x = time_origin, y = a11e), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp1
  geom_hline(yintercept = a11e)
#gc2_1a

p_gc2 <- plot_grid(
  plotlist = c(gc2_1,gc2_1a,gc1_1,gc1_1a),
  ncol = 4,
  byrow = TRUE
)
p_gc2#1048X1048
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig3b2.pdf", width=10.00, height=2.00, units = "in")


#------- figure 4 -------
#full paired data: log_combined_pair2new.csv
time_data_paired_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_combined_pair2new.csv", header = T))
time_data_paired <- subset(time_data_paired_full, rep == 1)

#only subset where replicates went extinct: log_combined_extinct.csv
time_data_extinct <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_combined_extinct.csv", header = T))

#make combined data frame to facet simulation runs by mutation kernel parameters?
mks <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutation_kernel_extra.csv", header = T))
mks_actual <- mks[which(mks$mutr > 0.002),]

mks_act_list <- unique(mks_actual$mutation_kernel)
mks_extinct <- unique(time_data_extinct$mutation_kernel)

mks_rest <- setdiff(mks_act_list, mks_extinct)

time_data_paired_sub <- subset(time_data_paired, mutation_kernel %in% mks_rest)
time_data_paired_new <- rbind(time_data_paired_sub, time_data_extinct) #data frame combining the data from extinct parameters

tdp_large <- merge(time_data_paired_new, mks_actual, by = "mutation_kernel")
tmax_p <- max(time_data_paired_new$tick)

min(time_data_paired$num_individuals_species1)
min(time_data_paired$num_individuals_species2)

#calculating extinction for all the data
tmax <- max(time_data_paired_new$tick)
exts <- which(time_data_paired_new$num_individuals_species1 < 1 | time_data_paired_new$num_individuals_species2 < 1)
time_data_sub <- time_data_paired_new[exts,]
temp_ext <- unique(time_data_sub[,c(1,2,3)])
print(nrow(temp_ext))
extinct <- rep(T, nrow(temp_ext))
temp2 <- cbind(temp_ext, time_data_sub[rownames(temp_ext),]$tick, extinct)
colnames(temp2) <- c("mutation_kernel","rep","mes","tick_extinct","extinct")

temp_all <- unique(time_data_paired_new[,c(1,2,3)])
temp_cox <- anti_join(temp_all, temp_ext)
print(nrow(temp_cox))
tick_extinct <- rep(tmax, nrow(temp_cox))
extinct <- rep(F, nrow(temp_cox))
temp3 <- cbind(temp_cox, tick_extinct, extinct)

time_ext <- rbind(temp2, temp3)
#write.table(time_ext, file = "~/Documents/SLiM/Rstuff/time_extinct_fig4.csv", append = F, sep = ",")
#read file in the future
time_ext <- read.csv(file = "~/Documents/SLiM/Rstuff/time_extinct_fig4.csv", header = T)

tt_ext_large <- merge(time_ext, mks, by = "mutation_kernel")
tt_ext_large$mut <- as.factor(tt_ext_large$mut)

#tt_ext_large1 <- subset(tt_ext_large, mutr > 0.002)

p_tt2 <- ggplot(data = tt_ext_large, aes(x = mut, y = tick_extinct, colour = mut, fill = mut)) +
  theme_bw()+
  #geom_boxplot(alpha = 0.6)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6)+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(-mutr~m3)+
  theme(legend.position = "none")
p_tt2
#ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig4a_horz.pdf", width=7.5, height=5, units = "in")

#side analysis for fig 6 (fig 5 in paper): calculate fraction of pops that went extinct
frac_extinct <- data_frame(mutation_kernel=integer(),
                           extinct = integer())
mks_list <- mks_actual$mutation_kernel
counter <- 1
for(mk in mks_list){
  ts <- subset(time_ext, mutation_kernel == mk)
  num_extinct <- length(which(ts$tick_extinct < 50000))
  frac_extinct[counter,1] <- mk
  frac_extinct[counter,2] <- num_extinct
  counter <- counter+1
}
#mutation kernel 74 has the most extinct (55%), use this for simulations in fig.6

library(cowplot)

#read fixed mutations
fixed_muts_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_pair2new.csv", header = T))
colnames_temp <- colnames(fixed_muts_full)
colnames_temp <- c("rep","mes","mutation_kernel", "species","time_origin", "time_fixed","effect_size","mutation_type","blank") 
colnames(fixed_muts_full) <- colnames_temp
fixed_muts <-  subset(fixed_muts_full, rep == 1)
fixed_muts$mutation_type <- as.factor(fixed_muts$mutation_type)

#1- compare pops that didn't go extinct, same mut rate (mk 64 vs 65)
#alternative mks - 44, 17 subset 1 65,70,80 subset 2

fm_subset1 <- subset(fixed_muts, mutation_kernel == 64)
fm_subset2 <- subset(fixed_muts, mutation_kernel == 65)

gcsub1 <- subset(time_data_paired_new, mutation_kernel == 64)
gcsub2 <- subset(time_data_paired_new, mutation_kernel == 65)

cmut <- c('1' = "#6ece58",
          '2' = "#6ece58",
          '3' = "#6ece58",
          '4' = "#6ece58",
          '5' = "#26828e",
          '6' = "#26828e",
          '7' = "#26828e",
          '8' = "#26828e") #colour names for the different species

pl1 <- ggplot(data = fm_subset1, aes(x = time_origin, y = effect_size, colour = mutation_type)) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = cmut) +
  facet_wrap(~ mutation_type, ncol = 4) +
  theme(legend.position = "none")

pl1

pl2 <- ggplot(data = fm_subset2, aes(x = time_origin, y = effect_size, colour = mutation_type)) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = cmut) +
  facet_wrap(~ mutation_type, nrow = 2) +
  theme(legend.position = "none")

pl2

#growth curve plots also highlighting the growth curves fron specific measurement reps
gc1full <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub1, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") +
  facet_wrap(~ mes, nrow = 10)
gc1full

gc2full <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub2, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub2, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") + 
  facet_wrap(~ mes, nrow = 10)
gc2full

#pick rep 63 for subset 1 and 73 for subset 2
#alt: 88 for subset 1 (44) 54 for subset2 (65) 75 for (80)
fm_subset1_mes <- subset(fm_subset1, mes == 73)
fm_subset2_mes <- subset(fm_subset2, mes == 63)

#growth curve subsets
gcsub1_mes <- subset(gcsub1, mes == 73)
gcsub2_mes <- subset(gcsub2, mes == 63)

pl1_spec <- pl1 +
  geom_point(data = subset(fm_subset1_mes, mutation_type %in% c(3,4,7,8)), aes(x = time_origin, y = effect_size), shape = 2, size = 2, stroke = 0.75, colour = '#000000') +
  geom_point(data = subset(fm_subset1_mes, mutation_type %in% c(1,2,5,6)), aes(x = time_origin, y = effect_size), shape = 6, size = 2, stroke = 0.75,  colour = '#B40000')
pl1_spec

pl2_spec <- pl2 +
  geom_point(data = subset(fm_subset2_mes, mutation_type %in% c(3,4,7,8)), aes(x = time_origin, y = effect_size), shape = 2, size = 2, stroke = 0.75, colour = '#000000') +
  geom_point(data = subset(fm_subset2_mes, mutation_type %in% c(1,2,5,6)), aes(x = time_origin, y = effect_size), shape = 6, size = 2, stroke = 0.75,  colour = '#B40000')

pl2_spec

p_all <- plot_grid(
  plotlist = c(pl2_spec,pl1_spec),
  ncol = 1,
  byrow = TRUE
)
p_all #1048X1048
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig4_supa.pdf", width=10.48, height=10.48, units = "in")

s <- 1
#growth curve plots
gc1_1 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 7), aes(x = time_origin, y = 0.245), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp1 
  geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 5), aes(x = time_origin, y = 0.245), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp1 
  geom_hline(yintercept = r1e)
#gc1_1

gc1_2 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") + 
  geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = r2, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 3), aes(x = time_origin, y = 0.20), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp2 
  geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 1), aes(x = time_origin, y = 0.20), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp2
  geom_hline(yintercept = r2e)
#gc1_2

gc2_1 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub2, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 7), aes(x = time_origin, y = 0.245), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp1 
  geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 5), aes(x = time_origin, y = 0.245), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp1 
  geom_hline(yintercept = r1e)
#gc2_1

gc2_2 <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub2, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") + 
  geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = r2, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 3), aes(x = time_origin, y = 0.20), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben r sp2 
  geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 1), aes(x = time_origin, y = 0.20), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del r sp2 
  geom_hline(yintercept = r2e)
#gc2_2

#plot the same graphs for a11 and a22
gc1_1a <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = a11, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = a11, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 8), aes(x = time_origin, y = -0.0001), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp1 
  geom_point(data = subset(fm_subset1_mes, species == 1 & mutation_type == 6), aes(x = time_origin, y = -0.0001), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp1 
  geom_hline(yintercept = a11e)
#gc1_1a

gc1_2a <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub1, mapping = aes(x = tick,y = a22, group = mes), alpha = 0.5, colour = "#6ece58") + 
  geom_line(data = gcsub1_mes, mapping = aes(x = tick,y = a22, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 4), aes(x = time_origin, y = -0.000035), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp2 
  geom_point(data = subset(fm_subset1_mes, species == 2 & mutation_type == 2), aes(x = time_origin, y = -0.000035), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp2 
  geom_hline(yintercept = a22e)
#gc1_2a

gc2_1a <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub2, mapping = aes(x = tick,y = a11, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = a11, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 8), aes(x = time_origin, y = -0.00008), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp1 
  geom_point(data = subset(fm_subset2_mes, species == 1 & mutation_type == 6), aes(x = time_origin, y = -0.00008), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp1
  geom_hline(yintercept = a11e)
#gc2_1a

gc2_2a <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub2, mapping = aes(x = tick,y = a22, group = mes), alpha = 0.5, colour = "#6ece58") + 
  geom_line(data = gcsub2_mes, mapping = aes(x = tick,y = a22, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 4), aes(x = time_origin, y = -0.00003), shape = 2, size = s, stroke = 0.75, colour = '#000000')+ #ben a sp2 
  geom_point(data = subset(fm_subset2_mes, species == 2 & mutation_type == 2), aes(x = time_origin, y = -0.00003), shape = 6, size = s, stroke = 0.75, colour = '#B40000')+ #del a sp2
  geom_hline(yintercept = a22e)
#gc2_2a

p_gc2 <- plot_grid(
  plotlist = c(gc2_1,gc2_1a,gc2_2,gc2_2a,gc1_1,gc1_1a,gc1_2,gc1_2a),
  ncol = 4,
  byrow = TRUE
)
p_gc2#1048X1048
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig4b.pdf", width=10.00, height=4.00, units = "in")

#expanded graphs over 500000 generations for subset of mutation kernels (supplementary)
time_data_paired_new <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_combined_cox_ext.csv", header = T))
fixed_muts_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_cox_ext.csv", header = T))
tdp_large <- merge(time_data_paired_new, mks_actual, by = "mutation_kernel")
#run the above code to a) check which populations went extinct
# + starting from 'analysis of fixed mutations and trajectory of evolution'

#-------- figure 6 in paper --------

#for mutation rate. did similar analysis for mutation effect size, different filenames. not included in paper
time_data_15 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_exp_15_ext2.csv", header = T))
mk15 <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutkern_15_ext2.csv", header = T))

#for mutation effect size-
#time_data_15 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_me_new.csv", header = T))
#mk15 <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutkern_me_new.csv", header = T))

#calculating extinction for all the data
tmax <- max(time_data_15$tick)

#check for which ones did not end simulation
length(which(time_data_15$tick == tmax))
temp_not_ext <- time_data_15[which(time_data_15$tick == tmax),]

#remove the ones that did not reach end of simulation from our analysis?
temp_rem <- unique(temp_not_ext[,c(1,2)])

exts <- which(time_data_15$num_individuals_species1 < 1 | time_data_15$num_individuals_species2 < 1)
time_data_sub <- time_data_15[exts,]
temp_ext <- unique(time_data_sub[,c(1,2)])
print(nrow(temp_ext))
temp2 <- cbind(temp_ext, time_data_sub[rownames(temp_ext),]$tick)
colnames(temp2) <- c("rep","mes","tick_extinct")

temp_all <- unique(time_data_15[,c(1,2)])
temp_cox2 <- anti_join(temp_all, temp_ext)
print(nrow(temp_cox2))
temp_cox <- anti_join(temp_cox2, temp_rem)
tick_extinct <- rep(500002, nrow(temp_cox))
temp3 <- cbind(temp_cox, tick_extinct)

time_ext_temp <- rbind(temp2, temp3)
time_ext <-  merge(time_ext_temp, mk15, by = c('rep', 'mes'))

#make bins for heatmaps
b <- 20 #binsize

time_ext_binned <- time_ext %>% mutate(mk1_bin = ntile(mut_sp1, n=b), mk2_bin = ntile(mut_sp2, n=b))

#write.table(time_ext_binned, file = "~/Documents/SLiM/Rstuff/time_extinct_16.csv", append = F, sep = ",")
#read previously saved file rather than having to do all the analysis again
time_ext_binned <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_extinct_16.csv", header = T))

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

#write.table(prop_ext, file = "~/Documents/SLiM/Rstuff/proportion_extinct_2.csv", append = F, sep = ",")
#prop_ext <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/proportion_extinct_new.csv", header = T))
frac_extinct <- prop_ext$extinct/prop_ext$total
frac_end <- prop_ext$end/prop_ext$total
prop_ext <- cbind(prop_ext, frac_extinct, frac_end)

p_h1 <- ggplot(data = prop_ext, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = frac_extinct)) +
  theme_bw()+
  geom_tile() +
  #xlim(10,20)+
  #ylim(10,20)+
  #scale_fill_viridis_c(option = "viridis", direction = -1, begin = 0.4856479, end = 1) #begin = 0.4856479 for mk15
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h1

p_h1me <- ggplot(data = prop_ext, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = frac_end)) +
  theme_bw()+
  geom_tile() +
  #xlim(10,20)+
  #ylim(10,20)+
  #scale_fill_viridis_c(option = "viridis", direction = -1, begin = 0.4856479, end = 1) #begin = 0.4856479 for mk15
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h1me

#ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig5a_exts15.pdf", width=4.5, height=3.75, units = "in")

#analysing if time to ext/ESS is mutation kernel dependant (see below for this same data for mk 16)

uniques <- unique(time_data_15[,c(1,2)], fromLast = T)
rns <- as.numeric(rownames(uniques))
rns <- rns - 9 #since simulation ends 1000 generations after conditions met
time_data_end_temp <- time_data_15[rns,]

time_data_end <- merge(time_data_end_temp, time_ext_binned[,c(1,2,12,13)], by = c('rep', 'mes'))

#proportion of populations extinct
prop_tick <- data_frame(mutation_kernel_1=numeric(),
                        mutation_kernel_2=numeric(),
                        total = integer(),
                        tick = integer())

mk_1 <- 1
mk_2 <- 1
counter <- 1

for (mk_1 in 1:b) {
  for (mk_2 in 1:b){
    ext_sub <-  subset(time_data_end, mk1_bin == mk_1 & mk2_bin == mk_2)
    nt <- nrow(ext_sub)
    if(nt > 0){
      n_ext <- median(ext_sub$tick)
      prop_tick[counter,1] <- mk_1
      prop_tick[counter,2] <- mk_2
      prop_tick[counter,3] <- nt
      prop_tick[counter,4] <- n_ext
      counter <- counter + 1
    }
  }
}
#write.table(prop_tick, file = "~/Documents/SLiM/Rstuff/median_tick_2.csv", append = F, sep = ",")
p_h2 <- ggplot(data = prop_tick, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = tick)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h2
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig5a_t_ext16.pdf", width=4.5, height=3.75, units = "in")


#do the same for lower prop of ben.muts

time_data_16 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_exp_16_ext2.csv", header = T))
mk16 <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutkern_16_ext2.csv", header = T))

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

#continue with previous code from here

uniques <- unique(time_data_16[,c(1,2)], fromLast = T)
rns <- as.numeric(rownames(uniques))
rns <- rns - 9
time_data_end_temp <- time_data_16[rns,]

#dealing with mutation effect size simulations (most sims did not end)

#calculating extinction for all the data
tmax <- max(time_data_15$tick)

#check for which ones did not end simulation
length(which(time_data_15$tick == tmax))
temp_not_ext <- time_data_15[which(time_data_15$tick == tmax),]

#remove the ones that did not reach end of simulation from our analysis?
temp_rem <- unique(temp_not_ext[,c(1,2)])
end <- rep(F, nrow(temp_rem))
tick_extinct <- rep(500002, nrow(temp_rem))
temp_rem2 <- cbind(temp_rem, tick_extinct, end) #most sims did not end, so add column denoting the ones that did

exts <- which(time_data_15$num_individuals_species1 < 1 | time_data_15$num_individuals_species2 < 1)
time_data_sub <- time_data_15[exts,]
temp_ext <- unique(time_data_sub[,c(1,2)])
print(nrow(temp_ext))
temp2 <- cbind(temp_ext, time_data_sub[rownames(temp_ext),]$tick, rep(T, nrow(temp_ext)))
colnames(temp2) <- c("rep","mes","tick_extinct", "end")

temp_all <- unique(time_data_15[,c(1,2)])
temp_cox2 <- anti_join(temp_all, temp_ext)
print(nrow(temp_cox2))
temp_cox <- anti_join(temp_cox2, temp_rem)
tick_extinct <- rep(500002, nrow(temp_cox))
end <- rep(T, nrow(temp_cox))
temp3 <- cbind(temp_cox, tick_extinct, end)

time_ext_temp <- rbind(temp2, temp3, temp_rem2)
time_ext <-  merge(time_ext_temp, mk15, by = c('rep', 'mes'))

write.table(time_ext, file = "~/Documents/SLiM/Rstuff/time_extinct_me_new.csv", append = F, sep = ",")
#continue previous code from here

#supplementary figure 9
time_data_paired_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_expand.csv", header = T)) #deleting this file (log_comb_init0.csv) because this analysis did not work, combined files still on cluster for reference

mks <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutation_kernel_expanded.csv", header = T))
mks_actual <- mks
mks_list <- mks_actual$mutation_kernel

#redo the extinction and growth rate plots from above
#calculating extinction for all the data
tmax <- max(time_data_paired_full$tick)
exts <- which(time_data_paired_full$num_individuals_species1 < 1 | time_data_paired_full$num_individuals_species2 < 1)
time_data_sub <- time_data_paired_full[exts,]
temp_ext <- unique(time_data_sub[,c(1,2,3,4)])
print(nrow(temp_ext))
temp2 <- cbind(temp_ext, time_data_sub[rownames(temp_ext),]$tick)
colnames(temp2) <- c("mutation_kernel_1","mutation_kernel_2","rep","mes","tick_extinct")

temp_all <- unique(time_data_paired_full[,c(1,2,3,4)])
temp_cox <- anti_join(temp_all, temp_ext)
print(nrow(temp_cox))
tick_extinct <- rep(tmax, nrow(temp_cox))
temp3 <- cbind(temp_cox, tick_extinct)

time_ext <- rbind(temp2, temp3)

write.table(time_ext, file = "~/Documents/SLiM/Rstuff/time_extinct_fig5_sup.csv", append = F, sep = ",")
time_ext <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_extinct_fig5_sup.csv", header = T))
time_ext_large <- merge(time_ext, mks[,c(1,2,6,7,8)], by.x = 'mutation_kernel_1', by.y = 'mutation_kernel')
time_ext_large <- merge(time_ext_large, mks[,c(1,2,6,7,8)], by.x = 'mutation_kernel_2', by.y = 'mutation_kernel')

#proportion of populations extinct (calculate separately for each )
prop_ext <- time_ext_large %>%
  group_by(mutation_kernel_1, mutation_kernel_2) %>%
  summarise(
    total = n(),
    extinct = sum(tick_extinct < 50000),
    .groups = "drop")

#group data by mutation kernel parameters
m1s <- unique(mks$m1) #supply of ben.muts
mes <- unique(mks$mutr) #effect size
muts <- unique(mks$mut) #mutation rates

prop_ext_m1 <- time_ext_large %>%
  group_by(m1.x, m1.y) %>%
  summarise(
    total = n(),
    extinct = sum(tick_extinct < 50000),
    .groups = "drop")

prop_ext_me <- time_ext_large %>%
  group_by(mutr.x, mutr.y) %>%
  summarise(
    total = n(),
    extinct = sum(tick_extinct < 50000),
    .groups = "drop")

prop_ext_mut <- time_ext_large %>%
  group_by(mut.x, mut.y) %>%
  summarise(
    total = n(),
    extinct = sum(tick_extinct < 50000),
    .groups = "drop")

p_h1 <- ggplot(data = prop_ext_m1, aes(x = as.factor(m1.x), y = as.factor(m1.y), fill = extinct/total)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h1
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig5_sup_m1.pdf", width=3.5, height=2.5, units = "in")

p_h2 <- ggplot(data = prop_ext_me, aes(x = as.factor(mutr.x), y = as.factor(mutr.y), fill = extinct/total)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h2
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig5_sup_me.pdf", width=3.5, height=2.5, units = "in")

p_h3 <- ggplot(data = prop_ext_mut, aes(x = as.factor(mut.x), y = as.factor(mut.y), fill = extinct/total)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h3
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig5_sup_mut.pdf", width=3.5, height=2.5, units = "in")

p_hfull <- ggplot(data = prop_ext, aes(x = as.factor(mutation_kernel_1), y = as.factor(mutation_kernel_2), fill = extinct/total)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_hfull
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig5_sup_all.pdf", width=8.75, height=8.00, units = "in")

#------- figure 5 in paper -------
temp_pred_coexist <- read.table(file = "/Users/kasturilele/Documents/SLiM/final_params/pairparams_new.csv", sep = ",", header = TRUE)

#tradeoff functions
fun_exp1 <- function(b) {-0.0000206 * exp(3.2*b)} 
fun_exp2 <- function(b) {-0.0000108 * exp(3.5*b)}

#1sp ESS
r1e <- 0.3125 
r2e <- 0.28571429 

#2spESS
r1_2sp <- 0.498344
r2_2sp <- 0.362690 

#alphas from tradeoff function
a11e <- fun_exp1(r1e)
a22e <- fun_exp2(r2e)
a11_2sp <- fun_exp1(r1_2sp)
a22_2sp <- fun_exp2(r2_2sp)

#distance function
fun_minim <- function(x, y, px, py){
  sqrt((x - px)*(x - px)/0.01 + (y - py)*(y - py)/1e-10)
}


#niche overlap and fitness difference for all combinations of initial and ESS parameters
coexist_all_combs <- data.frame(IDs = integer(),
                                rho_00=double(),
                                f2_f1_00=double(),
                                rho_01=double(),
                                f2_f1_01=double(),
                                rho_10=double(),
                                f2_f1_10=double())

i <- 1

for(i in 1:nrow(temp_pred_coexist)){
  
  r1_temp <- temp_pred_coexist$r1[i]
  r2_temp <- temp_pred_coexist$r2[i]
  a11_temp <- temp_pred_coexist$a11[i]
  a22_temp <- temp_pred_coexist$a22[i]
  coexist_all_combs[i,1] <- i-1
  
  #calculate predictions of coexistence
  temp_rho <- (-2.5e-05/a11_temp)*(-2.5e-05/a22_temp)
  temp1_f2f1 <- (a11_temp/a22_temp)*(-2.5e-05/-2.5e-05)
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(r2_temp/r1_temp)
  rho <- sqrt(abs(temp_rho))
  f2f1 <- (1/temp_f2f1)
  
  coexist_all_combs[i,2] <- rho 
  coexist_all_combs[i,3] <- f2f1
  
  #coexistence between species 1 ESS and species 2 params
  temp_rho <- (-2.5e-05/a11_2sp)*(-2.5e-05/a22_temp)
  temp1_f2f1 <- (a11_2sp/a22_temp)*(-2.5e-05/-2.5e-05)
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(r2_temp/r1_2sp)
  rho <- sqrt(abs(temp_rho))
  f2f1 <- (1/temp_f2f1)
  
  coexist_all_combs[i,4] <- rho 
  coexist_all_combs[i,5] <- f2f1
  
  #coexistence between species 1 params and species 2 ESS
  temp_rho <- (-2.5e-05/a11_temp)*(-2.5e-05/a22_2sp)
  temp1_f2f1 <- (a11_temp/a22_2sp)*(-2.5e-05/-2.5e-05)
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(r2_2sp/r1_temp)
  rho <- sqrt(abs(temp_rho))
  f2f1 <- (1/temp_f2f1)
  
  coexist_all_combs[i,6] <- rho 
  coexist_all_combs[i,7] <- f2f1
}

#read the simulation data for these populations
time_data_00 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_init_00.csv", header = T))
time_data <- subset(time_data_00, time_data_00$tick < 50102) #subset to make comparable with finite time data from figure 4

#did not use this one in paper
#time_data_12 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_init_12.csv", header = T))
#time_data <- time_data_12

#calculating extinction for all the data
tmax <- max(time_data$tick)

exts <- which(time_data$num_individuals_species1 < 1 | time_data$num_individuals_species2 < 1)
time_data_sub <- time_data[exts,]
temp_ext <- unique(time_data_sub[,c(1,2)])
print(nrow(temp_ext))
temp2 <- cbind(temp_ext, time_data_sub[rownames(temp_ext),]$tick)
colnames(temp2) <- c("rep","mes","tick_extinct")

temp_all <- unique(time_data[,c(1,2)])
temp_cox2 <- anti_join(temp_all, temp_ext)
print(nrow(temp_cox2))
tick_extinct <- rep(50002, nrow(temp_cox2))
temp3 <- cbind(temp_cox2, tick_extinct)

time_ext_temp <- rbind(temp2, temp3)
#write.table(time_ext_temp, file = "~/Documents/SLiM/Rstuff/time_extinct_fig6_12.csv", append = F, sep = ",")

time_ext_temp <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_extinct_fig6_00_50000.csv", header = T))

prop_ext <- data_frame(IDs=numeric(),
                       total = integer(),
                       extinct=integer())
counter <- 1

for (r in 0:99) {
  ext_sub <-  subset(time_ext_temp, rep == r)
  nt <- nrow(ext_sub)
  if(nt > 0){
    n_ext <- length(which(ext_sub$tick_extinct < 50000))
    prop_ext[counter,1] <- r
    prop_ext[counter,2] <- nt
    prop_ext[counter,3] <- n_ext
    counter <- counter + 1
  }
}

frac_extinct <- prop_ext$extinct/prop_ext$total

ext_none <- prop_ext$extinct == 0

prop_ext <- cbind(prop_ext, frac_extinct, ext_none)


coexist_00s <- merge(coexist_all_combs, prop_ext, by = 'IDs')

#plot these after because the colours are based on extinction fraction from previous data
fun_default <- function(x) {x}
fun_inverse <- function(x) {1/x}
fun_vert <- function(x,y) {y <- 1}

#coexistence between both ESS (for plot)
temp_rho <- (-2.5e-05/a11_2sp)*(-2.5e-05/a22_2sp)
temp1_f2f1 <- (a11_2sp/a22_2sp)*(-2.5e-05/-2.5e-05)
temp_f2f1 <- sqrt(abs(temp1_f2f1))*(r2_2sp/r1_2sp)
rho_2sp <- sqrt(abs(temp_rho))
f2f1_2sp <- (1/temp_f2f1)

pctemp <- ggplot(data = coexist_00s) +
  theme_bw() +
  geom_segment(aes(x = rho_00, y = f2_f1_00, xend = rho_10, yend = f2_f1_10, colour = frac_extinct), alpha = 0.5) +
  geom_point(aes(x = rho_00, y = f2_f1_00, colour = frac_extinct), size = 0.8) +
  geom_point(aes(x = rho_10, y = f2_f1_10, colour = frac_extinct, shape = ext_none), size = 2.5, stroke = 1) +
  annotate("point",x = rho_2sp, y = f2f1_2sp, colour = "#111111", size = 2.5) +
  scale_colour_viridis_c()+
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank()) #removing the facet labels - they are not necessary
pctemp
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/fig6a_1sup_new.pdf", width=7.50, height=5.00, units = "in")

pctemp2 <- ggplot(data = coexist_00s) +
  theme_bw() +
  #geom_segment(aes(x = rho_00, y = f2_f1_00, xend = rho_10, yend = f2_f1_10, colour = frac_extinct), alpha = 0.5) +
  geom_point(aes(x = rho_00, y = f2_f1_00, colour = frac_extinct), size = 0.8) +
  geom_point(aes(x = rho_10, y = f2_f1_10, colour = frac_extinct, shape = ext_none), size = 2.5, stroke = 1) +
  annotate("point",x = rho_2sp, y = f2f1_2sp, colour = "#111111", size = 2.5) +
  #xlim(0.3,0.4)+
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_y_log10() +
  scale_colour_viridis_c()+
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank()) #removing the facet labels - they are not necessary
pctemp2
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/fig6a_1_new.pdf", width=7.50, height=5.00, units = "in")

#part 1 of figure: extinction tick for different reps
time_data_long <- merge(time_ext_temp, prop_ext, by.x="rep", by.y="IDs")

p_tt6 <- ggplot(data = time_data_long, aes(x = 1, y = tick_extinct, colour = frac_extinct, fill = frac_extinct)) +
  theme_bw()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6)+
  xlim(0.5,1.5)+
  scale_colour_viridis_c()+
  scale_fill_viridis_c()+
  facet_wrap(~rep, nrow = 10)+
  theme(legend.position = "none")
p_tt6

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/fig6b.pdf", width=7.5, height=7.5, units = "in")

#-------- supplementary figure 3:distribution of parameter estimates for ri and aii  --------
single_data <-  read.table(file = "/Users/kasturilele/Documents/community/est_single_all_6-3.txt", sep = ",", header = TRUE)
scal <- 1000
single_data$a11 <- single_data$a11*1000*scal #(the second scaling factor is because the density function y axis is weird)

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum","C.paralimentarius","S.cerevisiae", "W.anomalus","K.humilis","K.servazzii")
strain_order <- c("17B2","0092a","232","460","550","253","163","228","177")
domain <- c("bacteria","bacteria","bacteria","bacteria","bacteria","yeast","yeast","yeast","yeast")
temp1 <- cbind(strain_order, domain)
colnames(temp1) <- c("Strains", "domain")
single_data_large <- merge(single_data,temp1, by="Strains")

colour_domain <- c("yeast" = "#46327e",
                   "bacteria" = "#fde725")
#read simulation param data
simparams <-  read.table(file = "/Users/kasturilele/Documents/SLiM/final_params/pairparams_new.csv", sep = ",", header = TRUE)
simparams$a11 <- simparams$a11*scal
simparams$a22 <- simparams$a22*scal
# #order data by correct order
single_data$Strains <- factor(single_data$Strains, levels=strain_order)
r1e <- 0.498344
r2e <- 0.362690
a11e <- bt1(r1)*scal
a22e <- bt2(r2)*scal

r1p <- simparams$r1[101]
r2p <- simparams$r2[101]
a11p <- simparams$a11[101]
a22p <- simparams$a22[101]
rdist <- ggplot(data = single_data_large) +
  theme_bw() +
  geom_density(aes(r,color = domain, fill = domain), adjust = 0.25, alpha = 0.1) +
  scale_color_manual(values = colour_domain) +
  scale_fill_manual(values = colour_domain) +
  xlim(0.1,0.5)+
  labs(x = "growth rate", y = "count") + 
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")
rdist
ggsave(file="/Users/kasturilele/Documents/SLiM/plotdump/sims_mar26/sup_ris_dist.pdf",plot=rdist, width=4.5, height=2.75, units = "in")
rbox <- ggplot() +
  theme_minimal()+
  annotate("point",x = r1e, y = -1, colour = "#111111", size = 2.5) +
  annotate("point",x = r2e, y = -2, colour = "#111111", size = 2.5, shape = 17) +
  annotate("point",x = r1p, y = -1, colour = "#26828e", size = 2.5) +
  annotate("point",x = r2p, y = -2, colour = "#6ece58", size = 2.5, shape = 17) +
  geom_boxplot(data = simparams, mapping = aes(x = r1, y= -1), color = "#26828e", fill = "#26828e", alpha = 0.5)+
  geom_boxplot(data = simparams, mapping = aes(x = r2, y= -2), color = "#6ece58", fill = "#6ece58", alpha = 0.5)+
  xlim(0.1,0.5)+
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")
rbox
ggsave(file="/Users/kasturilele/Documents/SLiM/plotdump/sims_mar26/sup_ris_box.pdf",plot=rbox, width=4.5, height=1.25, units = "in")

adist <- ggplot(data = single_data_large) +
  theme_bw() +
  geom_density(aes(a11,color = domain, fill = domain),adjust = 0.25, alpha = 0.1) +
  scale_color_manual(values = colour_domain) +
  scale_fill_manual(values = colour_domain) +
  xlim(-0.16, 0)+
  labs(x = "interaction coefficient", y = "count") + 
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")
adist
ggsave(file="/Users/kasturilele/Documents/SLiM/plotdump/sims_mar26/sup_aiis_dist.pdf",plot=adist, width=4.5, height=2.75, units = "in")
abox <- ggplot()+
  theme_minimal()+
  annotate("point",x = a11e, y = -1, colour = "#111111", size = 2.5) +
  annotate("point",x = a22e, y = -2, colour = "#111111", size = 2.5, shape = 17) +
  annotate("point",x = a11p, y = -1, colour = "#26828e", size = 2.5) +
  annotate("point",x = a22p, y = -2, colour = "#6ece58", size = 2.5, shape = 17) +
  geom_boxplot(data = simparams, mapping = aes(x = a11, y= -1), color = "#26828e", fill = "#26828e", alpha = 0.5)+
  geom_boxplot(data = simparams, mapping = aes(x = a22, y= -2), color = "#6ece58", fill = "#6ece58", alpha = 0.5)+
  xlim(-0.16, 0)+
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")
abox
ggsave(file="/Users/kasturilele/Documents/SLiM/plotdump/sims_mar26/sup_aiis_box.pdf",plot=abox, width=4.5, height=1.25, units = "in")


comb <- plot_grid(a, b, labels = c('A', 'B'), label_size = 12)
comb

#--------- supplementary figure 6 - can we recpature the tradeoff function from single-species simulation?---------

#open the old single species data file
fun_exp <- function(b, C1, C2) {C1 * exp(C2*b)} #make the function based on parameters C1 and C2
fun_exp_true <- function(b) {-0.0000206 * exp(3.2*b)}
mks_actual <- mks[which(mks$mutr > 0.002),]
time_data_sub <- subset(time_data_sub, mutation_kernel %in% mks_actual$mutation_kernel)

#subset of time_data at endpoint = time_data_sub
#temp: plot of r vs a, faceted by mutation_kernel
plot <- ggplot()+
  theme_bw()+
  geom_point(data = time_data_sub, mapping = aes(x = r1,y = a11)) +
  scale_colour_viridis_c(begin = 0, end = 0.8, direction = -1) +
  facet_wrap(~mutation_kernel, nrow = 9) +
  xlim(0.1, 0.6)+
  labs(title = "growth rate vs interaction coefficient, endpoint", x = "growth rate (r)", y = "interaction coefficient (a)") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot

mks_list <- mks_actual$mutation_kernel

fits <- data_frame(mk = integer(),
                   C1 = double(),
                   C1_p = double(),
                   C2 = double(),
                   C2_p = double())

for(m in mks_list){
  time_data_sub_mk <- subset(time_data_sub, mutation_kernel == m)
  n <- nrow(fits)
  y <- time_data_sub_mk$a11
  x <- time_data_sub_mk$r1
  newdata <- data.frame(x,y)
  tryCatch({
    m1 <- nls(y ~ fun_exp(x, C1, C2), newdata, start = list(C1 = -0.01, C2 = 1), trace = T) #getting values of C1 and C2 by fitting the functional form of exponential function
    m1_sum <- summary(m1)
    
    C1t <- round(m1_sum$coefficients[1], 8)
    C2t <- round(m1_sum$coefficients[2], 2)
    fits[n+1, 1] <- m
    fits[n+1, 2] <- C1t
    fits[n+1, 3] <- m1_sum$coefficients[7]
    fits[n+1, 4] <- C2t
    fits[n+1, 5] <- m1_sum$coefficients[8]
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
mks_new <- fits$mk
#reordering list of mks
mks_ordered <- c(79,78,77,76,75,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,84,83,82,81,80,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,89,88,87,86,85,74,73,72,71,70,69,68,67,66,65,64,63,62,61,60)

pList <- list()
r1s <- seq(0.1, 0.6, length.out = 500)
#fits_long <- data.frame(mk = integer(), r1s = double(), a11s = double())
for(m in mks_ordered){
  time_data_sub_mk <- subset(time_data_sub, mutation_kernel == m)
  n = which(fits$mk == m)
  C1t <- fits$C1[n]
  C2t <- fits$C2[n]
  if(length(C1t) > 0) {
    a11s <- fun_exp(r1s, C1t, C2t)
  } else {
    a11s <- rep(0, nrow(time_data_sub_mk)) }
  
  mk <- rep(m, nrow(time_data_sub_mk))
  temp1 <- cbind(mk, r1s, a11s)
  #fits_long <- rbind(fits_long, temp1)
  ptemp <- ggplot()+
    theme_bw()+
    geom_point(data = temp1, mapping = aes(x = r1s, y = a11s), colour = "#999999")+
    stat_function(fun = fun_exp_true, inherit.aes = F, linetype = 2, colour = "#26828E") +
    geom_point(data = time_data_sub_mk, mapping = aes(x = r1,y = a11)) +
    xlim(0.1, 0.6)+
    ylim(-1e-04, -5e-05)+
    labs(x = NULL, y = NULL, title = as.character(m)) +
    theme(#axis.line = element_blank(),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      legend.position="none")
  #filen <- paste("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/tempplots/", m, "_main.pdf", sep = "")
  #ggsave(filen,plot=ptemp, width=2.20, height=2.10, units = "in")
  pList[[length(pList)+1]] <- ptemp
  
}

p_all2 <- plot_grid(
  plotlist = pList,
  ncol = 5,
  byrow = T
)
p_all2
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/tempplots/all.pdf",plot=p_all2, width=7.5, height=18, units = "in")

#looking at the r1 and a11 values
fits_large <- merge(fits, mks_actual, by.x = "mk", by.y = "mutation_kernel")
C1_actual <- rep(-0.0000206, 59)
C2_actual <- rep(3.2, 59)
fits_large <- cbind(fits_large, C1_actual, C2_actual)
fits_large$mut <- as.character(fits_large$mut)

p_fits1 <- ggplot(data = fits_large, aes(x = mut, y = C1)) +
  theme_bw()+
  geom_point(aes(x = mut, y = C1_actual), colour="#000000")+
  geom_point(colour = "#26828E", shape = 1, stroke = 1.5)+
  facet_grid(-mutr~m3)+
  theme(legend.position = "none")
p_fits1

p_fits2 <- ggplot(data = fits_large, aes(x = mut, y = C2)) +
  theme_bw()+
  geom_point(aes(x = mut, y = C2_actual), colour="#000000")+
  geom_point(colour = "#26828E", shape = 1, stroke = 1.5)+
  facet_grid(-mutr~m3)+
  theme(legend.position = "none")
p_fits2

fits_long_large <- merge(fits_long, mks_actual, by.x = "mk", by.y = "mutation_kernel")
fits_long_large$mut <- as.character(fits_long_large$mut)

p_fits3 <- ggplot(data = fits_long_large, aes(x = r1s, y = a11s, group = mut, colour=mut)) +
  theme_bw()+
  geom_line(alpha = 0.6, linewidth = 2)+
  scale_colour_manual(values = c("#26828E","#2A5B6F","#203446"))+
  facet_grid(-mutr~m3)+
  theme(legend.position = "none")
p_fits3

p_fits4 <-  ggplot() +
  theme_bw()+
  geom_line(data = fits_long_large, aes(x = r1s, y = a11s, group = mut, colour=mut), colour = "#999999")+
  #stat_function(fun = fun_exp_true, inherit.aes = F, linetype = 2, colour = "#26828E") +
  xlim(0.1, 0.6)+
  ylim(-1e-04, -5e-05)+
  facet_wrap(~mk, nrow = 12)+
  theme(legend.position = "none")
p_fits4


