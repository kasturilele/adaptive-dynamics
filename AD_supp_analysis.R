#predictions of coexistence for parameters used in adaptive dynamics model
library(ggplot2)
library(dplyr)
setwd("/Users/kasturilele/Documents/SLiM/pair_coexist")

fun_default <- function(x) {x}
fun_inverse <- function(x) {1/x}
fun_vert <- function(x,y) {y <- 1}


#------- data from sourdough microbes ------ 
paired_data <- read.table(file = "/Users/kasturilele/Documents/SLiM/params_temp.csv", sep = ",", header = TRUE)
#paired_data <- subset(paired_data, Replicate < 20)
strain_order <- c("17B2","0092a","232","460","550","253","163","228","177")

#change data frame such that alphas are multiplied by 1000
paired_data$a11 <- paired_data$a11*1000
paired_data$a12 <- paired_data$a12*1000
paired_data$a21 <- paired_data$a21*1000
paired_data$a22 <- paired_data$a22*1000

coexist <- data.frame(Strain_1=character(),
                      Strain_2=character(),
                      rho=double(),
                      f2_f1=double())

for(r in 1:20){
  S1 <- "17B2"
  S2 <- "253"
  #calculate rho and f2/f1 for all pairs
  temp_rho <- (paired_data[r,"a12"]/paired_data[r,"a11"])*(paired_data[r,"a21"]/paired_data[r,"a22"])
  temp1_f2f1 <- (paired_data[r,"a11"]/paired_data[r,"a22"])*(paired_data[r,"a12"]/paired_data[r,"a21"])
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(paired_data[r,"r2"]/paired_data[r,"r1"])
  
  coexist[r,1] <- paired_data[r, "Strain_2"]
  coexist[r,2] <- paired_data[r, "Strain_1"]
  coexist[r,3] <- sqrt(abs(temp_rho))
  coexist[r,4] <- (1/temp_f2f1)
  
}

#plot 1 - separate bacteria-bacteria and yeast-yeast competitions, superimposed on OG chesson plot
plot_chesson <- ggplot(data = coexist, mapping = aes(x = rho, y = f2_f1)) +
  theme_bw() +
  geom_point(shape = 16, stroke = 1.0, alpha = 0.5) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
plot_chesson

#r-alpha on tradeoff function
fun_exp <- function(b) {-0.000005 * exp(b)}
plot_tradeoff <- ggplot(data= paired_data, aes(x = r1, y = a11))+
  theme_bw()+
  geom_point()+
  geom_point(data= paired_data, aes(x = r2, y = a22))+
  stat_function(fun = fun_exp, inherit.aes = F, linewidth = 2, colour = "#999999")+
  xlim(0.1,0.5)
plot_tradeoff

paired_data_evolpred <- paired_data[,c(1,2,3,4,5,10,11)]
paired_data_evolpred$r1 <- rep(1, 20)
paired_data_evolpred$a11 <- rep(0.0001359141, 20)
paired_data_evolpred$r2 <- rep(1, 20)
paired_data_evolpred$a22 <- rep(0.0001359141, 20)

coexist2 <- data.frame(Strain_1=character(),
                       Strain_2=character(),
                       rho=double(),
                       f2_f1=double())

for(r in 1:20){
  S1 <- "17B2"
  S2 <- "253"
  #calculate rho and f2/f1 for all pairs
  temp_rho <- (paired_data_evolpred[r,"a12"]/paired_data_evolpred[r,"a11"])*(paired_data_evolpred[r,"a21"]/paired_data_evolpred[r,"a22"])
  temp1_f2f1 <- (paired_data_evolpred[r,"a11"]/paired_data_evolpred[r,"a22"])*(paired_data_evolpred[r,"a12"]/paired_data_evolpred[r,"a21"])
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(paired_data_evolpred[r,"r2"]/paired_data_evolpred[r,"r1"])
  
  coexist2[r,1] <- paired_data_evolpred[r, "Strain_2"]
  coexist2[r,2] <- paired_data_evolpred[r, "Strain_1"]
  coexist2[r,3] <- sqrt(abs(temp_rho))
  coexist2[r,4] <- (1/temp_f2f1)
  
}

plot_chesson2 <- ggplot(data = coexist2, mapping = aes(x = rho, y = f2_f1)) +
  theme_bw() +
  geom_point(shape = 16, stroke = 1.0, alpha = 0.5) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
plot_chesson2

#-------- trying out toy parameters that are still in the range --------
# r range: 0.1 to 0.5
# a range: -1e-04 to -1e-06

r1 <- seq(0.1, 0.5, by = 0.1)
r2 <- seq(0.1, 0.5, by = 0.1)
a11 <- -1e-06*1.071519^(seq(10,50,by=10))
a12 <- -1e-06*1.071519^(seq(10,50,by=10))
a21 <- -1e-06*1.071519^(seq(10,50,by=10))
a22 <- -1e-06*1.071519^(seq(10,50,by=10))
toy_params <- expand.grid(r1 = r1, r2 = r2, a11 = a11, a12 = a12, a21 = a21, a22 = a22)
toy_params$IDs <- rownames(toy_params)

coexist3 <- data.frame(IDs = character(),
                       rho=double(),
                       f2_f1=double(),
                       coexist=logical())
for(r in 1:15625){
  S1 <- "a"
  S2 <- "b"
  #calculate rho and f2/f1 for all pairs
  temp_rho <- (toy_params[r,"a12"]/toy_params[r,"a11"])*(toy_params[r,"a21"]/toy_params[r,"a22"])
  temp1_f2f1 <- (toy_params[r,"a11"]/toy_params[r,"a22"])*(toy_params[r,"a12"]/toy_params[r,"a21"])
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(toy_params[r,"r2"]/toy_params[r,"r1"])
  rho <- sqrt(abs(temp_rho))
  f2f1 <- (1/temp_f2f1)
  
  coexist3[r,1] <- r
  coexist3[r,2] <- rho
  coexist3[r,3] <- (1/temp_f2f1)
  if(f2f1 > rho & f2f1 < 1/rho){
    coexist3[r,4] <- T
  } else{coexist3[r,4] <- F}
}

plot_chesson3 <- ggplot(data = coexist3, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
  theme_bw() +
  geom_point(shape = 16, stroke = 1.0, alpha = 0.5) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
plot_chesson3

toy_params_coexist <- merge(toy_params, coexist3, by = 'IDs')

plot_ra <- ggplot(data = toy_params_coexist, mapping = aes(x = r1, y = a11, colour = coexist)) +
  theme_bw() +
  geom_jitter(shape = 16, stroke = 1.0, alpha = 0.5)+
  facet_grid(r2~a22)
plot_ra

tp_coexist <- toy_params_coexist[which(toy_params_coexist$coexist == T),]

#trying to understand patterns in parameters that coexist
plot_ra2 <- ggplot(data = tp_coexist, mapping = aes(x = r1, y = a11, colour = a21)) +
  theme_bw() +
  geom_jitter(shape = 16, stroke = 1.0, alpha = 0.5)+
  scale_colour_viridis_c()+
  facet_grid(r2~a22)
plot_ra2

tp_subset <- subset(toy_params_coexist, a12 == a12[1] & a21 == a21[1])
plot_ra3 <- ggplot(data = tp_subset, mapping = aes(x = r1, y = a11, colour = coexist)) +
  theme_bw() +
  geom_point(shape = 16, stroke = 5.0)+
  facet_grid(r2~a22)
plot_ra3

#proportion of params that coexist that lie above and below tradeoff function
fun_exp <- function(b) {-0.00005 * exp(b)}
plot_tradeoff <- ggplot()+
  theme_bw()+
  geom_jitter(data= tp_coexist, aes(x = r1, y = a11, colour = a21), alpha = 0.5)+
  #geom_jitter(data= tp_coexist, aes(x = r2, y = a22, colour = a12))+
  scale_colour_viridis_c() +
  stat_function(fun = fun_exp, inherit.aes = F, linewidth = 2, colour = "#999999")+
  xlim(0,0.6)
plot_tradeoff

#-------do the same as before, except limit the values of a21 and a12------
#use -7.943214e-06 for a12 and a21 as it was the second value that gave a decent amount of coexistence
#also use linear sequence for making aii and ajj, old sequence- -1e-06*1.071519^(seq(5,75,by=5))
#multiply everything by 10 because pop sizes are too large

r1 <- seq(0.2, 0.5, by = 0.05)
r2 <- seq(0.2, 0.5, by = 0.05)
a11 <- seq(-1e-5, -1.25e-4, by = -5e-6)
a12 <- -2.5e-05
a21 <- -2.5e-05
a22 <-  seq(-1e-5, -1.25e-4, by = -5e-6)
toy_params <- expand.grid(r1 = r1, r2 = r2, a11 = a11, a12 = a12, a21 = a21, a22 = a22)
toy_params$IDs <- rownames(toy_params)

toy_params <- read.table(file = "/Users/kasturilele/Documents/SLiM/final_params/paired_params_2new.csv", sep = ",", header = TRUE)
coexist4 <- data.frame(IDs = character(),
                       rho=double(),
                       f2_f1=double(),
                       coexist=logical())
for(r in 1:nrow(toy_params)){
  S1 <- "a"
  S2 <- "b"
  #calculate rho and f2/f1 for all pairs
  temp_rho <- (toy_params[r,"a12"]/toy_params[r,"a11"])*(toy_params[r,"a21"]/toy_params[r,"a22"])
  temp1_f2f1 <- (toy_params[r,"a11"]/toy_params[r,"a22"])*(toy_params[r,"a12"]/toy_params[r,"a21"])
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(toy_params[r,"r2"]/toy_params[r,"r1"])
  rho <- sqrt(abs(temp_rho))
  f2f1 <- (1/temp_f2f1)
  
  coexist4[r,1] <- r
  coexist4[r,2] <- rho
  coexist4[r,3] <- (1/temp_f2f1)
  if(f2f1 > rho & f2f1 < 1/rho){
    coexist4[r,4] <- T
  } else{coexist4[r,4] <- F}
}

plot_chesson4 <- ggplot(data = coexist4, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
  theme_bw() +
  geom_point(shape = 16, stroke = 1.0, alpha = 0.5) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
plot_chesson4

toy_params_coexist <- merge(toy_params, coexist4, by = 'IDs')

#proportion of params that coexist that lie above and below tradeoff function
fun_exp <- function(b) {-0.0000173 * exp(3.2*b)}
plot_tradeoff <- ggplot()+
  theme_bw()+
  geom_tile(data= toy_params_coexist, aes(x = r1, y = a11, fill = coexist), alpha = 0.75)+
  facet_grid(r2~-a22)+ #-a22 to change order of facets to match scale in a11
  scale_fill_viridis_d(begin = 0.3, end = 0.95, option = 'G') +
  stat_function(fun = fun_exp, inherit.aes = F, linewidth = 2, colour = "#999999")+
  xlim(0.2,0.5)
plot_tradeoff

#toy_params_subset <- subset(toy_params_coexist, a22 == -3e-05 & r2 == 0.55)
#write.csv(toy_params_subset, file = 'test_params.csv', row.names = FALSE) 

#make another subset of param values, where coexistence is true
# toy_params_2 <- toy_params_coexist[which(toy_params_coexist$coexist == TRUE),]
# rand50 <- sample(1:nrow(toy_params_2), 50)
# params_coexist <- toy_params_2[rand50,]
# 
# ggplot(data = params_coexist, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
#   theme_bw() +
#   geom_point(shape = 16, stroke = 1.0, alpha = 0.5) +
#   scale_y_log10() +
#   stat_function(fun = fun_inverse, inherit.aes = F) +
#   stat_function(fun = fun_default, inherit.aes = F) +
#   theme(axis.text = element_text(size = 16),
#         legend.position = "none",
#         axis.title = element_blank()) #removing the facet labels - they are not necessary
# ggplot()+
#   theme_bw()+
#   geom_point(data= params_coexist, aes(x = r1, y = a11, colour = r2), alpha = 0.75)+
#   #facet_grid(r2~a22)+
#   scale_colour_viridis_c() +
#   stat_function(fun = fun_exp, inherit.aes = F, linewidth = 2, colour = "#999999")+
#   xlim(0.2,0.8)
# 
# write.csv(params_coexist, file = 'test_params_coexist.csv', row.names = FALSE) 

#subset of parameter values part 2, where 
#a) coexistence is true and 
#b) distance from tradeoff function greater than some given value
#to make sure that r-a values are not causing bias
toy_params_2 <- toy_params_coexist[which(toy_params_coexist$coexist == TRUE),]

#minimization function
fun_minim <- function(x, f, px, py){
  (x - px)*(x - px)/0.01 + (f(x) - py)*(f(x) - py)/1e-10
}

dist_1 <- rep(0, nrow(toy_params_2))
dist_2 <- rep(0, nrow(toy_params_2))
below <- rep(F, nrow(toy_params_2))
for(i in 1:nrow(toy_params_2)){
  ir1 <- toy_params_2$r1[i]
  ir2 <- toy_params_2$r2[i]
  ia11 <- toy_params_2$a11[i]
  ia22 <- toy_params_2$a22[i]
  distance1 <-  optimize(
    f = function(x) fun_minim(x, f = fun_exp, px = ir1, py = ia11),
    interval = c(0, 1)
  )$objective
  distance2 <-  optimize(
    f = function(x) fun_minim(x, f = fun_exp, px = ir2, py = ia22),
    interval = c(0, 1)
  )$objective
  dist_1[i] <- distance1
  dist_2[i] <- distance2
  
  if(ia11 < fun_exp(ir1)){
    if(ia22 < fun_exp(ir2)){
      below[i] <- T
    }
  }
}

toy_params_2_sub <- toy_params_2[which(dist_1 > 1 & dist_2 > 1 & dist_1 < 3 & dist_2 < 3 & below == T),]

rand50 <- sample(1:nrow(toy_params_2_sub), 50)
params_coexist <- toy_params_2_sub[rand50,]

ggplot(data = params_coexist, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
  theme_bw() +
  geom_point(shape = 16, stroke = 1.0, alpha = 0.5) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
ggplot()+
  theme_bw()+
  geom_text(data= params_coexist, aes(x = r1, y = a11, label = IDs, colour = r2), alpha = 0.75)+
  scale_colour_viridis_c() +
  stat_function(fun = fun_exp, inherit.aes = F, linewidth = 2, colour = "#999999")+
  xlim(0.2,0.5)
ggplot()+
  theme_bw()+
  geom_text(data= params_coexist, aes(x = r2, y = a22, label = IDs, colour = r1), alpha = 0.75)+
  scale_colour_viridis_c() +
  stat_function(fun = fun_exp, inherit.aes = F, linewidth = 2, colour = "#999999")+
  xlim(0.2,0.5)

write.csv(params_coexist, file = 'test_params_coexist.csv', row.names = FALSE) 

#------- checking if predictions of coexistence match with pres/abs of coexistence in slim model (sanity check) -------
pred_coexist <- read.table(file = "/Users/kasturilele/Documents/SLiM/pair_coexist/output_test_coexist.csv", sep = ",", header = TRUE)

pred_coexist_larger <- merge(pred_coexist, toy_params_subset, by.x = "ID", by.y = "IDs")

plot_sancheck <- ggplot(data = pred_coexist_larger, mapping = aes(x = coexist.x, y = coexist.y)) +
  theme_bw()+
  geom_jitter()
plot_sancheck

#look at the individual values that didn't fall in the correct bins
pred_coexist_larger[which(pred_coexist_larger$coexist.x == "0" & pred_coexist_larger$coexist.y == TRUE),]
pred_coexist_larger[which(pred_coexist_larger$coexist.x == "1" & pred_coexist_larger$coexist.y == FALSE),]

#checking the parameters
params_coexist <- read.table(file = "/Users/kasturilele/Documents/SLiM/pair_coexist/test_params_coexist.csv", sep = ",", header = TRUE)[1:20,]

#--------- repeat the above, but check the distance from the ESS, not the tradeoff function --------
#note: jan 2026: used the paired param bit of the code to generate new parameters for testing fraction of extinctions
# feb 2026: used the single param bit to generate new parameters on the tradeoff function to test time to ESS (new simulations)
#tradeoff functions
fun_exp1 <- function(b) {-0.0000206 * exp(3.2*b)} #function(b) {-0.0000173 * exp(3.2*b)} #function(b) {-0.0000206 * exp(3.2*b)}
fun_exp2 <- function(b) {-0.0000108 * exp(3.5*b)} #function(b) {-0.0000108 * exp(4.4*b)} #function(b) {-0.0000108 * exp(3.5*b)}

#1sp ESS
r1e <- 0.3125 #0.3125 #0.3125
r2e <- 0.28571429 #0.22727273 #0.28571429

#2spESS
r1_2sp <- 0.498344 #0.429727996 #0.498343958
r2_2sp <- 0.362690 #0.341218115 #0.362689891

#alphas from tradeoff function
a11e <- fun_exp1(r1e)
a22e <- fun_exp2(r2e)
a11_2sp <- fun_exp1(r1_2sp)
a22_2sp <- fun_exp2(r2_2sp)

#distance function
fun_minim <- function(x, y, px, py){
  sqrt((x - px)*(x - px)/0.01 + (y - py)*(y - py)/1e-10)
}

#plot tradeoff functions, ESSes
ggplot()+
  theme_bw()+
  stat_function(fun = fun_exp1, inherit.aes = F, linewidth = 2, colour = "#A0BBBB")+
  stat_function(fun = fun_exp2, inherit.aes = F, linewidth = 2, colour = "#BBA0AD")+
  geom_point(aes(x = r1e, y = a11e), colour = "#1E2A2A")+ 
  geom_point(aes(x = r2e, y = a22e), colour = "#45313B")+ 
  geom_point(aes(x = r1_2sp, y = a11_2sp), colour = "#1E2A2A", shape = 17, size = 2)+
  geom_point(aes(x = r2_2sp, y = a22_2sp), colour = "#45313B", shape = 17, size = 2)+ 
  
  xlim(0.1,0.6)

#distance between 1sp and 2sp ESS
dist1 <- fun_minim(r1e, a11e, r1_2sp, a11_2sp)
dist2 <- fun_minim(r2e, a22e, r2_2sp, a22_2sp)

#constraints for parameter choice: between dist 1 and 5, r < re
temp_params_single1 <- data.frame(IDs = integer(),
                                  N10 = integer(),
                                  r1 = double(),
                                  a11 = double(),
                                  dist = double())
temp_params_single2 <- data.frame(IDs = integer(),
                                  N10 = integer(),
                                  r1 = double(),
                                  a11 = double(),
                                  dist = double())

#for picking the single params
i <- 1
while (i < 8) {
  r1_temp <- runif(1, 0.1, r1e)
  r2_temp <- runif(1, 0.1, r2e)
  
  #a11_temp <- -runif(1, -a11e, 0.0001)
  #a22_temp <- -runif(1, -a22e, 0.0001)
  a11_temp <- fun_exp1(r1_temp)
  a22_temp <- fun_exp2(r2_temp)
  
  d1t <- fun_minim(r1e, a11e, r1_temp, a11_temp)
  d2t <- fun_minim(r2e, a22e, r2_temp, a22_temp)
  
  if(d1t > 1 & d2t > 1 & d1t < 5 & d2t < 5){
    temp_params_single1[i,1] <- i-1
    temp_params_single2[i,1] <- i-1
    temp_params_single1[i,2] <- 50
    temp_params_single2[i,2] <- 50
    temp_params_single1[i,3] <- r1_temp
    temp_params_single2[i,3] <- r2_temp
    temp_params_single1[i,4] <- a11_temp
    temp_params_single2[i,4] <- a22_temp
    temp_params_single1[i,5] <- d1t
    temp_params_single2[i,5] <- d2t
    
    i <- i+1
    
    print(paste(r1_temp, r2_temp, d1t, d2t, sep = ","))
  }
}

ggplot()+
  theme_bw()+
  stat_function(fun = fun_exp1, inherit.aes = F, linewidth = 2, colour = "#A0BBBB")+
  stat_function(fun = fun_exp2, inherit.aes = F, linewidth = 2, colour = "#BBA0AD")+
  geom_point(aes(x = r1e, y = a11e), colour = "#1E2A2A")+ 
  geom_point(aes(x = r2e, y = a22e), colour = "#45313B")+ 
  geom_point(aes(x = r1_2sp, y = a11_2sp), colour = "#1E2A2A", shape = 17, size = 2)+
  geom_point(aes(x = r2_2sp, y = a22_2sp), colour = "#45313B", shape = 17, size = 2)+ 
  
  geom_text(data= temp_params_single1, aes(x = r1, y = a11, label = IDs), alpha = 0.75, colour = "#109393")+
  geom_text(data= temp_params_single2, aes(x = r1, y = a11, label = IDs), alpha = 0.75, colour = "#EF6CAE")+
  scale_colour_viridis_c() +
  
  xlim(0.1,0.6)

#save this graph as record of parameter choice (800 px width)

#save params
write.csv(temp_params_single1, file = 'singleparams3_new.csv', row.names = FALSE) 
write.csv(temp_params_single2, file = 'singleparams4_new.csv', row.names = FALSE)

#repeat for paired parameters
temp_params_paired <- data.frame(IDs = integer(),
                                 N10 = integer(),
                                 N20 = integer(),
                                 r1 = double(),
                                 r2 = double(),
                                 a11 = double(),
                                 a22 = double(),
                                 a12 = double(),
                                 a21 = double(),
                                 dist1 = double(),
                                 dist2 = double(),
                                 coexist = logical(),
                                 coexist1 = logical(),
                                 coexist2 = logical())
#for picking the paired params
i <- 1
while (i < 101) {
  coexist_temp <- F
  coexist_temp1 <- F
  coexist_temp2 <- F
  
  r1_temp <- runif(1, 0.2, r1_2sp)
  r2_temp <- runif(1, 0.15, r2_2sp)
  
  a11_temp <- -runif(1, -a11_2sp, 0.0002)
  a22_temp <- -runif(1, -a22_2sp, 0.0002)
  
  d1t <- fun_minim(r1_2sp, a11_2sp, r1_temp, a11_temp)
  d2t <- fun_minim(r2_2sp, a22_2sp, r2_temp, a22_temp)
  
  if(d1t > 1 & d2t > 1 & d1t < 5 & d2t < 5){
    
    #calculate predictions of coexistence
    temp_rho <- (-2.5e-05/a11_temp)*(-2.5e-05/a22_temp)
    temp1_f2f1 <- (a11_temp/a22_temp)*(-2.5e-05/-2.5e-05)
    temp_f2f1 <- sqrt(abs(temp1_f2f1))*(r2_temp/r1_temp)
    rho <- sqrt(abs(temp_rho))
    f2f1 <- (1/temp_f2f1)
    
    if(f2f1 > rho & f2f1 < 1/rho){
      coexist_temp <- T
    } else {coexist_temp <- F}
    
    #coexistence between species 1 ESS and species 2 params
    temp_rho <- (-2.5e-05/a11_2sp)*(-2.5e-05/a22_temp)
    temp1_f2f1 <- (a11_2sp/a22_temp)*(-2.5e-05/-2.5e-05)
    temp_f2f1 <- sqrt(abs(temp1_f2f1))*(r2_temp/r1_2sp)
    rho <- sqrt(abs(temp_rho))
    f2f1 <- (1/temp_f2f1)
    
    if(f2f1 > rho & f2f1 < 1/rho){
      coexist_temp1 <- T
    } else {coexist_temp1 <- F}
    
    #coexistence between species 1 params and species 2 ESS
    temp_rho <- (-2.5e-05/a11_temp)*(-2.5e-05/a22_2sp)
    temp1_f2f1 <- (a11_temp/a22_2sp)*(-2.5e-05/-2.5e-05)
    temp_f2f1 <- sqrt(abs(temp1_f2f1))*(r2_2sp/r1_temp)
    rho <- sqrt(abs(temp_rho))
    f2f1 <- (1/temp_f2f1)
    
    if(f2f1 > rho & f2f1 < 1/rho){
      coexist_temp2 <- T
    } else {coexist_temp2 <- F}
    
    temp_params_paired[i,1] <- i-1
    temp_params_paired[i,2] <- 50
    temp_params_paired[i,3] <- 50
    temp_params_paired[i,4] <- r1_temp
    temp_params_paired[i,5] <- r2_temp
    temp_params_paired[i,6] <- a11_temp
    temp_params_paired[i,7] <- a22_temp
    temp_params_paired[i,8] <- -2.5e-05
    temp_params_paired[i,9] <- -2.5e-05
    temp_params_paired[i,10] <-  d1t
    temp_params_paired[i,11] <-  d2t
    temp_params_paired[i,12] <-  coexist_temp
    temp_params_paired[i,13] <-  coexist_temp1
    temp_params_paired[i,14] <-  coexist_temp2
    
    i <- i+1
    
    print(paste(r1_temp, r2_temp, d1t, d2t, sep = ","))
  }
}

ggplot()+
  theme_bw()+
  stat_function(fun = fun_exp1, inherit.aes = F, linewidth = 2, colour = "#A0BBBB")+
  stat_function(fun = fun_exp2, inherit.aes = F, linewidth = 2, colour = "#BBA0AD")+
  geom_point(aes(x = r1e, y = a11e), colour = "#1E2A2A")+ 
  geom_point(aes(x = r2e, y = a22e), colour = "#45313B")+ 
  geom_point(aes(x = r1_2sp, y = a11_2sp), colour = "#1E2A2A", shape = 17, size = 2)+
  geom_point(aes(x = r2_2sp, y = a22_2sp), colour = "#45313B", shape = 17, size = 2)+ 
  
  geom_text(data= temp_params_paired, aes(x = r1, y = a11, label = IDs), alpha = 0.75, colour = "#109393")+
  geom_text(data= temp_params_paired, aes(x = r2, y = a22, label = IDs), alpha = 0.75, colour = "#EF6CAE")+
  scale_colour_viridis_c() +
  
  xlim(0.1,0.6)
#save this graph as record of parameter choice (800 px width)

write.csv(temp_params_paired, file = 'pairparams_new.csv', row.names = FALSE) 

#--------- plot showing how aij choice alters coexistence -----------

tdoffs_aijs <- read.table(file = "/Users/kasturilele/Documents/SLiM/tradeoff_fun_ESS_aijs.csv", sep = ",", header = TRUE)[35:49,]
r1s <- unique(tdoffs_aijs$r1_sing)
r2s <- unique(tdoffs_aijs$r2_sing)
a11s <- unique(tdoffs_aijs$a1_sing)
a22s <- unique(tdoffs_aijs$a2_sing)
r1p <- unique(tdoffs_aijs$r1)
r2p <- unique(tdoffs_aijs$r2)
a11p <- unique(tdoffs_aijs$a11)
a22p <- unique(tdoffs_aijs$a22)

r1vec <- seq(r1s, r1p, 0.01)
r2vec <- seq(r2s, r2p, 0.01)

fun_exp1 <- function(b) {-0.0000206 * exp(3.2*b)} 
fun_exp2 <- function(b) {-0.0000108 * exp(3.5*b)}

a11_vec <- fun_exp1(r1vec)
a22_vec <- fun_exp2(r2vec)

coexist6 <- data.frame(r1 = double(),
                       r2 = double(),
                       a11 = double(),
                       a22 = double(),
                       a12 = double(),
                       a21 = double(),
                       rho=double(),
                       f2_f1=double(),
                       coexist=logical())
r1 <- 1
r2 <- 1
a <- 1
for(a in 1:nrow(tdoffs_aijs)){
  
  #calculate rho and f2/f1 for all pairs
  temp_rho <- (toy_params[r,"a12"]/toy_params[r,"a11"])*(toy_params[r,"a21"]/toy_params[r,"a22"])
  temp1_f2f1 <- (toy_params[r,"a11"]/toy_params[r,"a22"])*(toy_params[r,"a12"]/toy_params[r,"a21"])
  temp_f2f1 <- sqrt(abs(temp1_f2f1))*(toy_params[r,"r2"]/toy_params[r,"r1"])
  rho <- sqrt(abs(temp_rho))
  f2f1 <- (1/temp_f2f1)
  
  coexist4[r,1] <- r
  coexist4[r,2] <- rho
  coexist4[r,3] <- (1/temp_f2f1)
  if(f2f1 > rho & f2f1 < 1/rho){
    coexist4[r,4] <- T
  } else{coexist4[r,4] <- F}
}

plot_chesson4 <- ggplot(data = coexist4, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
  theme_bw() +
  geom_point(shape = 16, stroke = 1.0, alpha = 0.5) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
plot_chesson4

#------- supplement plot showing how predictions of coexistence ------

temp_pred_coexist <- read.table(file = "/Users/kasturilele/Documents/SLiM/final_params/pairparams_new.csv", sep = ",", header = TRUE)

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

plot_chesson5 <- ggplot(data = coexist_all_combs, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
  theme_bw() +
  geom_text(aes(x = rho_00, y = f2_f1_00, label = IDs), shape = 16, stroke = 1.0, alpha = 0.5, colour = "red") +
  geom_text(aes(x = rho_01, y = f2_f1_01, label = IDs), shape = 16, stroke = 1.0, alpha = 0.5, colour = "blue") +
  geom_text(aes(x = rho_10, y = f2_f1_10, label = IDs), shape = 16, stroke = 1.0, alpha = 0.5, colour = "green") +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
plot_chesson5

pctemp <- ggplot(data = coexist_all_combs, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
  theme_bw() +
  geom_segment(aes(x = rho_00, y = f2_f1_00, xend = rho_01, yend = f2_f1_01, colour = IDs)) +
  geom_point(aes(x = rho_00, y = f2_f1_00, colour = IDs)) +
  geom_point(aes(x = rho_01, y = f2_f1_01, colour = IDs)) +
  scale_colour_viridis_c()+
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
pctemp

pctemp2 <- ggplot(data = coexist_all_combs, mapping = aes(x = rho, y = f2_f1, colour = coexist)) +
  theme_bw() +
  geom_segment(aes(x = rho_00, y = f2_f1_00, xend = rho_10, yend = f2_f1_10),colour = "red") +
  geom_point(aes(x = rho_00, y = f2_f1_00), shape = 16, stroke = 1.0, alpha = 0.5, colour = "red") +
  geom_point(aes(x = rho_10, y = f2_f1_10), shape = 16, stroke = 1.0, alpha = 0.5, colour = "black") +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
pctemp2


#collecting useful analyses for both single species and paired species
#(mut_info.R and pred_coexist.R)

setwd("/Users/kasturilele/Documents/SLiM")
library(ggplot2)
library(dplyr)
library(viridis)
library(ggnewscale)
library(cowplot)

#testing r vs a SLiM code
time_data_full <- as.data.frame(read.csv("~/Documents/SLiM/output2/log_combined_single3.csv", header = T)) #use output2 for the most recent files
mks <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutation_kernel_extra.csv", header = T))
reps <- as.data.frame(read.csv("~/Documents/SLiM/final_params/single_params_3.csv", header = T)) #don't use this for the subset with no parameters

#functions: #function(b) {-0.0000173 * exp(3.2*b)} #function(b) {-0.0000108 * exp(4.4*b)} #function(b) {-0.0000206 * exp(3.2*b)} #function(b) {-0.0000108 * exp(3.5*b)}
fun_exp <- function(b) {-0.0000206 * exp(3.2*b)}
fun_K <- function(r) {-r/fun_exp(r)}
#distance function for thresholding
fun_minim <- function(x, y, px, py){
  sqrt((x - px)*(x - px)/0.01 + (y - py)*(y - py)/1e-10)
}

#equilibrium values of r and alpha (from python)
r1e <- 0.3125  #0.3125 #0.22727273 #0.3125 #0.28571429
a11e <- -5.59966057e-05 #-4.7116885e-05 #-2.93574437e-05 #-5.59966057e-05 #-2.93574437e-05

time_data <- subset(time_data_full, rep == 2)
tmax <- max(time_data$tick)
time_data_sub <- subset(time_data, tick == tmax)

#-------- growth curve plots ---------
plot5 <- ggplot()+
  theme_bw()+
  stat_function(fun = fun_exp, inherit.aes = F, linewidth = 2, colour = "#999999") +
  geom_line(data = time_data, mapping = aes(x = r1,y = a11, colour = tick, group = mes), alpha = 0.5) + #without correlation
  geom_point(data = time_data_sub, mapping = aes(x = r1,y = a11, colour = tick, group = mes)) +
  scale_colour_viridis_c(begin = 0, end = 0.8, direction = -1) +
  facet_wrap(~mutation_kernel, nrow = 5) +
  xlim(0.1, 0.6)+
  labs(title = "growth rate vs interaction coefficient over time", subtitle = "(for different starting conditions)", x = "growth rate (r)", y = "interaction coefficient (a)") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot5 
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_nov25/rvsa_4_5.pdf", width=12.00, height=8.00, units = "in")

plot6 <- ggplot()+
  theme_bw()+
  stat_function(fun = fun_K, inherit.aes = F, linewidth = 2, colour = "#999999") +
  geom_line(data = time_data, mapping = aes(x = r1, y = num_individuals_species1,colour = tick, group = mes), alpha = 0.1) + #without correlation
  geom_point(data = time_data_sub, mapping = aes(x = r1,y =num_individuals_species1,colour = tick, group = mes)) +
  scale_colour_viridis_c(begin = 0, end = 0.8, direction = -1) +
  facet_wrap(~mutation_kernel, nrow = 3) +
  xlim(0.1, 0.6)+
  labs(title = "growth rate vs population size over time", subtitle = "(for different starting conditions)", x = "growth rate (r)", y = "population size (K)") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot6 
#pdf: 10X11 in

#making growth curve plot for 
plotgr1 <- ggplot()+
  theme_bw()+
  geom_line(data = tdp_large1, mapping = aes(x = tick,y = num_individuals_species1, group = mes),colour = "#92d2bf",  alpha = 0.5) + #without correlation
  geom_line(data = tdp_large1, mapping = aes(x = tick,y = num_individuals_species2, group = mes),colour = "#313668", alpha = 0.5) + #without correlation
  #geom_point(data = time_data_subp, mapping = aes(x = r1,y = num_individuals_species1, colour = tick, group = mes)) +
  facet_grid(m1~mutr) +
  labs(title = "population size over time (mut = 1e-07)", x = "time", y = "population size") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotgr1

#--------- distance plots ----------
#what timepoint alpha reaches epsilon neighbourhood of ESS
#also calculate the distance from the ESS

#minimization function #we need to find x (on curve) that minimizes (x-px)^2 + (f(x)-py)^2
#fun_minim <- function(x, f, px, py){
#  (x - px)*(x - px)/0.01 + (f(x) - py)*(f(x) - py)/1e-10
#}

#(in the loop below, for distance)
# distance <- optimize(
#   f = function(x) fun_minim(x, f = fun_exp, px = r2_end, py = a22_end),
#   interval = c(0, 1)
# )$objective #save the minimum distance between the tradeoff function and the r/a values at the last timepoint

time_threshold <- data_frame(mutation_kernel=character(),
                             rep=integer(),
                             mes=integer(),
                             time=integer(),
                             distance=double())
mk <- 0
r <- 2
m <- 1
counter <- 1

for (mk in 0:74) {
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

tt_combine <- merge(time_threshold, mks, by = "mutation_kernel")
#tt_combine_large <- merge(tt_combine, reps, by.x = "rep", by.y = "IDs")
tt_combine_large <- tt_combine
tt_combine_large$mut <- as.factor(tt_combine_large$mut)

#plot time taken to reach threshold
p_tt <- ggplot(data = tt_combine_large, aes(x = mut, y = time, colour = mut, fill = mut)) +
  theme_bw()+
  geom_boxplot(alpha = 0.6)+
  #geom_point()+
  #scale_y_log10()+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(m1~mutr)
p_tt

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_nov25/single_tt_4_3.pdf", width=7.42, height=6.2, units = "in")

#plot distance between endpoint and tradeoff curve
p_dist <-ggplot(data = tt_combine_large, aes(x = mut, y = distance, colour = mut, fill = mut)) +
  theme_bw()+
  geom_boxplot(alpha = 0.6)+
  scale_y_log10()+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(m1~mutr)
p_dist
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_nov25/single_dist_4_3.pdf", width=7.42, height=6.2, units = "in")

#--------new analysis for paper -------

#redoing the distance and time analyses for parameter combos that started on the tradeoff curve
time_data_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_combined_single3_extra.csv", header = T)) #log_combined_single3_new.csv has old data, use log_combined_single3_extra.csv
mks <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutation_kernel_single3.csv", header = T))
reps <- as.data.frame(read.csv("~/Documents/SLiM/final_params/singleparams3_new.csv", header = T)) 

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

tmax <- max(time_data_full)
mindist_all <-  fun_minim(r1e, a11e, time_data_full$r1, time_data_full$a11)
mins <- which(mindist_all < 0.5)
time_data_mins <- time_data_full[mins,]
temp_mins <- unique(time_data_mins[,c(1,2,3)])
temp_mins_2 <- cbind(temp_mins, time_data_mins[rownames(temp_mins),]$tick,time_data_mins[rownames(temp_mins),]$r1,time_data_mins[rownames(temp_mins),]$a11)
colnames(temp_mins_2) <- c("mutation_kernel","rep","mes","tick_ESS","r1","a11") #tick at which populations reached ESS, or mindist < 0.5

temp_all <- unique(time_data_full[,c(1,2,3)])
temp_other <- anti_join(temp_all, temp_mins)
tick_ESS <- rep(tmax, nrow(temp_other))
r1 <- time_data_full[rownames(temp_other),]$r1
a11 <- time_data_full[rownames(temp_other),]$a11
temp_other_2 <- cbind(temp_other, tick_ESS, r1, a11)

#all of the data
mins_data <- rbind(temp_mins_2,temp_other_2)
mins_combine <- merge(mins_data, mks, by = "mutation_kernel")
mins_combine$rep <- as.factor(mins_combine$rep)
mins_combine$mut <- as.factor(mins_combine$mut)

write.table(mins_combine, file = "~/Documents/SLiM/Rstuff/time_ESS_single3_add.csv", append = F, sep = ",")

# p_tt <- list()
# 
# for(r in 0:6){
#   mins_subset <- subset(mins_combine, rep == r)
#   p_ttsub <- ggplot(data = mins_subset, aes(x = mut, y = tick_ESS, colour = mut, fill = mut)) +
#     theme_bw()+
#     geom_boxplot(alpha = 0.6)+
#     #geom_point()+
#     #scale_y_log10()+
#     scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
#     scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
#     facet_grid(-mutr~m3)
#   p_tt[[r+1]] <- p_ttsub
# }

# p_ttfull <- plot_grid(
#   plotlist = p_tt,
#   ncol = 4,
#   byrow = TRUE
# )
# p_ttfull 

#plot time to ESS for all reps, specific values of mutation kernels
#c(40,60,80,85)
mutkerns <- c(40,45,80,60,65,70,90,91,92)
mins_subset <- subset(mins_combine, mutation_kernel %in% mutkerns)
mins_subset$mutation_kernel <- factor(mins_subset$mutation_kernel, levels = mutkerns)
p_all <- ggplot(data = mins_subset, aes(x = rep, y = tick_ESS, colour = rep)) +
  theme_bw() +
  geom_boxplot(outlier.color = NA, whisker.color = NA, box.color = NA)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6)+
  stat_summary(fun="mean", geom="text", aes(label=after_stat(y)), vjust=-1)+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  #scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_wrap(~mutation_kernel, nrow = 3)
p_all

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/fig3_add1_new.pdf", width=9.0, height=9.0, units = "in")

p_all2 <- ggplot(data = mins_subset, aes(x = rep, y = r1, colour = rep)) +
  theme_bw() +
  geom_boxplot(outlier.color = NA, whisker.color = NA, box.color = NA)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6)+
  stat_summary(fun="mean", geom="text", aes(label=round(after_stat(y),5)), vjust=-2)+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  #scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_wrap(~mutation_kernel, nrow = 3)
p_all2

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/fig3_add2_new.pdf", width=9.0, height=9.0, units = "in")

p_all3 <- ggplot(data = mins_subset, aes(x = rep, y = a11, colour = rep)) +
  theme_bw() +
  geom_boxplot(outlier.color = NA, whisker.color = NA, box.color = NA)+
  geom_jitter(width = 0.2, height = 0, alpha = 0.6)+
  stat_summary(fun="mean", geom="text", aes(label=round(after_stat(y),8)), vjust=-2)+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  #scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_wrap(~mutation_kernel, nrow = 2)
p_all3

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_jan26/fig3_add3.pdf", width=9.0, height=6.0, units = "in")

#means for time to ESS + ri and aii value at time to ESS
mins_means <- data_frame(mutation_kernel=integer(),
                         rep=integer(),
                         tt_ESS_mean=integer(),
                         r1_mean=double(),
                         a11_mean=double())
mks_l <- mutkerns #unique(mins_data$mutation_kernel)
reps_l <- c(0,1,2,3,4,5,6)

for (mk in mks_l) {
  for(r in reps_l){
    mins_sub <- subset(mins_data, mutation_kernel == mk & rep == r)
    n <- nrow(mins_means)
    mins_means[n+1,1] <- mk
    mins_means[n+1,2] <- r
    mins_means[n+1,3] <- mean(mins_sub$tick_ESS)
    mins_means[n+1,4] <- mean(mins_sub$r1)
    mins_means[n+1,5] <- mean(mins_sub$a11)
  }
}

mins_means_combine <- merge(mins_means, mks, by = "mutation_kernel")
mins_means_combine <- merge(mins_means_combine, reps, by.x = "rep", by.y = "IDs")
write.table(mins_means_combine, file = "~/Documents/SLiM/Rstuff/time_ESS_single3_add_means.csv", append = F, sep = ",")

#--------- new analysis part 2 - can we recpature the tradeoff function from single-species simulation?---------

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



#---------- repeat above analyses on paired simulations ---------

#testing r vs a SLiM code
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

time_data_subp <- subset(time_data_paired, tick == tmax_p)
fun_exp1 <- function(b) {-0.0000206 * exp(3.2*b)} #functions: #function(b) {-0.0000173 * exp(3.2*b)} function(b) {-0.0000206 * exp(3.2*b)} 
fun_exp2 <- function(b) {-0.0000108 * exp(3.5*b)} #functions: #function(b) {-0.0000108 * exp(4.4*b)} function(b) {-0.0000108 * exp(3.5*b)}
fun_K1 <- function(r) {-r/fun_exp1(r)}
fun_K2 <- function(r) {-r/fun_exp2(r)}
fun_minim <- function(x, y, px, py){
  sqrt((x - px)*(x - px)/0.01 + (y - py)*(y - py)/1e-10)
}

plot5p <- ggplot()+
  theme_bw()+
  stat_function(fun = fun_exp1, inherit.aes = F, linewidth = 1.5, colour = "#cccccc") +
  stat_function(fun = fun_exp2, inherit.aes = F, linewidth = 1.5, colour = "#999999") +
  geom_line(data = time_data_paired, mapping = aes(x = r1,y = a11, colour = tick, group = mes), alpha = 0.5) + #without correlation
  geom_point(data = time_data_subp, mapping = aes(x = r1,y = a11, colour = tick, group = mes)) +
  scale_colour_viridis_c(begin = 0, end = 0.8, direction = -1) +
  new_scale_color()+
  geom_line(data = time_data_paired, mapping = aes(x = r2,y = a22, colour = tick, group = mes), alpha = 0.5) + #without correlation
  geom_point(data = time_data_subp, mapping = aes(x = r2,y = a22, colour = tick, group = mes)) +
  scale_colour_viridis_c(begin = 0, end = 0.8, direction = -1, option = "A") +
  facet_wrap(~mutation_kernel, nrow = 5) +
  xlim(0.1, 0.6)+
  ylim(-0.0005, 0)+
  labs(title = "growth rate vs interaction coefficient over time", subtitle = "(for different starting conditions)", x = "growth rate (r)", y = "interaction coefficient (a)") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot5p 
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_nov25/rvsa_paired_2_2.pdf", width=12.00, height=8.00, units = "in")

plot6p <- ggplot()+
  theme_bw()+
  stat_function(fun = fun_K, inherit.aes = F, linewidth = 2, colour = "#999999") +
  geom_line(data = time_data_paired, mapping = aes(x = r1,y = num_individuals_species1, colour = tick, group = mes), alpha = 0.5) + #without correlation
  geom_point(data = time_data_subp, mapping = aes(x = r1,y = num_individuals_species1, colour = tick, group = mes)) +
  scale_colour_viridis_c(begin = 0, end = 0.8, direction = -1) +
  new_scale_color()+
  geom_line(data = time_data_paired, mapping = aes(x = r2,y = num_individuals_species2, colour = tick, group = mes), alpha = 0.5) + #without correlation
  geom_point(data = time_data_subp, mapping = aes(x = r2,y = num_individuals_species2, colour = tick, group = mes)) +
  scale_colour_viridis_c(begin = 0, end = 0.8, direction = -1, option = "A") +
  facet_wrap(~mutation_kernel, nrow = 3) +
  xlim(0.1, 0.6)+
  labs(title = "growth rate vs population size over time", subtitle = "(for different starting conditions)", x = "growth rate (r)", y = "population size (K)") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot6p 
#pdf: 10X11 in

#plotting growth of both species over time
#tdp_large1 <- subset(tdp_large, mut < 2.5e-08)
#tdp_large2 <- subset(tdp_large, mut < 2.5e-07 & mut > 2.5e-08)
#tdp_large3 <- subset(tdp_large, mut > 2.5e-07)


# plotgr1 <- ggplot()+
#   theme_bw()+
#   geom_line(data = tdp_large1, mapping = aes(x = tick,y = num_individuals_species1, group = mes),colour = "#c36122",  alpha = 0.5) + #without correlation
#   geom_line(data = tdp_large1, mapping = aes(x = tick,y = num_individuals_species2, group = mes),colour = "#512c1d", alpha = 0.5) + #without correlation
#   #geom_point(data = time_data_subp, mapping = aes(x = r1,y = num_individuals_species1, colour = tick, group = mes)) +
#   facet_grid(m1~mutr) +
#   labs(title = "population size over time (mut = 1e-07)", x = "time", y = "population size") + 
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5))
# 
# plotgr1
# 
# plotgr2 <- ggplot()+
#   theme_bw()+
#   geom_line(data = tdp_large2, mapping = aes(x = tick,y = num_individuals_species1, group = mes),colour = "#c36122",  alpha = 0.5) + #without correlation
#   geom_line(data = tdp_large2, mapping = aes(x = tick,y = num_individuals_species2, group = mes),colour = "#512c1d", alpha = 0.5) + #without correlation
#   #geom_point(data = time_data_subp, mapping = aes(x = r1,y = num_individuals_species1, colour = tick, group = mes)) +
#   facet_grid(m1~mutr) +
#   labs(title = "population size over time (mut = 5e-06)", x = "time", y = "population size") + 
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5))
# 
# plotgr2

#temp extinction time data
time_ext <- data_frame(mutation_kernel=character(),
                       rep=integer(),
                       mes=integer(),
                       time=integer())
mk <- 0
r <- 1
m <- 1
counter <- 1
mks_list <- unique(tdp_large$mutation_kernel)

# for (mk in mks_list) {
#   for (m in 1:100) {
#     time_data_sub <- subset(time_data_paired_new, rep == r & mes == m & mutation_kernel == mk)
#     tmax <- 1
#     text1 <- which(time_data_sub$num_individuals_species1 < 1)
#     text2 <- which(time_data_sub$num_individuals_species2 < 1)
#       if(length(text1) > 0){
#         tmax <- min(text1)
#       }else if(length(text2) > 0){
#         tmax <- min(text2)
#       }else{
#         tmax <- nrow(time_data_sub)
#       }
#       tt <- time_data_sub$tick[tmax]
#       #save all the values
#       time_ext[counter,1] <- as.character(time_data_sub[1,1])
#       time_ext[counter,2] <- time_data_sub[1,2]
#       time_ext[counter,3] <- time_data_sub[1,3]
#       time_ext[counter,4] <- tt
#       counter <- counter + 1
#   }
# }

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
#ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/fig4a_horz.pdf", width=7.5, height=5.0, units = "in")

#side analysis for sup fig 6: calculate fraction of pops that went extinct
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
#mutation kernel 74 has the most extinct (55%)

#checking for the pops that went extinct

#time and distance analyses on paired species 
r1e <- 0.498344
r2e <- 0.362690
a11e <- -0.000101493
a22e <- -0.0000384
time_threshold <- data_frame(mutation_kernel=character(),
                             rep=integer(),
                             mes=integer(),
                             time1=integer(),
                             distance1=double(),
                             time2=integer(),
                             distance2=double(),
                             extinct=logical())
mk <- 0
r <- 1
m <- 1
counter <- 1
mks_list <- unique(tdp_large$mutation_kernel)

for (mk in mks_list) {
  #for (r in 0:5) {
  for (m in 1:100) {
    time_data_sub <- subset(time_data_paired_new, rep == r & mes == m & mutation_kernel == mk)
    sp_ext <- F
    nr <- nrow(time_data_sub)
    if(nrow(time_data_sub > 0)){
      #calculate minimum distance from ESS (to minimize oscillations around the ESS?)
      distance_all1 <- fun_minim(r1e, a11e, time_data_sub$r1, time_data_sub$a11)
      d1min <- min(distance_all1, na.rm = T)
      distance1 <- distance_all1[which(distance_all1 == d1min)][1]
      
      tmax1 <- max(time_data_sub$tick)
      tt1 <- 0
      if(length(which(distance_all1 < 0.5)) > 0){
        tt1 <- min(time_data_sub$tick[which(distance_all1 < 0.5)])
      } else {
        tt1 <- tmax1
      }
      
      #calculate minimum distance from ESS (to minimize oscillations around the ESS?)
      distance_all2 <- fun_minim(r2e, a22e, time_data_sub$r2, time_data_sub$a22)
      d2min <- min(distance_all2, na.rm = T)
      distance2 <- distance_all2[which(distance_all2 == d2min)][1]
      
      tmax2 <- max(time_data_sub$tick)
      tt2 <- 0
      if(length(which(distance_all2 < 0.5)) > 0){
        tt2 <- min(time_data_sub$tick[which(distance_all2 < 0.5)])
      } else {
        tt2 <- tmax2
      }
      
      if(time_data_sub$num_individuals_species1[nr] < 1 || time_data_sub$num_individuals_species2[nr] < 1 ){
        sp_ext <- T
      }
      
      #save all the values
      time_threshold[counter,1] <- as.character(time_data_sub[1,1])
      time_threshold[counter,2] <- time_data_sub[1,2]
      time_threshold[counter,3] <- time_data_sub[1,3]
      time_threshold[counter,4] <- tt1
      time_threshold[counter,5] <- distance1
      time_threshold[counter,6] <- tt2
      time_threshold[counter,7] <- distance2
      time_threshold[counter,8] <- sp_ext
      counter <- counter + 1
    }
  }
  #}
}

tt_combine <- merge(time_threshold, mks, by = "mutation_kernel")
#tt_combine_large <- merge(tt_combine, reps, by.x = "rep", by.y = "IDs")
tt_combine_large <- tt_combine
tt_combine_large$mut <- as.factor(tt_combine_large$mut)
#tt_combine_large1 <- subset(tt_combine_large, mutr > 0.002)

#plot distance between endpoint and tradeoff curve
p_dist <-ggplot(data = tt_combine_large, aes(x = mut, y = distance1, colour = extinct)) +
  theme_bw()+
  #geom_boxplot(alpha = 0.6)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.5)+
  scale_y_log10()+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(m1~mutr)
p_dist
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_nov25/paired_dist1_2_2_extra.pdf", width=7.42, height=6.2, units = "in")

p_dist2 <-ggplot(data = tt_combine_large, aes(x = mut, y = distance2, colour = extinct)) +
  theme_bw()+
  #geom_boxplot(alpha = 0.6)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.5)+
  scale_y_log10()+
  scale_colour_viridis_d(begin = 0, end = 0.8, direction = -1) +
  scale_fill_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(m1~mutr)
p_dist2 
ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_nov25/paired_dist2_2_2.pdf", width=7.42, height=6.2, units = "in")

#-------- analysis of fixed mutations and trajectory of mutation fixation for populations that went extinct vs those that didn't --------
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

#expanded graphs over 500000 generations for subset of mutation kernels
time_data_paired_new <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_combined_cox_ext.csv", header = T))
fixed_muts_full <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_cox_ext.csv", header = T))
tdp_large <- merge(time_data_paired_new, mks_actual, by = "mutation_kernel")
#run the above code to a) check which populations went extinct
# + starting from 'analysis of fixed mutations and trajectory of evolution'

#--------- 2 - compare pops that went extinct vs didn't go extinct, same mutkern --------

cmut_ext <- c('1' = "#026133",
              '2' = "#026133",
              '3' = "#026133",
              '4' = "#026133",
              '5' = "#1D394C",
              '6' = "#1D394C",
              '7' = "#1D394C",
              '8' = "#1D394C") #colour names for the different species

fixed_muts_ext <- as.data.frame(read.csv("~/Documents/SLiM/outputs/mutdata/fixed_combined_extinct.csv", header = T))
colnames_temp <- c("rep","mes","mutation_kernel", "species","time_origin", "time_fixed","effect_size","mutation_type","blank") 
colnames(fixed_muts_ext) <- colnames_temp
fixed_muts_ext$mutation_type <- as.factor(fixed_muts_ext$mutation_type)
fm_subset3 <- subset(fixed_muts_ext, mutation_kernel == 48)
time_ext_subset <- subset(time_ext, mutation_kernel == 48)
gcsub3 <- subset(time_data_paired_new, mutation_kernel == 48)

#growth curves
gc3full <- ggplot()+
  theme_bw()+
  geom_line(data = gcsub3, mapping = aes(x = tick,y = r1, group = mes), colour = "#26828e") + 
  geom_line(data = gcsub3, mapping = aes(x = tick,y = r2, group = mes), colour = "#6ece58") + 
  facet_wrap(~ mes, nrow = 10)
gc3full

mes_extinct <- time_ext_subset$mes[which(time_ext_subset$time < 50000)]
mes_coexist <- time_ext_subset$mes[which(time_ext_subset$time > 50000)]

gc3_ext <- subset(gcsub3, mes %in% mes_extinct)
gc3_cox <- subset(gcsub3, mes %in% mes_coexist)
fm3_ext <- subset(fm_subset3, mes %in% mes_extinct)
fm3_cox <- subset(fm_subset3, mes %in% mes_coexist)

pcomb <- ggplot() +
  theme_bw() +
  geom_point(data = fm3_cox, aes(x = time_origin, y = effect_size, colour = mutation_type), alpha = 0.75) +
  scale_colour_manual(values = cmut) +
  new_scale_colour()+
  geom_point(data = fm3_ext, aes(x = time_origin, y = effect_size, colour = mutation_type), alpha = 0.75) +
  scale_colour_manual(values = cmut_ext) +
  facet_wrap(~ mutation_type, ncol = 4) +
  theme(legend.position = "none")
pcomb

#pick measurement reps for one that went extinct and one that didn't
fm3_mes_ext <- subset(fm_subset3, mes == 56)
fm3_mes_cox <- subset(fm_subset3, mes == 93)

#growth curve subsets
gc3_mes_ext <- subset(gcsub3, mes == 56)
gc3_mes_cox <- subset(gcsub3, mes == 93)

#growth curve plots
gc3_1 <- ggplot()+
  theme_bw()+
  geom_line(data = gc3_cox, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = gc3_ext, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#1D394C") + 
  #geom_line(data = gc3_mes_ext, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
  geom_line(data = gc3_mes_cox, mapping = aes(x = tick,y = r1, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm3_mes_cox, species == 1 & mutation_type == 7), aes(x = time_origin, y = 0.245), shape = 2, size = 2, stroke = 0.75, colour = '#000000')+ #ben r sp1 
  geom_point(data = subset(fm3_mes_cox, species == 1 & mutation_type == 5), aes(x = time_origin, y = 0.245), shape = 6, size = 2, stroke = 0.75, colour = '#B40000')+ #del r sp1 
  geom_hline(yintercept = r1e)
gc3_1

gc3_2 <- ggplot()+
  theme_bw()+
  geom_line(data = gc3_cox, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") + 
  geom_line(data = gc3_ext, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#026133") + 
  geom_line(data = gc3_mes_cox, mapping = aes(x = tick,y = r2, group = mes), colour = "#000000") + 
  geom_point(data = subset(fm3_mes_cox, species == 2 & mutation_type == 3), aes(x = time_origin, y = 0.20), shape = 2, size = 2, stroke = 0.75, colour = '#000000')+ #ben r sp2 
  geom_point(data = subset(fm3_mes_cox, species == 2 & mutation_type == 1), aes(x = time_origin, y = 0.20), shape = 6, size = 2, stroke = 0.75, colour = '#B40000')+ #del r sp2 
  geom_hline(yintercept = r2e)
gc3_2

#---------- analysis on varying mutation parameters between populations ------------
#doing this on cluster now

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

#make subsets of this graph
mks_m1 <- mks_actual[which(mks_actual$m1 < 0.3),]$mutation_kernel #supply of ben. muts

prop_ext_m1 <- prop_ext[which(prop_ext$mutation_kernel_1 %in% mks_m1 & prop_ext$mutation_kernel_2 %in% mks_m1),]

p_h2 <- ggplot(data = prop_ext_m1, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = frac_extinct)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h2

mks_mut <- mks_actual[which(mks_actual$mut < 5e-7 & mks_actual$mut > 5e-8),]$mutation_kernel #mutation rate

prop_ext_mut <- prop_ext[which(prop_ext$mutation_kernel_1 %in% mks_mut & prop_ext$mutation_kernel_2 %in% mks_mut),]

p_h3 <- ggplot(data = prop_ext_mut, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = frac_extinct)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h3

mks_ef <- mks_actual[which(mks_actual$mutr < 0.1 & mks_actual$mutr > 0.005),]$mutation_kernel #effect size  & mks_actual$mutr > 0.005

prop_ext_ef <- prop_ext[which(prop_ext$mutation_kernel_1 %in% mks_ef & prop_ext$mutation_kernel_2 %in% mks_ef),]

p_h4 <- ggplot(data = prop_ext_ef, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = frac_extinct)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_h4

ggsave("/Users/kasturilele/Documents/SLiM/plotdump/sims_dec25/prop_ext_subset.pdf", width=4.00, height=4.00, units = "in")

time_data_ef <- time_data_paired_full[which(time_data_paired_full$mutation_kernel_1 %in% mks_ef & time_data_paired_full$mutation_kernel_2 %in% mks_ef),]

gcfull <- ggplot()+
  theme_bw()+
  geom_line(data = time_data_ef, mapping = aes(x = tick,y = r1, group = mes), alpha = 0.5, colour = "#26828e") + 
  geom_line(data = time_data_ef, mapping = aes(x = tick,y = r2, group = mes), alpha = 0.5, colour = "#6ece58") +
  facet_grid(mutation_kernel_1 ~ mutation_kernel_2)
gcfull

#time and distance analyses on paired species 
r1e <- 0.498344
r2e <- 0.362690
a11e <- -0.0001015
a22e <- -0.0000384

distance_all1 <- fun_minim(r1e, a11e, time_data_paired_full$r1, time_data_paired_full$a11)

#calculate minimum distance from ESS (to minimize oscillations around the ESS?)
distance_all2 <- fun_minim(r2e, a22e, time_data_paired_full$r2, time_data_paired_full$a22)

#time_data_ef_large <- cbind(time_data_ef, distance_all1, distance_all2)

mindists <- which(distance_all1 < 0.5 | distance_all2 < 0.5)
time_data_mindists <- time_data_paired_full[mindists,]
mindist_IDs <- rownames(unique(time_data_mindists[,c(1,2,3,4)]))
mindists_unique <- time_data_paired_full[mindist_IDs,]



mindists_count <- mindists_unique %>%
  count(mutation_kernel_1, mutation_kernel_2) %>%
  ungroup()

p_ht <- ggplot(data = mindists_count, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = n)) +
  theme_bw()+
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_ht

#more restrictive distances
mindists1 <- which(distance_all1 < 0.5 & distance_all2 < 0.5)
time_data_mindists1 <- time_data_paired_full[mindists1,]
mindist_IDs1 <- rownames(unique(time_data_mindists1[,c(1,2,3,4)]))
mindists_unique1 <- time_data_paired_full[mindist_IDs1,]



mindists_count1 <- mindists_unique1 %>%
  count(mutation_kernel_1, mutation_kernel_2) %>%
  ungroup()

p_ht1 <- ggplot(data = mindists_count1, aes(x = mutation_kernel_1, y = mutation_kernel_2, fill = n)) +
  theme_bw()+
  geom_tile() +
  xlim(0,26)+
  ylim(0,26)+
  scale_fill_viridis_c(option = "viridis", direction = -1)
p_ht1

both <- plot_grid(
  plotlist = c(p_ht, p_ht1)
)
both  

#------------ expanded versions of varying across mutation rate ----------

time_data_15 <- as.data.frame(read.csv("~/Documents/SLiM/outputs/logfiles/log_comb_exp_15_ext2.csv", header = T))
mk15 <- as.data.frame(read.csv("~/Documents/SLiM/final_params/mutkern_15_ext2.csv", header = T))

#time_data_15_full <- time_data_15
#time_data_15 <- subset(time_data_15_full, mes < 200) #temporary subset for checking with previous graph
#time_data_15 <- time_data_15_full

#temp - histogram of reps that did not go extinct
# exts_hist <- merge(temp_rem, mk15, by = c('rep', 'mes'))
# 
# hist1 <- ggplot(data = exts_hist, aes(x = mut_sp1))+
#   theme_bw()+
#   xlim(1e-08, 1e-07) +
#   geom_histogram(bins = 100)
# hist1
# 
# 
# hist2 <- ggplot(data = exts_hist, aes(x = mut_sp2))+
#   theme_bw()+
#   xlim(1e-08, 1e-07) +
#   geom_histogram(bins = 100)
# hist2

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
#time_ext_binned <- time_ext %>% mutate(mk1_bin = ntile(mutr_sp1, n=b), mk2_bin = ntile(mutr_sp2, n=b))
#write.table(time_ext_binned, file = "~/Documents/SLiM/Rstuff/time_extinct_16.csv", append = F, sep = ",")
time_ext_binned15 <- as.data.frame(read.csv("~/Documents/SLiM/Rstuff/time_extinct_15.csv", header = T))

#figure out the boundaries of the bin
binlims <- data_frame(species=integer(),
                      bin=integer(),
                      min=numeric(),
                      max=numeric())
mk <- 1
counter <- 1

#rerun code for species 2 bins without resetting counter

for (mk in 1:b){
  bins <- subset(time_ext_binned15, mk2_bin == mk)$mut_sp2
  
  binlims[counter,1]<- 2
  binlims[counter,2]<- mk
  binlims[counter,3]<- min(bins)
  binlims[counter,4]<- max(bins)
  counter <- counter + 1
  
}
write.table(binlims, file = "~/Documents/SLiM/Rstuff/fig5_bins_15.csv", append = F, sep = ",")


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
#------- do the same for lower prop of ben.muts ----------

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


#------- new supp - distribution of parameter estimates for ri and aii ------
#figuring out whether my species are more 'bacteria-like' or 'yeast-like'
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
