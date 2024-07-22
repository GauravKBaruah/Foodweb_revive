
rm(list=ls())

require(statmod)
require(ggplot2)
require(tidyr)
require(dplyr)
require(deSolve)
require(viridis)
library(cowplot)
library(igraph)
library(network)
require(flux)
require(akima)
library(ggnet)
source("01_functions_updated.R")
source("food_web_generator.R")
source("plot_foodweb.R")
source("phase_plane_fweb.R")
fact<-expand.grid( C=seq(0.09,0.35,0.025),
                   balance=seq(0.2,0.9,0.1),
                   S=c(24), 
                   delta=c(0.1, 1, 5),
                   random_seed=4327+(1:1)*100) %>% as_tibble() %>%
  mutate(foodweb_richness=0,
         foodweb_biomass=0,
         intermediate_richness=0,
         top_inter_richness=0,
         top_richness=0,
         intermediate_biomass=0,
         top_biomass=0,
         theory_biomass=0,
         difference_biomass=0,
         simulation_biomass=0)
#?expand.grid
plot_phase<-list()
data_t<-NULL
webs<-list()

for(i in 1:nrow(fact)){
  
  
  
  web<-pyramidal.food(S = fact$S[i],C = fact$C[i],prop=c(0.5,0.3,0.2), 
                      balance=fact$balance[i])
  webs[i]<- list(web$web)
  prop<-web$index
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  allee_b<-0
  allee_c<-0
  allee_p<-0
  h<-0.25 #
  feed_eff_c<-0.5
  feed_eff_p<-0.5
  alpha.c<-matrix((runif( S_C*S_C, 0.1,0.1)), nrow=S_C, ncol= S_C)
  alpha.n<-matrix((runif( S_N*S_N, 0.5,0.5)), nrow=S_N, ncol= S_N)
  
  attack_rate_c<-1
  attack_rate_p<-1
  gr_n<- 0.1 #growth rate respurce
  gr_c<- -0.1 # death rate consumer
  gr_p<- -0.05 #death rate pradtor
  mortality_n<-0.5
  tmax<-500
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 500 #fact$forcing_duration[r]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  noise_c_d<-replicate(S_C,rnorm(d,0,0.1))
  noise_p_d<-replicate(S_P,rnorm(d,0,0.1))
  #noise_consumer<- as.data.frame(cbind(seq(0,tmax,1),noise_c_d))
  #noise_top<-as.data.frame(cbind(seq(0,tmax,1), noise_p_d))
  
  duration_mat_N<-(replicate(S_N,d))
  duration_mat_C<-(replicate(S_C,d))
  duration_mat_P<-(replicate(S_P,d))
  times<-seq(0, tmax, 1)
  t1_N<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_C<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_N$import<-duration_mat_N[,1]
  t1_P$import<-duration_mat_P[,1]
  t1_C$import<-duration_mat_C[,1]
  
  #noise_c<-as.data.frame(list(times=times,import=rep(0,length(times))))
  #noise_p<-as.data.frame(list(times=times, import=rep(0,length(times))))
  
  #noise_c$import <- noise_c_d
  #noise_p$import <- noise_p_d
  
  degree<-degree_N<-degree_C<-degree_P<-numeric()
  for(r in 1:(S_N+S_P+S_C) ){
    degree[r]<-sum(web$web[r,])
  }
  species_index<-c(1,2,3,4,5,6,7,8,9,10,11,12) #which((degree) == max(degree))
  forcing_delta<-fact$delta[i]
  forcing_strength <- rep(0, (S_C+S_N+S_P)) #forcing strength,0 means 0 percent,1 means 100 percent
  parms<-list(S_N=S_N,S_C=S_C,Adjacency=web$web,prop=prop,S=(S_N+S_C+S_P),
              S_P=S_P,alpha.c=alpha.c,alpha.n=alpha.n,allee_b=allee_b,
              forcing_delta=forcing_delta,
              allee_c=allee_c,h=h,feed_eff_c=feed_eff_c,feed_eff_p=feed_eff_p,
              attack_rate_c=attack_rate_c,attack_rate_p=attack_rate_p,
              gr_n=gr_n,gr_c=gr_c,gr_p=gr_p,mortality_n=mortality_n,allee_p=allee_p,
              forcing_strength=forcing_strength, t1_N=t1_N, t1_C=t1_C, t1_P=t1_P,
              duration_mat_C=duration_mat_C, duration_mat_N=duration_mat_N,
              duration_mat_P=duration_mat_P, species_index=species_index)
  na<-rep(0,S_N) #initial N
  na[species_index]<- parms$forcing_delta
  nc<-rep(0.000001,S_C) # inital C
  np<-rep(0.000001,S_P) # initla P

  ic<-c(na,nc,np)
  sol<-ode(func=fwebdyn_perturb_constant  , y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
    organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe
  
  cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
  (timeseries1<- sol %>% filter(type=="C" | type =="P") %>% 
      ggplot(aes(x =time,y =density, color=factor(species)))+
      geom_line(size=2, alpha=0.6)+
      xlim(c(0,499))+
      scale_colour_manual(values = cols)+
      geom_hline(yintercept = 0.01, linetype ="dashed")+
      ylab("C and P density")+
      annotate("text", x= 100,y=5, 
               label="~Delta == 5",parse=TRUE, size=4)+
      # labs(color="Consumers \n and Predators")+
      theme_classic()+theme(legend.position = "none"))
  
  
  parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                   at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                   delta=forcing_delta, bi=gr_c,alpha_c=0.5)
  parameters$method <- "degree"
  parameters$delta <- forcing_delta
  ic_t<-c(0.000001,0.000001)
  theory_dyn<-ode(func=dynamics_reduced_model, y=ic_t, parms=parameters, 
                  times=seq(0, tmax, by=1)) 
  
  mean_density_c<-sum(mean(theory_dyn[2:tmax,2]))
  mean_density_p<-sum(mean(theory_dyn[2:tmax,3]))
  mean_density<- (mean_density_c+mean_density_p)/2
  
  
  C_eq<-((sol %>% filter(time > 100, type =="C")) %>% group_by(species) %>% summarise(mean_density_C=mean(density)))$mean_density_C
  P_eq<-((sol %>% filter(time > 100, type =="P")) %>% group_by(species) %>% summarise(mean_density_P=mean(density)))$mean_density_P
  
  inter_mean_density<-mean(C_eq)
  top_mean_density<-mean(P_eq)
  inter_richness<-length(which(((sol %>% filter(time > 100, type =="C")) %>% group_by(species) %>% summarise(mean_density_C=mean(density)))$mean_density_C > 0.05))
  top_richness<-length(which(((sol %>% filter(time > 100, type =="P")) %>% group_by(species) %>% summarise(mean_density_P=mean(density)))$mean_density_P > 0.05))
  
  
  fact$intermediate_biomass[i]<-inter_mean_density
  fact$top_biomass[i]<-top_mean_density
  fact$intermediate_richness[i]<-inter_richness
  fact$top_richness[i]<-top_richness
  fact$top_inter_richness[i]<-(top_richness+inter_richness)/(S_C+S_P)
  fact$simulation_biomass[i]<- mean(c(C_eq,P_eq))
  fact$theory_biomass[i]<- mean_density
  fact$difference_biomass_normalised[i]<- (mean_density-(fact$simulation_biomass[i]))/mean_density
  
  fact$difference_biomass[i]<- abs(mean_density-(fact$simulation_biomass[i]))
  print(i)
  
 
 
}



delta_0<-fact %>% filter(delta ==.1)
plot_delta_0<-ggplot(delta_0 , aes(x=balance,y=C,z=(top_inter_richness)))+
  geom_raster(aes(fill=top_inter_richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0,1), high = 'darkgreen', low = 'red')+
  xlab("Proportion of top predator links")+ylab("Connectance") + 
  ggtitle("A")+
  theme_cowplot()+
  labs(fill="Predator and \n consumer richness")+
  annotate("segment",x=0.25, xend=0.6,y=0.15,
           yend=0.39, colour="white", size=3,linetype="dashed")+
  annotate("text", x= 0.25,y=0.38, 
           label="~Delta == 0.1",parse=TRUE, size=5)+
  annotate("text", x= 0.65,y=0.1, 
           label="Not~recoverable",colour="grey", parse=TRUE, size=5.5)


plot_delta_0
####

delta_1<-fact %>% filter(delta ==1)
plot_delta_1<-ggplot(delta_1 , aes(x=balance,y=C,z=(top_inter_richness)))+
  geom_raster(aes(fill=top_inter_richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0,1), high = 'darkgreen', low = 'red')+
  xlab("Proportion of top predator links")+ylab("Connectance") + 
  ggtitle("B")+
  theme_cowplot()+
  labs(fill="Predator and \n consumer richness")+
  annotate("text", x= 0.25,y=0.38, 
           label="~Delta == 1",parse=TRUE, size=5)+
  annotate("segment",x=0.2, xend=0.7,y=0.1,
           yend=0.39, colour="white", size=4,linetype="dashed")+
 # annotate("text", x= 0.5,y=0.3, 
  #         label="Recoverable",colour="grey",
   #        parse=TRUE, size=5)+
  annotate("text", x= 0.66,y=0.1, 
           label="Not~recoverable",colour="grey",
           parse=TRUE, size=5)

#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
plot_delta_1


delta_5<-fact %>% filter(delta ==5)


plot_delta_5<-ggplot(delta_5 , aes(x=balance,y=C,z=(top_inter_richness)))+
  geom_raster(aes(fill=top_inter_richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0,1),high = 'darkgreen', low = 'red')+
  xlab("Proportion of top predator links")+ylab("Connectance") + 
  ggtitle("C")+
  theme_cowplot()+
  labs(fill="Predator and \n consumer richness")+
  annotate("text", x= 0.25,y=0.39, 
           label="~Delta == 5",parse=TRUE, size=5)
 # annotate("segment",x=0.3, xend=0.5,y=0.1,
  #         yend=0.39, colour="white", size=4,linetype="dashed")+
  #annotate("text", x= 0.25,y=0.3, 
   #        label="Recoverable",colour="grey",
    #       parse=TRUE, size=5)+
 # annotate("text", x= 0.76,y=0.15, 
  #         label="Not~recoverable",colour="grey", parse=TRUE, size=5)

plot_delta_5






# 
# TP<-c(rep(0,each=S_N), rep(1,each=S_C),rep(2,each=S_P))
# 
# web1<-pyramidal.food(S = 18,C = 0.1,prop=c(0.4,0.4,0.2), 
#                      balance=0.2)
# web2<-pyramidal.food(S = 18,C = 0.1,prop=c(0.4,0.4,0.2), 
#                      balance=1)
# web3<-pyramidal.food(S = 18,C = 0.3,prop=c(0.4,0.4,0.2), 
#                      balance=0.2)
# web4<-pyramidal.food(S = 18,C = 0.3,prop=c(0.4,0.4,0.2), 
#                      balance=1)
# 
# 
# 
# d1<- plotweb(web = web1[[1]],TP = TP, nspecies = 18)
# d2<- plotweb(web = web2[[1]],TP = TP, nspecies = 18)
# d3<-plotweb(web = web3[[1]],TP = TP, nspecies = 18)
# d4<-plotweb(web = web4[[1]],TP = TP, nspecies = 18)


grid_arrange_shared_legend(plot_delta_0,plot_delta_1,
                           plot_delta_5,nrow =1,ncol=3,
                           position = "bottom")



delta_1<-fact %>% filter(delta ==0.1)
plot_diff_1<-ggplot(delta_1 , aes(x=balance,y=C,z=difference_biomass))+
  geom_raster(aes(fill=(difference_biomass),show.legend =TRUE))+ 
  scale_fill_gradient(limits=range(abs(delta_1$difference_biomass)), high = 'darkgreen', low = 'red')+
  xlab("Proportion of top predator links")+ylab("Connectance") + 
  ggtitle("A")+
  theme_cowplot()+
  labs(fill="Difference in \n average density \n between models")+
  #annotate("segment",x=0.2, xend=0.5,y=0.2,
   #        yend=0.39, colour="white", size=3,linetype="dashed")+
  annotate("text", x= 0.25,y=0.39, 
           label="~Delta == 0.1",parse=TRUE, size=5)
  #annotate("text", x= 0.25,y=0.3, 
   #        label="Poor~ \n prediction",
    #       colour="grey", parse=TRUE, size=4)

plot_diff_1


delta_2<-fact %>% filter(delta ==1)
plot_diff_2<-ggplot(delta_2 , aes(x=balance,y=C,z=(difference_biomass)))+
  geom_raster(aes(fill=(difference_biomass),show.legend =TRUE))+ 
  scale_fill_gradient(limits=range((delta_2$difference_biomass)), high = 'darkgreen', low = 'red')+
  xlab("Proportion of top predator links")+ylab("Connectance") + 
  ggtitle("B")+
  theme_cowplot()+
  labs(fill="Difference in \n average density \n between models")+
  #annotate("segment",x=0.25, xend=0.4,y=0.27,
   #        yend=0.39, colour="white", size=3,linetype="dashed")+
  annotate("text", x= 0.25,y=0.39, 
           label="~Delta == 1",parse=TRUE, size=5)
 # annotate("text", x= 0.3,y=0.35, 
  #         label="Poor~prediction",
   #        colour="grey", parse=TRUE, size=4)


plot_diff_2

delta_5<-fact %>% filter(delta ==5)
plot_diff_5<-ggplot(delta_5 , aes(x=balance,y=C,z=(difference_biomass)))+
  geom_raster(aes(fill=difference_biomass),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(delta_5$difference_biomass), high = 'darkgreen', low = 'red')+
  xlab("Proportion of top predator links")+ylab("Connectance") + 
  ggtitle("C")+
  theme_cowplot()+
  labs(fill="Difference in \n average density \n between models")+
  annotate("segment",x=0.2, xend=0.85,y=0.1,
          yend=0.35, colour="white", size=3,linetype="dashed")+
  annotate("text", x= 0.25,y=0.39, 
           label="~Delta == 5",parse=TRUE, size=5)+
  annotate("text", x= 0.4,y=0.25, 
           label="Poor~prediction",
         colour="grey", parse=TRUE, size=5)



plot_diff_5




ggpubr::ggarrange(plot_diff_1,plot_diff_2,plot_diff_5,
                  nrow = 1)
grid_arrange_shared_legend(plot_diff_5,plot_diff_2,
                           plot_diff_1,nrow =1,ncol=3,
                           position = "right")


grid_arrange_shared_legend(plot_delta_0,plot_delta_1,
                           plot_delta_5,nrow =1,ncol=3,
                           position = "bottom")


#######

fact_1<-expand.grid( C=c(0.09, 0.15, 0.25),
                   balance=0.3,
                   S=c(16), delta=seq(0, 5, 0.5),
                   forcing_species=c(1,2,3,4,5,6,7,8),
                   random_seed=4327+(1:1)*100) %>% as_tibble() %>%
  mutate(foodweb_richness=0,
         foodweb_biomass=0,
         intermediate_richness=0,
         top_richness=0,
         intermediate_biomass=0,
         top_biomass=0)
#?expand.grid
plot_phase<-list()
data_t<-NULL
webs<-list()

for(i in 1:nrow(fact_1)){
  
  
  
  web<-pyramidal.food(S = fact_1$S[i],C = fact_1$C[i],prop=c(0.5,0.3,0.2), 
                      balance=fact_1$balance[i])
  webs[i]<- list(web$web)
  prop<-web$index
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  alpha.c<-matrix((runif( S_C*S_C, 0.1,0.1)), nrow=S_C, ncol= S_C)
  alpha.n<-matrix((runif( S_N*S_N, 0.5,0.5)), nrow=S_N, ncol= S_N)
  allee_b<-0
  allee_c<-0
  allee_p<-0
  h<-0.25 #
  feed_eff_c<-0.5
  feed_eff_p<-0.5
  attack_rate_c<-1.25
  attack_rate_p<-1.25
  gr_n<- 0.1 #growth rate respurce
  gr_c<- -0.1 # death rate consumer
  gr_p<- -0.05  #death rate pradtor
  mortality_n<-0.5
  tmax<-500
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 500 #fact_1$forcing_duration[r]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  noise_c_d<-replicate(S_C,rnorm(d,0,0.1))
  noise_p_d<-replicate(S_P,rnorm(d,0,0.1))
  #noise_consumer<- as.data.frame(cbind(seq(0,tmax,1),noise_c_d))
  #noise_top<-as.data.frame(cbind(seq(0,tmax,1), noise_p_d))
  
  duration_mat_N<-(replicate(S_N,d))
  duration_mat_C<-(replicate(S_C,d))
  duration_mat_P<-(replicate(S_P,d))
  times<-seq(0, tmax, 1)
  t1_N<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_C<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_N$import<-duration_mat_N[,1]
  t1_P$import<-duration_mat_P[,1]
  t1_C$import<-duration_mat_C[,1]
  
  #noise_c<-as.data.frame(list(times=times,import=rep(0,length(times))))
  #noise_p<-as.data.frame(list(times=times, import=rep(0,length(times))))
  
  #  noise_c$import <- noise_c_d
  #noise_p$import <- noise_p_d
  
  degree<-degree_N<-degree_C<-degree_P<-numeric()
  for(r in 1:(S_N+S_P+S_C) ){
    degree[r]<-sum(web$web[r,])
  }
  species_index<-seq(1,fact_1$forcing_species[i],1) #which((degree) == max(degree))
  forcing_delta<-fact_1$delta[i]
  forcing_strength <- rep(0, (S_C+S_N+S_P)) #forcing strength,0 means 0 percent,1 means 100 percent
  parms<-list(S_N=S_N,S_C=S_C,Adjacency=web$web,prop=prop,S=(S_N+S_C+S_P),
              S_P=S_P,alpha.c=alpha.c,alpha.n=alpha.n,allee_b=allee_b,
              forcing_delta=forcing_delta,
              allee_c=allee_c,h=h,feed_eff_c=feed_eff_c,feed_eff_p=feed_eff_p,
              attack_rate_c=attack_rate_c,attack_rate_p=attack_rate_p,
              gr_n=gr_n,gr_c=gr_c,gr_p=gr_p,mortality_n=mortality_n,allee_p=allee_p,
              forcing_strength=forcing_strength, t1_N=t1_N, t1_C=t1_C, t1_P=t1_P,
              duration_mat_C=duration_mat_C, duration_mat_N=duration_mat_N,
              duration_mat_P=duration_mat_P, species_index=species_index)
  na<-rep(0,S_N) #initial N
  na[species_index]<- parms$forcing_delta
  nc<-rep(0.00001,S_C) # inital C
  np<-rep(0.00001,S_P) # initla P
  ic<-c(na,nc,np)
  sol<-ode(func=fwebdyn_perturb_constant, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
    organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe
  
  (timeseries1<- sol %>% filter(type=="C" | type =="P") %>% 
      ggplot(aes(x =time,y =density, color=factor(species)))+
      geom_line(size=2, alpha=0.6)+
      xlim(c(0,499))+
      scale_colour_manual(values = cols)+
      geom_hline(yintercept = 0.01, linetype ="dashed", size =1.5)+
      ylab("C and P density")+
      annotate("text", x= 100,y=5, 
               label="~Delta == 5",parse=TRUE, size=4)+
      # labs(color="Consumers \n and Predators")+
      theme_classic()+theme(legend.position = "none"))
  (dens1<-sol %>% filter(type=="C" | type =="P", time == (20)) %>%
      ggplot(aes(x=species,y=density, fill = factor(species)))+
      geom_col()+
      labs(fill="Species")+
      annotate("text", x= 7,y=4, 
               label="~Delta == 5",parse=TRUE, size=4)+
      geom_hline(yintercept = 0.01, linetype ="dashed", size =1.5)+
      scale_fill_manual(values = cols)+
      theme_classic()+theme(legend.position = "none"))
  
  
  C_eq<-((sol %>% filter(time > 100, type =="C")) 
         %>% group_by(species) %>% summarise(mean_density_C=mean(density)))$mean_density_C
  P_eq<-((sol %>% filter(time > 100, type =="P")) 
         %>% group_by(species) %>% summarise(mean_density_P=mean(density)))$mean_density_P
  
  inter_mean_density<-mean(C_eq)
  top_mean_density<-mean(P_eq)
  inter_richness<-length(which(((sol %>% filter(time > 100, type =="C")) %>% group_by(species) %>% summarise(mean_density_C=mean(density)))$mean_density_C > 0.05))
  top_richness<-length(which(((sol %>% filter(time > 100, type =="P")) %>% group_by(species) %>% summarise(mean_density_P=mean(density)))$mean_density_P > 0.05))
  
  
  fact_1$intermediate_biomass[i]<-inter_mean_density
  fact_1$top_biomass[i]<-top_mean_density
  fact_1$intermediate_richness[i]<-inter_richness
  fact_1$top_richness[i]<-top_richness
  fact_1$top_inter_richness[i]<- (top_richness+inter_richness)
  fact_1$simulation_biomass[i]<- mean(c(C_eq,P_eq))
  print(i)
  

  
}


plotweb(web = web1[[1]],TP = TP, nspecies = 16)

delta_4<-fact_1 %>% filter(C ==0.09, delta > 0)
(plot_delta_4<-ggplot(delta_4 , aes(x=forcing_species,y=delta,z=(top_inter_richness)/8 ))+
  geom_raster(aes(fill=(top_inter_richness)/8),show.legend =TRUE)+ 
  scale_fill_gradient(limits=c(0,1), high = 'darkgreen', low = 'red')+
  xlab("No. of basal species forced")+ylab(expression(Delta)) + 
  ggtitle("D")+
  theme_cowplot()+
  labs(fill="Predator and \n consumer richness")+
  #annotate("segment",x=2, xend=7,y=5,
   #        yend=3, colour="white", size=3,linetype="dashed")+
  annotate("text", x= 2,y=5.5, 
           label="C == 0.09",parse=TRUE, size=5))
 # annotate("text", x= 3,y=5, 
  #         label="Not~recoverable",colour="grey", parse=TRUE, size=6)+
  #annotate("text", x= 5,y=5, 
   #        label="Recoverable",colour="grey", parse=TRUE, size=6)


delta_7<-fact_1 %>% filter(C ==0.15, delta > 0)
(plot_delta_7<-ggplot(delta_7 , aes(x=forcing_species,y=delta,z=(top_inter_richness/8)))+
  geom_raster(aes(fill=(top_inter_richness)/8,show.legend =TRUE))+ 
  scale_fill_gradient(limits=c(0,1), high = 'darkgreen', low = 'red')+
  xlab("No. of basal species forced")+ylab(expression(Delta)) + 
  ggtitle("E")+
  theme_cowplot()+
  labs(fill="Predator and \n consumer richness")+
#  annotate("segment",x=1, xend=7,y=5,
 #          yend=1, colour="white", size=3,linetype="dashed")+
  annotate("text", x= 2,y=5.5, 
           label="C == 0.15",parse=TRUE, size=5))
#  annotate("text", x= 3,y=5, 
 #          label="Not~recoverable",colour="grey", parse=TRUE, size=6)+
  #annotate("text", x= 6,y=5, 
   #        label="Recoverable",colour="grey", parse=TRUE, size=6)

delta_6<-fact_1 %>% filter(C ==0.25, delta > 0)
(plot_delta_6<-ggplot(delta_6 , aes(x=forcing_species,y=delta,z=(top_inter_richness/8)))+
  geom_raster(aes(fill=(top_inter_richness)/8,show.legend =TRUE))+ 
  scale_fill_gradient(limits=c(0,1), high = 'darkgreen', low = 'red')+
  xlab("No. of basal species forced")+ylab(expression(Delta)) + 
  ggtitle("F")+
  theme_cowplot()+
  labs(fill="Predator and \n consumer richness")+
  #annotate("segment",x=1, xend=6,y=5,
   #        yend=1, colour="white", size=3,linetype="dashed")+
  annotate("text", x= 2,y=5.5, 
           label="C == 0.25",parse=TRUE, size=5))
 # annotate("text", x= 3,y=5, 
  #         label="Not~recoverable",colour="grey", parse=TRUE, size=6)+
  #annotate("text", x= 6,y=5, 
   #        label="Recoverable",colour="grey", parse=TRUE, size=6)






grid_arrange_shared_legend(plot_delta_0,plot_delta_1,
                           plot_delta_5,
                           plot_delta_4,plot_delta_7,
                           plot_delta_6,nrow =2,
                           ncol=3,position = "right")


