
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
library(ggnet)
library(GGally)
library(network)
library(ggnet2)
library(GGally)
library(sna)
source("01_functions_updated.R")
source("food_web_generator.R")
source("plot_foodweb.R")
source("phase_plane_fweb.R")

#12 species food-web



fact<-expand.grid(delta = seq(0,5,0.25),
                  balance=c(0.2),
                  model=c("D-reduction-degree", "Simulations"),
                  reps=(1:20)) %>% as_tibble %>% mutate( consumer_density=0,top_predator_density=0,
                                                         Avg_density=0)

web<-pyramidal.food(S = 16,C = 0.1, prop=c(0.5,0.3,0.2), 
                    balance=0.1)


webs<- list(web$web)
prop<-web$index
S_N<- length(prop$Bottom)    # of resource species, N
S_C<- length(prop$Intermediate)  # of consumer species, C
S_P<- length(prop$Top)   # of predator species, P
TP<-c(rep(0,each=S_N), rep(1,each=S_C),rep(2,each=S_P))
d2<-plotweb(web = web$web,TP = TP,nspecies = 16)
d2

for(r in 1:nrow(fact)){
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
  gr_p<- -0.1 #death rate pradtor
  mortality_n<-0.5
  tmax<-500
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 500 #fact$forcing_duration[r]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  noise_c_d<-replicate(S_C,rnorm(d,0,0.1))
  noise_p_d<-replicate(S_P,rnorm(d,0,0.1))
  
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
  
  
  noise_c<-as.data.frame(list(times=times,import=rep(0,length(times))))
  noise_p<-as.data.frame(list(times=times, import=rep(0,length(times))))
  
  noise_c$import <- noise_c_d
  noise_p$import <- noise_p_d
  
  degree<-degree_N<-degree_C<-degree_P<-numeric()
  for(j in 1:(S_N+S_P+S_C) ){
    degree[j]<-sum(web$web[j,])
  }
  
  
  #the simulations over the data frame all factorial combinations
  
  noise_c_d<-replicate(S_C,rnorm(d,0,0.1))
  noise_p_d<-replicate(S_P,rnorm(d,0,0.1))
  
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
  
  
  noise_c<-as.data.frame(list(times=times,import=rep(0,length(times))))
  noise_p<-as.data.frame(list(times=times, import=rep(0,length(times))))
  
  noise_c$import <- noise_c_d
  noise_p$import <- noise_p_d
  
  
  species_index<-c(1,2,3,4,5,6,7,8) #which((degree) == max(degree))
  forcing_delta<-fact$delta[r]
  forcing_strength <- rep(0, (S_C+S_N+S_P)) #forcing strength,0 means 0 percent,1 means 100 percent
  parms<-list(S_N=S_N,S_C=S_C,Adjacency=web$web,prop=prop,S=(S_N+S_C+S_P),
              S_P=S_P,alpha.c=alpha.c,alpha.n=alpha.n,allee_b=allee_b,
              forcing_delta=forcing_delta,
              allee_c=allee_c,h=h,feed_eff_c=feed_eff_c,feed_eff_p=feed_eff_p,
              attack_rate_c=attack_rate_c,attack_rate_p=attack_rate_p,
              gr_n=gr_n,gr_c=gr_c,gr_p=gr_p,mortality_n=mortality_n,allee_p=allee_p,
              forcing_strength=forcing_strength, t1_N=t1_N, t1_C=t1_C, t1_P=t1_P,
              duration_mat_C=duration_mat_C, duration_mat_N=duration_mat_N,
              duration_mat_P=duration_mat_P, species_index=species_index, noise_c=noise_c,noise_p=noise_p)
  na<-rep(0,S_N) #initial N
  na[species_index]<- 0
  nc<-rep(0.000001,S_C) # inital C
  np<-rep(0.000001,S_P) # initla P
  
  if(fact$model[r] == "Simulations"){
    ic<-c(na,nc,np)
    sol<-ode(func=fwebdyn_perturb_continuous, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
      organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe
    cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
    
    (timeseries1<- sol %>% filter(type=="C" | type =="P") %>% 
        ggplot(aes(x =time,y =density, color=factor(species)))+
        geom_line(size=2, alpha=0.6)+
        xlim(c(0,499))+
        scale_colour_manual(values = cols)+
        geom_hline(yintercept = 0.01, 
                   linetype ="dashed", size =1)+
        ylab("C and P density")+
        annotate("text", x= 100,y=5, 
                 label="~Delta == 5",parse=TRUE, size=4)+
        # labs(color="Consumers \n and Predators")+
        theme_classic()+theme(legend.position = "none"))
    
    mean_density <- mean((sol %>% filter(type=="C" | type =="P", time == (tmax-5)))$density)
    mean_density_c <- mean((sol %>% filter(type=="C", time == (tmax-5)))$density)
    mean_density_p <- mean((sol %>% filter(type =="P", time == (tmax-5)))$density)
  }else if(fact$model[r] =="D-reduction-degree"){
    
    parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                     at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                     delta=forcing_delta, bi=gr_c,alpha_c=0.1)
    parameters$method <- "degree"
    parameters$delta <- forcing_delta
    ic_t<-c(0.000001,0.000001)
    theory_dyn<-ode(func=dynamics_reduced_model, y=ic_t, parms=parameters, times=seq(0, tmax, by=1)) 
    
    mean_density_c<-sum(theory_dyn[tmax-5,2])
    mean_density_p<-sum(theory_dyn[tmax-5,3])
    mean_density<- (mean_density_c+mean_density_p)/2
    
  }else{
    parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                     at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                     delta=forcing_delta, bi=gr_c,alpha_c=0.1)
    parameters$method <- "unweighted"
    parameters$delta <- forcing_delta
    ic_t<-c(0.000001,0.000001)
    theory_dyn<-ode(func=dynamics_reduced_model, y=ic_t, parms=parameters, times=seq(0, tmax, by=1)) 
    
    mean_density_c<-sum(theory_dyn[tmax-5,2])
    mean_density_p<-sum(theory_dyn[tmax-5,3])
    mean_density<- (mean_density_c+mean_density_p)/2
    
    
    
    
  }
  
  fact$consumer_density[r] <-  mean_density_c
  fact$top_predator_density[r] <- mean_density_p
  fact$Avg_density[r] <- mean_density
  
  print(r)
}

(a1<-fact %>% ggplot(aes(x=delta,y=Avg_density, color = model))+
    geom_point(size=5,alpha=0.6, aes(shape=model))+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_color_manual(values=c('#E69F00', '#56B4E9'))+
    xlab(expression(Delta))+ylab("Average density"))

species_index<-c(1,2,3,4,5,6) #which((degree) == max(degree))
forcing_delta<-5
forcing_strength <- rep(0, (S_C+S_N+S_P)) #forcing strength,0 means 0 percent,1 means 100 percent
parms<-list(S_N=S_N,S_C=S_C,Adjacency=web$web,prop=prop,S=(S_N+S_C+S_P),
            S_P=S_P,alpha.c=alpha.c,alpha.n=alpha.n,allee_b=allee_b,
            forcing_delta=forcing_delta,
            allee_c=allee_c,h=h,feed_eff_c=feed_eff_c,feed_eff_p=feed_eff_p,
            attack_rate_c=attack_rate_c,attack_rate_p=attack_rate_p,
            gr_n=gr_n,gr_c=gr_c,gr_p=gr_p,mortality_n=mortality_n,allee_p=allee_p,
            forcing_strength=forcing_strength, t1_N=t1_N, t1_C=t1_C, t1_P=t1_P,
            duration_mat_C=duration_mat_C, duration_mat_N=duration_mat_N,
            duration_mat_P=duration_mat_P, species_index=species_index, noise_c=noise_c,noise_p=noise_p)
na<-rep(0,S_N) #initial N
na[species_index]<- forcing_delta
nc<-rep(0.00001,S_C) # inital C
np<-rep(0.00001,S_P) # initla P
ic<-c(na,nc,np)
sol<-ode(func=fwebdyn_perturb_constant_noise, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
  organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe
cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))

(timeseries1<- sol %>% filter(type=="C" | type =="P") %>% 
    ggplot(aes(x =time,y =density, color=factor(species)))+
    geom_line(size=2, alpha=0.6)+
    xlim(c(0,499))+
    scale_colour_manual(values = cols)+
    geom_hline(yintercept = 0.01, 
               linetype ="dashed", size =1)+
    ylab("C and P density")+
    annotate("text", x= 100,y=5, 
             label="~Delta == 5",parse=TRUE, size=4)+
    # labs(color="Consumers \n and Predators")+
    theme_classic()+theme(legend.position = "none"))
(dens1<-sol %>% filter(type=="C" | type =="P", time == (tmax-5)) %>%
    ggplot(aes(x=species,y=density, fill = factor(species)))+
    geom_col()+
    labs(fill="Species")+
    annotate("text", x= 7,y=4, 
             label="~Delta == 5",parse=TRUE, size=4)+
    geom_hline(yintercept = 0.01, linetype ="dashed", 
               size =1)+
    scale_fill_manual(values = cols)+
    theme_classic()+theme(legend.position = "none"))



##########################################################
######################## Two predators ###################






fact1<-expand.grid(delta = seq(0,5,0.25),
                  balance=c(0.1),
                  model=c("D-reduction-degree", "Simulations"),
                  reps=(1:20)) %>% as_tibble %>% mutate( consumer_density=0,top_predator_density=0,
                                                         Avg_density=0)

web<-pyramidal.food(S = 16,C = 0.1, prop=c(0.5,0.4,.1), 
                    balance=0.05)


webs<- list(web$web)
prop<-web$index
S_N<- length(prop$Bottom)    # of resource species, N
S_C<- length(prop$Intermediate)  # of consumer species, C
S_P<- length(prop$Top)   # of predator species, P
TP<-c(rep(0,each=S_N), rep(1,each=S_C),rep(2,each=S_P))
d2<-plotweb(web = web$web,TP = TP,nspecies = 16)
d2

for(r in 1:nrow(fact1)){
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
  gr_p<- -0.1 #death rate pradtor
  mortality_n<-0.5
  tmax<-500
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 500 #fact1$forcing_duration[r]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  noise_c_d<-replicate(S_C,rnorm(d,0,0.1))
  noise_p_d<-replicate(S_P,rnorm(d,0,0.1))
  
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
  
  
  noise_c<-as.data.frame(list(times=times,import=rep(0,length(times))))
  noise_p<-as.data.frame(list(times=times, import=rep(0,length(times))))
  
  noise_c$import <- noise_c_d
  noise_p$import <- noise_p_d
  
  degree<-degree_N<-degree_C<-degree_P<-numeric()
  for(j in 1:(S_N+S_P+S_C) ){
    degree[j]<-sum(web$web[j,])
  }
  
  
  #the simulations over the data frame all fact1orial combinations
  
  noise_c_d<-replicate(S_C,rnorm(d,0,0.1))
  noise_p_d<-replicate(S_P,rnorm(d,0,0.1))
  
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
  
  
  noise_c<-as.data.frame(list(times=times,import=rep(0,length(times))))
  noise_p<-as.data.frame(list(times=times, import=rep(0,length(times))))
  
  noise_c$import <- noise_c_d
  noise_p$import <- noise_p_d
  
  
  species_index<-c(1,2,3,4,5,6,7,8) #which((degree) == max(degree))
  forcing_delta<-fact1$delta[r]
  forcing_strength <- rep(0, (S_C+S_N+S_P)) #forcing strength,0 means 0 percent,1 means 100 percent
  parms<-list(S_N=S_N,S_C=S_C,Adjacency=web$web,prop=prop,S=(S_N+S_C+S_P),
              S_P=S_P,alpha.c=alpha.c,alpha.n=alpha.n,allee_b=allee_b,
              forcing_delta=forcing_delta,
              allee_c=allee_c,h=h,feed_eff_c=feed_eff_c,feed_eff_p=feed_eff_p,
              attack_rate_c=attack_rate_c,attack_rate_p=attack_rate_p,
              gr_n=gr_n,gr_c=gr_c,gr_p=gr_p,mortality_n=mortality_n,allee_p=allee_p,
              forcing_strength=forcing_strength, t1_N=t1_N, t1_C=t1_C, t1_P=t1_P,
              duration_mat_C=duration_mat_C, duration_mat_N=duration_mat_N,
              duration_mat_P=duration_mat_P, species_index=species_index, noise_c=noise_c,noise_p=noise_p)
  na<-rep(0,S_N) #initial N
  na[species_index]<- 0
  nc<-rep(0.000001,S_C) # inital C
  np<-rep(0.000001,S_P) # initla P
  
  if(fact1$model[r] == "Simulations"){
    ic<-c(na,nc,np)
    sol<-ode(func=fwebdyn_perturb_continuous, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
      organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe
    cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
    
    (timeseries1<- sol %>% filter(type=="C" | type =="P") %>% 
        ggplot(aes(x =time,y =density, color=fact1or(species)))+
        geom_line(size=2, alpha=0.6)+
        xlim(c(0,499))+
        scale_colour_manual(values = cols)+
        geom_hline(yintercept = 0.01, 
                   linetype ="dashed", size =1)+
        ylab("C and P density")+
        annotate("text", x= 100,y=5, 
                 label="~Delta == 5",parse=TRUE, size=4)+
        # labs(color="Consumers \n and Predators")+
        theme_classic()+theme(legend.position = "none"))
    
    mean_density <- mean((sol %>% filter(type=="C" | type =="P", time == (tmax-5)))$density)
    mean_density_c <- mean((sol %>% filter(type=="C", time == (tmax-5)))$density)
    mean_density_p <- mean((sol %>% filter(type =="P", time == (tmax-5)))$density)
  }else if(fact1$model[r] =="D-reduction-degree"){
    
    parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                     at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                     delta=forcing_delta, bi=gr_c,alpha_c=0.1)
    parameters$method <- "degree"
    parameters$delta <- forcing_delta
    ic_t<-c(0.000001,0.000001)
    theory_dyn<-ode(func=dynamics_reduced_model, y=ic_t, parms=parameters, times=seq(0, tmax, by=1)) 
    
    mean_density_c<-sum(theory_dyn[tmax-5,2])
    mean_density_p<-sum(theory_dyn[tmax-5,3])
    mean_density<- (mean_density_c+mean_density_p)/2
    
  }else{
    parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                     at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                     delta=forcing_delta, bi=gr_c,alpha_c=0.1)
    parameters$method <- "unweighted"
    parameters$delta <- forcing_delta
    ic_t<-c(0.000001,0.000001)
    theory_dyn<-ode(func=dynamics_reduced_model, y=ic_t, parms=parameters, times=seq(0, tmax, by=1)) 
    
    mean_density_c<-sum(theory_dyn[tmax-5,2])
    mean_density_p<-sum(theory_dyn[tmax-5,3])
    mean_density<- (mean_density_c+mean_density_p)/2
    
    
    
    
  }
  
  fact1$consumer_density[r] <-  mean_density_c
  fact1$top_predator_density[r] <- mean_density_p
  fact1$Avg_density[r] <- mean_density
  
  print(r)
}



(a1<-fact1 %>% ggplot(aes(x=delta,y=Avg_density, color = model))+
    geom_point(size=5,alpha=0.6, aes(shape=model))+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_color_manual(values=c('#E69F00', '#56B4E9'))+
    xlab(expression(Delta))+ylab("Average density"))

