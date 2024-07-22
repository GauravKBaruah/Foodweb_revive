
rm(list=ls())

library(gridBase)
library(ggplotify)
require(statmod)
require(ggplot2)
require(tidyr)
require(dplyr)
require(deSolve)
require(viridis)
library(cowplot)
library(igraph)
library(network)
library(ggnet2)
library(phaseR)
library(GGally)
source("01_functions_updated.R")
source("food_web_generator.R")
source("plot_foodweb.R")
source("phase_plane_fweb.R")

#12 species food-web
  
  
  web<-pyramidal.food(S = 12,C = 0.1, prop=c(0.5,0.4,0.1), 
                      balance=0.5)
  
  webs<- list(web$web)
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
  duration <- 500 #fact$forcing_duration[r]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  
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
  

  degree<-degree_N<-degree_C<-degree_P<-numeric()
  for(r in 1:(S_N+S_P+S_C) ){
    degree[r]<-sum(web$web[r,])
  }
  species_index<-c(1,2,3,4,5,6) #which((degree) == max(degree))
  forcing_delta<-0.1
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
  
  cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
  (timeseries1<- sol %>% filter(type=="C" | type =="P") %>% 
    ggplot(aes(x =time,y =density, color=factor(species)))+
    geom_line(size=2, alpha=0.6)+
      xlim(c(0,498))+
    scale_colour_manual(values = cols)+
      geom_hline(yintercept = 0.05, linetype ="dashed", size =1.5)+
    ylab("C and P density")+
      annotate("text", x= 100,y=1.5, 
               label="~Delta == 0.1",parse=TRUE, size=4)+
   # labs(color="Consumers \n and Predators")+
    theme_classic()+theme(legend.position = "none"))
  
  (dens1<-sol %>% filter(type=="C" | type =="P", time == (tmax-1)) %>%
    ggplot(aes(x=species,y=density, fill = factor(species)))+
    geom_col()+
      labs(fill="Species")+
      annotate("text", x= 7,y=1.5, 
               label="~Delta == 0.1",parse=TRUE, size=4)+
    geom_hline(yintercept = 0.05, linetype ="dashed", size =1.5)+
    scale_fill_manual(values = cols)+
    theme_classic()+theme(legend.position = "none"))
  
  
  
  
  total_density <- sum((sol %>% filter(type=="C" | type =="P", time == (tmax-5)))$density)
  total_density_c <- sum((sol %>% filter(type=="C", time == (tmax-5)))$density)
  total_density_p <- sum((sol %>% filter(type =="P", time == (tmax-5)))$density)
  
  parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                   at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                   delta=forcing_delta, bi=gr_c,alpha_c=0.5)
  parameters$method <- "degree"
  parameters$delta <- 0.1
  
  
  flowField(dynamics_reduced_model, xlim = c(-0.5, 3.5), 
            ylim = c(-0.5, 6.5),
            main=(expression(Delta == 0.1)),
            parameters = parameters, 
            ylab = expression(P[eff]),
            xlab= expression(C[eff]),
            points = 20, add = FALSE)  
  
  nullclines_2(dynamics_reduced_model,xlim = c(-0.5,3.5), 
               ylim = c(-0.5, 6.5),
               parameters = parameters, 
               points = 500, lwd = 3,
               col=c("#E69F00","#009E73"),
               add.legend = TRUE) 
  d<-trajectory(dynamics_reduced_model, y0 = c(0.0001,0.0001), tlim =c(0,500), 
                parameters = parameters,lwd=2 )
  
  
  a1<- recordPlot()  
  
  p1<-ggpubr::ggarrange(a1, timeseries1,dens1,d2, nrow=1,ncol = 4,
                        labels=c("A","B","C","D"))
  
#### saame network delta =5 


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
            duration_mat_P=duration_mat_P, species_index=species_index)
na<-rep(0,S_N) #initial N
na[species_index]<- parms$forcing_delta
nc<-rep(0.00001,S_C) # inital C
np<-rep(0.00001,S_P) # initla P
ic<-c(na,nc,np)
sol1<-ode(func=fwebdyn_perturb_constant, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
  organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe

cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
(timeseries2<- sol1 %>% filter(type=="C" | type =="P") %>% 
    ggplot(aes(x =time,y =density, color=factor(species)))+
    xlim(c(0,499))+
    geom_line(size=2, alpha=0.6)+
    scale_colour_manual(values = cols)+
    geom_hline(yintercept = 0.05, linetype ="dashed", size =2)+
    ylab("C and P density")+
    annotate("text", x= 110,y=2.5, 
             label="~Delta == 5",parse=TRUE, size=4)+
   # labs(color="Consumers \n and Predators")+
    theme_classic()+theme(legend.position = "none"))

(dens2<-sol1 %>% filter(type=="C" | type =="P", time == (tmax-4)) %>%
    ggplot(aes(x=species,y=density, fill = factor(species)))+
    geom_col()+
   # labs(fill="Species")+
    annotate("text", x= 7,y=2.5, 
             label="~Delta ==5",parse=TRUE, size=4)+
    geom_hline(yintercept = 0.05, linetype ="dashed", size =1.5)+
    scale_fill_manual(values = cols)+
    theme_classic()+theme(legend.position = "none"))

total_density <- sum((sol1 %>% filter(type=="C" | type =="P", time == (tmax-5)))$density)
total_density_c <- sum((sol1 %>% filter(type=="C", time == (tmax-5)))$density)
total_density_p <- sum((sol1 %>% filter(type =="P", time == (tmax-5)))$density)


# Plotting the phase plane for each of the parameters in the fact dataframe

parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                 at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                 delta=forcing_delta, bi=gr_c,alpha_c=0.5)
parameters$method <- "degree"
parameters$delta <- 5


flowField(dynamics_reduced_model, xlim = c(-0.5, 3.5), 
          ylim = c(-0.5, 6.5),
          main=(expression(Delta == 5)),
          parameters = parameters, 
          ylab = expression(P[eff]),
          xlab= expression(C[eff]),
          points = 20, add = FALSE)  

nullclines_2(dynamics_reduced_model,xlim = c(-0.5,3.5), 
             ylim = c(-0.5, 6.5),
             parameters = parameters, 
             points = 500, lwd = 3,
             col=c("#E69F00","#009E73"),
             add.legend = TRUE) 
d<-trajectory(dynamics_reduced_model, y0 = c(0.0001,0.0001), tlim =c(0,500), 
              parameters = parameters,lwd=2 )


a2<- recordPlot()  

p2<-ggpubr::ggarrange(a2, timeseries2,dens2, nrow=1,ncol = 4,
                      labels=c("E","F","G"))


TP<-c(rep(0,each=S_N), rep(1,each=S_C),rep(2,each=S_P))
d2<-plotweb(web = webs[[1]],TP = TP,nspecies = 12)


plot_grid(p1,p2,nrow=2,scale=0.85)



#############################
web<-pyramidal.food(S = 14,C = 0.09, prop=c(0.5,0.3,0.2), 
                    balance=0.2)
TP<-c(rep(0,each=7), rep(1,each=4),rep(2,each=3))
d4<-plotweb(web = web$web,TP = TP,nspecies = 14)
d4
#web<-adjmat

prop<-web$index
S_N<- length(prop$Bottom)    # of resource species, N
S_C<- length(prop$Intermediate)  # of consumer species, C
S_P<- length(prop$Top)   # of predator species, P
TP<-c(rep(0,each=S_N), rep(1,each=S_C),rep(2,each=S_P))
d4<-plotweb(web = web[[1]],TP = TP,nspecies = 15)

#web<-list(web=web, index=index)

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
duration <- 500 #fact$forcing_duration[r]
d<- c(rep(1,duration),rep(0,(deltat-duration)))

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


species_index<-c(1,2,3,4,5,6,7) #which((degree) == max(degree))
forcing_delta<-0.1
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
na[species_index]<-forcing_delta
nc<-rep(0.000001,S_C) # inital C
np<-rep(0.000001,S_P) # initla P
ic<-c(na,nc,np)
sol3<-ode(func=fwebdyn_perturb_constant, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
  organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe



cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
(timeseries3<- sol3 %>% filter(type=="C" | type =="P") %>% 
  ggplot(aes(x =time,y =density, color=factor(species)))+
  geom_line(size=2, alpha=0.6)+
  scale_colour_manual(values = cols)+
    geom_hline(yintercept = 0.05, linetype ="dashed", size =2)+
    ylab("C and P density")+
    annotate("text", x= 110,y=1, 
             label="~Delta == 0.1",parse=TRUE, size=4)+
    #labs(color="Consumers \n and Predators")+
    theme_classic()+theme(legend.position = "none"))
    
(dens3<-sol3 %>% filter(type=="C" | type =="P", time == (tmax-2)) %>%
    ggplot(aes(x=species,y=density, fill = factor(species)))+
    geom_col()+
   # labs(fill="Species")+
    annotate("text", x= 8,y=2, 
             label="~Delta == 0.1",parse=TRUE, size=4)+
    geom_hline(yintercept = 0.05, linetype ="dashed", size =2)+
    scale_fill_manual(values = cols)+
    theme_classic()+theme(legend.position = "none"))


total_density <- sum((sol3 %>% filter(type=="C" | type =="P", time == (tmax-5)))$density)
total_density_c <- sum((sol3 %>% filter(type=="C", time == (tmax-5)))$density)
total_density_p <- sum((sol3 %>% filter(type =="P", time == (tmax-5)))$density)

parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                 at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                 delta=forcing_delta, bi=gr_c,alpha_c=0.5)
parameters$method <- "degree"



flowField(dynamics_reduced_model, xlim = c(-0.5, 2), 
          ylim = c(-0.5, 2),
          main=(expression(Delta == 0.1)),
          parameters = parameters, 
          ylab = expression(P[eff]),
          xlab= expression(C[eff]),
          points = 20, add = FALSE)  

nullclines_2(dynamics_reduced_model,xlim = c(-0.5, 2), 
             ylim = c(-0.5, 2),
             parameters = parameters, 
             points = 200, lwd = 3,
             col=c("#E69F00","#009E73"),
             add.legend = TRUE) 
d<-trajectory(dynamics_reduced_model, y0 = c(0.00001,0.00001), tlim =c(0,1000), 
              parameters = parameters )


a3<-recordPlot()

p3<-ggpubr::ggarrange(a3, timeseries3,dens3,d4, nrow=1,ncol = 4 , labels=c("I","J","K","L"))

p3
#findEquilibrium(dynamics_reduced_model,max.iter = 50, y0 = NULL,parameters = parameters)



#### delta 10 for lagoons




species_index<-c(1,2,3,4,5,6,7) #which((degree) == max(degree))
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
            duration_mat_P=duration_mat_P, species_index=species_index)
na<-rep(0,S_N) #initial N
na[species_index]<-forcing_delta
nc<-rep(0.00001,S_C) # inital C
np<-rep(0.00001,S_P) # initla P
ic<-c(na,nc,np)
sol4<-ode(func=fwebdyn_perturb_constant, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
  organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe


cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
(timeseries4<- sol4 %>% filter(type=="C" | type =="P", time < (tmax-1)) %>% 
    ggplot(aes(x =time,y =density, color=factor(species)))+
    geom_line(size=2, alpha=0.6)+
    scale_colour_manual(values = cols)+
    geom_hline(yintercept = 0.05, linetype ="dashed", size =2)+
    ylab("C and P density")+
    annotate("text", x= 110,y=3, 
             label="~Delta == 5",parse=TRUE, size=4)+
  #  labs(color="Consumers \n and Predators")+
    theme_classic()+theme(legend.position = "none"))

(dens4<-sol4 %>% filter(type=="C" | type =="P", time == (tmax-5)) %>%
    ggplot(aes(x=species,y=density, fill = factor(species)))+
    geom_col()+
    #labs(fill="Species")+
    annotate("text", x= 8,y=2, 
             label="~Delta == 5",parse=TRUE, size=4)+
    geom_hline(yintercept = 0.05, linetype ="dashed", size =2)+
    scale_fill_manual(values = cols)+
    theme_classic()+theme(legend.position = "none"))


parameters<-list(h=0.25, eff_c = feed_eff_c, eff_p  =feed_eff_p,web=web,
                 at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                 delta=forcing_delta, bi=gr_c,alpha_c=0.5)
parameters$method <- "degree"



flowField(dynamics_reduced_model, xlim = c(-0.5, 4), 
          ylim = c(-0.5, 4),
          main=(expression(Delta == 5)),
          parameters = parameters, 
          ylab = expression(P[eff]),
          xlab= expression(C[eff]),
          points = 20, add = FALSE)  

nullclines_2(dynamics_reduced_model,xlim = c(-0.5, 4), 
             ylim = c(-0.5, 4),
             parameters = parameters, 
             points = 200, lwd = 3,
             col=c("#E69F00","#009E73"),
             add.legend = TRUE) 
d<-trajectory(dynamics_reduced_model, y0 = c(0.0001,0.0001), tlim =c(0,1000), 
              parameters = parameters )


a4<-recordPlot()

p4<-ggpubr::ggarrange(a4, timeseries4,dens4, nrow=1,ncol = 4, labels=c("M","N","O"))

p4




pdf(file = "figure_2_new1.pdf",   # The directory you want to save the file in
   width = 16, # The width of the plot in inches
   height = 14) # The height of the plot in inches
#ggpubr::ggarrange(p1,p2,p3,p4,nrow = 4)
plot_grid(p1,p2,p3,p4,nrow=4,scale=0.86)

# Step 3: Run dev.off() to create the file!
dev.off()
