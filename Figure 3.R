
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
source("01_functions_updated.R")
source("food_web_generator.R")
source("plot_foodweb.R")
source("phase_plane_fweb.R")
fact<-expand.grid( C=0.09,
                   balance=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                   S=c(18), delta=c(3),
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
for(i in 1:nrow(fact)){
  

  
  web<-pyramidal.food(S = fact$S[i],C = fact$C[i],prop=c(0.4,0.4,0.2), 
                      balance=fact$balance[i])
  webs[i]<- list(web$web)
  prop<-web$index
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  alpha.c<-matrix((runif( S_C*S_C, 0.01,0.01)), nrow=S_C, ncol= S_C)
  alpha.n<-matrix((runif( S_N*S_N, 0.005,0.005)), nrow=S_N, ncol= S_N)
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
  gr_p<- -0.1  #death rate pradtor
  mortality_n<-0.5
  tmax<-1000
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 1000 #fact$forcing_duration[r]
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
  species_index<-c(1,2,3,4,5,6,7) #which((degree) == max(degree))
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
  nc<-rep(0.0001,S_C) # inital C
  np<-rep(0.0001,S_P) # initla P
  ic<-c(na,nc,np)
  sol<-ode(func=fwebdyn_perturb_constant, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
    organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe
 
  
  
 fact$foodweb_biomass[i]<-sum( (sol %>% filter(time == tmax))$density)
  fact$foodweb_richness[i]<- length(which((sol %>% filter(time == tmax))$density > 0.1))
  fact$intermediate_biomass[i]<-sum( (sol %>% filter(time == tmax, type == "C"))$density)
  fact$top_biomass[i]<-sum( (sol %>% filter(time == tmax, type == "P"))$density)
  fact$intermediate_richness[i]<-length(which((sol %>% filter(time == tmax, type =="C"))$density > 0.1))
  fact$top_richness[i]<-length(which((sol %>% filter(time == tmax, type =="P"))$density > 0.1))
  fact$top_inter_richness[i]<-length(which((sol %>% filter(time == tmax, type=="C"| type =="P"))$density > 0.1))
  
  print(i)
  
  # Plotting the phase plane for each of the parameters in the fact dataframe

  
  parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,
                   at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                   delta=fact$delta[i], bi=gr_c,alpha_c=0.02)
  
  dat<-phase_plot(parameters =parameters, web=web, delta = fact$delta[i] )
  
  temp<-rbind(cbind(dat$dna,dat$dnp,
                    dat$na,
                    dat$np,
                    rep(fact$random_seed[i],each=length(dat$dna)),
                    rep(fact$delta[i],each=length(dat$dna)),
                    rep(fact$balance[i], each = length(dat$dna)),
                    rep(fact$S[i],each=length(dat$dna)),
                    rep(fact$C[i],each=length(dat$dna))
  ))
  
  
  data_t<-rbind(data_t,temp)
}


(r1<-fact %>% 
  ggplot( aes(x=balance, y = (top_inter_richness)/11, color=factor(random_seed)))+
    geom_point(size =6, alpha =0.7, color = "firebrick")+
    theme_classic()+
    labs(color="")+
    ylim(c(0,1))+
    theme(legend.position = "none")+
    geom_smooth(method = "lm", formula = y~x, se=FALSE,size=1.25, color = "grey")+
    xlab("Proportion of top predator links")+
   ylab("Recoverability (fraction richness)"))#+
  #facet_wrap(.~Connecrance))

(r2<-fact %>% 
    ggplot( aes(x=balance, y = (intermediate_biomass+top_biomass), color=factor(random_seed)))+
    geom_point(size =6, alpha =0.7, color = "firebrick")+
    theme_classic()+
    labs(color="")+
    theme(legend.position = "none")+
    geom_smooth(method = "lm", formula = y~x, se=FALSE,size=1.25, color = "grey")+
    xlab("Proportion of top predator links")+
    #geom_point(x = 0.2, y = 157, color = "darkgreen", size=6)+
   # geom_point(x = 0.6, y = 22.6, color = "darkgreen", size=6)+
   # geom_point(x = 1, y = 8, color = "darkgreen", size=6)+
    ylab("Recoverability (biomass)"))#+


data_t<-as.data.frame(data_t)
na.omit(data_t)
colnames(data_t)<-c("dnC", "dnP","Nc","Np","reps", "Delta","Balance", "S", "Connectance")
str(data_t)

data_t$Delta<-as.factor(data_t$Delta)
(p1<-data_t %>%filter(Balance==0.2) %>% ggplot(aes(x=Nc,y=dnC, color=(Delta)))+
  geom_line(aes(y=Np,x=dnP), size=1.5, color="#009E73")+
    labs(color=expression(Delta))+
    geom_line(aes(col=factor(Delta)),size =1.5, color ="#E69F00")+
 
    theme_classic()+
  geom_hline(yintercept = 0, linetype ="dashed", size=1.2)+
  xlab(expression(P))+ylab(expression(C))+
    xlim(c(0,80))+ylim(c(-20,20))+scale_color_brewer(palette = 'Dark2')+
  geom_point(x = 10, y = 1.65, color = "black", size=4, pch=7)+
  annotate("text", x= 80,y=50, label="Prop. of top predator links = 0.2", size=3))

(p2<-data_t %>%filter( Balance==0.5) %>% ggplot(aes(x=Nc,y=dnC, color=Delta))+
    geom_line(aes(y=Np,x=dnP), size=1.5, color="#009E73")+
    geom_line(aes(col=factor(Delta)),size =1.5,color ="#E69F00")+
    labs(color=expression(Delta))+
    theme_classic()+
    geom_hline(yintercept = 0, linetype ="dashed", size=1.2)+
    xlab(expression(P))+ylab(expression(C))+
    xlim(c(0,80))+ylim(c(-20,20))+scale_color_brewer(palette = 'Dark2')+
    #geom_point(x = 2.75, y = 0.7, color = "black", size=4, pch=7)+
   annotate("text",  x= 80,y=50, label="Prop. of top predator links = 0.5", size=3))

(p3<-data_t %>%filter( Balance==1) %>% ggplot(aes(x=Nc,y=dnC, color=Delta))+
    geom_line(aes(y=Np,x=dnP), size=1.5, color="#009E73")+
    geom_line(aes(col=factor(Delta)),size =1.5,color ="#E69F00")+
    labs(color=expression(Delta))+
    theme_classic()+
    geom_hline(yintercept = 0, linetype ="dashed", size=1.2)+
    xlab(expression(P))+ylab(expression(C))+
    xlim(c(0,80))+ylim(c(-20,20))+scale_color_brewer(palette = 'Dark2')+
   # geom_point(x = 1, y = 0.05, color = "black", size=4, pch=7)+
    annotate("text",  x= 80,y=50, label="Prop. of top predator links = 1", size=3))

### plotting foodwebs
TP<-c(rep(0,each=S_N), rep(1,each=S_C),rep(2,each=S_P))
d1<-plotweb(web = webs[[1]],TP = TP,nspecies = 18)
d2<-plotweb(web = webs[[5]],TP = TP,nspecies = 18)
d3<-plotweb(web = webs[[9]],TP = TP,nspecies = 18)


gg<-ggpubr::ggarrange(ggarrange(d1,d2,d3,ncol = 3, 
                    labels=c("A","B","C")),
          ggarrange(r2,r1,labels ="D"),nrow=3,
          ggarrange(p1,p2,p3,ncol=3,labels = c("E","F","G")), 
          heights = c(1.5, 2,2)
          )

gg



## empirical foodweb rocky shore california foodweb

adjmat<-read.csv("rocky_shore_cal.csv",sep=",", header = FALSE)
adjmat<-as.matrix(adjmat)
adjmat<-unname(adjmat)
TP<-c(rep(0,each=4), rep(1,each=7),rep(2,each=3))
d4<-plotweb(web = adjmat,TP = TP,nspecies = 14)
web<-adjmat

S_N<- 4  # of resource species, N
S_C<- 7  # of consumer species, C
S_P<- 3   # of predator species, P
prop$Bottom<- seq(1, S_N )
prop$Intermediate <- seq(5, ((5+S_C)-1))
prop$Top<- seq(12, ((S_P+12)-1))

index<-NULL
index$Bottom <- seq(1, S_N )
index$Intermediate<- seq(5, ((5+S_C)-1))
index$Top<-seq(12, ((S_P+12)-1))
web<-list(web=web, index=index)

alpha.c<-matrix((runif( S_C*S_C, 0.01,0.01)), nrow=S_C, ncol= S_C)
alpha.n<-matrix((runif( S_N*S_N, 0.005,0.005)), nrow=S_N, ncol= S_N)
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
gr_p<- -0.1  #death rate pradtor
mortality_n<-0.5
tmax<-1000
time_range<-c(0,tmax)
deltat<- (time_range[2]-time_range[1])/1 + 1
duration <- 1000 #fact$forcing_duration[r]
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


species_index<-c(1,2,3,4) #which((degree) == max(degree))
forcing_delta<-fact$delta[1]
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
nc<-rep(0.0001,S_C) # inital C
np<-rep(0.0001,S_P) # initla P
ic<-c(na,nc,np)
sol<-ode(func=fwebdyn_perturb_constant, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
  organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe


cols<- c(rep("#E69F00" , each= S_C), rep("#009E73", each=S_P))
timeseries<- sol %>% filter(type=="C" | type =="P", density > 0.5 ) %>% 
  ggplot(aes(x =time,y =density, color=factor(species)))+
    geom_line(size=2, alpha=0.6)+
   scale_colour_manual(values = cols)+
  ylab("C and P density")+
  # labs(color="species")+
    theme_classic()+
 theme(legend.position = "none")+
  annotate("text", x= 600,y=10, 
           label="Prop. of top predator links = 0.52", size=3)
  

 
 parameters<-list(h=h, eff_c = feed_eff_c, eff_p  =feed_eff_p,
                  at_c=attack_rate_c, at_p=attack_rate_p, allee_c=allee_c, bp=gr_p,
                  delta=forcing_delta, bi=gr_c,alpha_c=0.02)
 
 dat<-phase_plot(parameters =parameters, web=web, delta = forcing_delta )
 
 
 prop.predator.links<- 2*sum(web$web[web$index$Intermediate,web$index$Top])/sum(web$web)
 temp<-rbind(cbind(dat$dna,dat$dnp,
                   dat$na,
                   dat$np,
                   rep(forcing_delta,each=length(dat$dna)),
                   rep(prop.predator.links, each=length(dat$dna)),
                   rep(fact$C[i],each=length(dat$dna))
 ))
 
 
na.omit(temp)
temp<-as.data.frame(temp)
 colnames(temp)<-c("dnC", "dnP","Nc","Np","Delta","Prop_predator_links", "Connectance")
 str(temp)
 
 
 
 (p4<-temp %>%ggplot(aes(x=Nc,y=dnC, color=(Delta)))+
     geom_line(aes(y=Np,x=dnP), size=1.5, color="#009E73")+
     labs(color=expression(Delta))+
     geom_line(aes(col=factor(Delta)),size =1.5, color ="#E69F00")+
     
     theme_classic()+
     geom_hline(yintercept = 0, linetype ="dashed", size=1.2)+
     xlab(expression(P))+ylab(expression(C))+
     xlim(c(0,100))+ylim(c(-20,30))+
     #geom_point(x = 3, y = 0.8, color = "black", size=4, pch=7)+
     annotate("text", x= 80,y=50, label="Prop. of top predator links = 0.52", size=3))
 
 
 d4
 
 gg<-ggpubr::ggarrange(ggarrange(d1,d2,d3,ncol = 3, 
                                 labels=c("A","B","C")),
                       ggarrange(r2,labels ="D"),nrow=4,
                       ggarrange(p1,p2,p3,ncol=3,labels = c("E","F","G")), 
                       ggarrange(d4,p4,timeseries, ncol=3,labels = c("H", "I", "J")),
                       heights = c(2, 2,2)
 )
 
 gg
 
 annotate_figure(gg,
                 bottom = text_grob("Food-web data: \n Rocky Shore, California", color = "black",
                                    hjust = 1, x = 0.23,face="italic", size = 11),
                 
 )
 gg
 