#code for Ecology Letters paper: 
#" The impact of individual variation on abrupt collapses in mutualistic networks" 2021. Gaurav Baruah
#email: gbaruahecoevo@gmail.com


cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)



network_structure<-function(Na,Np, g ){
  new_g<- matrix(0, nrow= length(Np), ncol =length(Na) )
  
  Na[which(Na < 5e-1)]<-0
  Np[which(Np < 5e-1)]<-0
  
  for(i in 1:length(Np)){
    for(j in 1:length(Na)){
      new_g[i,j]<-g[i,j]*Na[j]*Np[i]    
      
    }
  }
  new_g[which(new_g > 0)]<-1
  
  return(new_g) 
}


#conversion to a matrix
adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(unname(d))
  dat[dat > 0] = 1
  dat[dat < 0] = 0
  dat<-apply(dat,2,as.numeric)
  return(dat)}



#time: time is the vector at which the ode is solved
#y : the states of the system
# pars: list of parameters that are needed for the ODE solver.
# b: is the growth rate vector of resources N, and consumers C
# alpha_n : intra and interspecific competition matrix for resources N
# alpha_c: intra and interspecific competition matrix for consumers C
# mu: strong allee threshold for consumers C.
# A : adjacency matrix of interactions for the food-web: only has 0s and 1s.
# h : handling time which is fixed at 0.2
# e: efficiency in coversion of resources to energy for consumers and predators, 0.5
# cr: consumption rate of consumers and predators, 0.5
# a_p: density-dependence/intraspecific competition in top predator which is 0.01
#d : death rate of top predators
fwebdyn <- function(time, y, pars){
  A<-pars$Adjacency
  prop<-pars$prop # this function takes the web and gives you the
  # the indexes of the basal/resource species, consumer intermediate, and top predator indexes in the matrix
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  
  S<-S_N+S_P + S_C #total species in the web
  na <- y[1:S_N] ## species densities of resources
  nc <- y[(S_N+1): (S_C+S_N)] ## species densities of plants
  np <- y[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  #competition matrices for all the species i.e., inter and intraspecific competition terms
  
  alpha.n<- pars$alpha.n 
  diag(alpha.n)<- 0.1 #intraspecific competition
  alpha.c<-pars$alpha.c # pars$Cmatrix # competition matrix for consumers
  diag(alpha.c)<- 0.05 #intraspecific competition
  bi<-bc<-bp<-numeric()
  allee_b<- pars$allee_b # allee threshold basal, 0
  allee_c<- pars$allee_c #allee threshold consumers, 0.25  
  allee_p<- pars$allee_p
  h<-pars$h # handling time, 0.15
  a_r_c<- pars$attack_rate_c # 1.5
  a_r_p<- pars$attack_rate_p  # 1.5
  f_eff_c<- pars$feed_eff_c # 0.2
  f_eff_p<- pars$feed_eff_p # 0.25
  bi<- pars$gr_n # 1
  bc<-pars$gr_c # 0.05
  bp<-pars$gr_p # -0.001
  u_transition_p<-pars$mortality_n # 0
  a1 <- matrix(0, nrow=S_N, ncol = S_C)
  a2 <- matrix(0, nrow=S_C, ncol = S_N)
  a3 <- matrix(0, nrow=S_C,ncol = S_P)
  a4 <- matrix(0, nrow=S_P, ncol=S_C)
  summed_a1<-summed_a2<-summed_a3<-summed_a4<-numeric()
  
#for( t in 1:100){  
  ##predation of resources by consumers
  for(i in 1:S_N){
    for(j in 1:S_C){
      a1[i,j]<- pars$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j]/(1+h*pars$attack_rate_c*sum(A[prop$Bottom[i],prop$Intermediate]*nc))
    }
    summed_a1[i]  <- sum(a1[i,])
  }
  
  for(r in 1:S_C){
    for(l in 1:S_N){
      a2[r,l]<- pars$attack_rate_c* pars$feed_eff_c*A[prop$Bottom[l],prop$Intermediate[r]]*na[l]/(1+h*pars$attack_rate_c*A[prop$Bottom[l],prop$Intermediate[r]]*nc[r])
    }
      summed_a2[r]<- sum(a2[r,])
}
  for(m in 1:S_C){
    for(n in 1:S_P){
      a3[m,n]<-pars$attack_rate_p*A[prop$Intermediate[m],prop$Top[n]]*np[n]/(1+h*pars$attack_rate_p*A[prop$Intermediate[m],prop$Top[n]]*np[n])
      summed_a3[m] <-   sum(a3[m,]) 
      }
  }
  for(o in 1:S_P){
    for(p in 1:S_C){
      a4[o,p]<-pars$attack_rate_p*pars$feed_eff_p*A[prop$Intermediate[p],prop$Top[o]]*nc[p]/(1+h*pars$attack_rate_p*A[prop$Intermediate[p],prop$Top[o]]*np[o])
    }
    summed_a4[o]<- sum(a4[o,])
  }
  
  
  #predation of consumers by top predators
  
  
  #equations for population dynamics
  dndt <- na*(bi - alpha.n%*%(na) - summed_a1 -pars$mortality_n*na )*cutoff(na/(1e-4)) # + perturb
  dcdt <- nc*(bc - alpha.c%*%(nc) + summed_a2 - summed_a3)*cutoff(nc/(1e-4)) #
  dpdt <- np*(bp -0.001*np + summed_a4)*cutoff(np/1e-4)
  
 # na<-na+dndt
 # nc <- nc+dcdt
 # np<-np+dpdt
 # 
#}
  return(list(c(dndt, dcdt,dpdt)))
}



#time: time is the vector at which the ode is solved
#y : the states of the system
# pars: list of parameters that are needed for the ODE solver.
# b: is the growth rate vector of resources N, and consumers C
# alpha_n : intra and interspecific competition matrix for resources N
# alpha_c: intra and interspecific competition matrix for consumers C
# mu: strong allee threshold for consumers C.
# A : adjacency matrix of interactions for the food-web: only has 0s and 1s.
# h : handling time which is fixed at 0.2
# e: efficiency in coversion of resources to energy for consumers and predators, 0.5
# cr: consumption rate of consumers and predators, 0.5
# a_p: density-dependence/intraspecific competition in top predator which is 0.01
#d : death rate of top predators
fwebdyn_perturb <- function(time, y, parms){
  A<-parms$Adjacency
  prop<-parms$prop # this function takes the web and gives you the
  # the indexes of the basal/resource species, consumer intermediate, and top predator indexes in the matrix
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  
  S<-S_N+S_P + S_C #total species in the web
  na <- y[1:S_N] ## species densities of resources
  nc <- y[(S_N+1): (S_C+S_N)] ## species densities of plants
  np <- y[(S_N+S_C+1):(S_C+S_N+S_P)]
  #
  #competition matrices for all the species i.e., inter and intraspecific competition terms
  
  alpha.n<- parms$alpha.n 
  diag(alpha.n)<- 1 #intraspecific competition
  alpha.c<-parms$alpha.c # parms$Cmatrix # competition matrix for consumers
  diag(alpha.c)<- 0.01 #intraspecific competition
  bi<-bc<-bp<-numeric()
  allee_b<- parms$allee_b # allee threshold basal, 0
  allee_c<- parms$allee_c #allee threshold consumers, 0.25  
  allee_p<-parms$allee_p
  h<-parms$h # handling time, 0.15
  a_r_c<- parms$attack_rate_c # 1.5
  a_r_p<- parms$attack_rate_p  # 1.5
  f_eff_c<- parms$feed_eff_c # 0.2
  f_eff_p<- parms$feed_eff_p # 0.25
  bi<- parms$gr_n # 1
  bc<-parms$gr_c # 0.05
  bp<-parms$gr_p # -0.001
  u_transition_p<-parms$mortality_n # 0
  a1 <- matrix(0, nrow=S_N, ncol = S_C)
  a2 <- matrix(0, nrow=S_C, ncol = S_N)
  a3 <- matrix(0, nrow=S_C,ncol = S_P)
  a4 <- matrix(0, nrow=S_P, ncol=S_C)
  summed_a1<-summed_a2<-summed_a3<-summed_a4<-numeric()
  
  time_func_N<-approxfun(parms$t1_N,method="linear", rule =2)
  time_func_C<-approxfun(parms$t1_C,method="linear", rule =2)
  time_func_P<-approxfun(parms$t1_P,method="linear", rule =2)
  
  
  
  e_N<-replicate(S_N,time_func_N(time))
  e_C<-replicate(S_C,time_func_C(time))
  e_P<-replicate(S_P,time_func_P(time))
  
  forcing_strength_N<-parms$forcing_strength[1:S_N]
  forcing_strength_C<-parms$forcing_strength[(S_N+1): (S_C+S_N)]
  forcing_strength_P<-parms$forcing_strength[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  forcing_index<-rep(0,(S_N+S_C+S_P))
  forcing_index[parms$species_index]<-1
  forcing_index_sp_N<-forcing_index[1:S_N]
  forcing_index_sp_C<-forcing_index[(S_N+1): (S_C+S_N)]
  forcing_index_sp_P<-forcing_index[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  ##predation of resources by consumers
  for(i in 1:S_N){
    for(j in 1:S_C){
      a1[i,j]<- parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j]/(1+h*parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j])
    }
    summed_a1[i]  <- sum(a1[i,])
  }
  
  for(r in 1:S_C){
    for(l in 1:S_N){
      a2[r,l]<- parms$attack_rate_c* parms$feed_eff_c*A[prop$Bottom[l],prop$Intermediate[r]]*na[l]/(1+h*parms$attack_rate_c*A[prop$Bottom[l],prop$Intermediate[r]]*nc[r])
    }
    summed_a2[r]<- sum(a2[r,])
  }
  for(m in 1:S_C){
    for(n in 1:S_P){
      a3[m,n]<-parms$attack_rate_p*A[prop$Intermediate[m],prop$Top[n]]*np[n]/(1+h*parms$attack_rate_p*A[prop$Intermediate[m],prop$Top[n]]*np[n])
      summed_a3[m] <-   sum(a3[m,]) 
    }
  }
  for(o in 1:S_P){
    for(p in 1:S_C){
      a4[o,p]<-parms$attack_rate_p*parms$feed_eff_p*A[prop$Intermediate[p],prop$Top[o]]*nc[p]/(1+h*parms$attack_rate_p*A[prop$Intermediate[p],prop$Top[o]]*np[o])
    }
    summed_a4[o]<- sum(a4[o,])
  }
  
  
  #predation of consumers by top predators
  
  #equations for population dynamics
 
  dndt <- na*(bi - alpha.n%*%na  - summed_a1 -parms$mortality_n)*cutoff(na/(1e-5)) +  forcing_index_sp_N*forcing_strength_N*e_N    # + perturb
  dcdt <- nc*(bc- alpha.c%*%(nc) + summed_a2 - summed_a3 +  forcing_index_sp_C*forcing_strength_C*e_C)*cutoff(nc/(1e-5))+rnorm(5,0,0.1) #
  dpdt <- np*(bp  + summed_a4 +  forcing_index_sp_P*forcing_strength_P*e_P)*cutoff(np/1e-5)+rnorm(1,0,0.1)
  

  return(list(c(dndt, dcdt,dpdt)))
}






#time: time is the vector at which the ode is solved
#y : the states of the system
# pars: list of parameters that are needed for the ODE solver.
# b: is the growth rate vector of resources N, and consumers C
# alpha_n : intra and interspecific competition matrix for resources N
# alpha_c: intra and interspecific competition matrix for consumers C
# mu: strong allee threshold for consumers C.
# A : adjacency matrix of interactions for the food-web: only has 0s and 1s.
# h : handling time which is fixed at 0.2
# e: efficiency in coversion of resources to energy for consumers and predators, 0.5
# cr: consumption rate of consumers and predators, 0.5
# a_p: density-dependence/intraspecific competition in top predator which is 0.01
#d : death rate of top predators
fwebdyn_perturb_constant <- function(time, y, parms){
  A<-parms$Adjacency
  prop<-parms$prop # this function takes the web and gives you the
  # the indexes of the basal/resource species, consumer intermediate, and top predator indexes in the matrix
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  
  S<-S_N+S_P + S_C #total species in the web
  na <- y[1:S_N] ## species densities of resources
  nc <- y[(S_N+1): (S_C+S_N)] ## species densities of plants
  np <- y[(S_N+S_C+1):(S_C+S_N+S_P)]
  #
  #na[parms$species_index]<- parms$forcing_delta
  
  #competition matrices for all the species i.e., inter and intraspecific competition terms
  
  alpha.n<- parms$alpha.n 
  diag(alpha.n)<- 1 #intraspecific competition
  alpha.c<-parms$alpha.c # parms$Cmatrix # competition matrix for consumers
  diag(alpha.c)<- 0.5 #intraspecific competition
  bi<-bc<-bp<-numeric()
  allee_b<- parms$allee_b # allee threshold basal, 0
  allee_c<- parms$allee_c #allee threshold consumers, 0.25  
  allee_p<-parms$allee_p
  h<-parms$h # handling time, 0.15
  a_r_c<- parms$attack_rate_c # 1.5
  a_r_p<- parms$attack_rate_p  # 1.5
  f_eff_c<- parms$feed_eff_c # 0.2
  f_eff_p<- parms$feed_eff_p # 0.25
  bi<- parms$gr_n # 1
  bc<-parms$gr_c # 0.05
  bp<-parms$gr_p # -0.001
  u_transition_p<-parms$mortality_n # 0
  a1 <- matrix(0, nrow=S_N, ncol = S_C)
  a2 <- matrix(0, nrow=S_C, ncol = S_N)
  a3 <- matrix(0, nrow=S_C,ncol = S_P)
  a4 <- matrix(0, nrow=S_P, ncol=S_C)
  summed_a1<-summed_a2<-summed_a3<-summed_a4<-numeric()
  
  time_func_N<-approxfun(parms$t1_N,method="linear", rule =2)
  time_func_C<-approxfun(parms$t1_C,method="linear", rule =2)
  time_func_P<-approxfun(parms$t1_P,method="linear", rule =2)
  
  
  
  e_N<-replicate(S_N,time_func_N(time))
  e_C<-replicate(S_C,time_func_C(time))
  e_P<-replicate(S_P,time_func_P(time))
  na<- e_N*na
  y[1:S_N]<-na
  forcing_strength_N<-parms$forcing_strength[1:S_N]
  forcing_strength_C<-parms$forcing_strength[(S_N+1): (S_C+S_N)]
  forcing_strength_P<-parms$forcing_strength[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  forcing_index<-rep(0,(S_N+S_C+S_P))
  forcing_index[parms$species_index]<-1
  forcing_index_sp_N<-forcing_index[1:S_N]
  forcing_index_sp_C<-forcing_index[(S_N+1): (S_C+S_N)]
  forcing_index_sp_P<-forcing_index[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  ##predation of resources by consumers
  for(i in 1:S_N){
    for(j in 1:S_C){
      a1[i,j]<- parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j]/(1+h*parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j])
    }
    summed_a1[i]  <- sum(a1[i,])
  }
  
  for(r in 1:S_C){
    for(l in 1:S_N){
      a2[r,l]<- parms$attack_rate_c* parms$feed_eff_c*A[prop$Bottom[l],prop$Intermediate[r]]*na[l]/(1+h*parms$attack_rate_c*sum(A[prop$Bottom,prop$Intermediate[r]]*na))
    }
    summed_a2[r]<- sum(a2[r,])
  }
  for(m in 1:S_C){
    for(n in 1:S_P){
      a3[m,n]<-parms$attack_rate_p*A[prop$Intermediate[m],prop$Top[n]]*np[n]/(1+h*parms$attack_rate_p*sum(A[prop$Intermediate,prop$Top[n]]*nc))
      summed_a3[m] <-   sum(a3[m,]) 
    }
  }
  for(o in 1:S_P){
    for(p in 1:S_C){
      a4[o,p]<-parms$attack_rate_p*parms$feed_eff_p*A[prop$Intermediate[p],prop$Top[o]]*nc[p]/(1+h*parms$attack_rate_p*sum(A[prop$Intermediate,prop$Top[o]]*nc))
    }
    summed_a4[o]<- sum(a4[o,])
  }
  
  #predation of consumers by top predators
  
  #equations for population dynamics
 # nc_d_noise<- rnorm(S_C, 0,0.01)
  ##=
 # nc_d_noise<- rnorm(S_C, 0,0.01)
  #np_d_noise<- rnorm(S_P,0,0.01)
  
  dndt <- rep(0, S_N) # na*(bi - alpha.n%*%na  - summed_a1)*cutoff(na/(1e-5))     # + perturb
  dcdt <- nc*(bc- alpha.c%*%(nc) + summed_a2 - summed_a3 +  forcing_index_sp_C*forcing_strength_C*e_C )*cutoff(nc/(1e-6))
  dpdt <- np*(bp -0.1*np + summed_a4 + forcing_index_sp_P*forcing_strength_P*e_P)*cutoff(np/1e-6)

  
  return(list(c(dndt, dcdt,dpdt)))
}



fwebdyn_perturb_constant_noise <- function(time, y, parms){
  A<-parms$Adjacency
  prop<-parms$prop # this function takes the web and gives you the
  # the indexes of the basal/resource species, consumer intermediate, and top predator indexes in the matrix
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  
  S<-S_N+S_P + S_C #total species in the web
  na <- y[1:S_N] ## species densities of resources
  nc <- y[(S_N+1): (S_C+S_N)] ## species densities of plants
  np <- y[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  #na[parms$species_index]<- parms$forcing_delta
  
  #competition matrices for all the species i.e., inter and intraspecific competition terms
  
  alpha.n<- parms$alpha.n 
  diag(alpha.n)<- 1 #intraspecific competition
  alpha.c<-parms$alpha.c # parms$Cmatrix # competition matrix for consumers
  diag(alpha.c)<- 0.5 #intraspecific competition
  bi<-bc<-bp<-numeric()
  allee_b<- parms$allee_b # allee threshold basal, 0
  allee_c<- parms$allee_c #allee threshold consumers, 0.25  
  allee_p<-parms$allee_p
  h<-parms$h # handling time, 0.15
  a_r_c<- parms$attack_rate_c # 1.5
  a_r_p<- parms$attack_rate_p  # 1.5
  f_eff_c<- parms$feed_eff_c # 0.2
  f_eff_p<- parms$feed_eff_p # 0.25
  bi<- parms$gr_n # 1
  bc<-parms$gr_c # 0.05
  bp<-parms$gr_p # -0.001
  u_transition_p<-parms$mortality_n # 0
  a1 <- matrix(0, nrow=S_N, ncol = S_C)
  a2 <- matrix(0, nrow=S_C, ncol = S_N)
  a3 <- matrix(0, nrow=S_C,ncol = S_P)
  a4 <- matrix(0, nrow=S_P, ncol=S_C)
  summed_a1<-summed_a2<-summed_a3<-summed_a4<-numeric()
  
  time_func_N<-approxfun(parms$t1_N,method="linear", rule =2)
  time_func_C<-approxfun(parms$t1_C,method="linear", rule =2)
  time_func_P<-approxfun(parms$t1_P,method="linear", rule =2)
  
  #E_c<-as.numeric(noise_of_consumer(time,parms$noise_consumer,S_C))
  #E_p<-as.numeric(noise_of_top(time,parms$noise_top,S_P))
  noise_f_C0<-approxfun(parms$noise_c[,2],method="linear",rule=2)
  #noise_f_C1<-approxfun(parms$noise_c[,2][,2],method="linear",rule=2)
  #noise_f_C2<-approxfun(parms$noise_c[,2][,3],method="linear",rule=2)
  #noise_f_C3<-approxfun(parms$noise_c[,2][,4],method="linear",rule=2)
  
  noise_f_P1<-approxfun(parms$noise_p[,2],method="linear", rule=2)
  #noise_f_P2<-approxfun(parms$noise_p[,2][,2],method="linear", rule=2)
  #noise_f_P3<-approxfun(parms$noise_p[,2][,3],method="linear", rule=2)
  
  e_N<-replicate(S_N,time_func_N(time))
  
  e_C_1<-replicate(S_C,time_func_C(time))
  e_P<-replicate(S_P,time_func_P(time))
  
  E_c<-replicate(S_C, noise_f_C0(time))
 # E_c2<-replicate(S_C, noise_f_C1(time))
  #E_c3<-replicate(S_C, noise_f_C2(time))
  #E_c4<-replicate(S_C, noise_f_C3(time))
  
  E_p1<-replicate(S_P,noise_f_P1(time))
  #E_p2<-replicate(S_P,noise_f_P2(time))
  #E_p3<-replicate(S_P,noise_f_P3(time))
  
  E_p<-c(E_p1)#,E_p2[1], E_p3[1])
  E_c<-c(E_c)#,E_c2[1],E_c3[1],E_c4[1])
  forcing_strength_N<-parms$forcing_strength[1:S_N]
  forcing_strength_C<-parms$forcing_strength[(S_N+1): (S_C+S_N)]
  forcing_strength_P<-parms$forcing_strength[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  forcing_index<-rep(0,(S_N+S_C+S_P))
  forcing_index[parms$species_index]<-1
  forcing_index_sp_N<-forcing_index[1:S_N]
  forcing_index_sp_C<-forcing_index[(S_N+1): (S_C+S_N)]
  forcing_index_sp_P<-forcing_index[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  ##predation of resources by consumers
  ##predation of resources by consumers
  ##predation of resources by consumers
  for(i in 1:S_N){
    for(j in 1:S_C){
      a1[i,j]<- parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j]/(1+h*parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j])
    }
    summed_a1[i]  <- sum(a1[i,])
  }
  
  for(r in 1:S_C){
    for(l in 1:S_N){
      a2[r,l]<- parms$attack_rate_c* parms$feed_eff_c*A[prop$Bottom[l],prop$Intermediate[r]]*na[l]/(1+h*parms$attack_rate_c*sum(A[prop$Bottom,prop$Intermediate[r]]*na))
    }
    summed_a2[r]<- sum(a2[r,])
  }
  for(m in 1:S_C){
    for(n in 1:S_P){
      a3[m,n]<-parms$attack_rate_p*A[prop$Intermediate[m],prop$Top[n]]*np[n]/(1+h*parms$attack_rate_p*sum(A[prop$Intermediate,prop$Top[n]]*nc))
      summed_a3[m] <-   sum(a3[m,]) 
    }
  }
  for(o in 1:S_P){
    for(p in 1:S_C){
      a4[o,p]<-parms$attack_rate_p*parms$feed_eff_p*A[prop$Intermediate[p],prop$Top[o]]*nc[p]/(1+h*parms$attack_rate_p*sum(A[prop$Intermediate,prop$Top[o]]*nc))
    }
    summed_a4[o]<- sum(a4[o,])
  }  
  #predation of consumers by top predators
  
  #equations for population dynamics
 # nc_d_noise<- rnorm(S_C, 0,0.001)
 # np_d_noise<- rnorm(S_P,0,0.001)
  dndt <- rep(0, S_N) #na*(bi - alpha.n%*%na  - summed_a1 -parms$mortality_n)*cutoff(na/(1e-5)) +  forcing_index_sp_N*forcing_strength_N*e_N    # + perturb
  dcdt <- nc*(bc- alpha.c%*%(nc) +E_c+ summed_a2 - summed_a3 +  forcing_index_sp_C*forcing_strength_C*e_C_1)*cutoff(nc/(1e-6))
  dpdt <- np*(bp  -0.1*np + summed_a4 +E_p + forcing_index_sp_P*forcing_strength_P*e_P)*cutoff(np/1e-6)

  
  return(list(c(dndt, dcdt,dpdt)))
}



#time: time is the vector at which the ode is solved
#y : the states of the system
# pars: list of parameters that are needed for the ODE solver.
# b: is the growth rate vector of resources N, and consumers C
# alpha_n : intra and interspecific competition matrix for resources N
# alpha_c: intra and interspecific competition matrix for consumers C
# mu: strong allee threshold for consumers C.
# A : adjacency matrix of interactions for the food-web: only has 0s and 1s.
# h : handling time which is fixed at 0.2
# e: efficiency in coversion of resources to energy for consumers and predators, 0.5
# cr: consumption rate of consumers and predators, 0.5
# a_p: density-dependence/intraspecific competition in top predator which is 0.01
#d : death rate of top predators
fwebdyn_perturb_continuous <- function(time, y, parms){
  A<-parms$Adjacency
  prop<-parms$prop # this function takes the web and gives you the
  # the indexes of the basal/resource species, consumer intermediate, and top predator indexes in the matrix
  S_N<- length(prop$Bottom)    # of resource species, N
  S_C<- length(prop$Intermediate)  # of consumer species, C
  S_P<- length(prop$Top)   # of predator species, P
  
  S<-S_N+S_P + S_C #total species in the web
  na <- y[1:S_N] ## species densities of resources
  nc <- y[(S_N+1): (S_C+S_N)] ## species densities of plants
  np <- y[(S_N+S_C+1):(S_C+S_N+S_P)]
  #
  #competition matrices for all the species i.e., inter and intraspecific competition terms
  
  alpha.n<- parms$alpha.n 
  diag(alpha.n)<- 1 #intraspecific competition
  alpha.c<-parms$alpha.c # parms$Cmatrix # competition matrix for consumers
  diag(alpha.c)<- 0.5 #intraspecific competition
  bi<-bc<-bp<-numeric()
  allee_b<- parms$allee_b # allee threshold basal, 0
  allee_c<- parms$allee_c #allee threshold consumers, 0.25  
  allee_p<-parms$allee_p
  h<-parms$h # handling time, 0.15
  a_r_c<- parms$attack_rate_c # 1.5
  a_r_p<- parms$attack_rate_p  # 1.5
  f_eff_c<- parms$feed_eff_c # 0.2
  f_eff_p<- parms$feed_eff_p # 0.25
  bi<- parms$gr_n # 1
  bc<-parms$gr_c # 0.05
  bp<-parms$gr_p # -0.001
  u_transition_p<-parms$mortality_n # 0
  a1 <- matrix(0, nrow=S_N, ncol = S_C)
  a2 <- matrix(0, nrow=S_C, ncol = S_N)
  a3 <- matrix(0, nrow=S_C,ncol = S_P)
  a4 <- matrix(0, nrow=S_P, ncol=S_C)
  summed_a1<-summed_a2<-summed_a3<-summed_a4<-numeric()

  time_func_N<-approxfun(parms$t1_N,method="linear", rule =2)
  time_func_C<-approxfun(parms$t1_C,method="linear", rule =2)
  time_func_P<-approxfun(parms$t1_P,method="linear", rule =2)
  
  #E_c<-as.numeric(noise_of_consumer(time,parms$noise_consumer,S_C))
  #E_p<-as.numeric(noise_of_top(time,parms$noise_top,S_P))
  noise_f_C<-approxfun(parms$noise_c$import[,1],method="linear",rule=2)
  noise_f_P<-approxfun(parms$noise_p$import[,1],method="linear", rule=2)
  
  e_N<-replicate(S_N,time_func_N(time))
  e_C<-replicate(S_C,time_func_C(time))
  e_P<-replicate(S_P,time_func_P(time))
  
  E_c<-replicate(S_C, noise_f_C(time))
  E_p<-replicate(S_P,noise_f_P(time))
  
  forcing_strength_N<-parms$forcing_strength[1:S_N]
  forcing_strength_C<-parms$forcing_strength[(S_N+1): (S_C+S_N)]
  forcing_strength_P<-parms$forcing_strength[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  forcing_index<-rep(0,(S_N+S_C+S_P))
  forcing_index[parms$species_index]<-1
  forcing_index_sp_N<-forcing_index[1:S_N]
  forcing_index_sp_C<-forcing_index[(S_N+1): (S_C+S_N)]
  forcing_index_sp_P<-forcing_index[(S_N+S_C+1):(S_C+S_N+S_P)]
  
  ##predation of resources by consumers
  for(i in 1:S_N){
    for(j in 1:S_C){
      a1[i,j]<- parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j]/(1+h*parms$attack_rate_c*A[prop$Bottom[i],prop$Intermediate[j]]*nc[j])
    }
    summed_a1[i]  <- sum(a1[i,])
  }
  
  for(r in 1:S_C){
    for(l in 1:S_N){
      a2[r,l]<- parms$attack_rate_c* parms$feed_eff_c*A[prop$Bottom[l],prop$Intermediate[r]]*na[l]/(1+h*parms$attack_rate_c*sum(A[prop$Bottom,prop$Intermediate[r]]*na))
    }
    summed_a2[r]<- sum(a2[r,])
  }
  for(m in 1:S_C){
    for(n in 1:S_P){
      a3[m,n]<-parms$attack_rate_p*A[prop$Intermediate[m],prop$Top[n]]*np[n]/(1+h*parms$attack_rate_p*sum(A[prop$Intermediate,prop$Top[n]]*nc))
      summed_a3[m] <-   sum(a3[m,]) 
    }
  }
  for(o in 1:S_P){
    for(p in 1:S_C){
      a4[o,p]<-parms$attack_rate_p*parms$feed_eff_p*A[prop$Intermediate[p],prop$Top[o]]*nc[p]/(1+h*parms$attack_rate_p*sum(A[prop$Intermediate,prop$Top[o]]*nc))
    }
    summed_a4[o]<- sum(a4[o,])
  }   
  #predation of consumers by top predators
  
  #equations for population dynamics
  
  dndt <- na*(bi - alpha.n%*%na  - summed_a1)*cutoff(na/(1e-5)) + rep(parms$forcing_delta,S_N)   # + perturb
  dcdt <- nc*(bc- alpha.c%*%(nc) + E_c+  summed_a2 - summed_a3 +  forcing_index_sp_C*forcing_strength_C*e_C)*cutoff(nc/(1e-5))
  dpdt <- np*(bp   -0.1*np + summed_a4 + E_p+  forcing_index_sp_P*forcing_strength_P*e_P)*cutoff(np/1e-5)
  
  
  return(list(c(dndt, dcdt,dpdt)))
}


sim_webs_function<-function(web,S_N,S_C,S_P, prop,delta){
  
  alpha.c<-matrix((runif( S_C*S_C, 0.1,0.1)), nrow=S_C, ncol= S_C)
  alpha.n<-matrix((runif( S_N*S_N, 0.5,0.5)), nrow=S_N, ncol= S_N)
  allee_b<-0
  allee_c<-0
  allee_p<-0
  h<-0.25 #
  feed_eff_c<-0.5
  feed_eff_p<-0.5
  attack_rate_c<-1.2
  attack_rate_p<-1.2
  gr_n<- 0.1 #growth rate respurce
  gr_c<- -0.1 # death rate consumer
  gr_p<- -0.05  #death rate pradtor
  mortality_n<-0.5
  tmax<-500
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 500 #fact$forcing_duration[r]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  
  noise_c_d<-replicate(S_C,rnorm(d,0,0.2))
  noise_p_d<-replicate(S_P,rnorm(d,0,0.2))
  
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
  
  
  species_index<-prop$Bottom #which((degree) == max(degree))
  forcing_delta<-delta
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
  na[species_index]<-forcing_delta
  nc<-rep(0.00001,S_C) # inital C
  np<-rep(0.00001,S_P) # initla P
  ic<-c(na,nc,np)
  sol3<-ode(func=fwebdyn_perturb_constant_noise, y=ic, parms=parms, times=seq(0, tmax, by=1)) %>%
    organize_results(parms)  # function that calculatesy dynamics and stores in a dataframe
  
  C_eq<-((sol3 %>% filter(time > 400, type =="C")) %>% group_by(species) %>% summarise(mean_density_C=mean(density)))$mean_density_C
  P_eq<-((sol3 %>% filter(time > 400, type =="P")) %>% group_by(species) %>% summarise(mean_density_P=mean(density)))$mean_density_P
  
  richness<-length(which(c(C_eq,P_eq)> 0.05))/(S_C+S_P)
  #mean_density_c<-mean(C_eq)
  #mean_density_p<-mean(P_eq)
  #inter_richness<-length(which(((sol3 %>% filter(time > 100, type =="C")) %>% group_by(species) %>% summarise(mean_density_C=mean(density)))$mean_density_C > 0.05))
  #top_richness<-length(which(((sol3 %>% filter(time > 100, type =="P")) %>% group_by(species) %>% summarise(mean_density_P=mean(density)))$mean_density_P > 0.05))
  
  mean_density <- mean((sol3 %>% filter(type=="C" | type =="P", time == (tmax-5)))$density)
  mean_density_c <- mean((sol3 %>% filter(type=="C", time == (tmax-5)))$density)
  mean_density_p <- mean((sol3 %>% filter(type =="P", time == (tmax-5)))$density)
  
  #mean_density <- mean(c(C_eq,P_eq))
  #mean_density_c <- mean_density_c #mean((sol3 %>% filter(type=="C", time == (tmax-5)))$density)
  #mean_density_p <- mean_density_p #mean((sol3 %>% filter(type =="P", time == (tmax-5)))$density)
  
  return(list(mean_density=mean_density, mean_density_c=mean_density_c,mean_density_p=mean_density_p,richness=richness,
              parms=parms))
}

nullclines_2<-function (deriv, xlim, ylim, parameters = NULL, system = "two.dim", 
                        points = 101, col = c("blue", "cyan"), add = TRUE, add.legend = TRUE, 
                        state.names = if (system == "two.dim") c("x", "y") else "y", 
                        ...) 
{
  if (any(!is.vector(xlim), length(xlim) != 2)) {
    stop("xlim is not a vector of length 2, as is required")
  }
  if (xlim[2] <= xlim[1]) {
    stop("xlim[2] is less than or equal to xlim[1]")
  }
  if (any(!is.vector(ylim), length(ylim) != 2)) {
    stop("ylim is not a vector of length 2, as is required")
  }
  if (ylim[2] <= ylim[1]) {
    stop("ylim[2] is less than or equal to ylim[1]")
  }
  if (points <= 0) {
    stop("points is less than or equal to zero")
  }
  if (!(system %in% c("one.dim", "two.dim"))) {
    stop("system must be set to either \"one.dim\" or \"two.dim\"")
  }
  if (!is.vector(col)) {
    stop("col is not a vector as required")
  }
  if (length(col) != 2) {
    if (length(col) == 1) {
      col <- rep(col, 2)
    }
    else if (length(col) > 2) {
      col <- col[1:2]
    }
    message("Note: col has been reset as required")
  }
  if (!is.logical(add)) {
    stop("add must be logical")
  }
  if (!is.logical(add.legend)) {
    stop("add.legend must be logical")
  }
  x <- seq(xlim[1], xlim[2], length.out = points)
  y <- seq(ylim[1], ylim[2], length.out = points)
  dx <- dy <- matrix(0, ncol = points, nrow = points)
  if (system == "one.dim") {
    for (i in 1:points) {
      dy[1, i] <- deriv(0, stats::setNames(c(y[i]), state.names), 
                        parameters)[[1]]
    }
    for (i in 2:points) {
      dy[i, ] <- dy[1, ]
    }
    graphics::contour(x, y, dy, levels = 0, add = add, col = col[1], 
                      drawlabels = F, ...)
    if (add.legend) {
      graphics::legend("bottomright", paste0("d", state.names, 
                                             "/dt = 0 for all t"), lty = 1, lwd = 1, col = col[1])
    }
    return(list(add = add, add.legend = add.legend, col = col, 
                deriv = deriv, dy = dy, parameters = parameters, 
                points = points, system = system, x = x, xlim = xlim, 
                y = y, ylim = ylim))
  }
  else {
    for (i in 1:points) {
      for (j in 1:points) {
        df <- deriv(0, stats::setNames(c(x[i], y[j]), 
                                       state.names), parameters)
        dx[i, j] <- df[[1]][1]
        dy[i, j] <- df[[1]][2]
      }
    }
    graphics::contour(x, y, dx, levels = 0, add = add, col = col[1], 
                      drawlabels = F, ...)
    graphics::contour(x, y, dy, levels = 0, add = T, col = col[2], 
                      drawlabels = F, ...)
    if (add.legend) {
      legend("top", paste(c("",""), ""), bty="n",inset = c(-0.45, 0),
             lty = 1, lwd = 2, col = col, cex  =1.1)
    }
    return(list(add = add, add.legend = add.legend, col = col, 
                deriv = deriv, dx = dx, dy = dy, parameters = parameters, 
                points = points, system = system, x = x, xlim = xlim, 
                y = y, ylim = ylim))
  }
}



#%>% 
 # organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs

## Organize simulation results into tidy table (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - sol: output produced by the function ode()
## - pars: list of parameters, with the following elements:
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma
organize_results <- function(sol, pars) {
  S <- length(pars$S) ## number of species
  N<-pars$S_N # no. of animals
  C<-pars$S_C
  P<-pars$S_P # no. of plants
  temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  ## name the first column "time"
  # temp<- temp %>% filter(time >= pars$cutoff.time)
  names(temp)[2:(N+1)] <- paste0("N_", 1:(N)) ## name abundance columns (n_k)
  names(temp)[1] <- "time"
  names(temp)[(N+2):(N+C+1)] <- paste0("C_",  N+1:C) ## name trait mean columns
  names(temp)[(N+C+2):(N+C+P+1)] <- paste0("P_",N+C+1:P)
  temp <- temp %>%
    gather("variable", "density", 2:ncol(temp)) %>% ## normalize the data
    separate(variable, c("type", "species"), sep="_") %>%
    # spread(type, v) %>% ## separate columns for animal densities n and plant densities m
    dplyr::select(time, type, species, density) %>% ## rearrange columns
    mutate(species=as.integer(species)) ## add params
  return(as_tibble(temp))
}


## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
  plot_grid(plot_density(dat), ncol=1, align="hv") %>%
    return
}



## code adapted from Barabas and D'Andrea 2016 Eco.Letts paper
plot_density<- function(dat) {
  dat %>%
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    theme(legend.position="none") + facet_wrap(.~type) %>%
    return
}


#lay out function for multiple plots
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 



#multiplot of ggplot2 figures with a common shared legend. Code taken from :https://rpubs.com/sjackman/grid_arrange_shared_legend
grid_arrange_shared_legend <- function(..., ncol, nrow, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol =ncol , nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


# plots the density distribution  at a particular timepoint. This function was used to produce figure 1.
# Na: Abundance of animals at equilibrium
# Np: Abundance of plants at equilibrium
# m: mean traits at equilibrium
# sigma: variance of traits
# moment: mean
# limits: limits of the mean trait axis which in the study are -1,1

plot_snapshot <- function(Na, Np, m, sigma, moment=0, limits=c(-1, 1), res=1001) {
  Sa <- length(Na) ## number of species
  Sp <- length(Np)
  ma<- m[1:(Sa)]
  mp<- m[Sa+1:Sp]
  sigma_a <-sigma[1:(Sa)]
  sigma_p <- sigma[Sa+1:Sp]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:Sa, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=Sa+1:Sp, trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  
  for (i in 1:Sa) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:Sp) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==(Sa+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==(Sa+i))], 
                                                                mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+scale_fill_viridis_d()+
    facet_wrap(.~species_group, nrow = 2)+
    theme(legend.title = element_text(size = 14, face = "bold"), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14, face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14, face ="bold"))+
    #geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
    #         colour="darkred", alpha=0.5, na.rm=TRUE) +
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    theme(legend.position="none") %>%
    return }




#computes the raw NODF taken from Song et al 2017 J. Animal Ecology
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}



# measures connectance of a web network
# web: interaction network
Connectance<-function(web){
  return(sum(web)/(ncol(web)*nrow(web)))}



# function for sampling competitive coefficients from random uniform distribution 
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength

mat.comp<-function(matrix){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-matrix(runif(Aspecies^2, 0.0001, 0.001), nrow=Aspecies, ncol = Aspecies)/Aspecies #scaled by number of competitors within a guild
  diag(Amatrix)<-1 #intraspecific competition for animals
  Pmatrix<-matrix(runif(Plantspecies^2, 0.0001, 0.001), nrow=Plantspecies, ncol = Plantspecies)/Plantspecies ##scaled by number of competitors within a guild
  diag(Pmatrix)<-1 #intraspecific competion for plants
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}


cubic_solution<-function(params,web,delta){
  Degrees <- colSums(web$web[web$index$Bottom,web$index$Intermediate]) #degree distribution
   
  web$web[web$index$Intermediate, web$index$Top] #web of interaction of top and consumer
  
  if( length(web$index$Top) > 1){
    degree_top<- colSums(web$web[web$index$Intermediate, web$index$Top])
    degree_int <- rowSums(web$web[web$index$Intermediate, web$index$Top])#degree of intermediate
    w_cn<-sum(colSums(web$web[web$index$Bottom,web$index$Intermediate])*Degrees)/sum(Degrees)*parameters$at_c*parameters$eff_c #reduced effective interaction strength of consumers/intermeiate
    w_pc <- sum( colSums(web$web[web$index$Intermediate, web$index$Top])*degree_top)/sum(degree_top)*parameters$eff_p*parameters$at_p #reduced effective strength predataor
    w_cp <- sum( rowSums(web$web[web$index$Intermediate, web$index$Top])*degree_int)/(sum(degree_int))*parameters$eff_p*parameters$at_p
    
  }else{
    degree_top<-(web$web[web$index$Intermediate, web$index$Top])
    degree_int <- (web$web[web$index$Intermediate, web$index$Top])#degree of intermediate
    w_cn<-sum(colSums(web$web[web$index$Bottom,web$index$Intermediate])*Degrees)/sum(Degrees)*parameters$at_c*parameters$eff_c #reduced effective interaction strength of consumers/intermeiate
    w_pc <- sum( (web$web[web$index$Intermediate, web$index$Top])*degree_top)/sum(degree_top)*parameters$eff_p*parameters$at_p #reduced effective strength predataor
    w_cp <- sum( (web$web[web$index$Intermediate, web$index$Top])*degree_int)/(sum(degree_int))*parameters$eff_p*parameters$at_p
  } 

  B<-params$b2 - params$f
  a<- -params$h^2*params$alpha_N*params$alpha_c
  b <- params$F_c*params$h^2*params$alpha_c+ params$h^2*params$b1*params$alpha_c -2*params$h*params$alpha_N*params$alpha_c
  c <- 2*params$h*params$b1*params$alpha_c - w_cn*w_nc*params$eff_c^2 - w_cn*params$eff_c*params$h*B
  
  d <-2*params$h*params$F_c*params$alpha_c -params$alpha_c*params$alpha_N
  
  d0<- round( (b^2 - 3*a*c),2)
  d1<- round( (2*b^3 - 9*a*b*c + 27*a^2*d),2)
  
  if(d0 ==0 | d1 ==0 ){
    ep<- (-1+ sqrt(as.complex(-3)))/2
    C<- ((d1 - sqrt( as.complex(d1^2 - 4*d0^3)))/2)^1/3
    x1<- -1/(3*a)*(b + ep^0*C )
    x2<- -1/(3*a)*(b + ep^1*C )
    x3<- -1/(3*a)*(b + ep^2*C )
    
  }else{
  
  ep<- (-1+ sqrt(as.complex(-3)))/2
  C<- ((d1 -sqrt( as.complex(d1^2 - 4*d0^3)))/2)^1/3
  x1<- -1/(3*a)*(b + ep^0*C + d0/(ep^0*C))
  x2<- -1/(3*a)*(b + ep^1*C + d0/(ep^1*C))
  x3<- -1/(3*a)*(b + ep^2*C + d0/(ep^2*C))
  }
  return(list(c(x1=x1,x2=x2,x3=x3)))
}

