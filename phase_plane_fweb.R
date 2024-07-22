

phase_plot<-function(parameters,web,delta, method){
  
  Degrees <- colSums(web$web[web$index$Bottom,web$index$Intermediate]) #degree distribution
  web$web[web$index$Intermediate, web$index$Top] #web of interaction of top and consumer
  if( length(web$index$Top) > 1){
    degree_top<- colSums(web$web[web$index$Intermediate, web$index$Top])
    degree_int <- colSums(web$web[web$index$Top, web$index$Intermediate])#degree of intermediate
    
  }else{
    degree_top<-(web$web[web$index$Intermediate, web$index$Top])
    degree_int <- (web$web[web$index$Top, web$index$Intermediate])#degree of intermediate
    
    }
 #degree_int <- colSums(web$web[web$index$Top, web$index$Intermediate])#degree of intermediate
  if(method == "unweighted"){
  w_cn<-sum(colSums(web$web[web$index$Bottom,web$index$Intermediate]))/sum(1)*parameters$at_c*parameters$eff_c #reduced effective interaction strength of consumers/intermeiate
  w_pc <- sum(web$web[web$index$Intermediate, web$index$Top]*1)/sum(1)*parameters$eff_p*parameters$at_p #reduced effective strength predataor
  w_cp <- sum(web$web[web$index$Top, web$index$Intermediate]*1)/(sum(1))*parameters$eff_p*parameters$at_p
  }else{
  w_cn<-sum(colSums(web$web[web$index$Bottom,web$index$Intermediate])*Degrees)/sum(Degrees)*parameters$at_c*parameters$eff_c #reduced effective interaction strength of consumers/intermeiate
  w_pc <- sum(web$web[web$index$Intermediate, web$index$Top]*degree_top)/sum(degree_top)*parameters$eff_p*parameters$at_p #reduced effective strength predataor
  w_cp <- sum(web$web[web$index$Top, web$index$Intermediate]*degree_int)/(sum(degree_int))*parameters$eff_p*parameters$at_p
  }
    #initial densities of the reduced model
  
  na<-seq(-5,20,0.1)
  np<-seq(-5,20,0.1)
  
  
  dna <- (parameters$bi) -parameters$alpha_c*na + w_cn*(parameters$delta/(1+parameters$h*w_cn*parameters$delta)) - 
    w_cp*np/(1+ parameters$h*w_cp*na)
  dnp <- parameters$bp -0.05*np  + w_pc*na/(1+parameters$h*w_pc*na)
  
  return(list(dna=dna, dnp=dnp,na=na,np=np))
}



dynamics_reduced_model<-function(t,y,parameters){
  
  na<-y[1]
  np<-y[2]
  
  web<-parameters$web
  
  Degrees <- colSums(web$web[web$index$Bottom,web$index$Intermediate]) #degree of intermediates
  web$web[web$index$Intermediate, web$index$Top] #web of interaction of top and consumer
  if( length(web$index$Top) > 1){
    degree_top<- colSums(web$web[web$index$Intermediate, web$index$Top])
    degree_int <- rowSums(web$web[web$index$Intermediate, web$index$Top])#degree of intermediate
    w_cn<-sum(colSums(web$web[web$index$Bottom,web$index$Intermediate])*Degrees)/sum(Degrees)*parameters$at_c #reduced effective interaction strength of consumers/intermeiate
    w_pc <- sum( colSums(web$web[web$index$Intermediate, web$index$Top])*degree_top)/sum(degree_top)*parameters$at_p #reduced effective strength predataor
    w_cp <- sum( rowSums(web$web[web$index$Intermediate, web$index$Top])*degree_int)/(sum(degree_int))*parameters$at_p
    w_cp_p<- sum( rowSums(web$web[web$index$Intermediate, web$index$Top])*degree_top)/(sum(degree_top))*parameters$at_p
  }else{
    degree_top<-(web$web[web$index$Intermediate, web$index$Top])
    degree_int <- (web$web[web$index$Intermediate, web$index$Top])#degree of intermediate
    w_cn<-sum(colSums(web$web[web$index$Bottom,web$index$Intermediate])*Degrees)/sum(Degrees)*parameters$at_c #reduced effective interaction strength of consumers/intermeiate
    w_pc <- sum( (web$web[web$index$Intermediate, web$index$Top])*degree_top)/sum(degree_top)*parameters$at_p #reduced effective strength predataor
    w_cp <- sum( (web$web[web$index$Intermediate, web$index$Top])*degree_int)/(sum(degree_int))*parameters$at_p
    w_cp_p<-sum( (web$web[web$index$Intermediate, web$index$Top])*degree_top)/(sum(degree_int))*parameters$at_p
  }
  
  dna <- (parameters$bi)*na -parameters$alpha_c*na^2 + w_cn*parameters$eff_c*(parameters$delta/(1+parameters$h*w_cn*parameters$delta))*na - 
    w_cp*np/(1+ parameters$h*na*w_cp_p)*na
  dnp <- parameters$bp*np -0.1*np^2  + w_pc*parameters$eff_p*na/(1+parameters$h*w_pc*na)*np
  
  
    return(list(c(dna, dnp)))
}




# parameters$web <- web
# ic<-c(0.00001,0.00001)
# parameters$delta<-1
# 
# 
# soln<-ode(func=dynamics_reduced_model, y=ic, parms=parameters, times=seq(0, 1000, by=1)) 
# soln[1000,2]
# soln[1000,3]
# 
#  dd<-expand.grid(delta = c(0, 0.01,0.1,0.5,5,10,20,100))
# # 
#  for(i in 1: nrow(dd)){
#    
#    ic<-c(1,1)
#    parameters$delta <- dd$delta[i]
#    soln<-ode(func=dynamics_reduced_model, y=ic, parms=parameters, times=seq(0, 1000, by=1)) 
#    c<-soln[970,2]
#    p<-soln[970,3]
#    dd$c[i]<-c
#    dd$p[i] <- p
#    
#  }
#  dd
# # 
# phase<-function(t, y, parameters){
#   
#   web <- parameters$web
#   
#   Degrees <- colSums(web$web[web$index$Bottom,web$index$Intermediate]) #degree distribution
#   if( length(web$index$Top) > 1){
#     degree_top<- colSums(web$web[web$index$Intermediate, web$index$Top])
#     degree_int <- colSums(web$web[web$index$Top, web$index$Intermediate])#degree of intermediate
#     
#   }else{
#     degree_top<-(web$web[web$index$Intermediate, web$index$Top])
#     degree_int <- (web$web[web$index$Top, web$index$Intermediate])#degree of intermediate
#     
#   }    
#   
#     w_cn<-sum(colSums(web$web[web$index$Bottom,web$index$Intermediate])*Degrees)/sum(Degrees)*parameters$at_c*parameters$eff_c #reduced effective interaction strength of consumers/intermeiate
#     w_pc <- sum(web$web[web$index$Intermediate, web$index$Top]*degree_top)/sum(degree_top)*parameters$eff_p*parameters$at_p #reduced effective strength predataor
#     w_cp <- sum(web$web[web$index$Top, web$index$Intermediate]*degree_int)/(sum(degree_int))*parameters$eff_p*parameters$at_p
#   
#   #initial densities of the reduced modelhttp://127.0.0.1:33879/graphics/plot_zoom_png?width=1693&height=861
#   
#   na<-seq(0, 500,0.1)
#   np<-seq(0, 500,0.1)
#   
#   #
#   #dnp <- (1+parameters$h*parameters$at_p*na)/(w_cp*na)*((parameters$bi)*na - parameters$alpha_c*na^2 +  w_cn*(parameters$delta/(1+parameters$h*w_cn*parameters$delta)) )
#   #dna <-  (0.001*np - parameters$bp)/(w_pc*np -w_pc*0.001*np*parameters$h+parameters$bp*w_pc*parameters$h)
#  
#   
  #denominator_n<-(parameters$bi*parameters$h +parameters$bi*parameters$h^2*w_cn*parameters$delta -parameters$h*w_cn*parameters$delta)
  #c<- (w_cp*np +w_cn*parameters$delta*w_cp*np - parameters$bi - parameters$bi*w_cn*parameters$delta*parameters$h + w_cn*parameters$delta)/denominator_n

#  alpha_c<-parameters$alpha_c
#  p<- (parameters$bp/1 + 1/1*(w_pc*na/(1+w_pc*parameters$h*na)))
#
#  h<-parameters$h
#  delta<-parameters$delta
#  bi<-parameters$bi
#  a1 = -alpha_c*h*w_cp -alpha_c*w_cn*h^2*delta*w_cp
#  b1 = bi*h*w_cp +w_cn*bi*delta*h^2*w_cp - alpha_c- alpha_c*w_cn*h*delta + h*w_cn*delta*w_cp
#  c1 = w_cn*delta + bi + bi*w_cn*h*delta - w_cp*np - w_cn*h*delta*w_cp*np
#
# x1_c = (- b1 - sqrt(b1^2 - 4*a1*c1))/(2*a1)
# x2_c = (-b1 + sqrt(b1^2 - 4*a1*c1))/(2*a1)
# #   
# # plot(na,p)
# # lines(x2_c,na)
# # lines(x1_c,na)
# # 
# # plot(na,x2_c) 
# # 
# # lines(np,x2_c)
# 
#   return(list(dna=dna, dnp=dnp,na=na,np=np))
# }
# 




# 
# parameters$delta <-1
# parameters$method <- "degree"
# 
# 
# 
# flowField(dynamics_reduced_model, xlim = c(-1, 5), 
#           ylim = c(-1, 5),
#           main=(expression(Delta == 20)),
#           parameters = parameters, 
#           ylab = expression(P[eff]),
#           xlab= expression(C[eff]),
#           points = 20, add = FALSE)  
# 
# nullclines_2(dynamics_reduced_model,xlim = c(-1,5), 
#            ylim = c(-1, 5),
#            parameters = parameters, 
#            points = 500, lwd = 3,
#            col=c("#E69F00","#009E73"),
#            add.legend = TRUE) 
# d<-trajectory(dynamics_reduced_model, y0 = c(0.0001,0.0001), tlim =c(0,1000), 
#            parameters = parameters )
# 
# 
# findEquilibrium(dynamics_reduced_model,max.iter = 50, y0 = NULL,parameters = parameters)
