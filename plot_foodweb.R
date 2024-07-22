#plotting the foodweb in igraph


plotweb<-function(web,TP, nspecies){
  lay<-matrix(nrow=nspecies,ncol=2)
  
  basal <- length(which(TP==0))
  intermediate<- length(which(TP==1))
  top<- length(which(TP==2))
  
  lay[,1][1:basal]<-seq(0,0.3,(0.3-0)/(basal-1)) 
  lay[,1][(basal+1):(basal+intermediate)]<-seq(0.0,0.3,(0.3-0)/(intermediate-1)) 
  lay[,1][(basal+intermediate+1):(basal+intermediate+top)]<-seq(0.05,0.3, (0.3-0)/top)
  
  lay[,2]<-TP
  web<-unname(web)
  net<-network(web, directed=F)
  
  net %v% "groups" = c(rep("#0072B2", each=basal), rep("#E69F00" , each= intermediate), rep("#009E73", each=top))
  g<- ggnet2(net, mode=lay, label=1: (basal+intermediate+top), label.size = 3, color = "groups",  edge.size = 1.1,
             node.size = 5,max_size =3,edge.alpha = 1) 
  return(g)
}

#plotting the foodweb in igraph


plotweb_PMN<-function(web,TP,prop, nspecies){
  lay<-matrix(nrow=nspecies,ncol=2)
  
  basal <- length(which(TP==0))
  intermediate<- length(which(TP==1))
  top<- length(which(TP==2))
  
  lay[,1][1:basal]<-seq(0,0.3,(0.3-0)/(basal-1)) 
  lay[,1][(basal+1):(basal+intermediate)]<-seq(0.0,0.3,(0.3-0)/(intermediate-1)) 
  lay[,1][(basal+intermediate+1):(basal+intermediate+top)]<-seq(0.05,0.3, (0.3-0)/top)
  
  lay[prop$Bottom,2]<-0
  lay[as.numeric(prop$Intermediate),2]<-1
  lay[prop$Top,2]<-2
  
  web<-unname(web)
  net<-network(web, directed=F)
  
  
 # rep("#0072B2", each=basal), rep("#E69F00" , each= intermediate), 
#  rep("#009E73", each=top)
  colors<-c("#0072B2","#0072B2",
            "#E69F00","#0072B2",
            "#E69F00","#E69F00",
            "#E69F00","#0072B2",
            "#E69F00","#0072B2",
            "#0072B2","#E69F00",
            "#E69F00","#E69F00",
            "#E69F00","#0072B2",
            "#E69F00","#E69F00",
            "#009E73","#009E73",
            "#009E73","#009E73",
            "#009E73","#009E73")
  

  
  net %v% "groups" = colors
  g<- ggnet2(net, mode=lay, label=1: (basal+intermediate+top), color ="groups", label.size = 3,  edge.size = 1.1,
             node.size = 5,max_size =3,edge.alpha = 1) 
  return(g)
}

