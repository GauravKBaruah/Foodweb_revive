
Niche.model <- function(S,C,N=1){
  L<- C*S^2 
  if(N==1){
    n <- sort(runif(S))
    beta <- (1 - 2 * C) / (2 * C)
    r <- n*(1 - (1 - runif(S))^(1/beta))
    c <- r/2 + runif(S) * (n - r/2)
    web <- matrix(0,S,S)
    min.n <- c-r/2
    max.n <- c+r/2
    for(i in 1:S){
      diet <- c(1:S)[c(which(n>min.n[i]), which(n<max.n[i]))[duplicated(c(which(n>min.n[i]), which(n<max.n[i])))]]
      web[diet,i] <- 1
    }
    dimnames(web) <- list(1:length(web[,1]), 1:length(web[,1]))
  }
  if(N>1){
    web <- list()
    for(j in 1:N){
      n <- sort(runif(S))
      beta <- (1 - 2 * C) / (2 * C)
      r <- n*(1 - (1 - runif(S))^(1/beta))
      c <- r/2 + runif(S) * (n - r/2)
      web[[j]] <- matrix(0,S,S)
      min.n <- c-r/2
      max.n <- c+r/2
      for(i in 1:S){
        diet <- c(1:S)[c(which(n>min.n[i]), which(n<max.n[i]))[duplicated(c(which(n>min.n[i]), which(n<max.n[i])))]]
        web[[j]][diet,i] <- 1
      }
      dimnames(web[[j]]) <- list(1:length(web[[j]][,1]), 1:length(web[[j]][,1]))
    }
  }
  web
}


## Make N cascade food webs with S species and L links
## The output is a matrix if N=1 or a list of matrices if N>1
## The cascade model is that of Cohen et al
Cascade.model <- function(S, L, N=1){
  if(N==1){
    web <- matrix(0, S, S)
    web[upper.tri(web)] <- c(rep(1, L), rep(0, (S^2-S)/2 - L))[order(runif((S^2-S)/2))]
    dimnames(web) <- list(1:length(web[,1]), 1:length(web[,1]))
  }
  if(N>1){
    web <- list()
    for(i in 1:N){
      web[[i]] <- matrix(0, S, S)
      web[[i]][upper.tri(web[[i]])] <- c(rep(1, L), rep(0, (S^2-S)/2-L))[order(runif((S^2-S)/2))]
      dimnames(web[[i]]) <- list(1:length(web[[i]][,1]), 1:length(web[[i]][,1]))
    }
  }
  web
}





## return the proportion of bottom, intermediate, top,
## unconnected and purely cannibalistic specie
Bottom.Intermediate.Top <- function(web, proportion=TRUE, check.web=FALSE){
  
  ## find the names and numbers of BIT, unconnected and pure cannibals
  names.all <- 1:length(web[1,])
  dimnames(web) <- list(names.all, names.all)
  names.Bottom <-  names.all[apply(web, 2, sum)==0 & apply(web, 1, sum)!=0]
  number.Bottom <- length(names.Bottom)
  names.Top <-   names.all[apply(web, 2, sum)!=0 & apply(web, 1, sum)==0]
  number.Top <- length(names.Top)
  names.Unconnected <- names.all[apply(web, 2, sum)==0 & apply(web, 1, sum)==0]
  number.Unconnected <- length(names.Unconnected)
  number.Pure.cannibals <- 0
  names.Pure.cannibals <- character(0)
  for(i in 1:length(web[,1]))
    if(web[i,i]!=0 && sum(web[-i,i]+web[i,-i])==0){
      if(number.Pure.cannibals==0){
        names.Pure.cannibals <- dimnames(web)[[1]][i]
        number.Pure.cannibals <- 1
      }
      else{
        names.Pure.cannibals <- c(names.Pure.cannibals, dimnames(web)[[1]][i])
        number.Pure.cannibals <- number.Pure.cannibals + 1
      }
    }
  names.Intermediate <- dimnames(web)[[1]][is.na(match(dimnames(web)[[1]],
                                                       c(names.Bottom, names.Top, names.Unconnected, names.Pure.cannibals)))]
  number.Intermediate <- length(web[1,])-number.Bottom-number.Top-number.Unconnected-number.Pure.cannibals
  
  if(proportion)
    result <- list(Bottom=names.Bottom,
                   Intermediate=names.Intermediate,
                   Top=names.Top,
                   Unconnected=names.Unconnected,
                   Pure.cannibals=names.Pure.cannibals,
                   Proportions.of.each=c(number.Bottom, number.Intermediate, number.Top, number.Unconnected, number.Pure.cannibals)/length(web[1,]))
  if(!proportion)
    result <- list(Bottom=names.Bottom,
                   Intermediate=names.Intermediate,
                   Top=names.Top,
                   Unconnected=names.Unconnected,
                   Pure.cannibals=names.Pure.cannibals,
                   Proportions.of.each=c(number.Bottom, number.Intermediate, number.Top, number.Unconnected, number.Pure.cannibals))
  result
}




## Make N cascade food webs with S species and L links
## The output is a matrix if N=1 or a list of matrices if N>1
## The cascade model is that of Cohen et al
Cascade.model <- function(S, L, N=1){
  if(N==1){
    web <- matrix(0, S, S)
    web[upper.tri(web)] <- c(rep(1, L), rep(0, (S^2-S)/2-L))[order(runif((S^2-S)/2))]
    dimnames(web) <- list(1:length(web[,1]), 1:length(web[,1]))
  }
  if(N>1){
    web <- list()
    for(i in 1:N){
      web[[i]] <- matrix(0, S, S)
      web[[i]][upper.tri(web[[i]])] <- c(rep(1, L), rep(0, (S^2-S)/2-L))[order(runif((S^2-S)/2))]
      dimnames(web[[i]]) <- list(1:length(web[[i]][,1]), 1:length(web[[i]][,1]))
    }
  }
  web
}


pyramidal.food<-function(S,C, prop, balance){
  index<-list()
  L <-S^2*C
  top.predators<- round(prop[3]*S,0)
  intermediate<-round(prop[2]*S,0)
  basal<-round(prop[1]*S,0)
  #first trophic level
  Links_1<-round( L -balance*top.predators*intermediate,0)
  #2nd trophic level
  Links_2<-round((L-Links_1),0)
  bassal.matr<-matrix(0,nrow=basal,ncol=basal)
  top.pred.matr<-matrix(0, nrow=top.predators,ncol=top.predators)
  consumer.matr<-matrix(0,nrow=intermediate,ncol=intermediate)
  basal.consumer.matr<-matrix(c(rep(1, Links_1/2), rep(0,abs(basal*intermediate-Links_1/2)))[order(runif(basal*(intermediate)))] ,
                              nrow=basal, ncol=(intermediate))
  consumer.to.pred.matrix<-matrix(c(rep(1,Links_2),
                                    rep(0, abs(top.predators*intermediate-Links_2/2)))[order(runif(intermediate*top.predators))],
                                  nrow=intermediate, ncol = top.predators)
  basal.to.top.pred.matr<-matrix(rep(0, basal*top.predators), nrow=basal,ncol=top.predators)
  transpose.basal.top.pred.matr<-t(basal.to.top.pred.matr)
  transpose.basal.consumer.matr<-t(basal.consumer.matr)
  transpose.consumer.top.pred.matr<-t(consumer.to.pred.matrix)
  web<-rbind(cbind(bassal.matr,basal.consumer.matr,basal.to.top.pred.matr),
             cbind(transpose.basal.consumer.matr,consumer.matr,consumer.to.pred.matrix),
             cbind(transpose.basal.top.pred.matr,transpose.consumer.top.pred.matr,top.pred.matr))
  web
  index<-NULL
  index$Bottom <- seq(1:basal)
  index$Intermediate<-seq((basal+1),(basal+intermediate),1)
  index$Top<-seq( (basal+intermediate+1),(basal+intermediate+top.predators),1)
  web[is.na(web)]<-0
  tp1<-web[index$Bottom,index$Intermediate]
  for (i in 1:ncol(web[index$Bottom,index$Intermediate])){
    if (all(web[index$Bottom,index$Intermediate][,i]==0) ==TRUE){
      xval<-round(runif(1,1,length(web[index$Bottom,index$Intermediate][,i])),0)
      web[index$Bottom,index$Intermediate][xval,i]<-1
      web[index$Intermediate,index$Bottom][i,xval]<-1}else {
        web<-web
      }
  }
  for (i in 1:nrow(web[index$Bottom,index$Intermediate])){
    if (all(web[index$Bottom,index$Intermediate][i,]==0) ==TRUE){
      xval<-round(runif(1,1,length(web[index$Bottom,index$Intermediate][i,])),0)
      web[index$Bottom,index$Intermediate][i,xval]<-1
      web[index$Intermediate,index$Bottom][xval,i]<-1}else {
        web<-web
      }
  }
  tp<-web[index$Top,index$Intermediate]
  if(length(dim(tp)) ==0 ){
    for (i in 1:length(web[index$Top,index$Intermediate])){
      if (all(web[index$Top,index$Intermediate]==0) == FALSE){
        web<-web }else{
          xval<-round(runif(1,1,length(web[index$Top,index$Intermediate][i,])) )
          web[index$Top,index$Intermediate][xval]<-1
        }
    }
  }else{
    for (i in 1:nrow(web[index$Top,index$Intermediate])){
      if (all(web[index$Top,index$Intermediate][i,]==0)==TRUE){
        xval<-round(runif(1,1,length(web[index$Top,index$Intermediate][i,])),0)
        web[index$Top,index$Intermediate][i,xval]<-1
        web[index$Intermediate,index$Top][xval,i]<-1
      }else {
        web<-web
      }
    }
    for (i in 1:ncol(web[index$Top,index$Intermediate])){
      if (all(web[index$Top,index$Intermediate][,i]==0)==TRUE){
        xval<-round(runif(1,1,length(web[index$Top,index$Intermediate][,i])),0)
        web[index$Top,index$Intermediate][xval,i]<-1
        web[index$Intermediate,index$Top][i,xval]<-1
      }else {
        web<-web
      }
    }
  }
  return(list(web=web,index=index))
}




restructuring_niche_web<-function(S,C,N=1){
 s = 0
  while(s < 1){
   web<-Niche.model(S = S,C=C, N = 1)
   index<-Bottom.Intermediate.Top(web = web)
   if( length(index$Unconnected)==0 & length(index$Pure.cannibals) ==0 & length(index$Top)>=1 ){
     s = 1
     web<- list(web=web, index=index)
   }else{ s = 0}
   }
   
 for (i in 1:ncol(web$web[web$index$Bottom,web$index$Intermediate])){
    
    if (all(web$web[web$index$Bottom,web$index$Intermediate][,i]==0) ==TRUE){
      xval<-round(runif(1,1,length(web$web[web$index$Bottom,web$index$Intermediate][,i])),0)
      web$web[web$index$Bottom,web$index$Intermediate][xval,i]<-1
      web$web[web$index$Intermediate,web$index$Bottom][i,xval]<-1}else {
        web<-web
      }
  }
  tp<-web$web[web$index$Top,web$index$Intermediate]
  if(length(dim(tp)) ==0 ){
    for (i in 1:length(web$web[web$index$Top,web$index$Intermediate])){
      if (all(web$web[web$index$Top,web$index$Intermediate]==0) == FALSE){
        web<-web }else{
          xval<-round(runif(1,1,length(web[web$index$Top,web$index$Intermediate][i,])) )
          web$web[web$index$Top,web$index$Intermediate][xval]<-1
        }
    }
  }else{
    for (i in 1:ncol((web$web[web$index$Intermediate,web$index$Top]))){
      if (all(web$web[web$index$Intermediate,web$index$Top][,i]==0)==TRUE){
        xval<-round(runif(1,1,length(web$web[web$index$Intermediate,web$index$Top][,i])),0)
        
        web$web[web$index$Intermediate,web$index$Top][xval,i]<-1
      }else {
        web<-web
      }
    }
    
    
  }
  return(web)
}
