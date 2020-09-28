args = commandArgs(trailingOnly=TRUE)
library(deSolve) # for ODE simulation
# library(matlib)
setwd("Data/Raw/Serial_extinction_sub_model")

## Functions ######
compModel <- function(Time, B, Pars) {
  with(as.list(Pars), {
    dB <- as.vector((r + C %*% B)*B)
    return(list(dB))
  })
}

## Gives number of species in community
richness <- function(){
  length(pars$r)
}

## Add a species to a community
invade <- function(){
  S <- richness()
  i <- S+1
  pars$r[i] <- growth
  pars$C <- pars$C[c(1:S,S),,drop=F][,c(1:S,S),drop=F]
  pars$C[i,] <- 0.2*rbinom(n = i, size = 1,prob = 0.2)
  pars$C[,i] <- 0.2*rbinom(n = i, size = 1,prob = 0.2)
  pars$C[i,i] <- 1
  pars$i[i] <- inv
  BRelaxed[i] <<- inv
  return(pars)
}

## Simulate community dynamics
relax <- function(pars,BB,tRelax = 200){
  times <- c(0,tRelax)
  sol <-
    ode(y = c(B=BB),times = times,maxsteps = 5000,
        func = compModel,parms = pars )
  return(initfunc = sol[2,-1])
  # return(BB<--inv(pars$C)%*%pars$r)
}

## Remove extinct species
cleanup <- function(){
  extinct <- BRelaxed < 100*inv
  if(any(extinct)){
    pars$r <- pars$r[!extinct]
    pars$i <- pars$i[!extinct]
    pars$C <- pars$C[!extinct,!extinct,drop=FALSE]
  }
  BRelaxed <<- BRelaxed[!extinct]
  return(pars)
}

## Initialise community
inv <- 1e-6 # invasion pressure for all species

par(mfrow=c(2,2))

ser_ext_temp_num<-function(pars,BRelaxed){
  dead<-which(BRelaxed<(inv*10))
  if(any(dead)){
    pars$C<-pars$C[-dead,-dead,drop=F]
    pars$r<-pars$r[-dead]
    pars$i<-pars$i[-dead]
  }
  return(pars)
}

serial_extinction_t<-function(Hi){
  # Simulates the serial extinction process at equilibrium
  repeat{
    sum_Hi <- sum(Hi)
    sum_H2i <- sum(Hi^2)
    max_Hi <- max(Hi)
    out<-Hi
    if(sum_H2i >=  max_Hi*(sum_Hi-1)) break
    deadd<-which.max(Hi)
    Hi <- Hi[-deadd]
  }
  return(out)
}

check_extinction<-function(HH){
  A<-resp*HH/(K*epsilon)
  C_P<-diag(rep(growth/K,dim(HH)[1]))
  C_hat<-epsilon*(t(A)%*%chol2inv(C_P)%*%A)
  s<-matrix(rep(growth,dim(HH)[1]),nrow=dim(HH)[1])
  s_hat<-epsilon*t(A)%*%chol2inv(C_P)%*%s-matrix(rep(resp,dim(HH)[2]),nrow=dim(HH)[2])
  cChat<-chol(C_hat)
  b_ccc<-chol2inv(cChat)%*%s_hat
  b_R<-(s-A%*%b_ccc)/K
  return(list(R=b_R,C=b_ccc))
}

ser_ext_step<-function(sigma,agg,links){
  
  links<- rbind(links,exp(rnorm(1,0,sigma))*exp(agg))
  
  repeat{
    if(length(which(check_extinction(links)$R<Mmin))>0){
      links<- matrix(links[-which.max(links),],ncol=1)
    }else{break}
    
  }
  return(links)
  
}


### Parameters #####

sigma<-4
numb_sr<-300
Mmin<-10^-5
inv <- 1e-7 # invasion pressure for all species
epsilon<-0.1
resp<-0.1
growth<-0.1
# K<-5
K<-1
# K<-1
intra<-growth/K
factorr<-K*epsilon/resp
number_adds<-500
nreps<-100
rlx_t2<-1000
rlx_t<-500
# agg_1<-exp(-11.75);agg_2<-exp(-11);agg_3<-exp(-10.25)
agg_1<-10^-5.5;agg_2<-10^-4.75;agg_3<-10^-4
### Sub-model #####

for (k in 275){
  
  #####################
  # Initiate matrices #
  #####################
  
  links_1_DD_bb<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_2_DD_bb<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_3_DD_bb<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  
  links_1_ccc<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_2_ccc<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_3_ccc<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  
  links_1_cc<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_2_cc<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_3_cc<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_1_cc2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_2_cc2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_3_cc2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_1_cc2_2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_1_cc_2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_2_cc2_2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_2_cc_2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_3_cc2_2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  links_3_cc_2<-matrix(rep(0,number_adds*nreps),ncol=nreps)
  
  numb_sr<-k
  
  #####################
  numb_sr<-k
  
  for (j in 1:nreps){
    
    cat(k," ",signif(j/nreps*100,3),"%","\r")
    
    ##############################################
    # Find invassively fit consumers both models #
    ##############################################
    
    repeat{
      
      repeat{
        links_1<-matrix(exp(rnorm(numb_sr,0,sigma))*agg_1,nrow=1)
        if(sum(links_1)>resp) break
      }
      
      C_1<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
      C_1[1,]<-c(0,epsilon*links_1)
      
      for (i in 1:numb_sr){
        C_1[(i+1),1]<--links_1[i]
        C_1[(i+1),(i+1)]<--growth/K
      }
      
      pars_1  <- list(
        C = C_1,
        r = c(-resp,rep(growth,numb_sr)),
        i = rep(inv,numb_sr+1)
      )
      
      BRelaxed_1 <- rep(inv,numb_sr+1)
      BRelaxed_1 <- relax(pars_1,BRelaxed_1,rlx_t2)
      
      biom_1<-check_extinction(matrix(links_1,ncol=1))
      if(!is.na(BRelaxed_1[1])){
        if(BRelaxed_1[1]>(inv*10) & biom_1$C>0) break
      }
    }
    
    repeat{
      
      repeat{
        links_2<-matrix(exp(rnorm(numb_sr,0,sigma))*agg_2,nrow=1)
        if(sum(links_2)>resp) break
      }
      
      C_2<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
      C_2[1,]<-c(0,epsilon*links_2)
      
      for (i in 1:numb_sr){
        C_2[(i+1),1]<--links_2[i]
        C_2[(i+1),(i+1)]<--growth/K
      }
      
      pars_2  <- list(
        C = C_2,
        r = c(-resp,rep(growth,numb_sr)),
        i = rep(inv,numb_sr+1)
      )
      
      BRelaxed_2 <- rep(inv,numb_sr+1)
      BRelaxed_2 <- relax(pars_2,BRelaxed_2,rlx_t2)
      biom_2<-check_extinction(matrix(links_2,ncol=1))
      if(!is.na(BRelaxed_2[1])){
        if(BRelaxed_2[1]>(inv*10) & biom_2$C>0) break
      }
      
    }
    
    repeat{
      
      repeat{
        links_3<-matrix(exp(rnorm(numb_sr,0,sigma))*agg_3,nrow=1)
        if(sum(links_3)>resp) break
      }
      
      C_3<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
      C_3[1,]<-c(0,epsilon*links_3)
      
      for (i in 1:numb_sr){
        C_3[(i+1),1]<--links_3[i]
        C_3[(i+1),(i+1)]<--growth/K
      }
      
      pars_3  <- list(
        C = C_3,
        r = c(-resp,rep(growth,numb_sr)),
        i = rep(inv,numb_sr+1)
      )
      
      BRelaxed_3 <- rep(inv,numb_sr+1)
      BRelaxed_3 <- relax(pars_3,BRelaxed_3,rlx_t2)
      biom_3<-check_extinction(matrix(links_3,ncol=1))
      if(!is.na(BRelaxed_3[1])){
        if(BRelaxed_3[1]>(inv*10) & biom_3$C>0) break
      }
      
    }
    
    ################################################################
    # Create different variables for fulll and deconstructed model #
    ################################################################
    
    links_11<-links_1;links_22<-links_2;links_33<-links_3
    links_1<-matrix(links_1,ncol=1);links_2<-matrix(links_2,ncol=1);links_3<-matrix(links_3,ncol=1)
    
    ########################
    # Find pseudo invader #
    ########################
    
    repeat{
      links_r<-matrix(exp(rnorm(length(links_11),0,sigma))*mean(agg_1,agg_2,agg_3),nrow=1)
      if(sum(links_r)>1) break
    }
    
    ###################################
    # Save values for first time-step #
    ###################################
    {
      links_1_ccc[1,j]<-sum(pars_1$C[1,]/epsilon)
      links_2_ccc[1,j]<-sum(pars_2$C[1,]/epsilon)
      links_3_ccc[1,j]<-sum(pars_3$C[1,]/epsilon)
      
      links_1_DD_bb[1,j]<-1-sum(links_11^2)/(sum(links_11)^2)
      links_2_DD_bb[1,j]<-1-sum(links_22^2)/(sum(links_22)^2)
      links_3_DD_bb[1,j]<-1-sum(links_33^2)/(sum(links_33)^2)
      
      links_1_cc2[1,j]<-sum(links_11*links_r)
      links_1_cc2_2[1,j]<-sum(links_r^2)
      links_1_cc_2[1,j]<-sum(links_r)
      
      links_2_cc2[1,j]<-sum(links_22*links_r)
      links_2_cc2_2[1,j]<-sum(links_r^2)
      links_2_cc_2[1,j]<-sum(links_r)
      
      links_3_cc2[1,j]<-sum(links_33*links_r)
      links_3_cc2_2[1,j]<-sum(links_r^2)
      links_3_cc_2[1,j]<-sum(links_r)
      
      links_1_cc[1,j]<-sum(links_11)
      links_2_cc[1,j]<-sum(links_22)
      links_3_cc[1,j]<-sum(links_33)
    }
    ###################################
    
    for (i in 2:number_adds){
      
      # cat(k," ",signif(i/number_adds*100,3),"%","\r")
      
      {
        
        pars_1<- ser_ext_temp_num(pars_1,BRelaxed_1)
        links_1_ccc[i,j]<-sum(pars_1$C[1,]/epsilon)
        if(any(BRelaxed_1<(inv*10))){
          BRelaxed_1<-BRelaxed_1[-which(BRelaxed_1<(inv*10))]
        }
        links_1<-pars_1$C[1,]
        if(length(links_1)<k){
          new_link<- exp(rnorm(k-length(links_1),0,sigma))*agg_1
          for (ik in 1:(k-length(links_1))){
            pars_1$C<-rbind(pars_1$C,rep(0,dim(pars_1$C)[1]))
            pars_1$C<-cbind(pars_1$C,rep(0,dim(pars_1$C)[1]))
            pars_1$C[1,dim(pars_1$C)[1]]<-epsilon*new_link[ik]
            pars_1$C[dim(pars_1$C)[1],1]<--new_link[ik]
            pars_1$C[dim(pars_1$C)[1],dim(pars_1$C)[1]]<--growth/K
            pars_1$r<-c(pars_1$r,growth);pars_1$i<-c(pars_1$i,inv)
            BRelaxed_1<-c(BRelaxed_1,inv)
          }
        }else{
          repeat{
            temp<-links_1/epsilon
            new_link<-exp(rnorm(1,0,sigma))*agg_1
            replacee<-sample(1:length(temp),1)
            temp[replacee]<-new_link
            if(sum(factorr*temp)>1) break
          }
          pars_1$C[1,replacee]<-epsilon*new_link
          pars_1$C[replacee,1]<--new_link
          BRelaxed_1[replacee]<-inv
        }
        # BRelaxed_1<-c(BRelaxed_1,inv)
        BRelaxed_1<-as.numeric(relax(pars_1,BRelaxed_1,rlx_t))
        
        pars_2<- ser_ext_temp_num(pars_2,BRelaxed_2)
        links_2_ccc[i,j]<-sum(pars_2$C[1,]/epsilon)
        if(any(BRelaxed_2<(inv*10))){
          BRelaxed_2<-BRelaxed_2[-which(BRelaxed_2<(inv*10))]
        }
        links_2<-pars_2$C[1,]
        if(length(links_2)<k){
          new_link<- exp(rnorm(k-length(links_2),0,sigma))*agg_2
          for (ik in 1:(k-length(links_2))){
            pars_2$C<-rbind(pars_2$C,rep(0,dim(pars_2$C)[1]))
            pars_2$C<-cbind(pars_2$C,rep(0,dim(pars_2$C)[1]))
            pars_2$C[1,dim(pars_2$C)[1]]<-epsilon*new_link[ik]
            pars_2$C[dim(pars_2$C)[1],1]<--new_link[ik]
            pars_2$C[dim(pars_2$C)[1],dim(pars_2$C)[1]]<--growth/K
            pars_2$r<-c(pars_2$r,growth);pars_2$i<-c(pars_2$i,inv)
            BRelaxed_2<-c(BRelaxed_2,inv)
          }
        }else{
          repeat{
            temp<-links_2/epsilon
            new_link<-exp(rnorm(1,0,sigma))*agg_2
            replacee<-sample(1:length(temp),1)
            temp[replacee]<-new_link
            if(sum(factorr*temp)>1) break
          }
          pars_2$C[1,replacee]<-epsilon*new_link
          pars_2$C[replacee,1]<--new_link
          BRelaxed_2[replacee]<-inv
        }
        # BRelaxed_2<-c(BRelaxed_2,inv)
        BRelaxed_2<-as.numeric(relax(pars_2,BRelaxed_2,rlx_t))
        
        
        pars_3<- ser_ext_temp_num(pars_3,BRelaxed_3)
        links_3_ccc[i,j]<-sum(pars_3$C[1,]/epsilon)
        if(any(BRelaxed_3<(inv*10))){
          BRelaxed_3<-BRelaxed_3[-which(BRelaxed_3<(inv*10))]
        }
        links_3<-pars_3$C[1,]
        if(length(links_3)<k){
          new_link<- exp(rnorm(k-length(links_3),0,sigma))*agg_3
          for (ik in 1:(k-length(links_3))){
            pars_3$C<-rbind(pars_3$C,rep(0,dim(pars_3$C)[1]))
            pars_3$C<-cbind(pars_3$C,rep(0,dim(pars_3$C)[1]))
            pars_3$C[1,dim(pars_3$C)[1]]<-epsilon*new_link[ik]
            pars_3$C[dim(pars_3$C)[1],1]<--new_link[ik]
            pars_3$C[dim(pars_3$C)[1],dim(pars_3$C)[1]]<--growth/K
            pars_3$r<-c(pars_3$r,growth);pars_3$i<-c(pars_3$i,inv)
            BRelaxed_3<-c(BRelaxed_3,inv)
          }
        }else{
          repeat{
            temp<-links_3/epsilon
            new_link<-exp(rnorm(1,0,sigma))*agg_3
            replacee<-sample(1:length(temp),1)
            temp[replacee]<-new_link
            if(sum(factorr*temp)>1) break
          }
          pars_3$C[1,replacee]<-epsilon*new_link
          pars_3$C[replacee,1]<--new_link
          BRelaxed_3[replacee]<-inv
        }
        # BRelaxed_3<-c(BRelaxed_3,inv)
        BRelaxed_3<-as.numeric(relax(pars_3,BRelaxed_3,rlx_t))
      }
      
      {
        # links_1_DD_bb[i,j]<-1-sum(sort(links_11)^2)/(sum(sort(links_11))^2)
        links_11<- serial_extinction_t(links_11)
        # links_11<- c(links_11,exp(rnorm(1,0,sigma))*agg_1)
        repeat{
          links_r<-matrix(exp(rnorm(length(links_22),0,sigma))*mean(agg_1,agg_2,agg_3),nrow=1)
          if(sum(links_r)>1) break
        }
        links_1_DD_bb[i,j]<-1-sum(sort(links_11)^2)/(sum(sort(links_11))^2)
        links_1_cc2[i,j]<-sum(links_11*links_r)
        links_1_cc[i,j]<-sum(links_11)
        links_1_cc2_2[i,j]<-sum(links_r^2)
        links_1_cc_2[i,j]<-sum(links_r)
        
        if(length(links_11)<k){
          links_11<- c(links_11,exp(rnorm(k-length(links_11),0,sigma))*agg_1)
        }else{
          repeat{
            temp<-links_11
            temp[sample(1:length(temp),1)]<- exp(rnorm(1,0,sigma))*agg_1
            if(sum(temp)>1) break
          }
          links_11<-temp
        }
        # links_11<- serial_extinction_t(links_11)
        # links_1_cc[i,j]<-sum(links_11)


        # links_2_DD_bb[i,j]<-1-sum(links_22^2)/(sum(links_22)^2)
        links_22<- serial_extinction_t(links_22)
        # links_22<- c(links_22,exp(rnorm(1,0,sigma))*agg_2)
        repeat{
          links_r<-matrix(exp(rnorm(length(links_22),0,sigma))*mean(agg_1,agg_2,agg_3),nrow=1)
          if(sum(links_r)>1) break
        }
        links_2_DD_bb[i,j]<-1-sum(sort(links_22)^2)/(sum(sort(links_22))^2)
        links_2_cc2[i,j]<-sum(links_22*links_r)
        links_2_cc[i,j]<-sum(links_22)
        links_2_cc2_2[i,j]<-sum(links_r^2)
        links_2_cc_2[i,j]<-sum(links_r)
      
        if(length(links_22)<k){
          links_22<- c(links_22,exp(rnorm(k-length(links_22),0,sigma))*agg_2)
        }else{
          repeat{
            temp<-links_22
            temp[sample(1:length(temp),1)]<- exp(rnorm(1,0,sigma))*agg_2
            if(sum(temp)>1) break
          }
          links_22<-temp
        }
        # links_22<- serial_extinction_t(links_22)
        # links_2_cc[i,j]<-sum(links_22)


        
        # links_3_DD_bb[i,j]<-1-sum(links_33^2)/(sum(links_33)^2)
        links_33<- serial_extinction_t(links_33)
        # links_33<- c(links_33,exp(rnorm(1,0,sigma))*agg_3)
        repeat{
          links_r<-matrix(exp(rnorm(length(links_33),0,sigma))*mean(agg_1,agg_2,agg_3),nrow=1)
          if(sum(links_r)>1) break
        }
        links_3_DD_bb[i,j]<-1-sum(sort(links_33)^2)/(sum(sort(links_33))^2)
        links_3_cc2[i,j]<-sum(links_33*links_r)
        links_3_cc[i,j]<-sum(links_33)
        links_3_cc2_2[i,j]<-sum(links_r^2)
        links_3_cc_2[i,j]<-sum(links_r)
        
        if(length(links_33)<k){
          links_33<- c(links_33,exp(rnorm(k-length(links_33),0,sigma))*agg_3)
        }else{
          repeat{
            temp<-links_33
            temp[sample(1:length(temp),1)]<- exp(rnorm(1,0,sigma))*agg_3
            if(sum(temp)>1) break
          }
          links_33<-temp
        }
        # links_33<- serial_extinction_t(links_33)
        # links_3_cc[i,j]<-sum(links_33)
  
   

      }
    }
    
    write.table(links_1_cc,file="links_1_cc.txt");write.table(links_2_cc,file="links_2_cc.txt");write.table(links_3_cc,file="links_3_cc.txt")
    write.table(links_1_cc_2,file="links_1_cc_2.txt");write.table(links_2_cc_2,file="links_2_cc_2.txt");write.table(links_3_cc_2,file="links_3_cc_2.txt")
    write.table(links_1_cc2,file="links_1_cc2.txt");write.table(links_2_cc2,file="links_2_cc2.txt");write.table(links_3_cc2,file="links_3_cc2.txt")
    write.table(links_1_cc2_2,file="links_1_cc2_2.txt");write.table(links_2_cc2_2,file="links_2_cc2_2.txt");write.table(links_3_cc2_2,file="links_3_cc2_2.txt")
    write.table(links_1_ccc,file="links_1_ccc.txt");write.table(links_2_ccc,file="links_2_ccc.txt");write.table(links_3_ccc,file="links_3_ccc.txt")
    write.table(links_1_DD_bb,file="links_1_DD_bb.txt");write.table(links_2_DD_bb,file="links_2_DD_bb.txt");write.table(links_3_DD_bb,file="links_3_DD_bb.txt");
    
  }
}

write.table(links_1_cc,file="links_1_cc.txt");write.table(links_2_cc,file="links_2_cc.txt");write.table(links_3_cc,file="links_3_cc.txt")
write.table(links_1_cc_2,file="links_1_cc_2.txt");write.table(links_2_cc_2,file="links_2_cc_2.txt");write.table(links_3_cc_2,file="links_3_cc_2.txt")
write.table(links_1_cc2,file="links_1_cc2.txt");write.table(links_2_cc2,file="links_2_cc2.txt");write.table(links_3_cc2,file="links_3_cc2.txt")
write.table(links_1_cc2_2,file="links_1_cc2_2.txt");write.table(links_2_cc2_2,file="links_2_cc2_2.txt");write.table(links_3_cc2_2,file="links_3_cc2_2.txt")
write.table(links_1_ccc,file="links_1_ccc.txt");write.table(links_2_ccc,file="links_2_ccc.txt");write.table(links_3_ccc,file="links_3_ccc.txt")
write.table(links_1_DD_bb,file="links_1_DD_bb.txt");write.table(links_2_DD_bb,file="links_2_DD_bb.txt");write.table(links_3_DD_bb,file="links_3_DD_bb.txt");

