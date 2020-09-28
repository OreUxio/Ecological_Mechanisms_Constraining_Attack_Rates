library(plyr)
library(deSolve)

###################
### Functions #####
###################

parColMax <- function(M){
  colMaxs(M,value=T,parallel=T)
}

parColWhichMax <- function(M){
  colMaxs(M,value=F)
}

parColWhichMins <- function(M){
  colMins(M,value=F)
}

pickRows <- function(M,indices){
  M[matrix(c(indices,seq(length(indices))),ncol=2)]
}

parColSums <- function(M){
  colsums(M,parallel=T)
}

parColWhichMax <- function(M){
  colMaxs(M,value=F)
}


sample_and_ser_ext<-function(agg){
  repeat{
    links<-matrix(exp(rnorm(numb_sr,0,sigma))*agg,nrow=1)*alpha_0
    if(sum(links)>1) break
  }
  links<- serial_extinction_t(links)
  return(links)
}

sample_links<-function(agg){
  repeat{
    links<-matrix(exp(rnorm(numb_sr,0,sigma))*agg,nrow=1)*alpha_0
    if(sum(links)>1) break
  }
  return(links)
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

###################
### Parameters ####
###################

a1<-(0.8^0.5)*1 #as.numeric(args[2])
a2 <- (1.3^0.5)*1
sigma<-4
epsilon<-0.1
resp<-0.1
growth<-0.1
K<-1
inv<-10^-6
rlx_t2<-2000
alpha_0<-epsilon*K/resp
amount_def<-0
amount_inc<-0.05

##############################################################################
### Experiment 1) Both consumer share the same resource:  max_C1 = max_C2 ####
##############################################################################

aggz<-10^seq(-5,-4,length.out = 10) # Granularity of aggressiveness values
agg_all<-rep(aggz,100) # Number of replicates per aggressiveness value
amount_list_diff<-c()
numb_sr<-300# Number of resources

# Generate list of links that can overcome respiration per aggressiveness values:
links_1<-llply(agg_all,function(x) sample_and_ser_ext(x)) 

for (j in 1:length(links_1)){
  cat(j/length(agg_all)*100,"%","\r")
  # Build interaction matrix with one consmer two resources
    repeat{
      # Sample new consumer and set condition for each consumer to rely on different resources
      repeat{
      inv_i <- sample(1:length(aggz),1)
      # Serial extinction of C1
      C1_temp<-sample_and_ser_ext(aggz[inv_i])
      C1_max<-which.max(C1_temp)
      inv_i <- sample(1:length(aggz),1)
      # Serial extinction of C2
      C2 <- sample_and_ser_ext(aggz[inv_i])
        if(length(C2)<length(C1_temp)){
          C1_temp<-C1_temp[-sample((1:length(C1_temp))[-C1_max],length(C1_temp)-length(C2))]
        } else if (length(C2)>length(C1_temp)){
          C2<-C2[-sample((1:length(C2))[-which.max(C2)],length(C2)-length(C1_temp))]
        }
      
      ###################
      ## Conditions #####
      ###################
      
        if(which.max(C2) == C1_max & max(C1_temp) < C2[C1_max] & sum(C1_temp)>1 & sum(C2)>1) break
      }
      C1 <- C1_temp
      dominance <- max(C1)> max(C2)
      
      numb_sr<-length(C1)
      C_matrix<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
      C_matrix[1,]<-c(0,epsilon*C1/alpha_0)
      for (i in 1:numb_sr){
        C_matrix[(i+1),1]<--C1[i]/alpha_0
        C_matrix[(i+1),(i+1)]<--growth/K
      }
      pars_1  <- list(
        C = C_matrix,
        r = c(-resp,rep(growth,numb_sr)),
        i = rep(inv,numb_sr+1)
      )
      C_matrix<-cbind(C_matrix,rep(0,numb_sr+1))
      C_matrix<-rbind(C_matrix,rep(0,numb_sr+2))
      
      C_matrix[numb_sr+2,] <- c(0,epsilon*C2/alpha_0,0)
      for (i in 1:numb_sr){
        C_matrix[(i+1),numb_sr+2] <- -C2[i]/alpha_0
      }
      
      pars_1  <- list(
        C = C_matrix,
        r = c(-resp,rep(growth,numb_sr),-resp),
        i = rep(inv,numb_sr+2)
      )
      
      BRelaxed_1<-c("no","no")
      try({BRelaxed_1 <- -solve(C_matrix)%*%as.matrix(pars_1$r,col=1)},silent = T)
      if(is.numeric(BRelaxed_1[1])) break
      
    }
  amount_list_diff<-rbind(amount_list_diff,list(d1=BRelaxed_1[1]<inv*10,d2=BRelaxed_1[1+which.max(C2)]<inv*10))
}

mean( unlist(amount_list_diff[,1]) )

# Fraction of times shared main resource goes extinct
mean( unlist(amount_list_diff[,2]) )

##############################################################################
### Reverse Experiment 1) Both consumer do not share the same resource:  max_C1 = max_C2 ####
##############################################################################

aggz<-10^seq(-5,-4,length.out = 10) # Granularity of aggressiveness values
agg_all<-rep(aggz,100) # Number of replicates per aggressiveness value
amount_list_diff<-c()
numb_sr<-300# Number of resources

# Generate list of links that can overcome respiration per aggressiveness values:
links_1<-llply(agg_all,function(x) sample_and_ser_ext(x)) 

for (j in 1:length(links_1)){
  cat(j/length(agg_all)*100,"%","\r")
  # Build interaction matrix with one consmer two resources
  repeat{
    # Sample new consumer and set condition for each consumer to rely on different resources
    
    repeat{
      inv_i <- sample(1:length(aggz),1)
      # Serial extinction of C1
      C1_temp<-sample_and_ser_ext(aggz[inv_i])
      C1_max<-which.max(C1_temp)
      inv_i <- sample(1:length(aggz),1)
      # Serial extinction of C2
      C2 <- sample_and_ser_ext(aggz[inv_i])
      if(length(C2)<length(C1_temp)){
        C1_temp<-C1_temp[-sample((1:length(C1_temp))[-C1_max],length(C1_temp)-length(C2))]
      } else if (length(C2)>length(C1_temp)){
        C2<-C2[-sample((1:length(C2))[-which.max(C2)],length(C2)-length(C1_temp))]
      }
      
      ###################
      ## Conditions #####
      ###################
      
      if(which.max(C2) != C1_max & sum(C1_temp)>1 & sum(C2)>1) break
    }
    C1 <- C1_temp
    dominance <- max(C1)> max(C2)
    
    numb_sr<-length(C1)
    C_matrix<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
    C_matrix[1,]<-c(0,epsilon*C1/alpha_0)
    for (i in 1:numb_sr){
      C_matrix[(i+1),1]<--C1[i]/alpha_0
      C_matrix[(i+1),(i+1)]<--growth/K
    }
    pars_1  <- list(
      C = C_matrix,
      r = c(-resp,rep(growth,numb_sr)),
      i = rep(inv,numb_sr+1)
    )
    C_matrix<-cbind(C_matrix,rep(0,numb_sr+1))
    C_matrix<-rbind(C_matrix,rep(0,numb_sr+2))
    
    C_matrix[numb_sr+2,] <- c(0,epsilon*C2/alpha_0,0)
    for (i in 1:numb_sr){
      C_matrix[(i+1),numb_sr+2] <- -C2[i]/alpha_0
    }
    
    pars_1  <- list(
      C = C_matrix,
      r = c(-resp,rep(growth,numb_sr),-resp),
      i = rep(inv,numb_sr+2)
    )
    
    BRelaxed_1<-c("no","no")
    try({BRelaxed_1 <- -solve(C_matrix)%*%as.matrix(pars_1$r,col=1)},silent = T)
    if(is.numeric(BRelaxed_1[1])) break
    
  }
  amount_list_diff<-rbind(amount_list_diff,list(d1=BRelaxed_1[1]<inv*10,d2=BRelaxed_1[1+which.max(C2)]<inv*10))
}

mean( unlist(amount_list_diff[,1]) )

# Fraction of times shared main resource goes extinct
mean( unlist(amount_list_diff[,2]) )

######################################################################################################################################
## MAIN Type I: Experiment 2.a) Condition that consumers DO NOT share the same resource AND C1's main trophic link is weaker than equivalent of C2 ##
######################################################################################################################################
### i.e : (a) A_{max_C1,1} < A_{max_C1,2} & (b) A_{max_C1,2} != A_{max_C2,2} (i.e.  A_{max_C1,2} < A_{max_C2,2} )

# Note accuracy improves as aggressiveness increases

# aggz<-10^seq(-5,-4,length.out = 10) # Granularity of aggressiveness values
agg_all<-10^seq(-4.5,-4.5,length.out = 500) # Number of replicates per aggressiveness value
amount_list_1_diff<-c()
 # Number of resources
# inv_i_all<-sample(1:length(aggz),length(agg_all),replace = T)
# Generate list of links that can overcome respiration per aggressiveness values:

args<-c(100,200,300,400,500)
# Be sure that no serial extinction occurs!
for (ll in args){
  numb_sr_og<-ll
  for (j in 1:length(agg_all)){
  # Build interaction matrix with one consmer two resources
  # cat(j/length(agg_all)*100,"%","\r")
  repeat{
  numb_sr<-numb_sr_og
  {
  # repeat{
  #   
  #   repeat{
  #     agg_1<-sample(agg_all,1)
  #     repeat{
  #       C1_temp <- sample_and_ser_ext(agg_1)
  #       if(length(C1_temp)==numb_sr) break
  #     }
  #     C1_max<-which.max(C1_temp)
  #     agg_2<-sample(agg_all,1)
  #     repeat{
  #       C2 <- sample_and_ser_ext(agg_2)
  #       if(length(C2)==numb_sr) break
  #     }
  #     C1_max<-which.max(C1_temp)
  #     if(sum(C1_temp)>1 & sum(C2)>1 & max(C1_temp) < max(C2) & length(serial_extinction_t(C2))==length(C2) &
  #        length(serial_extinction_t(C1_temp))==length(C1_temp)) {break}
  #   }
  #     
  #     C2_max<-which.max(C2)
  #     C1_max<-which.max(C1_temp)
  #     C2[C1_max]<-sample(seq(C1_temp[C1_max],C2[C2_max],length.out=100)[-1],1) # Ensures Phyrric competition
  #     t_HH<-sum(C1_temp*C2) 
  #     
  #     # No exploitative competition
  #     if (length(C2)>=3 &  
  #         ( ((sum(C1_temp)-1)/t_HH) / ( (sum(C2)-1)/sum(C2^2) ) > 1) &
  #         ( ((sum(C2)-1)/t_HH) / ( (sum(C1_temp)-1)/sum(C1_temp^2) ) > 1)  ){break}
  # }
  }
  # Sample new consumer and set condition for each consumer to rely on different resources
  {
    repeat{
    repeat{
      # Serial extinction of C1
      agg_1<-sample(agg_all,1)
      C1 <- sample_and_ser_ext(agg_1)
      C1_temp<-C1
      C1_max<-which.max(C1_temp)
      agg_2<-sample(agg_all,1)
  
      # Serial extinction of C2 #
      C2 <- sample_and_ser_ext(agg_2)
      if(length(C2)<length(C1_temp)){
        C1_temp<-C1_temp[-sample((1:length(C1_temp))[-C1_max],length(C1_temp)-length(C2))]
      } else if (length(C2)>length(C1_temp)){
        C2<-C2[-sample((1:length(C2))[-which.max(C2)],length(C2)-length(C1_temp))]
      }
      if(sum(C1_temp)>1 & sum(C2)>1 & max(C1_temp) < max(C2) & which.max(C2) != C1_max) break
    }
    C1_max<-which.max(C1_temp)
    # repeat{
      C2[C1_max]<-sample(seq(C1_temp[C1_max],max(C2),length.out=100)[-1],1) # Enforce Phyrric competition
      # C2[C1_max]<-max(C2)*0.9#sample(seq(C1_temp[C1_max],,length.out=100)[-1],1)
      t_HH<-sum(C1_temp*C2)
      if (length(C2)>=3  &
        C1_temp[C1_max] < C2[C1_max] &
        ( ((sum(C1_temp)-1)/t_HH) / ( (sum(C2)-1)/sum(C2^2) ) > 1) &
        ( ((sum(C2)-1)/t_HH) / ( (sum(C1_temp)-1)/sum(C1_temp^2) ) > 1) &
        length(serial_extinction_t(C2))==length(C2) ){break}
    # }
    # break
    }
  }
    
    C1<-C1_temp
    numb_sr<-length(C1)
    C_matrix<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
    C_matrix[1,]<-c(0,epsilon*C1/alpha_0)
    for (i in 1:numb_sr){
      C_matrix[(i+1),1]<--C1[i]/alpha_0
      C_matrix[(i+1),(i+1)]<--growth/K
    }
    pars_1  <- list(
      C = C_matrix,
      r = c(-resp,rep(growth,numb_sr)),
      i = rep(inv,numb_sr+1)
    )
    C_matrix<-cbind(C_matrix,rep(0,numb_sr+1))
    C_matrix<-rbind(C_matrix,rep(0,numb_sr+2))
    
    C_matrix[numb_sr+2,] <- c(0,epsilon*C2/alpha_0,0)
    for (i in 1:numb_sr){
      C_matrix[(i+1),numb_sr+2] <- -C2[i]/alpha_0
    }
    
    pars_1  <- list(
      C = C_matrix,
      r = c(-resp,rep(growth,numb_sr),-resp),
      i = rep(inv,numb_sr+2)
    )
    
    BRelaxed_1<-c("no","no")
    try({BRelaxed_1 <- -solve(C_matrix)%*%as.matrix(pars_1$r,col=1)},silent = T)
    if(is.numeric(BRelaxed_1[1]) & BRelaxed_1[1]>0 & tail(BRelaxed_1,1) > 0 ) break
  }
  # Check d1 = if C1 goes extinct and, d2 = if C2 main resource goes extinct
  amount_list_1_diff<-rbind(amount_list_1_diff,list(d1=BRelaxed_1[1+which.max(C1)]<inv*10,d2= C1[C1_max] < C2[C1_max]))
  }
  cat(mean( unlist(amount_list_1_diff[,1])),"\r")
  write.csv(cbind(amount_list_1_diff,agg_all),file=paste("/Data/Raw/Validation/",ll,"phyrric_tests.txt",sep=""))
}

######################################################################################################################################
## MAIN Type II: Experiment 2.a) Condition that consumers DO NOT share the same resource AND C1's main trophic link is stonger than equivalent of C2 ##
######################################################################################################################################
### i.e : (a) A_{max_C1,1} < A_{max_C1,2} & (b) A_{max_C1,2} != A_{max_C2,2} (i.e.  A_{max_C1,2} < A_{max_C2,2} )

# Note accuracy improves as aggressiveness increases

# aggz<-10^seq(-5,-4,length.out = 10) # Granularity of aggressiveness values
amount_list_2_diff<-c()
# Number of resources
# inv_i_all<-sample(1:length(aggz),length(agg_all),replace = T)
# Generate list of links that can overcome respiration per aggressiveness values:

args<-c(100,200,300,400,500)
for (ll in args){
  numb_sr_og<-ll
  numb_sr<-ll
  for (j in 1:length(agg_all)){
    # Build interaction matrix with one consmer two resources
    cat(j/length(agg_all)*100,"%","\r")
    repeat{
      numb_sr<-numb_sr_og
      {
        # repeat{
        #   
        #   repeat{
        #     agg_1<-sample(agg_all,1)
        #     repeat{
        #       C1_temp <- sample_and_ser_ext(agg_1)
        #       if(length(C1_temp)==numb_sr) break
        #     }
        #     C1_max<-which.max(C1_temp)
        #     agg_2<-sample(agg_all,1)
        #     repeat{
        #       C2 <- sample_and_ser_ext(agg_2)
        #       if(length(C2)==numb_sr) break
        #     }
        #     C1_max<-which.max(C1_temp)
        #     if(sum(C1_temp)>1 & sum(C2)>1 & max(C1_temp) < max(C2) & length(serial_extinction_t(C2))==length(C2) &
        #        length(serial_extinction_t(C1_temp))==length(C1_temp)) {break}
        #   }
        #     
        #     C2_max<-which.max(C2)
        #     C1_max<-which.max(C1_temp)
        #     C2[C1_max]<-sample(seq(C1_temp[C1_max],C2[C2_max],length.out=100)[-1],1) # Ensures Phyrric competition
        #     t_HH<-sum(C1_temp*C2) 
        #     
        #     # No exploitative competition
        #     if (length(C2)>=3 &  
        #         ( ((sum(C1_temp)-1)/t_HH) / ( (sum(C2)-1)/sum(C2^2) ) > 1) &
        #         ( ((sum(C2)-1)/t_HH) / ( (sum(C1_temp)-1)/sum(C1_temp^2) ) > 1)  ){break}
        # }
      }
      # Sample new consumer and set condition for each consumer to rely on different resources
      {
          repeat{
            # Serial extinction of C1
            agg_1<-sample(agg_all,1)
            C1 <- sample_and_ser_ext(agg_1)
            C1_temp<-C1
            C1_max<-which.max(C1_temp)
            agg_2<-sample(agg_all,1)
            
            # Serial extinction of C2 #
            C2 <- sample_and_ser_ext(agg_2)
            if(length(C2)<length(C1_temp)){
              C1_temp<-C1_temp[-sample((1:length(C1_temp))[-C1_max],length(C1_temp)-length(C2))]
            } else if (length(C2)>length(C1_temp)){
              C2<-C2[-sample((1:length(C2))[-which.max(C2)],length(C2)-length(C1_temp))]
            }
            C1_max<-which.max(C1_temp)
            t_HH<-sum(C1_temp*C2)
            if(sum(C1_temp)>1 & sum(C2)>1 & which.max(C2) != C1_max & length(C2)>=3  & C1_temp[C1_max] > C2[C1_max] &
               ( ((sum(C1_temp)-1)/t_HH) / ( (sum(C2)-1)/sum(C2^2) ) > 1) &
               ( ((sum(C2)-1)/t_HH) / ( (sum(C1_temp)-1)/sum(C1_temp^2) ) > 1)) {break}
          }
      }
      
      C1<-C1_temp
      numb_sr<-length(C1)
      C_matrix<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
      C_matrix[1,]<-c(0,epsilon*C1/alpha_0)
      for (i in 1:numb_sr){
        C_matrix[(i+1),1]<--C1[i]/alpha_0
        C_matrix[(i+1),(i+1)]<--growth/K
      }
      pars_1  <- list(
        C = C_matrix,
        r = c(-resp,rep(growth,numb_sr)),
        i = rep(inv,numb_sr+1)
      )
      C_matrix<-cbind(C_matrix,rep(0,numb_sr+1))
      C_matrix<-rbind(C_matrix,rep(0,numb_sr+2))
      
      C_matrix[numb_sr+2,] <- c(0,epsilon*C2/alpha_0,0)
      for (i in 1:numb_sr){
        C_matrix[(i+1),numb_sr+2] <- -C2[i]/alpha_0
      }
      
      pars_1  <- list(
        C = C_matrix,
        r = c(-resp,rep(growth,numb_sr),-resp),
        i = rep(inv,numb_sr+2)
      )
      
      BRelaxed_1<-c("no","no")
      try({BRelaxed_1 <- -solve(C_matrix)%*%as.matrix(pars_1$r,col=1)},silent = T)
      if(is.numeric(BRelaxed_1[1]) & BRelaxed_1[1]>0 & tail(BRelaxed_1,1) > 0 ) break
    }
    # Check d1 = if C1 goes extinct and, d2 = if C2 main resource goes extinct
    amount_list_2_diff<-rbind(amount_list_2_diff,list(d1=BRelaxed_1[1+which.max(C1)]>inv*10,d2= C1[C1_max] > C2[C1_max]))
  }
  # cat(mean( unlist(amount_list_2_diff[,1])),"\r")
  write.csv(cbind(amount_list_2_diff,agg_all),file=paste("/Data/Raw/Validation/",ll,"phyrric_2_tests.txt",sep=""))
}
# Fraction of times  C1 main resource goes extinct
mean( unlist(amount_list_2_diff[,1]))

# Fraction of times C1's main equivalent resource is larger
mean( unlist(amount_list_2_diff[,2]))