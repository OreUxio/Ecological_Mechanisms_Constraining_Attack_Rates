args = commandArgs(trailingOnly=TRUE)
library(quantreg)
library(Rfast)
library(plyr)

######## Functions ##########
findFitConsumer_tt <- function(agg){
  
repeat{
    H <- agg*exp(rnorm(Sr,sd=sigma))*K*epsilon/resp
    if(sum(H)>1) break
}
  return(list(H=H,agg=agg))
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

find_mu_sd<-function(s_r,sigmaa){
  #Function for estimating the expected value of the strongest link strength, i.e. E(max(exp(sigmaa*xi_{max,i})))
  
  n_samples<-1000#Number of samples
  sum_m<-rep(0,n_samples)
  for (i in 1:n_samples){
    sum_m[i]<-sum(exp(sigmaa*rnorm(s_r,0,1)))
  }
  mean_sum_m<-mean(log(sum_m))
  sd_sum_m<-(var(log(sum_m)))^0.5
  out<-list(mean_sum_m,sd_sum_m)
  return(out)
}

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

## Is the invader outcompeted?
outcompeted <- function(H,t_H_HH,HHSums,HH2Sums){
  directly <-
    any( (HHSums-1) * t_H_HH >
           (sum(H)-1) * HH2Sums )
  if(directly) return(T)
  return(directly) #directly |
}

## Does the invader outcompete residents?
## [ equivalent to apply(X = HH, MARGIN = 2, FUN = outcompeted,matrix(H,ncol=1)) ]
outcompetes_directly <- function(H,t_H_HH,HHSums){
  directly <- (sum(H)-1) * t_H_HH > (HHSums-1) * sum(H^2)
  return(directly) 
}

## Which resident (output row) outcompetes which other resident (output column)?
## [ equivalent to t(apply(X = HH, MARGIN = 2, FUN = outcompetes, HH)),
## except for diagonal]
outcompetition <- function(HH){
  mainResource <- parColWhichMax(HH)
  indirectly <-  t( HH[mainResource,] > pickRows(HH,mainResource) )
  deadd<-colSums(indirectly)>0
  if(any(deadd)){
    victim<-which(deadd)
    assassin<-llply(victim,function(x) which(indirectly[,x]==1))
    victim<-llply(1:length(victim),
                  function(x) victim[x]*(!(mainResource[victim[x]] %in% mainResource[unlist(assassin[x])])))
  }
  indirectly<-rep(F,dim(HH)[2])
  if(any(deadd)){indirectly[unlist(victim)]<-T}
  return(indirectly)
}

outcompetition_directly <- function(HH,sum_HH,sum_HH2){
  colSumsHHm1 <- sum_HH-1
  directly <-
    colSumsHHm1 * mat.mult(t(HH), HH) >
    outer(  sum_HH2, colSumsHHm1  )
  diag(directly) <- FALSE #should always be FALSE, apparently some numerical issues
  return( directly)
}

findFitConsumer <- function(aggs){
  repeat{
    parent<-sample(aggs,1)
    agg <- a1*a2^rnorm(1)*parent
    H <- agg*exp(rnorm(Sr,sd=sigma))*constanttt
    if(sum(H)>1) break
  }
  return(list(H=H,agg=agg,parent_number=which(aggs==parent)))
}

## One step in serial extinction for one consumer
resourceEx1 <- function(H){
  if((sum(H)-1)*max(H) > sum(H^2)){
    return(which.max(H));
  }else{
    return(0);
  }
}

## One step in serial extinction for many consumers
resourceEx <- function(HH){
  ex <-
    (parColSums(HH)-1)*parColMax(HH) > parColSums(HH^2)
  if(any(ex)){
    return(list(apply(HH[,ex,drop=F],MARGIN = 2,which.max),ex));
  }else{
    return(list(0,ex));
  }
}

check_extinction<-function(HH){
  A<-resp*HH/(K*epsilon)
  C_P<-diag(rep(K,dim(HH)[1]))
  C_hat<-epsilon*(t(A)%*%chol2inv(C_P)%*%A)
  s<-matrix(rep(growth,dim(HH)[1]),nrow=dim(HH)[1])
  s_hat<-epsilon*t(A)%*%chol2inv(C_P)%*%s-matrix(rep(resp,dim(HH)[2]),nrow=dim(HH)[2])
  b_C<-chol2inv(C_hat)%*%s_hat
  b_R<-(s-A%*%b_C)/K
  return(b_R)
}

## Parameters ############

doSerialEx_bin <- T
doSerialEx <- T # set to F for no serial extinction
specialist_check <- T
doOver <- T
rand_samp_animals <- T
choose_animals<-c(T,F)

analytic_check <- T
machine_learn_f1 <- F
doComp_dir <- T
Resp_ser_ext <- T
Kick_res <- T
Resp_over <- T
Resp_res <- T
do_stats <- F
doPyrrhic<-T

cut_gran<-100
cut_gran_com<-150
aim_HH<-150
steps <- 2e06
cut_gran_HH<-round(steps/50)
plot_gran<-20

sigma<-as.numeric(args[1]) 
sigmaa<-sigma
a1 <- (0.8^0.5)*as.numeric(args[2]) ##1.0125#1.05
a2 <- (1.3^0.5)*1 #1.01#1.3^0.5
a0 <- 10^as.numeric(args[3])
K<-1 #as.numeric(args[4])
epsilon<-0.1;growth<-1;resp<-0.1
constanttt<-epsilon*K/resp
if(constanttt==0){
  epsilon<-as.numeric(args[4])
  resp<-as.numeric(args[5])
  K<-as.numeric(args[6])
}
disp_prob<-0.2

print(c(sigma,a1,a0))
Mmin<-1e-7
#(sigma) #1e-02 #exp(-8)#initial_aggressivity
a1_first<-a1
a1_old<-a1
first<-0
perc<-0.01
gran<-0.001
range<-0
iii<-1
exp_burn_in <- 4000
Sr <- 20
Sc <- 10


{
  false_true_negative_positive_gen<-c(0,0,0,0) # TN, TP, FN, FP.
  false_true_negative_positive_spec<-c(0,0,0,0)
  false_true_negative_positive_over<-c(0,0,0,0)
  fish.test.check<-500
  # Variable parameters 
  
  # Make sure all consumers are reasonably fit ( sum(H) > 1 )
  aggs <- a0*abs(rnorm(Sc)*2)# adjust initial aggs here     # adjust initial aggs here
  HH <- matrix(exp(rnorm(Sr*Sc,sd = sigma))*aggs,nrow=Sr,ncol=Sc,byrow = T)*constanttt

  repeat{
    
    unfit <- colSums(HH) < 1;
    Sex <- sum(unfit)
    HH[,which(unfit)] <- matrix(exp(rnorm(Sr*Sex,sd = sigma))*aggs[unfit],
                                nrow=Sr,ncol=Sex,byrow = T)*constanttt
    # repeat{
    #   unfit <- which(outcompetition(HH))
    #   HH[,unfit] <- matrix(exp(rnorm(Sr*Sex,sd = sigma))*aggs[unfit],
    #                               nrow=Sr,ncol=Sex,byrow = T)*constanttt
    #   unfit <- which(outcompetition(HH))
    #   if(length(which(unfit))==0){break}
    # }
    if(!any(sum(colSums(HH) < 1))){break}
  }
  
  total_species<-cut_gran*2
  Sc<-length(aggs)
  Sr<-dim(HH)[1]
  {
    list_species<-as.data.frame(list(
      agg=c(rep(0,total_species)),#1
      id=c(rep(0,total_species)),#2
      outi=c(rep(-1,total_species)),#3
      extinction=c(rep(-1,total_species)),#4
      parent_id=c(rep(-1,total_species)), #5
      n_C=c(rep(-1,total_species)), #6
      n_P=c(rep(-1,total_species)), #7
      m_l_a=c(rep(-1,total_species)),#8
      v_l_a=c(rep(-1,total_species)),#9
      inv_counter=c(rep(-1,total_species))#10
      # A1_m_ja=rep(-1,total_species),#11
      # A2_m_ja=rep(-1,total_species),#12
      # B1_m_ja=rep(-1,total_species),#13
      # B2_m_ja=rep(-1,total_species),#14
      # ine_M_ja=rep(-1,total_species),#15
      # A1_m_jb=rep(-1,total_species),#16
      # A2_m_jb=rep(-1,total_species),#17
      # B1_m_jb=rep(-1,total_species),#18
      # B2_m_jb=rep(-1,total_species),#119
      # ine_M_jb=rep(-1,total_species),
      # max_ja=rep(-1,total_species),
      # max_jb=rep(-1,total_species),
      # DD_ja=rep(-1,total_species),
      # DD_jb=rep(-1,total_species)#20
    ))
    
    gran_table_ext<-as.data.frame(list(
      agg=c(rep(0,total_species)),
      id=c(rep(0,total_species)),
      A1=c(rep(0,total_species)),
      A2=c(rep(0,total_species)),
      B1=c(rep(0,total_species)),
      B2=c(rep(0,total_species)),
      A1_l_c=c(rep(0,total_species)),
      A2_l_c=c(rep(0,total_species)),
      B1_l_c=c(rep(0,total_species)),
      B2_l_c=c(rep(0,total_species)),
      ext_type=c(rep(0,total_species))
    ))
    
  }
  num_dead<-0
  n_C<-dim(HH)[2]
  n_P<-dim(HH)[1]
  num_inv<-0
  counter<-0
  only_ans<-list_species[which(list_species$outi==0),]
  n_plants_added<-n_P
  i<-1
  # try(load("deconstructed_PDMM.RData"))
  try(system("rm *.txt && rm *.csv"),silent=T)
  
  if(i>1){
    steps_start<-i
  }else{steps_start<-1}
  
  history <- c(exp(mean(log(aggs))),min(aggs),max(aggs),1-mean(parColSums(HH^2)/parColSums(HH)^2),Sc)
  
  par(mfrow=c(1,1))
  n_C<-dim(HH)[2]
  n_P<-dim(HH)[1]
  
  # plot(0,n_P,col="green",cex=0.1,ylim=c(0,150),xlim=c(0,steps),xlab="",ylab="Number of Species")
  # points(0,n_C,col="red",cex=0.1)
  
  only_ans<-as.data.frame(list(
    agg=aggs,#1
    id=rep(0,n_C),#2
    outi=rep(0,n_C),#3
    extinction=rep(-1,n_C),#4
    parent_id=rep(-1,n_C),
    n_C=rep(-1,n_C),
    n_P=rep(-1,n_C),
    m_l_a=rep(-1,n_C),
    v_l_a=rep(-1,n_C),
    inv_counter=rep(-1,n_C)
  ))
  
  par_col_sum_yo <- parColSums(HH)
  par_col_sum_yo2 <- parColSums(HH^2)
  t_inv_HHH <- mat.mult(t(HH),HH)
  ine_all <-  (t(((par_col_sum_yo -1) / t_inv_HHH))) / ( (par_col_sum_yo -1)/par_col_sum_yo2 )
  num_save <- 0
  gran_table_ext_i<-0
  HH_temp_g<-HH
}
## simulation ##########

for(i in 1:steps){
  
  
  {
    n_C<-dim(HH)[2]
    n_P<-dim(HH)[1]
    
    # Sample producer/consumer
    if(rand_samp_animals){
      animal<-sample(choose_animals,1)
    }else{
      if(runif(1)>n_C/(n_P+n_C)){
        animal<-F
      }else{
        animal<-T
      }
    }
  }
  
  
  if(animal){ 
    
  who<-runif(n_C)<=disp_prob
  if(sum(who)>0){
    for (w in 1:sum(who)){
      
    {
      par_col_sum_yo<-parColSums(HH)
      par_col_sum_yo2<-parColSums(HH^2)
      num_inv<-num_inv+1
      Sr<-dim(HH)[1]
      # Find new invader
      counter<-0
      repeat{
        counter<-counter+1
        Sr<-dim(HH)[1]
        invader<- findFitConsumer(only_ans[,1])# Caan it overcome respiration?
        t_inv_HH<-t(invader$H)%*%HH 
        d_a <- !outcompeted(invader$H,t_inv_HH, par_col_sum_yo, par_col_sum_yo2 ) # Is it outcompeted?
        
        if( d_a ) break;
      }
      
      ## Summary stats: The code now requires that id in line Sc_start+num_inv is exactly Sc_start+num_inv
      {
        new_inv<-c(
          invader$agg, #1
          num_inv,     #2
          0,           #3
          0,           #4
          only_ans[invader$parent_number,2],#5
          n_C, #6
          n_P, #7
          m_l_a=mean(log(only_ans[,1])), #8
          v_l_a=var(log(only_ans[,1])), #9
          inv_counter=counter #10
        )
      }
    }
    
    still<-1
    
    # Check to see if strongest link of new arrival causes overexploitation of main resource
    if(max(invader$H)<=log(K/Mmin)){
      
      # Do serial extinction of invader, pseudo invaders (used for summary stats only) and residents
      {
        if(doSerialEx_bin){
          if(doSerialEx){
            exR <- resourceEx1(invader$H)[[1]]
            while(sum(exR)){
              ex<-exR
              HH<-HH[-ex, , drop=F]
              invader$H<-invader$H[-ex]
              exR <- resourceEx1(invader$H)
            }  
          }else{
            temp_H<-invader$H
            exR <- resourceEx1(temp_H)[[1]]
            while(sum(exR)){
              ex_r <- sample(1:exR,1)
              ex<-exR
              HH<-HH[-ex_r, ,drop=F]
              invader$H<-invader$H[-ex_r]
              temp_H<-temp_H[-ex]
              HH_temp_g<-HH_temp_g[-ex,]
              exR <- resourceEx1(temp_H)
            }  
          }
        }
      }
      
      # Check which species cannot overcome respiration
      if(Resp_ser_ext){
        exC1 <- parColSums(HH) < 1  
        if(any(exC1)){
          only_ans[which(exC1),c(3,4)]<-matrix(rep(c(num_inv,2),sum(exC1)),ncol=2,byrow=T)
          list_species[(num_dead+1:sum(exC1)),]<-only_ans[which(exC1),]
          HH <- HH[,which(!exC1), drop=F];Sc <- length(only_ans[,1]); 
          only_ans<-only_ans[which(!exC1),]
          if(!doSerialEx){ HH_temp_g <- HH_temp_g[,which(!exC1), drop=F]}
          num_dead<-num_dead+sum(exC1)
        }
      } 
      
      over<-0
      
    }else if(doOver){  # Overexploitation # 
      # New invader overexploits most important resource.
      over<-1
      over_2<-T
      still<-1
      # Removal of overexploited resource from community.
      {
        if(doSerialEx){
          while(over_2){
            ex<-which.max(invader$H)
            HH<-HH[-ex, , drop=F]        
            invader$H<-invader$H[-ex]
            over_2<-max(invader$H)>=log(K/Mmin)
          }
        }else{
          temp_H<-invader$H
          while(over_2){
            ex<-which.max(temp_H)
            dd<-sample(1:dim(HH)[1],1)
            HH<-HH[-dd, , drop=F]
            
            invader$H<-invader$H[-dd]
            temp_H<-temp_H[-ex]
            HH_temp_g <- HH_temp_g[-ex, , drop=F]
            over_2<-max(temp_H)>=log(K/Mmin)
          }
        }
      }
      
      # Check to see if consumer can overcome respiration
      deaddd<-0
      if(sum(invader$H)<1){
        deaddd<-1
        par_col_sum_yo<-parColSums(HH)
        par_col_sum_yo2<-parColSums(HH^2)
        t_inv_HH<-t(invader$H)%*%HH
        new_inv[c(3,4)]<-c(num_inv,3)
        list_species[(num_dead+1),]<-new_inv
        num_dead<-num_dead+1
        still<-0
      }else{
        
        # Decide if any (random or normal) serial extinction will occur
        if(doSerialEx_bin){
          if(doSerialEx){
            #Perform serial extinction of resource test on new invader
            exR <- resourceEx1(invader$H)[[1]]
            while(sum(exR)){
              ex<-exR
              HH<-HH[-ex, , drop=F]
              invader$H<-invader$H[-ex]
              if(Resp_ser_ext){ # Remove consumers that cannot overcome respiration.
                exC1 <- parColSums(HH) < 1  
                if(any(exC1)){
                  only_ans[which(exC1),c(3,4)]<-matrix(rep(c(num_inv,2),sum(exC1)),ncol=2,byrow=T)
                  list_species[(num_dead+1:sum(exC1)),]<-only_ans[which(exC1),]
                  HH <- HH[,which(!exC1), drop=F];Sc <- length(only_ans[,1]); 
                  only_ans<-only_ans[which(!exC1),]
                  if(!doSerialEx){ HH_temp_g <- HH_temp_g[,which(!exC1), drop=F]}
                  num_dead<-num_dead+sum(exC1)
                }
              }     
              exR <- resourceEx1(invader$H)
            }  
            
          }else{
            # Perform serial extinction by randomisng which resources actually go extinct.
            temp_H<-invader$H
            exR <- resourceEx1(temp_H)[[1]]
            while(sum(exR)){
              ex_r <- sample(1:exR,1)
              ex<-exR
              HH<-HH[-ex_r, , drop=F]
              chosen_links <- chosen_links[-ex, , drop=F]
              invader$H<-invader$H[-ex_r]
              temp_H<-temp_H[-ex]
              HH_temp_g<-HH_temp_g[-ex,]
              
              if(Resp_ser_ext){
                exC1 <- parColSums(HH) < 1 # Remove consumers that cannot overcome respiration.
                if(any(exC1)){
                  only_ans[which(exC1),c(3,4)]<-matrix(rep(c(num_inv,2),sum(exC1)),ncol=2,byrow=T)
                  list_species[(num_dead+1:sum(exC1)),]<-only_ans[which(exC1),]
                  HH <- HH[,which(!exC1), drop=F];Sc <- length(only_ans[,1]); 
                  only_ans<-only_ans[which(!exC1),]
                  if(!doSerialEx){ HH_temp_g <- HH_temp_g[,which(!exC1), drop=F]}
                  num_dead<-num_dead+sum(exC1)
                }
              } 
              exR <- resourceEx1(temp_H)
            }  
          }
        }
        par_col_sum_yo<-parColSums(HH)
        par_col_sum_yo2<-parColSums(HH^2)
        t_inv_HH<-t(invader$H)%*%HH
      }
      # Exit overexploitation check module
    }
    
    # Remove consumers that cannot overcome respiration.
    if(Resp_ser_ext){
      exC1 <- parColSums(HH) < 1  
      if(any(exC1)){
        only_ans[which(exC1),c(3,4)]<-matrix(rep(c(num_inv,2),sum(exC1)),ncol=2,byrow=T)
        list_species[(num_dead+1:sum(exC1)),]<-only_ans[which(exC1),]
        HH <- HH[,which(!exC1), drop=F];Sc <- length(only_ans[,1]); 
        only_ans<-only_ans[which(!exC1),]
        if(!doSerialEx){ HH_temp_g <- HH_temp_g[,which(!exC1), drop=F]}
        num_dead<-num_dead+sum(exC1)
      }
    } 
    
    # Still represents whether a consumer has gone that has just invaded is still alive after
    # checking for overexploitation. Hence we will add it if to "only_ans" matrix.
    if(still){
      HH<-cbind(HH,invader$H)
      if(!doSerialEx){HH_temp_g<-cbind(HH_temp_g,invader$H)}
      only_ans<-rbind(only_ans,new_inv)
    }
  }
  }
  }else{       # Producer invasion #
    # cat("\r","prod",i,format(i/steps*100,width=5),sprintf("%#7g ",history,counter),a1_first,a1_old,a1)
    who<-runif(dim(HH)[1])<=disp_prob
    if(sum(who)>0){
        for (p in 1:sum(who)){
      
      n_plants_added<-n_plants_added+1
      Sc <- length(only_ans[,1])
      if(dim(HH)[2]!=Sc){stop("Sc")}
      
      #  Sample new trophic links for all consumers
      new<-exp(rnorm(Sc,sd = sigma))*only_ans[,1]*constanttt 
      HH<-rbind(HH,new)
      if(!doSerialEx){HH_temp_g<-rbind(HH_temp_g,new)}
      caused_Sr<-rep(0,Sc)
      
      if(doSerialEx_bin){ # Binary optionn for activating serial extinction.
        
        if(doSerialEx){
          # Perform normal serial extinction given the new trophic links.
          res_ext_P_index<-resourceEx(HH)
          while(sum(res_ext_P_index[[1]])){
            if(length(res_ext_P_index[[1]])>1){ex<-sample(res_ext_P_index[[1]],1)}else{ex<-res_ext_P_index[[1]]
            }
            HH<-HH[-ex, , drop=F]
            exC1 <- parColSums(HH) < 1  ## Delete consumers that cant overcome respiration
            if(any(exC1)){
              only_ans[which(exC1),c(3,4)]<-matrix(rep(c(num_inv,2),sum(exC1)),ncol=2,byrow=T)
              list_species[(num_dead+1:sum(exC1)),]<-only_ans[which(exC1),]
              HH <- HH[,which(!exC1), drop=F];Sc <- length(only_ans[,1]); 
              only_ans<-only_ans[which(!exC1),]
              if(!doSerialEx){ HH_temp_g <- HH_temp_g[,which(!exC1), drop=F]}
              num_dead<-num_dead+sum(exC1)
            }
            res_ext_P_index<-resourceEx(HH)
          }
          
        }else{ 
          # Perform random serial extinction given the new trophic links by randonly choosing extinct resources.
          res_ext_P_index<-resourceEx(HH_temp_g)
          while(sum(res_ext_P_index[[1]])){
            if(length(res_ext_P_index[[1]])>1){ex<-sample(res_ext_P_index[[1]],1)}else{ex<-res_ext_P_index[[1]]
            }
            ex_r<-sample(1:dim(HH)[1],1)
            HH<-HH[-ex_r, , drop=F]
            HH_temp_g<-HH_temp_g[-ex, , drop=F]
            exC1 <- parColSums(HH) < 1  ## Delete consumers that cant overcome respiration
            if(any(exC1)){
              only_ans[which(exC1),c(3,4)]<-matrix(rep(c(num_inv,2),sum(exC1)),ncol=2,byrow=T)
              list_species[(num_dead+1:sum(exC1)),]<-only_ans[which(exC1),]
              HH <- HH[,which(!exC1), drop=F];Sc <- length(only_ans[,1]); 
              only_ans<-only_ans[which(!exC1),]
              if(!doSerialEx){ HH_temp_g <- HH_temp_g[,which(!exC1), drop=F]}
              num_dead<-num_dead+sum(exC1)
            }
            res_ext_P_index<-resourceEx(HH_temp_g)
          }
        }
      }
      }
    }
  }
  
  # Summary stats  
  
  {
    par_col_sum_yo <- parColSums(HH)
    par_col_sum_yo2 <- parColSums(HH^2)
    t_inv_HHH <- mat.mult(t(HH), HH)
    ine_all <-  (t(((par_col_sum_yo -1) / t_inv_HHH))) / ( (par_col_sum_yo -1)/par_col_sum_yo2 )
    
    ine <- do.call(cbind,llply(1:dim(HH)[2],function(x){
      return(ine_all[-x,x])
    }))
    
    Sr_lost <- Sr-dim(HH)[1] + 1
    
    if(i%%cut_gran_com == 0 & i > exp_burn_in ){
      
      who_save <- sample(1:dim(HH)[2])
      repeat{
        # For each resident sample one corresponding comp. val. at random
        # Focal
        who_save_2 <- sample(1:dim(HH)[2])#(1:dim(HH)[2])[-who_save] 
        if(sum(who_save == who_save_2)==0){break}
      }
      
      
      write.table(matrix(c(unlist(only_ans[who_save,c(1,2)]),
                           par_col_sum_yo[who_save],
                           unlist(llply(1:length(who_save_2),function(x) t_inv_HHH[who_save_2[x],who_save[x]])),
                           par_col_sum_yo[who_save_2],par_col_sum_yo2[who_save_2]),
                         ncol=6,byrow = F),
                  file=paste("gran_table_com_",i,".txt",sep=""),col.names =F,row.names = F)
    }
  }
  
  # Check for exploitative competition  
  {
    # {
    #   survive_dir <- rep(T,dim(HH)[2])
    #   if(doComp_dir){
    #     # dead<-colAny(outcompetition_directly(HH,par_col_sum_yo,par_col_sum_yo2))
    #     ine_all_T_F <- ine_all < 1; diag(ine_all_T_F) <- F;
    #     dead<-colAny(ine_all_T_F)
    #     if(sum(dead)>0){
    #       survive_dir <- !dead
    #       only_ans[which(dead),c(3,4)]<-matrix(rep(c(num_inv,62),sum(dead)),ncol=2,byrow=T)
    #       list_species[(num_dead+1:sum(dead)),]<-only_ans[which(dead),]
    #       num_dead<-num_dead+sum(dead)
    #     }
    #   }
    #   
    #   dead_final3<-which(!survive_dir)
    #   if(length(dead_final3)>0){
    #     only_ans<-only_ans[-dead_final3,]
    #     HH <- HH[,-dead_final3, drop=F]
    #     if(!doSerialEx){ HH_temp_g <- HH_temp_g[,-dead_final3, drop=F]}
    #   }
    # }
    
    if(doComp_dir){
      par_col_sum_yo <- parColSums(HH)
      par_col_sum_yo2 <- parColSums(HH^2)
      dead<-colAny(outcompetition_directly(HH,par_col_sum_yo,par_col_sum_yo2))
      if(sum(dead)>0){
        only_ans[which(dead),c(3,4)]<-matrix(rep(c(num_inv,62),sum(dead)),ncol=2,byrow=T)
        list_species[(num_dead+1:sum(dead)),]<-only_ans[which(dead),]
        dead_final3<-which(dead)
        num_dead<-num_dead+sum(dead)
        only_ans<-only_ans[-dead_final3,]
        HH <- HH[,-dead_final3, drop=F]
        if(!doSerialEx){ HH_temp_g <- HH_temp_g[,-dead_final3, drop=F]}
      }
    }
  }
  
  # Check for Pyrrhic competition
  {
    if(doSerialEx){
      cons_class<-which(outcompetition(HH))
      if(specialist_check){
        if(length(cons_class)>1){
          mainResource_rep <- parColWhichMax(HH[,cons_class])
        }else{
          mainResource_rep <- which.max(HH[,cons_class])
        }
        if(length(mainResource_rep)>0 ){
          HH<-HH[-mainResource_rep,]
        }
      }
    }else{
      cons_class<-which(outcompetition(HH))
      if(specialist_check){
        if(length(cons_class)>1){
          mainResource_rep <- parColWhichMax(HH[,cons_class])
        }else{
          mainResource_rep <- which.max(HH[,cons_class])
        }
        if(length(mainResource_rep)>0 ){
          HH<-HH[-sample(1:dim(HH)[1],length(unique(mainResource_rep))),]
          HH_temp_g<-HH_temp_g[-mainResource_rep,]
        }
      }
    }
    
    exC1 <- parColSums(HH) < 1  ## consumers that fell below threshold
    
    if(any(exC1)){
      only_ans[which(exC1),c(3,4)]<-matrix(rep(c(num_inv,1),sum(exC1)),ncol=2,byrow=T)
      list_species[(num_dead+1:sum(exC1)),]<-only_ans[which(exC1),]
      HH <- HH[,which(!exC1), drop=F];Sc <- length(only_ans[,1]); 
      only_ans<-only_ans[which(!exC1),]
      if(!doSerialEx){ HH_temp_g <- HH_temp_g[,which(!exC1), drop=F]}
      num_dead<-num_dead+sum(exC1)
    }
    
    n_C<-dim(HH)[2]
    n_P<-dim(HH)[1]
  }
  
  ## The rest is book keeping...
  aggs<-only_ans[,1]
  Sc <- length(aggs)
  history <- c(exp(mean(log(aggs))),min(aggs),max(aggs),1-mean(parColSums(HH^2)/parColSums(HH)^2),Sc)
  
  # if(i %% plot_gran == 0){
  # cat("\r",i,format(i/steps*100,width=5),sprintf("%#7g ",history,counter),a1_first,a1_old,a1)
  # points(i,dim(HH)[2],col="red",cex=0.1,ylim=c(0,200),xlim=c(0,steps))
  # points(i,dim(HH)[1],col="green",cex=0.1,ylim=c(0,200),xlim=c(0,steps))
  # }
  
  if(i %% cut_gran == 0){
    write.csv(list_species[which(list_species$agg!=0),],file=paste("list_species",i,".csv",sep=""))
    save.image("deconstructed_PDMM_end.RData")
    # write.table(gran_table_ext,file=paste("gran_table_ext_",i+10000,".txt",sep=""),col.names =F,row.names = F)
    # 
    num_dead<-0
    {
      list_species<-as.data.frame(list(
        agg=c(rep(0,total_species)),#1
        id=c(rep(0,total_species)),#2
        outi=c(rep(-1,total_species)),#3
        extinction=c(rep(-1,total_species)),#4
        parent_id=c(rep(-1,total_species)), #5
        n_C=c(rep(-1,total_species)), #6
        n_P=c(rep(-1,total_species)), #7
        m_l_a=c(rep(-1,total_species)),#8
        v_l_a=c(rep(-1,total_species)),#9
        inv_counter=c(rep(-1,total_species))#10
        # A1_m_ja=rep(-1,total_species),#11
        # A2_m_ja=rep(-1,total_species),#12
        # B1_m_ja=rep(-1,total_species),#13
        # B2_m_ja=rep(-1,total_species),#14
        # ine_M_ja=rep(-1,total_species),#15
        # A1_m_jb=rep(-1,total_species),#16
        # A2_m_jb=rep(-1,total_species),#17
        # B1_m_jb=rep(-1,total_species),#18
        # B2_m_jb=rep(-1,total_species),#119
        # ine_M_jb=rep(-1,total_species),
        # max_ja=rep(-1,total_species),
        # max_jb=rep(-1,total_species),
        # DD_ja=rep(-1,total_species),
        # DD_jb=rep(-1,total_species)#20
      ))
    }
    
    if(i%% (cut_gran*2) == 0){
      system("cat list_species* >> all_species.txt && rm list_species*") #cat other_stats* >> all_stats.txt && rm other_stats* &&
      # system("cat i_l_RAW_min_* >> all_i_l_RAW_min.txt && rm  i_l_RAW_min_* ") #&&
      # system("cat i_l_RAW_* >> all_i_l_RAW.txt && rm  i_l_RAW_* ") #&&
      # system("cat l_m_MEANS_MINS_P_* >> all_l_m_MEANS_MINS_P.txt && rm  l_m_MEANS_MINS_P_* ") #&&
      # system("cat l_m_MEANS_MINS_C_* >> all_l_m_MEANS_MINS_C.txt && rm  l_m_MEANS_MINS_C_* ") #&&
      # system("cat gran_table_com_* >> all_gran_table_com.txt && rm  gran_table_com_* ") #&&
      # system("cat min_table_com_* >> all_min_table_com.txt && rm  min_table_com_* ") #&&
      # system("cat gran_table_ext_* >> all_gran_table_ext.txt && rm  gran_table_ext_*") #&&
      # cat i_k_RAW_* >> all_i_k_RAW.txt && rm  i_k_RAW_* && 
    }
  }
  
  # The following is important for subsequent analysis
  if(i%% (cut_gran_com*2) == 0 & i > exp_burn_in ){
    system("cat gran_table_com_* >> all_gran_table_com.txt && rm  gran_table_com_* ")
  }
  
  if(i%% (cut_gran_HH) == 0 & i > exp_burn_in ){
    write(only_ans[,1],paste("agg",i,".txt",sep=""),sep="\n")
    write.table(as.matrix(HH),file=paste("HH",i,".txt",sep=""),col.names = F,row.names = F)
  }
}

system("cat list_species* >> all_species.txt && rm list_species*") #cat other_stats* >> all_stats.txt && rm other_stats* &&
system("cat gran_table_com_* >> all_gran_table_com.txt && rm  gran_table_com_* ")
# save.image("deconstructed_PDMM_end.RData")




