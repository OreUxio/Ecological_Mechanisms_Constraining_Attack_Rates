# Metacommunity.R - Simulates the approximated metacommunity assembly model analysed in Chapter 4.
# 
# - Requires: ../../Data/mu_sd_num_dt.RData, ../../Data/Pre_Processed/Data_deconstructed.RData and ../../Data/95_optimised_K.RData to run
# 
# - Saves simulation data in ../Data/Metacommunity/New/ as *_*_meta.txt

argzzz = commandArgs(trailingOnly=TRUE)
library(Rfast)
library(plyr)
library(zoo)
library(survival)
library(maxLik)
library(plyr)
library(ggplot2)
library(deSolve)
library(RColorBrewer)

load("../../Data/Pre_Processed/Data_deconstructed.RData")
xx<-with(species_fitness_equilibrium,rollmean(sort(log10(agg)),k=bin))
yy<-with(species_fitness_equilibrium,1/rollmean(lifetime[order(agg)],k=bin))
nlmod <- nls(yy ~  Const + A * xx + B * xx^2,start=c(Const=1,A=1,B=1)) #+ C * xx^3 ,C=1
coefs_death<-coef(nlmod)
xx<-with(for_nls_birth,rollmean(sort(log10(agg)),k=bin))
yy<-with(for_nls_birth,rollmean(birth_rate[order(agg)],k=bin))
nlmod <- nls(yy ~  Const + A * xx + B * xx^2,start=c(Const=1,A=1,B=1)) #+ C * xx^3 ,C=1
coefs_birth<-coef(nlmod)
m_ratio<-mean(all_species_keep$n_C/all_species_keep$n_P)
rm(list=c(setdiff(ls(),c("coefs_death","coefs_birth","argzzz","m_ratio"))))

load("../../Data/Pre_Processed/mu_sd_num_dt.RData")
load("../../Data/Pre_Processed/95_optimised_K.RData")
args<-argzzz

args<-c(0,0,0.02,0.8)

# Functions ####

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

pInv <- function(l_agg){
  
  pass<-(log(1/(10^l_agg)/K_dt)-as.numeric(mu))/(as.numeric(sd))
  
  return(1-pnorm(pass,0,1))
}

pBirth <- function(l_agg){
  coefs_birth[1]+coefs_birth[2]*l_agg+coefs_birth[3]*l_agg^2
}

pDeath <- function(l_agg){
  coefs_death[1]+coefs_death[2]*l_agg+coefs_death[3]*(l_agg^2)
}

samplee <- function(x,mu,var){
  return(mu*(var^rnorm(1))*x)
} 

# Parameters ####

Do_intra_exp_comp<-T
Random_disp_struct<-T 
  
Sc<-250*m_ratio
min_Sr<-150
max_Sr<-500
mu_sd<-find_mu_sd(250,4)
mu<-unlist(mu_sd[1]);sd<-unlist(mu_sd[2]);
sqrt_sites <- 100
nSites <- sqrt_sites^2
sLogAgg <- log10(1.3^0.5);sAgg<-10^sLogAgg
dLogAgg <- log10(0.8^0.5);
dAgg<-10^dLogAgg
disp_prob <- 1/(Sc)*0.5
sLogDisp <- as.numeric(args[3]) 
sDisp<- (10^sLogDisp) 
dLogDisp <- dLogAgg 
dDisp <- as.numeric(args[4])
initialAgg <- -5
sim_time<-round(5e+05*0.06)
savee<-T
percz<-seq(0.5,10,length.out = 40)

percz_og<-seq(0.5,10,length.out = 20)
percz<-seq(percz_og[6],percz_og[7],length.out = 20)

all_filez<-list.files(pattern="TRUE_meta")
xxx<-unlist(llply(all_filez, function(x) {
  as.numeric(unlist(strsplit(x,"_"))[1])
} ))

# Simulation ####
for (disp_prob in xxx){
  start_point<-sqrt_sites*(sqrt_sites/2-1)+sqrt_sites/2
  how_many<- round((sqrt_sites^2)*0.1)
  occupied <- rep(F,nSites);occupied[sample(1:(sqrt_sites^2),how_many)]<-T
  agg <- rnorm(nSites,initialAgg,sLogAgg)
  D_ind <- rnorm(nSites,disp_prob,disp_prob*0.1)
  
  ind_val <- round(abs(agg)*occupied/3*100+1)-100+1
  ind_val[which(ind_val<0)]<-1
  ind_val_matrix<-matrix(ind_val,ncol=sqrt_sites)
  indexx<-0
  index_m<-0
  save_gran<-100
  occupancy_m<-rep(0,sim_time/save_gran)
  agg_m<-rep(0,sim_time/save_gran)
  diff<-c()
{
  cat("Disp P","% occ","\t","% sim comp","\t","mean log_10(Agg)","\t","mean(Disp)","                \n")
  while(indexx <= sim_time){
    if(indexx%%save_gran==0){
      index_m<-index_m+1
      occupancy_m[index_m]<-mean(occupied)*100
      agg_m[index_m]<-mean(agg[occupied])
    }
    indexx <- indexx+1
    
    # One type per patch may attempt disperal or go extinct 
    disp<-runif(length(which(occupied)),0,1) < disp_prob 
    
    if(any(disp)){
      from <- which(occupied)[which(disp)]
      # Local dispersal
      sum_invs<-0
      for (fromm in from){
        
        if(fromm %% sqrt_sites == 0 ){
          to <- sample(abs(fromm+c(-1,-sqrt_sites,sqrt_sites,-(sqrt_sites+1),(sqrt_sites-1))),1)
        }else if(fromm %% sqrt_sites == 1 ){
          to <- sample(abs(fromm+c(1,-sqrt_sites,sqrt_sites,-(sqrt_sites+1),-(sqrt_sites-1))),1)
        }else{
          to <- sample(abs(fromm+c(1,-1,-sqrt_sites,sqrt_sites,-(sqrt_sites+1),sqrt_sites+1,-(sqrt_sites-1),
                                   sqrt_sites-1)),1)
        }
        if(to>nSites){
          to<-sample(c(nSites-1,nSites-sqrt_sites,nSites-sqrt_sites-1),1)
        }else if(to==0){
          to<-sample(c(2,11,1),1)
        }
        to <- as.numeric(!Random_disp_struct) * to + as.numeric(Random_disp_struct) * round(runif(1,1,nSites))
        newLogAgg <- log10(samplee(10^agg[fromm],dAgg,sAgg)) 
        if(length(occupied[to])>1){stop("ey")}
        if(occupied[to]){
          if(newLogAgg > agg[to] & Do_intra_exp_comp){
            agg[to] <- newLogAgg
            sum_invs<-sum_invs+1
          }
        }else{
          if(runif(1) < pInv(newLogAgg)){
            occupied[to] <- TRUE
            agg[to] <- newLogAgg
            sum_invs<-sum_invs+1
          }
        }
      }
    }
    
    dead<-runif(length(which(occupied)),0,1)<pDeath(agg[which(occupied)])
    if(any(dead)){
      # It's a death
      occupied[which(occupied)[which(dead)]]<-FALSE
    }
    if(all(!occupied)){ 
      if(savee){
        write.table(cbind(occupancy_m,agg_m),paste("Data/Metacommunity/New/",disp_prob,"_",Do_intra_exp_comp,"_meta.txt",sep=""))
      }
      break
    }
    cat(disp_prob,"\t",signif(mean(occupied)*100,2),"\t",signif(indexx/sim_time*100,2),"\t","\t",signif(mean(agg[occupied]),4),"\t","\t",signif(mean(D_ind[occupied]),4),"                \r")
    }
  }
  if(savee){
    write.table(cbind(occupancy_m,agg_m),paste("Data/Metacommunity/New/",disp_prob,"_",Do_intra_exp_comp,"_meta.txt",sep=""))
  }
}



