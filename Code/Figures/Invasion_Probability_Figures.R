library(Rfast)
library(plyr)
library(zoo)
library(survival)
library(maxLik)
library(plyr)
library(ggplot2)
library(deSolve)
library(RColorBrewer)

################################
#        Functions             #
################################

outcompeted <- function(H,t_H_HH,HHSums,HH2Sums){
  directly <-
    any( (HHSums-1) * t_H_HH >
           (sum(H)-1) * HH2Sums )
  if(directly) return(T)
  return(directly) #directly |
}

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

gen_linnks <- function(aggs){
    parent<-sample(aggs,1)
    agg <- parent #(0.8^0.5)*(1.3^0.5)^rnorm(1)*
    H <- agg*exp(rnorm(Sr,0,sigma))*constanttt
  return(list(H=H,agg=agg,parent_number=which(aggs==parent)))
}

gen_linnks_mut <- function(aggs){
  parent<-sample(aggs,1)
  agg <- a1*(a2^rnorm(1,0,1))*parent
  H <- agg*exp(rnorm(Sr,0,sigma))*constanttt
  return(list(H=H,agg=agg,parent_number=which(aggs==parent)))
}

find_mu_sd<-function(s_r,sigmaa){
  #Function for estimating the expected value of the strongest link strength, i.e. E(max(exp(sigmaa*xi_{max,i})))
  
  n_samples<-10000#Number of samples
  sum_m<-rep(0,n_samples)
  for (i in 1:n_samples){
    sum_m[i]<-sum(exp(sigmaa*rnorm(s_r,0,1)))
  }
  mean_sum_m<-mean(log(sum_m))
  sd_sum_m<-(var(log(sum_m)))^0.5
  out<-list(mean_sum_m,sd_sum_m)
  return(out)
}

invasion_prob_c<-function(l_a_m,mu,sd){
  pass<-(log(resp/(epsilon*exp(l_a_m)*K*P_b_m))-as.numeric(mu))/(as.numeric(sd))
  return(1-pnorm(pass,0,1))
}

invasion_prob_K<-function(l_a_m,mu,sd,K){
  pass<-(log(resp/(epsilon*exp(l_a_m)*K))-as.numeric(mu))/(as.numeric(sd))
  return(1-pnorm(pass,0,1))
}


birth_prob_c<-function(l_a_rr,mu,sd){
  # return(invasion_prob_c(l_a_rr+m_a,mu,sd))
  return(integrate(joint_prob_c,lower=-Inf,upper=Inf,l_a_rr=l_a_rr,mu=mu,sd=sd)$value)
}

birth_prob_K<-function(l_a_rr,mu,sd,K){
  return(invasion_prob_K(l_a_rr+m_a,mu,sd,K))
  # return(integrate(joint_prob_c,lower=-10,upper=10,l_a_rr=l_a_rr,mu=mu,sd=sd)$value)
}

joint_prob_c<-function(l_a_mm,l_a_rr,mu,sd){
  return(mut_prob_given_ar(l_a_mm,l_a_rr)*invasion_prob_c(l_a_mm,mu,sd))
}

mut_prob_given_ar<-function(l_a_m,l_a_r){
  mu_a<-m_a+l_a_r
  return(dnorm(l_a_m,mu_a,s_a))
}


birth_prob_optim<-function(par){
  mu<-par[1];sd<-par[2]
  return(sum((invasion_prob_c(log(aggz)+m_a,mu,sd)-alll_mut_m)^2))
}

invasion_prob_optim<-function(par){
  mu<-par[1];sd<-par[2]
  return(sum((invasion_prob_c(log(aggz),mu,sd)-alll_m)^2))
}

#############################
#        Parameters         #
#############################

n_points<-15 # Number of fixed base attack rate values
aggz<-10^c(seq(-6.25,-4,length.out = n_points)) # Range of base attack rate values
a1<-0.8^0.5
a2<-1.3^0.5
m_a<-log(a1)
s_a<-log(a2)
sigma<-4; #Inverse niche width
resp<-0.1; # Respiration
K<-1;
epsilon<-0.1;
constanttt<-K*epsilon/resp
doEst<-F
lwd=lwdd<-2

##############################################################
#     Invasion probability of a mutant once selected         #
##############################################################

# Sr variable loaded with this command too
load("Data/Pre_Procossed/invasion_prob_full.RData") # outputed from invasion_prob_full.R
Sr<-250
mu_sd<-find_mu_sd(Sr,sigma) 
alll<-data.frame()
nreps<-10000
for (i in 1:nreps){
  H_list<-llply(aggz, function(x) gen_linnks(x)$H)
  inv_list<-llply(H_list, function(x) sum(x)>1)
  alll<-rbind(alll,unlist(inv_list))
}

alll_2<-data.frame()
nreps<-10000
correction<-median(read.table("Data/Raw/Full_assembly_model/list_of_plant_biomass_values.dat")[,3])
for (i in 1:nreps){
  H_list<-llply(aggz, function(x) gen_linnks(x)$H)
  inv_list<-llply(H_list, function(x) sum(x)>1/correction)
  alll_2<-rbind(alll_2,unlist(inv_list))
}
alll_m<-colMeans(alll)
alll_2_m<-colMeans(alll_2)

mult<-1.75
comm_cex_lab<-2
leg.cex<-1
# Figure 2.4 #
{
  
  pdf("/Images/Invasion_verification.pdf",height=10*0.5*mult,width=10*0.5*mult)
  par(mar=c(4,5,1,1))
  plot(log10(aggz),alll_m,cex=0.1,type="l",xlab=expression(log[10](a[m])),ylab=expression(I(a[m])),
       cex.lab=comm_cex_lab,cex.axis=comm_cex_lab)
  with(mean_invasion_fitness_per_web,lines(log10(agg),sim_prob,lty=2))
  lines(log10(aggz),alll_2_m,lty=1,col="blue")
  legend("bottomright",legend=c("Approximation","Approximation (corrected)","Numeric estimate"),
         lty=c(1,2),col=c("black","black"),bty="n",cex=leg.cex,lwd=c(2),seg.len = 4)
  dev.off()
}

# Figure 3.2 #
{
  pdf("/Images/Invasion_Probabillity.pdf",
      height=10*0.5*mult,width=10*0.5*mult)
  par(mar=c(4,5,1,1))
  plot(log10(aggz),alll_m,cex=0.1,type="l",xlab=expression(log[10](a[m])),ylab=expression(I(a[m])),
       cex.lab=comm_cex_lab,cex.axis=comm_cex_lab,col="blue")
  with(mean_invasion_fitness_per_web,lines(log10(agg),sim_prob,lty=2,col="blue"))
  
  legend("bottomright",legend=c("Numeric estimate","Deconstructed model","Full model",
                                "Numeric estimate of mean and variance",
                                "Analytic estimate of mean and variance","Optimised estimate of mean and variance"),
         lty=c(1,1,2,1,1,1),col=c("blue","black","black","red","orange","green"),bty="n",cex=leg.cex,lwd=c(2),seg.len = 4)
  # Prediction using numeric estimate of mu and sd
  P_b_m<-1
  prob_pred<-llply(aggz,function(x) invasion_prob_c(log(x),mu_sd[[1]],mu_sd[[2]]))
  lines(log10(aggz),unlist(prob_pred),cex=0.1,type="l",col="red")
  # Prediction using analytic estimate of mu and sd
  g_mean_star<-sigma*sqrt(2*log(Sr))
  g_sd_star<-sqrt(((sigma^2)/sqrt(sigma^2+g_mean_star^2)))#1/sqrt(1+2*log(dim(HH)[1]))
  prob_pred_anal<-llply(aggz,function(x) invasion_prob_c(log(x),g_mean_star,g_sd_star))
  lines(log10(aggz),unlist(prob_pred_anal),cex=0.1,type="l",col="orange")
  # Prediction using optimised mu and sd
  mu_sd_2<-optim(par=c(mu=mu_sd[[1]],sd=mu_sd[[2]]),fn=invasion_prob_optim)$par
  prob_pred_opt<-llply(aggz,function(x) invasion_prob_c(log(x),mu_sd_2[1],mu_sd_2[2]))
  lines(log10(aggz),unlist(prob_pred_opt),cex=0.1,type="l",col="green")
  P_b_m<-mean(read.table("../Full/list_of_plant_biomass_values.dat")[,3])
  prob_pred_opt<-llply(aggz,function(x) invasion_prob_c(log(x),mu_sd_2[1],mu_sd_2[2]))
  lines(log10(aggz),unlist(prob_pred_opt),cex=0.1,type="l",col="green",lty=2)
  prob_pred<-llply(aggz,function(x) invasion_prob_c(log(x),mu_sd[[1]],mu_sd[[2]]))
  lines(log10(aggz),unlist(prob_pred),cex=0.1,type="l",col="red",lty=2)
  dev.off()
}

#############################################################
#     Birth probability of a resident once selected         #
#############################################################

alll_mut<-data.frame()
for (i in 1:nreps){
  H_list<-llply(aggz, function(x) gen_linnks_mut(x)$H)
  inv_list<-llply(H_list, function(x) sum(x)>1)
  alll_mut<-rbind(alll_mut,unlist(inv_list))
}

birth_prob_c<-function(l_a_rr,mu,sd){
  return(invasion_prob_c(l_a_rr+m_a,mu,sd))
  # return(integrate(joint_prob_c,lower=-15,upper=-1,l_a_rr=l_a_rr,mu=mu,sd=sd)$value)
}

joint_prob_c<-function(l_a_mm,l_a_rr,mu,sd){
  return(mut_prob_given_ar(l_a_mm,l_a_rr)*invasion_prob_c(l_a_mm,mu,sd))
}

mut_prob_given_ar<-function(l_a_m,l_a_r){
  mu_a<-m_a+l_a_r
  return(dnorm(l_a_m,mu_a,s_a))
}



alll_mut_m<-colMeans(alll_mut)
#
pdf("/Images/Birth_rate_selected.pdf",height=10*0.5*mult,width=10*0.5*mult)
par(mar=c(4,5,1,1))
plot(log10(aggz),alll_mut_m,cex=0.1,type="l",xlab=expression(log[10](a[m])),ylab=expression(b[S](a[r])),
     cex.lab=comm_cex_lab,cex.axis=comm_cex_lab)

# Prediction using numeric estimate of mu and sd
mu_sd<-find_mu_sd(Sr,sigma)
prob_pred<-llply(aggz,function(x) birth_prob_c(log(x),mu_sd[[1]],mu_sd[[2]]))
lines(log10(aggz),unlist(prob_pred),cex=0.1,type="l",col="red")

# Prediction using analytic estimate of mu and sd
g_mean_star<-sigma*sqrt(2*log(Sr))
g_sd_star<-sqrt(((sigma^2)/sqrt(sigma^2+g_mean_star^2)))#1/sqrt(1+2*log(dim(HH)[1]))
prob_pred_anal<-llply(aggz,function(x) birth_prob_c(log(x),g_mean_star,g_sd_star))
lines(log10(aggz),unlist(prob_pred_anal),cex=0.1,type="l",col="blue")

# Prediction using optimised mu and sd
mu_sd_2<-optim(par=c(mu=mu_sd[[1]],sd=mu_sd[[2]]),fn=birth_prob_optim)$par
prob_pred_opt<-llply(aggz,function(x) birth_prob_c(log(x),mu_sd_2[1],mu_sd_2[2]))
lines(log10(aggz),unlist(prob_pred_opt),cex=0.1,type="l",col="green",lty=1)

legend("bottomright",legend=c("Numeric","Numeric estimate of mean and variance",
                              "Analytic estimate of mean and variance","Optimised estimate of mean and variance" ), pch=c(NA, NA), lty=c(1), 
       col=c("black","red","blue","green"),bty="n",cex=leg.cex,lwd=c(2),seg.len = 4)
dev.off()

############################################
#     Birth rate of a resident      #
############################################

# Figure 3.4 #
{
alll_mut<-data.frame()
for (i in 1:nreps){
  aggzz<-sample(aggz,10,replace = T)
  H_list<-llply(aggzz, function(x) gen_linnks(x)$H)
  inv_list<-llply(H_list, function(x) sum(x)>1)
  alll_mut<-rbind(alll_mut,unlist(inv_list)*aggzz)
}

pdf("Images/Birth_rate_no_comp.pdf",height=10*0.5*mult,width=10*0.5*mult)

mu<-unlist(mu_sd[1]);sd<-unlist(mu_sd[2])
here<-unlist(as.list(alll_mut));here<-here[which(here!=0)]
probss<-unlist(llply(aggz,function(x) length(which(here==x))/length(here)))
P_b_m<-1
probs_b<-unlist(llply(aggz, function(x) birth_prob_c(log(x),mu,sd)))
probs_f<-unlist(llply(1:length(aggz),function(x) probs_b[x]/sum(probs_b)))

g_mean_star<-sigma*sqrt(2*log(Sr))
g_sd_star<-sqrt(((sigma^2)/sqrt(sigma^2+g_mean_star^2)))#1/sqrt(1+2*log(dim(HH)[1]))
probs_b_anal<-unlist(llply(aggz, function(x) birth_prob_c(log(x),g_mean_star,g_sd_star)))
probs_anal<-unlist(llply(1:length(aggz),function(x) probs_b[x]/sum(probs_b_anal)))

lines(log10(aggz),unlist(prob_pred_anal),cex=0.1,type="l",col="blue")

probs_b_o<-unlist(llply(aggz, function(x) birth_prob_c(log(x),mu_sd_2[1],mu_sd_2[2])))
probs_f_o<-unlist(llply(1:length(aggz),function(x) probs_b_o[x]/sum(probs_b_o)))

par(mar=c(4,5,1,1))
plot(log10(aggz),log10(probss),xlab=expression(log[10](a[m])),ylab=expression(b(a[r])),
     cex.lab=comm_cex_lab,cex.axis=comm_cex_lab,type="l",lwd=lwdd)
lines(log10(aggz),log10(probs_f),col="red",lwd=lwdd)
# lines(log10(aggz),log10(probs_f_o),col="green")
lines(log10(aggz),log10(probs_b_anal),col="blue")
legend("bottomright",legend=c("Observed numeric estimate","Predicted"), 
       pch=c(NA), lty=c(1), 
       col=c("black","red"),bty="n",cex=leg.cex,lwd=c(2),seg.len = 4)
dev.off()
}


#############################################################
#     Birth probability of a resident with competition      #
#############################################################

# Figure 3.3 #
{
filez<-list.files(pattern = "HH",recursive = T,full.names = T,path="Data/Raw/Deconstructed_assembly_model")
part1<-strsplit(strsplit(filez[1],".txt")[[1]],"HH")[[1]]
HH<-as.matrix(read.table(paste(part1[[1]],"HH",as.integer(part1[[2]]),".txt",sep="")));
par_col_sum_yo<-parColSums(HH)
par_col_sum_yo2<-parColSums(HH^2)
Sr<-dim(HH)[1]
alll_mut_comp<-data.frame()
nreps<-100000
for (i in 1:nreps){
  H_list<-llply(aggz, function(x) gen_linnks_mut(x)$H)
  t_inv_HH_list<-llply( H_list, function(x) sum(t(x)%*%HH ))
  d_a_list <- llply(1:length(aggz), function(x) !outcompeted(unlist(H_list[x]),unlist(t_inv_HH_list[x]), par_col_sum_yo, par_col_sum_yo2) ) 
  inv_list<-llply(1:length(aggz), function(x) sum(unlist(H_list[x]))>1 & unlist(d_a_list[x]))
  alll_mut_comp<-rbind(alll_mut_comp,unlist(inv_list))
}

m_a<-log(0.8^0.5)
s_a<-log(1.3^0.5)
alll_mut_m_comp<-colMeans(alll_mut_comp)
mu_sd<-find_mu_sd(Sr,sigma)
mu<-unlist(mu_sd[1]);sd<-unlist(mu_sd[2])
probss<-alll_mut_m_comp
P_b_m<-3
probs_b<-unlist(llply(aggz, function(x) birth_prob_c(log(x),mu,sd)))
probs_f<-unlist(llply(1:length(aggz),function(x) probs_b[x]/sum(probs_b)))

g_mean_star<-sigma*sqrt(2*log(Sr))
g_sd_star<-sqrt(((sigma^2)/sqrt(sigma^2+g_mean_star^2)))#1/sqrt(1+2*log(dim(HH)[1]))
probs_b_anal<-unlist(llply(aggz, function(x) birth_prob_c(log(x),g_mean_star,g_sd_star)))
probs_anal<-unlist(llply(1:length(aggz),function(x) probs_b[x]/sum(probs_b_anal)))
probs_b_o<-unlist(llply(aggz, function(x) birth_prob_K(log(x),mu,sd,1.65)))
probs_f_o<-unlist(llply(1:length(aggz),function(x) probs_b_o[x]/sum(probs_b_o)))

pdf("Images/Birth_rate_comp.pdf",height=10*0.5*mult,width=10*0.5*mult)
comm_cex_lab<-1.3
par(mar=c(4,5,1,1))
plot(log10(aggz),log10(probs_f),col="red",xlab=expression(log[10](a[m])),ylab=expression(b(a[r])),
     cex.lab=comm_cex_lab,cex.axis=comm_cex_lab,type="l",lwd=lwdd)
# lines(log10(aggz),log10(probs_f_o),col="green")
lines(log10(aggz),log10(probss),lwd=lwdd)
legend("bottomright",legend=c("Observed numeric estimate","Predicted"), 
       pch=c(NA), lty=c(1), 
       col=c("black","red"),bty="n",cex=leg.cex,seg.len = 4,lwd=lwdd)
dev.off()
}

####################################
#     Optimisation procedures      #
####################################
birth_prob_optim<-function(par){
  mu<-par[1];sd<-par[2]
  return(sum((invasion_prob_c(log(aggz)+m_a,mu,sd)-alll_mut_m)^2))
}

if(doEst){
  #Problem is in optimising the parameters for the right number of resources.
  # Find optimal mu and sd for each set of resources.
  mu_sd<-find_mu_sd(275,sigma)
  mu_sd_dt<-data.frame()
  for (S in 50:800){
    cat(S,"\r")
    Sr<-S
    alll_mut<-data.frame()
    nreps<-1000
    for (i in 1:nreps){
      H_list<-llply(aggz, function(x) gen_linnks_mut(x)$H)
      inv_list<-llply(H_list, function(x) sum(x)>1)
      alll_mut<-rbind(alll_mut,unlist(inv_list))
    }
    alll_mut_m<-colMeans(alll_mut)
    mu_sd_dt<-rbind(mu_sd_dt,optim(par=c(mu=mu_sd[[1]],sd=mu_sd[[2]]),fn=birth_prob_optim)$par)
  }
  colnames(mu_sd_dt)<-c("mu","sd");mu_sd_dt$Sr<-50:800
  save.image("Data/Pre_Processed/optimised_mu_sd_2.RData")


birth_rate_optim<-function(par){
  K<-par[1]
  probs_b<-unlist(llply(aggz, function(x) birth_prob_K(log(x),mu,sd,K)))
  probs_f<-unlist(llply(1:length(aggz),function(x) probs_b[x]/sum(probs_b)))
  return(sum((probs_f-alll_mut_m_comp)^2))
}
  
  K_dt<-data.frame()
  filez<-list.files(pattern = "HH",recursive = T,full.names = T)
  for (j in 1:length(filez)){
    part1<-strsplit(strsplit(filez[j],".txt")[[1]],"HH")[[1]]
    HH<-as.matrix(read.table(paste(part1[[1]],"HH",as.integer(part1[[2]]),".txt",sep="")));
    par_col_sum_yo<-parColSums(HH)
    par_col_sum_yo2<-parColSums(HH^2)
    Sr<-dim(HH)[1]
    mu_sd<-find_mu_sd(Sr,sigma)
    mu<-unlist(mu_sd[1]);sd<-unlist(mu_sd[2])
    alll_mut_comp<-data.frame()
    nreps<-10000
    for (i in 1:nreps){
      H_list<-llply(aggz, function(x) gen_linnks_mut(x)$H)
      t_inv_HH_list<-llply( H_list, function(x) sum(t(x)%*%HH ))
      d_a_list <- llply(1:length(aggz), function(x) !outcompeted(unlist(H_list[x]),unlist(t_inv_HH_list[x]),par_col_sum_yo, par_col_sum_yo2) ) 
      inv_list<-llply(1:length(aggz), function(x) sum(unlist(H_list[x]))>1 & unlist(d_a_list[x]))
      alll_mut_comp<-rbind(alll_mut_comp,unlist(inv_list))
    }
    alll_mut_m_comp<-colMeans(alll_mut_comp)
    K_dt<-rbind(K_dt,optim(par=c(K=1),fn=birth_rate_optim,method = "Brent",lower=0,upper=1)$par)
    cat(j/length(filez),"% ",unlist(K_dt[,1][j]))
  }
  colnames(K_dt)<-c("K");
  mean(K_dt[,1])
  save.image("Data/Pre_Processed/optimised_K.RData")
}
