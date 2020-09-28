library(plyr)
load("/Data/Pre_Processed/Data_full.RData")
#Functions #######

invasion_prob<-function(l_a_m,mu,sd){
  pass<-(log(respiration/(the_conversion_efficiency_epsilon*exp(l_a_m)*K*0.5))-as.numeric(mu))/as.numeric(sd)
  return(1-pnorm(pass,0,1))
}

#Joint probability of yielding mutant (given resident) and that the mutant succesfuly invades
joint_prob<-function(l_a_mm,l_a_rr,mu,sd){
  return(dnorm(l_a_mm,m_a+l_a_rr,s_a)*invasion_prob(l_a_mm,mu,sd))
}

#Birth probability: intergrate joint probability over all possible mutants.
birth_prob<-function(l_a_rr,mu,sd){
  return(integrate(joint_prob,lower=-20,upper=0,l_a_rr=l_a_rr,mu=mu,sd=sd)$value)
}

find_mu_sd<-function(s_r,sigmaa){
  n_samples<-1000#Number of samples
  sum_m<-rep(0,n_samples)
  for (i in 1:n_samples){
    # sum_m[i]<-sum(head(sort(exp(sigmaa*rnorm(s_r,0,1))),-1))
    sum_m[i]<-sum(exp(sigmaa*rnorm(s_r,0,1)))
  }
  mean_sum_m<-mean(log(sum_m))
  sd_sum_m<-(var(log(sum_m)))^0.5
  out<-list(mean_sum_m,sd_sum_m)
  return(out)
}

##### options ####
run_with_singularity = F
invasion_prob_number_files<-100
replicates_per_web<-1000
agg_granularity<-20
####### Statistics for the invasion probability using speciate() function #### 
# Data preperation #####

webfiles<-sample(list.files(pattern = "web*[1-9][0-9][0-9]"),invasion_prob_number_files,path="Data/Raw/Full_assembly_model")
agg<-species_fitness_equilibrium$agg
web_numbs<-webfiles
test_agg<-10^seq(-6.25,-4,length.out=agg_granularity)
store<-data.frame(rep(test_agg,replicates_per_web*length(webfiles)));colnames(store)<-"agg";store<-data.frame(store[order(store$agg),]);colnames(store)<-"agg"
store$web<-rep(web_numbs,length(test_agg));store<-store[order(store$agg,store$web),];store$replicates<-rep(1:replicates_per_web,length(web_numbs)*length(test_agg));
store$fitness<-0;
seeds<-sample(1:10^4, size=replicates_per_web+2, replace = FALSE, prob = NULL)
total<-dim(store)[1]
system("rm before*")
if(run_with_singularity){
  command<-paste("singularity exec PDMM.img /Project_Code/C_++/NewWeb_new_new/build/Inspect -k random_seed=",seeds[(i%%replicates_per_web+1)],
                 " -W ",0.1," ",
                 store[1,]$web,sep="")
}else{
  command<-paste("/Project_Code/C_++/NewWeb_new_new/build/Inspect -k random_seed=",seeds[(i%%replicates_per_web+1)],
                 " -W ",0.1," ",
                 store[1,]$web,sep="")
}
system(command)
com_agg<-read.table("Data/Raw/Full_assembly_model/before_agg_for_matrix.dat")$V2
Sr<-sum(com_agg==com_agg[1])
com_agg<-com_agg[which(com_agg!=com_agg[1])]

system("rm fitness_only.dat", wait=T,ignore.stdout=F,ignore.stderr = F)

# Simulation ####
for (i in (1:total)){
  if (i%%replicates_per_web==0){seeds<-sample(1:10^4, size=replicates_per_web+1, replace = FALSE, prob = NULL)}
  cat(i/total*100,"%","\r")
  
  if(run_with_singularity){
     command<-paste("singularity exec PDMM.img /Project_Code/C_++/NewWeb_new_new/build/Inspect -k random_seed=",seeds[(i%%replicates_per_web+1)],
                    " -W ",store[i,]$agg," ",
                    store[i,]$web,sep="")
  }else{
    command<-paste("/Project_Code/C_++/NewWeb_new_new/build/Inspect -k random_seed=",
                   seeds[(i%%replicates_per_web+1)],
                   " -W ",store[i,]$agg," ",
                   store[i,]$web,sep="")
  }
  system("ulimit -c 0")
  system(command, wait=T,ignore.stdout=T,ignore.stderr = T)
  
  store[i,]$fitness<-read.table("fitness_only.dat",sep = " ")[,4];
  system("rm fitness_only.dat", wait=T,ignore.stdout=F,ignore.stderr = F)
}

#Plot simulated  invasion fitness of a mutant and ... ####

invasion_prob_replicates<-replicates_per_web*length(webfiles)

fitness_plyr<-function(x){
  return(sum(x>=0)/(invasion_prob_replicates))
}

mean_invasion_fitness_per_web<-ddply(store,.(agg,web),summarise,sim_prob=fitness_plyr(fitness))
with(mean_invasion_fitness_per_web,plot(log10(agg),sim_prob))

mu_sd<-find_mu_sd(Sr,sigma)
mu<-unlist(mu_sd[1]);sd<-unlist(mu_sd[2])


invasion_prob<-function(a_m,mu,sd){
  pass<-(log(respiration/(the_conversion_efficiency_epsilon*a_m*K))-as.numeric(mu))/as.numeric(sd)
  return(1-pnorm(pass,0,1))
}

#... compare to predicted ####
x_agg<-mean_invasion_fitness_per_web$agg
mean_invasion_fitness_per_web$pred_prob<-invasion_prob(x_agg,mu,sd)
with(mean_invasion_fitness_per_web,lines(log10(agg),pred_prob,col="red"))

rm(list=setdiff(ls(),c("mean_invasion_fitness_per_web","Sr")))
save.image("Data/Pre_Processed/invasion_prob_full.RData")