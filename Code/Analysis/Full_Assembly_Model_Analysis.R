####### Load Libraries ######

library(plyr)
library(zoo)
library(maxLik)
library(survival)
library(numDeriv)
library(Rfast)

# Functions ##############

{
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
    colSumsHHm1 <- parColSums(HH)-1
    
    # directly <-
    #   colSumsHHm1 * mat.mult(t(HH), HH) > 
    #   outer(  parColSums(HH^2), colSumsHHm1  )
    # diag(directly) <- FALSE #should always be FALSE, apparently some numerical issues
    
    mainResource <- parColWhichMax(HH)
    indirectly <- t( HH[mainResource,] > pickRows(HH,mainResource) )
    return(  indirectly ) #directly |
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
  predict_prob<-function(beta0,beta1,beta2,agg){#predict probability from parameter estimates
    x<-agg
    p<-quadratic(x,beta0,beta1,beta2)
    return(p)
  }
  quadratic<-function(x,c,b,a){c+b*log(x)+a*log(x)^2};
  aicc.loess <- function(fit) {
    # compute AIC_C for a LOESS fit, from:
    #
    # Hurvich, C.M., Simonoff, J.S., and Tsai, C. L. 1998. Smoothing
    # parameter selection in nonparametric regression using an improved
    # Akaike Information Criterion. Journal of the Royal Statistical
    # Society B 60: 271–293.
    #
    # @param fit loess fit
    # @return 'aicc' value
    stopifnot(inherits(fit, 'loess'))
    # parameters
    n <- fit$n
    trace <- fit$trace.hat
    sigma2 <- sum(resid(fit) ^ 2) / (n - 1)
    return(log(sigma2) + 1 + 2 * (2 * (trace + 1)) / (n - trace - 2))
  }
  autoloess <- function(fit, span=c(.1, .9)) {
    # compute loess fit which has span minimizes AIC_C
    #
    # @param fit loess fit; span parameter value doesn't matter
    # @param span a two-value vector representing the minimum and
    # maximum span values
    # @return loess fit with span minimizing the AIC_C function
    stopifnot(inherits(fit, 'loess'), length(span) == 2)
    # loss function in form to be used by optimize
    f <- function(span) aicc.loess(update(fit, span=span))
    # find best loess according to loss function
    o <- optimize(f, span)
    print(o$minimum)
    return(update(fit, span=o$minimum))
  }
  inc <- function(x,what){
    eval.parent(substitute(x <- x + what))
  }
  fitness_func<-function(p_id){
    return(length(which(for_plyr$parent_id==p_id)))
  }
}

# General Parameters #####
run_with_singularity = F

bin<-1000
retrieve_EC_comps<-F

system("cp Master/*.cfg R_params.R")
system("sed -i -e 's/!//g' R_params.R")
system("sed -i -e 's/false/0/g' R_params.R")
system("sed -i -e 's/true/1/g' R_params.R")
system("sed -i -e 's/random/NA/g' R_params.R")

source("R_params.R")
a1<-aggressivity_decay
typical_bodymass_ratio_d<-1
the_mean_bodymass_M<-typical_bodymass_ratio_d^rnorm(1,0,1)
allometric_base_unit<-1000
M<-the_mean_bodymass_M/allometric_base_unit

#Respiration 
respiration_rate_allometric_prefactor = r
respiration_rate_allometric_exponent = 0
respiration<-respiration_rate_allometric_prefactor

# Animal physiology version = 3
the_turnover_rate_r<-respiration_rate_allometric_prefactor*(M^respiration_rate_allometric_exponent);
respiration<-the_turnover_rate_r
production_over_respiration = -1
max_growth_rate<-production_over_respiration*the_turnover_rate_r;
if(max_growth_rate < 0){
  the_handling_time_T=0;
}else{
  the_handling_time_T=conversion_efficiency_epsilon/max_growth_rate;
}
conversion_efficiency_epsilon = epsilon
the_conversion_efficiency_epsilon<-conversion_efficiency_epsilon+the_turnover_rate_r*the_handling_time_T;

the_plant_growth_rate_sigma=growth_rate_allometric_prefactor*(M^growth_rate_allometric_exponent);
the_carrying_capacity=dominating_species_primary_production/the_plant_growth_rate_sigma;
K<-the_carrying_capacity
the_turnover_rate_r=the_plant_growth_rate_sigma;

# Data Loading and Aggregation #######

{
  par(mfrow=c(3,2))
  agg_data = read.table("Data/Raw/Full_assembly_model/list_of_aggressivity_values.dat")
  # with(agg_data,plot(V1,V2))
  species_fitness_data = na.omit(as.data.frame(read.table("Data/Raw/Full_assembly_model/list_of_inserted_deleted_species.dat",fill=F)))
  check<-dim(species_fitness_data)[2]
  
  # if(check >5){
  #   if(species_fitness_data[1,]$V4==species_fitness_data[1,]$V5){
  #     species_fitness_data<-species_fitness_data[c(1:4,6)]
  #   }
  # }
  
  # list_of_inserted_deleted_species << web.s(new_species).unique_id() << " " << web.s(new_species).parent_unique_id() 
  # << " " << old_n_additions << " " << web.s(new_species).aggressivity_g() << " "  << n_additions 
  # << " " << 1-mean_c_star  << " " << web.number_of_plants() << " " << web.number_of_animals() <<endl;
  
  colnames(species_fitness_data)<-c("unique_id","parent_id","sequential_id","agg","addition_number","Sr","Sc","DD") #c("unique_id","parent_id","sequential_id","agg","addition_number","DD")
  
  inc(species_fitness_data[which(species_fitness_data$addition_number!=0),]$addition_number,1)#Increase addition number of each species by one
  onez<-(which(species_fitness_data$parent_id==0.1)[1]-1)
  species_fitness_data$addition_number[1:onez]<-rep(1,onez)#Force addition numbers of first species

  #Clean data
  if (length(which(species_fitness_data$agg==2))==0){
    species_fitness_data[which(species_fitness_data$agg==0.2),]$addition_number<-0.2#Set marker for new mutatn arrival
    species_fitness_data[which(species_fitness_data$agg==0.1),]$addition_number<-0.1#set marker for species extinctions
  }else{
    species_fitness_data[which(species_fitness_data$agg==2),]$addition_number<-2#Set marker for new mutatn arrival
    species_fitness_data[which(species_fitness_data$agg==1),]$addition_number<-1#set marker for species extinctions
  }

  #Get rid of markers
  if (length(which(species_fitness_data$agg==2))==0){
    species_fitness_data<-species_fitness_data[which(species_fitness_data$agg!=0.1),]
    species_fitness_data<-species_fitness_data[which(species_fitness_data$agg!=0.2),]
  }else{
    species_fitness_data<-species_fitness_data[which(species_fitness_data$agg!=1),]
    species_fitness_data<-species_fitness_data[which(species_fitness_data$agg!=2),]
  }


  #Remove plants
  #All plants have the same aggressiveness, find most frequent aggressiveness and eliminate species that have it.
  #save table with frquency of aggressivness values
  agg_freq_table<-as.data.frame(table(species_fitness_data$agg))
  #Agressivenss values with highest frequency will be ones belonging to plants
  max(agg_freq_table[which(agg_freq_table$Freq>20),]$Freq)
  total_number_plants<-agg_freq_table[which(agg_freq_table$Freq>=max(agg_freq_table[which(agg_freq_table$Freq>20),]$Freq)),]$Freq
  plant_agg<-as.numeric(as.character(agg_freq_table[which(agg_freq_table$Freq>=max(agg_freq_table[which(agg_freq_table$Freq>20),]$Freq)),]$Var1))
  species_fitness_no_plants<-species_fitness_data[which(species_fitness_data$agg!=plant_agg),]
  species_fitness_plants<-species_fitness_data[which(species_fitness_data$agg==plant_agg),]

  dim(species_fitness_data)[1]-total_number_plants == dim(species_fitness_no_plants)[1] #Peform check

  species_fitness_data<-species_fitness_no_plants

  #reassign addition number due to deleted plants
  reassign<-species_fitness_data[which(species_fitness_data$addition_number!=0),]
  reassign$addition_number<-seq(1:dim(reassign)[1])
  species_fitness_data[which(species_fitness_data$addition_number!=0),]$addition_number<-reassign$addition_number

  #Remember who died
  species_fitness_data$dead<-0
  species_fitness_data[which(species_fitness_data$addition_number==0),]$dead<-1
  zero_addition<-which(species_fitness_data$addition_number==0)#Keep the index of species that died
  zero_addition_p1<-zero_addition[2:length(zero_addition)]#Create vector where each index in zero_addition is shifted up by one
  diff<-zero_addition_p1-zero_addition[1:(length(zero_addition)-1)]-1#Calculate the difference between each index: number of additions between deletion events
  #Add the cumulative sum of these differences of each index to the addition number of the first deletion event: Time of each deletion
  cum_sum_diff<-cumsum(diff);cum_sum_diff<-cum_sum_diff+species_fitness_data[(zero_addition[1]-1),]$addition_number
  species_fitness_data[zero_addition[1],]$addition_number<-species_fitness_data[(zero_addition[1]-1),]$addition_number
  species_fitness_data[zero_addition[2:length(zero_addition)],]$addition_number<-cum_sum_diff

  born<-species_fitness_data[which(species_fitness_data$dead==0),]
  dead<-species_fitness_data[which(species_fitness_data$dead!=0),]

  keep<-intersect(unique(born$unique_id),unique(dead$unique_id))
  species_fitness_data[which(species_fitness_data$unique_id==keep[1]),]
  keep2<-species_fitness_data[species_fitness_data$unique_id %in% keep,]
  dim(keep2)[1]==2*length(keep)
  species_fitness_data<-species_fitness_data[species_fitness_data$unique_id %in% keep,]

  tbalez<-as.data.frame(table(keep2$unique_id),ncol=2)
  reps<-tbalez[which(tbalez$Freq>2),]$Var1
  #For the moment remove repetition of id's
  keep3 = subset(keep2, !(unique_id %in% reps))
  tbalez<-as.data.frame(table(keep3$unique_id),ncol=2)

  born<-keep3[which(keep3$dead==0),]
  dead<-keep3[which(keep3$dead!=0),]

  born_dead<-merge(born,dead[,c("addition_number","unique_id")], by.x = "unique_id", by.y = "unique_id")
  born_dead$lifetime<-born_dead$addition_number.y-born_dead$addition_number.x
  
   
  colnames(born_dead)  
  born_dead<-born_dead[,-9]
  colnames(born_dead)<- c("unique_id","parent_id","sequential_id","agg","addition_number",
                         "DD","Sr","Sc","death","lifetime")

  #WARNING: Rest of analysis will only deal with species that died if the following is not implemented:
  miss<-length(unique(born$unique_id))-length(unique(dead$unique_id))# We are ommitting 33 species that did not die, that is 5% of total sample space
  missing<-setdiff(unique(born$unique_id),unique(dead$unique_id))
  missing_m<-species_fitness_data[which(species_fitness_data$unique_id==missing[1]),]
  for (i in 2:miss){
    missing_m<-rbind(missing_m,species_fitness_data[which(species_fitness_data$unique_id==missing[i]),])
  }
  missing_m$time_dead<-tail(missing_m,1)$addition_number
  
  missing_m$lifetime<-missing_m$time_dead-missing_m$addition_number
  
  born_dead<-rbind(born_dead,missing_m)
  unique_id<-unique(species_fitness_data$unique_id)
  for_plyr<-species_fitness_data[which(species_fitness_data$dead==0),]

  
  
  fitness<-data.frame(unlist(llply(unique_id,fitness_func)))
  colnames(fitness)<-"fitness"
  
  
  fitness$unique_id<-unique_id
  
  born_dead_fitness<-merge(born_dead,fitness, by.x = "unique_id", by.y = "unique_id")
  
  burn_perc <- 1/5;burn_perc2 <- 4/4
  burn_in <- round(burn_perc*max(species_fitness_data$addition_number))
  burn_in2 <- round(burn_perc2*max(species_fitness_data$addition_number))
  species_fitness_equilibrium<-born_dead_fitness[which(born_dead_fitness$addition_number>burn_in &
                                                       born_dead_fitness$addition_number<burn_in2),]
  

  with(species_fitness_equilibrium,plot(log(rollmean(sort(agg),k=bin)), 
                                                rollmean(lifetime[order(agg)], k=bin),type="l"))
  }

# Numeric Analysis: Lifetime, Birth rate and Fitness ######

{
  species_fitness_equilibrium$birth_rate<-species_fitness_equilibrium$fitness/(species_fitness_equilibrium$lifetime+1)
  species_fitness_equilibrium$l_agg<-log(species_fitness_equilibrium$agg)
  for_nls<-species_fitness_equilibrium[order(species_fitness_equilibrium$agg),]
  for_nls_birth<-species_fitness_equilibrium[which(species_fitness_equilibrium$lifetime>0),];
  for_nls_birth<-for_nls_birth[order(for_nls_birth$agg),]
  for_nls_birth$birth_rate<-with(for_nls_birth,fitness/lifetime)
  species_fitness_equilibrium$birth_rate<-species_fitness_equilibrium$fitness/(species_fitness_equilibrium$lifetime+1)
  lm_birth_agg<-lm(birth_rate ~ log(agg)+I(log(agg)^2),for_nls_birth);coef(lm_birth_agg)
  b_nls_quadratic<-nls(birth_rate ~ c+b*log(agg)+a*log(agg)^2,for_nls_birth,
                       start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3]))
  b_cor_quadratic<-cor(for_nls_birth$birth_rate,predict(b_nls_quadratic))
  x<-seq(min(species_fitness_equilibrium$agg),max(species_fitness_equilibrium$agg),length.out = 1000)
  paramz61<-coef(b_nls_quadratic)
  prob<-predict_prob(as.numeric(paramz61[1]),as.numeric(paramz61[2]),as.numeric(paramz61[3]),x)#Probabillity of yielding fit offspring
  lm_fitness_agg<-lm(fitness ~ log(agg)+I(log(agg)^2),species_fitness_equilibrium);coef(lm_fitness_agg)
  f_nls_quadratic<-nls(fitness ~ c+b*log(agg)+a*log(agg)^2,for_nls,
                       start=list(c=coef(lm_fitness_agg)[1],b=coef(lm_fitness_agg)[2],a=coef(lm_fitness_agg)[3]))
}

# Analysis: Retrieval of eploitative comp. components #####

if(retrieve_EC_comps){
  all_webs<-sample(list.files(pattern = "web*[0-9][0-9][0-9][0-9]",path="/Data/Raw/Full_assembly_model"),300)
  all<-data.frame();all_mins<-data.frame()
  counter<-0
  
  for (ww in all_webs){
    counter<-counter+1
    cat(signif(counter/length(all_webs)*100,3),"%","\r")
    seed<-sample(c(1:100000),1)
    system("rm before_*")
  
    if(run_with_singularity){
      command<-paste("singularity exec PDMM.img Code/Simulations/Full_Assembly_Model/NewWeb_new_new_new_new/build/Inspect -k random_seed=",
                      seed," -W ",1," ",ww,sep="")
    }else{
      command<-paste("Code/Simulations/Full_Assembly_Model/NewWeb_new_new_new_new/build/Inspect -k random_seed=",
                     seed," -W ",1," ",ww,sep="")
    }
  
    
    system(command,ignore.stdout = TRUE,wait=TRUE)
    before_agg<-read.table("/Data/Raw/Full_assembly_model/before_agg_for_matrix.dat");colnames(before_agg)<-c("id","agg","attack","B")
    before_test<-read.table("/Data/Raw/Full_assembly_model/before_interaction_matrix.dat",sep=",");
    plant_agg<-before_agg$agg[1]
    S_b<-dim(before_test)[1];Sr_b<-sum(before_agg$agg==plant_agg);Sc_b<-S_b-Sr_b
    before_agg<-before_agg[which(before_agg$agg!=plant_agg),]
    before<-t(data.frame(before_test[(Sr_b+1):S_b,1:Sr_b])*before_agg$agg) #Each row represents a consumer
    
    who_save <- 1:dim(before)[2]
    repeat{
      # For each resident sample one corresponding comp. val. at random
      # Focal
      who_save_2 <- sample(1:dim(before)[2]) #sample(1:dim(HH)[2],dim(HH)[2]-5)#sample(1:dim(HH)[2],5) # Secondary
      if(sum(who_save == who_save_2)==0){break}
    }
    
    par_col_sum_yo <- colSums(before);  par_col_sum_yo2 <- colSums(before^2)
    
    t_inv_HHHH<-llply(who_save,function(x){
      sum(before[x,]*before[who_save_2[x],])
    })
    
    t_inv_HHH <- mat.mult(t(before), before)
    ine_all <-  (t(((par_col_sum_yo -1) / t_inv_HHH))) / ( (par_col_sum_yo -1)/par_col_sum_yo2 )
    
    ine <- do.call(cbind,llply(1:dim(before)[2],function(x){
      return(ine_all[-x,x])
    }))
    
    t_inv_HHH_all <- do.call(cbind,llply(1:dim(before)[2],function(x){
      return(t_inv_HHH[-x,x])
    }))
    
    
    ine_mins<-parColWhichMins(ine) 
    
    temp_save<-matrix(c(before_agg$agg,
                        log(par_col_sum_yo),
                        unlist(t_inv_HHHH),
                        par_col_sum_yo[who_save_2],par_col_sum_yo2[who_save_2]),
                      ncol=5,byrow = F)
    
    all<-rbind(all,temp_save)
    
    
    who_save <- 1:dim(before)[2]
    who_save_min <- colMins(ine)
    
    
    temp_save_min<-matrix(c(unlist(before_agg$agg),
                            log(par_col_sum_yo[who_save]),
                            unlist(llply(1:length(who_save),function(x) t_inv_HHH_all[who_save_min[x],who_save[x]])),
                            par_col_sum_yo[who_save_min],par_col_sum_yo2[who_save_min],
                            matrix(colMaxs(before,value=T),ncol=1),
                            matrix(1-colSums(before^2)/(colSums(before)^2),ncol=1)),
                          ncol=7,byrow = F)
    
    all_mins<-rbind(all_mins,temp_save_min)
    
    
    
  }
  colnames(all)<-c("agg","A1_com","A2_com","B1_com","B2_com")
  colnames(all_mins)<-c("agg","A1_min","A2_min","B1_min","B2_min","maxx","DD")
  all_mins$l_agg<-log(all_mins$agg)
  all$ine_com<-with(all,((A1_com)/A2_com)/((B1_com)/B2_com))
  all$l_agg<-log(all$agg)
  # all_com<-which(all$ine>0)
  plot_x<-rollmean(sort(log(all$agg)),k=bin)
  cc<-(range(plot_x)[2]-range(plot_x)[1])/5
  range(plot_x)[1]+cc
  hist(log(all$agg))
  colz<-heat.colors(5)
}

# Analysis: Aanlytic Birth Rate #####

sigma<-4
resp<-0.1;epsilon<-0.1;K<-1

# choose agg value to test
agg<-species_fitness_equilibrium$agg
a1<-(0.8^0.5)
a2<-(1.3^0.5)
m_a<-log(a1)
s_a<-log(a2)


# Load webs, claculate invasion prob from these

# Pre-prepare data.frame
all<-data.frame()
load("Data/Pre_processed/1_optimised_K.RData")
P_b_m<-K_dt#mean(read.table("list_of_plant_biomass_values.dat")$V3)
{
  invasion_prob_c<-function(l_a_m,mu,sd){
    pass<-(log(resp/(epsilon*exp(l_a_m)*K*P_b_m))-as.numeric(mu))/(as.numeric(sd))
    return(1-pnorm(pass,0,1))
  }
  
  invasion_prob<-function(l_a_m,mu,sd){
    pass<-(log(resp/(epsilon*exp(l_a_m)*K))-as.numeric(mu))/(as.numeric(sd))
    return(1-pnorm(pass,0,1))
  }
  
  #loop of webs
  constanttt<-K*epsilon/resp
  findFitConsumer_without <- function(aggs,HH){
    Sr<-dim(HH)[1]
    repeat{
      parent<-sample(aggs,1)
      agg <- a1*a2^rnorm(1)*parent
      H <- agg*exp(rnorm(Sr,sd=sigma))*constanttt
      if(sum(H)>1) break
    }
    return(parent)
  }
  
  findFitConsumer_with <- function(aggs,HH){
    Sr<-dim(HH)[1];par_col_sum_yo<-parColSums(HH)
    par_col_sum_yo2<-parColSums(HH^2)
    repeat{
      parent<-sample(aggs,1)
      agg <- a1*a2^rnorm(1)*parent
      H <- agg*exp(rnorm(Sr,sd=sigma))*constanttt
      t_inv_HH<-t(H)%*%HH
      if(sum(H)>1 & !outcompeted(H,t_inv_HH, par_col_sum_yo, par_col_sum_yo2 )) break
    }
    return(parent)
  }
  
  #Probability of yielding mutant given resident
  mut_prob_given_ar<-function(l_a_m,l_a_r){
    mu_a<-m_a+l_a_r
    return(exp(-((l_a_m-mu_a)^2)/(2*(s_a^2)))/((2*pi*(s_a^2))^0.5))
  }
  
  #Joint probability of yielding mutant (given resident) and that the mutant succesfuly invades
  joint_prob<-function(l_a_mm,l_a_rr,mu,sd){
    return(mut_prob_given_ar(l_a_mm,l_a_rr)*invasion_prob(l_a_mm,mu,sd))
  }
  
  joint_prob_c<-function(l_a_mm,l_a_rr,mu,sd){
    return(mut_prob_given_ar(l_a_mm,l_a_rr)*invasion_prob_c(l_a_mm,mu,sd))
  }
  
  #Birth probability: intergrate joint probability over all possible mutants.
  birth_prob_c<-function(l_a_rr,mu,sd){
    return(invasion_prob_c(l_a_rr+m_a,mu,sd))
    # return(integrate(joint_prob,lower=-Inf,upper=Inf,l_a_rr=l_a_rr,mu=mu,sd=sd)$value)
  }
  
  birth_prob<-function(l_a_rr,mu,sd){
    return(invasion_prob(l_a_rr+m_a,mu,sd))
    # return(integrate(joint_prob,lower=-10,upper=10,l_a_rr=l_a_rr,mu=mu,sd=sd)$value)
  }
}

loop_webs<-400
all_webs<-list.files(pattern = "web*")


for (i in 1:loop_webs){
  
  try(system("rm *agg_for_matrix.dat"),silent=T);
  try(system("rm *interaction_matrix.dat"),silent=T);
  
  k<-all_webs[i]
  cat(i/loop_webs*100,"%","\r")
  web_no<-as.numeric(strsplit(strsplit(k,split="web")[[1]][2],split="x")[[1]][1])
  seed<-sample(c(1:100000),1)
  if(run_with_singularity){
    command<-paste("singularity exec PDMM.img Code/Simulations/Full_Assembly_Model/NewWeb_new_new_new_new/build/Inspect -k random_seed=",
    seed," -W ",1," ",k,sep="")
  }else{
    command<-paste("Code/Simulations/Full_Assembly_Model/NewWeb_new_new_new_new/build/Inspect -k random_seed=",
                   seed," -W ",1," ",k,sep="")
  }
 

  system(command,ignore.stdout = TRUE,wait=TRUE)
  
  # Original community before any invasion.  
  before_agg<-read.table("Data/Raw/Full_assembly_model/before_agg_for_matrix.dat");colnames(before_agg)<-c("id","agg","attack","B")
  before_test<-read.table("Data/Raw/Full_assembly_model/before_interaction_matrix.dat",sep=",");
  plant_agg<-before_agg$agg[1]
  S_b<-dim(before_test)[1];Sr_b<-sum(before_agg$agg==plant_agg);Sc_b<-S_b-Sr_b
  before<-data.frame(before_test[(Sr_b+1):S_b,1:Sr_b]) #Each row represents a consumer
  Srs<-data.frame(table(before_agg$agg));agg_Sr<-as.numeric(as.character(Srs[which(Srs$Freq>1),]$Var1));
  
  HH<-t(before);aggz<-before_agg$agg[(Sr_b+1):S_b]
  current<-data.frame(aggz);colnames(current)<-"agg";
  mu_sd<-find_mu_sd(round(mean(dim(HH)[1])),sigma)
  current$b_p_Axel<-0
  current$b_p_Axel_c<-0
  # current$b_p_Axel_a<-0
  # current$b_p_Axel_a_c<-0
  mu<-mu_sd[1];sd<-mu_sd[2]
  # mu_sd_2<-mu_sd_dt[which(mu_sd_dt$Sr==Sr_b),]
  # mu_c<-mu_sd_2[1];sd_c<-mu_sd_2[2]
  
  all_b_p<-unlist(llply(aggz, function(x) birth_prob(log(x),mu,sd)))
  all_b_p_c<-unlist(llply(aggz, function(x) birth_prob_c(log(x),mu,sd)))
  # all_b_p_2<-unlist(llply(aggz, function(x) birth_prob(log(x),mu_c,sd_c)))
  # all_b_p_2_c<-unlist(llply(aggz, function(x) birth_prob_c(log(x),mu_c,sd_c)))

  
  current[,2]<-unlist(llply(1:dim(HH)[2],function(x) all_b_p[x]/sum(all_b_p)))
  current[,3]<-unlist(llply(1:dim(HH)[2],function(x) all_b_p_c[x]/sum(all_b_p_c)))
  # current[,4]<-unlist(llply(1:dim(HH)[2],function(x) all_b_p_2[x]/sum(all_b_p_c_2)))
  # current[,5]<-unlist(llply(1:dim(HH)[2],function(x) all_b_p_2_c[x]/sum(all_b_p_c_3)))
 
  current$Sr<-dim(HH)[1]
  all<-rbind(all,current)
}

colnames(all)<-c("agg","b_p_Axel","b_p_Axel_c","Sr")
save.image("Data_full.RData")

