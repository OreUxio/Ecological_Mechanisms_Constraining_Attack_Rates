args = commandArgs(trailingOnly=TRUE)
library(Rfast)
library(plyr)
library(zoo)
library(survival)
library(maxLik)
library(plyr)
library(ggplot2)
library(deSolve)
library(RColorBrewer)
library(scales)
library(plotfunctions)
do_others<-0
setwd("/home/ore/Documents/Academic/Ecological_Mechanisms_Constraining_Attack_Rates/all/Data")

################ Start #################

########################################
            #  Chapter Two  #  
########################################

# Figure 2.1 #
{
mult<-0.65
load("Pre_Processed/Data_deconstructed.RData")
all_spec<-read.table("Pre_Processed/agg_time_series_aggregate.txt") # O
all_spec_time<-read.table("Pre_Processed/agg_time_series_aggregate_time.txt") # 0
all_spec<-all_spec[,-21]
all_spec_time<-all_spec_time[,-20]
colfunc <- colorRampPalette(c("blue", "red"))
ordered<-all_spec[,(1+order(all_spec[1,2:dim(all_spec)[2]]))]
try({pdf("../Images/spec_a_full.pdf",height=10*0.5*mult,width=10*mult)
  #par(mfrow=c(1,2))
  comm_cexx<-1.5*mult
  mtext_cex<-1*mult
  colz<-colfunc(dim(all_spec)[2])#rainbow(dim(all_spec)[2])
  opac<-0.6
  {
    layout.matrix <- matrix(c(1, 2,3,3), nrow = 2, ncol = 2,byrow = T)
    layout(mat = layout.matrix,
           widths = c(2, 2),
           heights = c(2,0.5)) # Widths of the two columns
    par(mar=c(0,4,2,1))
    
    animals<-read.table("Raw/Full_assembly_model/list_of_noofanimalspecies_values.dat") # O
    plants<-read.table("Raw/Full_assembly_model/list_of_noofplantspecies_values.dat") # O
    plants$V1<-plants$V1*1.3
    animals$V1<-plants$V1*1.3
    
    ylimm<-max(plants$V2)
    plot(plants$V1,plants$V2,col="green",type="l", #ylim=c(min(log(plants)),max(log(plants)))
         ylab="Species richness",xlab="",ylim=c(0,ylimm),main="",yaxs="i",xaxt="n",
         cex.axis=comm_cexx,cex.lab= comm_cexx)#xlim=c(0,500000))
    axis(1,cex.axis= comm_cexx,cex.lab= comm_cexx,at=c(0,100000,200000,300000,400000,500000),
         labels=c(0,1,2,3,4,5))
    lines(plants$V1,animals$V2,ylim=c(0,ylimm),col="red",type="l")
    mtext("(a) Species richness",side=3,line=-1,adj=0.01,cex=mtext_cex)
    legend(x="bottomright",legend = c("Consumers","Resources"),col=c("red","green"),lty=1,lwd=2,bty="n")
  }
  
  par(mar=c(0,4.5,2,1))
  all_spec_temp<-all_spec[,2]
  plot(all_spec_time[,1],ordered[,1],cex=0.1,type="l", #plot(all_spec_decon$id*2,all_spec_temp,cex=0.1,type="l",
       xlab="",cex.axis=comm_cexx,cex.lab=comm_cexx,xaxt="n",#xaxs="i",yaxs="i",
       ylab=expression(bar(log[10](a[i][,][t]))),xlim=c(0,500000),
       ylim=c(log10(exp(log(1*10^-5.5))),log10(exp(log(1*10^-3.25)))),col=alpha(colz[1],opac))
  axis(1,cex.axis=comm_cexx,cex.lab=comm_cexx,at=c(0,100000,200000,300000,400000,500000),labels=c(0,1,2,3,4,5))
  
  for (i in sample((2:dim(ordered)[2]))){
    # i<-i+1
    all_spec_temp<-ordered[,i]
    lines(all_spec_time[,(i)],all_spec_temp,cex=0.1,type="l", #all_spec_decon$id*2,all_spec_temp,
          xlab="",
          ylab=expression(bar(log[10](a[i][,][t]))),ylim=c(log(1*10^-6),log(1*10^-4.25)),col=alpha(colz[i],opac))
  }
  mtext("(b) Base attack rate",side=3,line=-1,adj=0.01,cex=mtext_cex)
  mtext("Initial base attack rate",side=4,line=-3.75,adj=0.95,cex=mtext_cex*0.75)
  gradientLegend(valRange=c(-6,-4), color = alpha(colz,opac), coords=TRUE, depth=0.1,pos=0.8,pos.num=2,
                 side=4,border.col = "black")
  
  plot.new()
  par(mar=c(3,4,2,1))
  mtext(expression(Time~(number~of~consumers~added)~~x~10^5),side=1,line=2,cex=mtext_cex)
dev.off()},silent=T)
}

# Figure 2.2 #
{
mult<-0.65
load("Pre_Processed/Data_deconstructed.RData")
load("Pre_Processed/Data_full.RData") # O
try({pdf("../Images/animal_plants.pdf",height=10*0.5*mult,width=10*mult)
  #par(mfrow=c(1,2))
  comm_cexx<-1.5*mult
  mtext_cex<-1*mult
  {
    layout.matrix <- matrix(c(1, 2,3,3), nrow = 2, ncol = 2,byrow = T)
    layout(mat = layout.matrix,
           widths = c(2, 2),
           heights = c(2,0.5)) # Widths of the two columns
    par(mar=c(0,4,2,1))
  
    animals<-read.table("Raw/Full_assembly_model/list_of_noofanimalspecies_values.dat") #O
    plants<-read.table("Raw/Full_assembly_model/list_of_noofplantspecies_values.dat") #O
    plants$V1<-plants$V1*1.3
    animals$V1<-plants$V1*1.3
    
    ylimm<-max(all_species_keep$n_P,plants$V2)
    plot(plants$V1,plants$V2,col="green",type="l", #ylim=c(min(log(plants)),max(log(plants)))
         ylab="Species richness",xlab="",ylim=c(0,ylimm),main="",yaxs="i",xaxt="n",
         cex.axis=comm_cexx,cex.lab= comm_cexx)#xlim=c(0,500000))
    axis(1,cex.axis= comm_cexx,cex.lab= comm_cexx,at=c(0,100000,200000,300000,400000,500000),
         labels=c(0,1,2,3,4,5))
    lines(plants$V1,animals$V2,ylim=c(0,ylimm),col="red",type="l")
    mtext("(a) Full model",side=3,line=-1,adj=0.01,cex=mtext_cex)
    
    
    with(all_species_keep[which((all_species_keep$id*2) %in% plants$V1),],
         plot(id*2,n_P,type="l",col="green",ylim=c(0,ylimm),xlab="",main="",
              cex.axis= comm_cexx,cex.lab= comm_cexx,yaxs="i",xaxt="n",
              ylab=""))
    axis(1,cex.axis= comm_cexx,cex.lab= comm_cexx,at=c(0,100000,200000,300000,400000,500000),
         labels=c(0,1,2,3,4,5))
    with(all_species_keep[which((all_species_keep$id*2) %in% plants$V1),],
         lines(id*2,n_C,type="l",col="red",ylim=c(0,max(n_P))))
    
    mtext("(b) Deconstructed model",side=3,line=-1,adj=0.01,cex=mtext_cex)
    legend("topright", legend=c("Resources","Consumers" ), pch=c(NA, NA), lty=c(1,1), 
           col=c("green","red"),cex=1.3*mult,bty="n")
    plot.new()
    par(mar=c(3,4,2,1))
    mtext(expression(Time~(number~of~consumers~added)~~x~10^5),side=1,line=2,cex=mtext_cex)
  }
dev.off()},silent=T)
all_spec_decon<-read.table("Pre_Processed/agg_time_series_deconstructed_aggregate.txt")
all_spec_decon_time<-read.table("Pre_Processed/agg_time_series_deconstructed_aggregate_time.txt")
all_spec<-read.table("Pre_Processed/agg_time_series_aggregate.txt")
all_spec_time<-read.table("Pre_Processed/agg_time_series_aggregate_time.txt")
all_spec<-all_spec[,-21]
all_spec_time<-all_spec_time[,-20]
mult<-0.65
colfunc <- colorRampPalette(c("blue", "red"))
}

# Figure 2.3 #
{
pdf("../Images/agg_time_series.pdf",height=10*0.5*mult,width=10*mult)
{
  mult<-0.65
  load("Pre_Processed/Data_deconstructed.RData")
  load("Pre_Processed/Data_full.RData") # O
  ordered<-all_spec[,(1+order(all_spec[1,2:dim(all_spec)[2]]))]
  comm_cexx<-1.5*mult
  mtext_cex<-1*mult
  colz<-colfunc(dim(all_spec_decon)[2])
  layout.matrix <- matrix(c(1, 2,3,3), nrow = 2, ncol = 2,byrow = T)
  layout(mat = layout.matrix,
         widths = c(2, 2),
         heights = c(2,0.5)) # Widths of the two columns
  par(mar=c(0,4.5,2,1))
  all_spec_temp<-all_spec[,2]
  plot(all_spec_time[,1],ordered[,1],cex=0.1,type="l", #plot(all_spec_decon$id*2,all_spec_temp,cex=0.1,type="l",
       xlab="",cex.axis=comm_cexx,cex.lab=comm_cexx,xaxt="n",#xaxs="i",yaxs="i",
       ylab=expression(bar(log[10](a[i][,][t]))),xlim=c(0,500000),
       ylim=c(log10(exp(log(1*10^-6.25))),log10(exp(log(1*10^-3.25)))),col=alpha(colz[1],opac))
  axis(1,cex.axis=comm_cexx,cex.lab=comm_cexx,at=c(0,100000,200000,300000,400000,500000),labels=c(0,1,2,3,4,5))
  
  for (i in sample((2:dim(ordered)[2]))){
    # i<-i+1
    all_spec_temp<-ordered[,i]
    lines(all_spec_time[,(i)],all_spec_temp,cex=0.1,type="l", #all_spec_decon$id*2,all_spec_temp,
          xlab="",
          ylab=expression(bar(log[10](a[i][,][t]))),ylim=c(log(1*10^-6),log(1*10^-4.25)),col=alpha(colz[i],opac))
  }
  mtext("(a) Full model",side=3,line=-2,adj=0.1,cex=mtext_cex)
  
  ordered<-all_spec_decon[,(1+order(all_spec_decon[1,2:dim(all_spec_decon)[2]]))]
  all_spec_decon_temp<-ordered[,1]
  plot(all_spec_decon_time[,1]*2,log10(exp(all_spec_decon_temp)),type="l",col=alpha(colz[1],opac),xlim=c(0,500000),
       ylim=c(log10(exp(log(1*10^-6.25))),log10(exp(log(1*10^-3.25)))),#,cex=0.1,,#yaxs="i",xaxs="i",
       xlab="",cex.axis=comm_cexx,cex.lab=comm_cexx,xaxt="n",ylab="")
  #      col=colz[1])
  axis(1,cex.axis=comm_cexx,cex.lab=comm_cexx,at=c(0,100000,200000,300000,400000,500000),labels=c(0,1,2,3,4,5))
       # labels=c(0,1,2,3,4,5))
  for (i in  sample((2:dim(ordered)[2]))){
    all_spec_decon_temp<-ordered[,i]
    lines(all_spec_decon_time[,i]*2,log10(exp(all_spec_decon_temp)),cex=0.1,type="l",
          xlab="Time (Number of consumers added)",
          ylab=expression(bar(log(a[i][,][t]))),ylim=c(log(1*10^-6),log(1*10^-4.25)),col=alpha(colz[i],opac))
  }
  mtext("(b) Deconstructed model",side=3,line=-2,adj=0.1,cex=mtext_cex)
  mtext("Initial base attack rate",side=4,line=-3.75,adj=0.95,cex=mtext_cex*0.75)
  gradientLegend(valRange=c(-6,-4), color = alpha(colz,opac), coords=TRUE, depth=0.1,pos=0.8,pos.num=2,
                  side=4,border.col = "black")

  plot.new()
  par(mar=c(3,4,2,1))
  # mtext("Time (number of consumers added)",side=1,cex=1*mult)
  mtext(expression(Time~(number~of~consumers~added)~~x~10^5),side=1,line=2,cex=mtext_cex)


}
dev.off()
}

# Figure 2.4: Please go to Invasion_Probability.R

# Figure 2.5 #
{
  load("Pre_Processed/Data_deconstructed.RData")
  growth<-0.1
  compModel <- function(Time, B, Pars) {
    with(as.list(Pars), {
      dB <- as.vector((r + C %*% B)*B) + inv*0.1
      return(list(dB))
    })
  }
  
  relax_all <- function(pars,BB,tRelax = 200){
    times <- c(0:tRelax)
    sol <-
      ode(y = c(B=BB),times = times,maxsteps = 5000,
          func = compModel,parms = pars )
    return(initfunc = sol)#
    # return(BB<--inv(pars$C)%*%pars$r)
  }
  
  numb_sr<-1
  
  
  check_min<-function(x,initial_K){
    
    links_1<-x*log(K/(inv*growth))*resp/(K*epsilon)
    C_1<-matrix(rep(0,(numb_sr+1)^2),ncol=numb_sr+1)
    C_1[1,]<-c(0,epsilon*links_1)
    
    for (i in 1:numb_sr){
      C_1[(i+1),1]<--links_1[i]
      C_1[(i+1),(i+1)]<--growth/K
    }
    
    pars <- list(
      C = C_1,
      r = c(-resp,rep(growth,numb_sr)),
      i = rep(inv,numb_sr+1)
    )
    
    BRelaxed_1<-c(inv,initial_K)
    BRelaxed_1 <- relax_all(pars,BRelaxed_1,200)
    if(min(BRelaxed_1[,3])<inv){ 
      out <- 1
    }else{
      out<-0
    }
    return(min(BRelaxed_1[,3]))
  }
  
  # percc_all <- data.frame(rep(seq(1.95,1.96,length.out = 10),100))
  percc_all <- data.frame(seq(1,2,length.out = 100))
  colnames(percc_all)<-"perc";
  ini_valz<- seq(1,0.25,length.out=8)
  for (i in ini_valz){
    percc_all<-cbind(percc_all,unlist(llply(percc_all$perc,function(x) check_min(x,i*K))))
  }
  colfunc <- colorRampPalette(c("green","blue"))
  colz<-colfunc(length(ini_valz))#rainbow(dim(all_spec)[2])
  library(plotfunctions)
  mult<-1
  cex.plot<-2
  comm_cexx<-1.5*mult
  mtext_cex<-2*mult
  pdf("Images/Over_exploitative.pdf",height=10*mult,width=10*mult)
  par(mar=(c(4,6,2,4)))
  with(percc_all,plot(perc,log10(percc_all[,2]),type="l",xlab=expression(omega),col=colz[1],cex.axis=comm_cexx,cex.lab=comm_cexx*2,cex=cex.plot,
                      ylab=expression(min[t]~B^R*(t),B[Min]),ylim=c(-6.5,-5),xaxs="i",yaxs="i"));
  for (i in 2:(length(ini_valz)-1)){
    with(percc_all,lines(perc,log10(percc_all[,1+i]),col=colz[i]));
  }
  
  abline(h=log10(inv),lty=2)
  colfunc <- colorRampPalette(c("blue","green"))
  colz<-colfunc(length(ini_valz))
  gradientLegend(valRange=c(min(ini_valz),max(ini_valz)), color = colz, coords=TRUE, depth=0.1,pos=0.8,pos.num=2,
                 side=4,border.col = "black")
  mtext(expression(B^R*(0)),side=3,line=-1.5,adj=0.99)#side=4,line=-3.75,adj=0.95,cex=mtext_cex*0.75)
  dev.off()
}

# Figure 2.6 #
{
  load("Pre_Processed/Data_deconstructed.RData")
  load("Pre_Processed/Numerical_serial_extinction_sub_model/Numerical_serial_extinction_sub_model.RData") #O
  library(scales)
  mult<-1.75
  
  try({pdf("../Images/serial_extinction.pdf",height=10*0.5*mult,width=10*0.5*mult)
    points_to_plot<-50
    cut_of<--10
    comm_cex_lab<-1.5
    par(mfrow=c(1,1))
    par(mar=c(4,6.5,1,1))
    {
      blue_red<-colorRampPalette(c("blue", "red"))
      colz<-blue_red(3)
      maxi<-max(log10((10^links_1_c_mm[1])+1),log10((10^links_2_c_mm[1])+1),log10((10^links_3_c_mm[1])+1))
      mini<-min(log10(2-links_1_DD_b_mm),log10(2-links_2_DD_b_mm),log10(2-links_3_DD_b_mm))
      number_adds<-length(links_1_c_mm)
      minnn<-which.min(c(agg_1,agg_2,agg_3))
      maxxx<-which.max(c(agg_1,agg_2,agg_3))
      sample_n<-100
      plot(log10(1:number_adds),log10((10^links_1_c_mm)+1),col=colz[1],cex=0.1,type="l",cex.axis=1.2,cex.lab=1.2,xaxs="i",yaxs="i",
           ylim=c(mini,maxi),xlab="Number of potential serial extinctions",ylab=expression(log[10](sum(H[k][i],k==1, S[R]))),xaxt="n",
           cex.lab=comm_cex_lab,cex.axis=comm_cex_lab)
      # abline(h=log10(2-links_1_DD_b_mm)[1],cex=0.1,col=colz[1],lty=3)
      axis(1,cex.axis=comm_cex_lab,cex.lab=comm_cex_lab,at=log10(c(1,2,5,10,50,100,250)),
           labels=c("Invasion",1,5,10,50,100,250))
      axis(1,cex.axis=comm_cex_lab,cex.lab=comm_cex_lab,at=log10(2),
           labels=c(1))
      # mtext(text = expression(log[10](C)),
      #       side = 2,line = 2.5,cex=1)
      
      lines(log10(1:number_adds),log10((10^links_2_c_mm)+1),cex=0.1,col=colz[2],type="l")
      # abline(h=log10(2-links_2_DD_b_mm)[1],cex=0.1,col=colz[2],lty=3)
      
      lines(log10(1:number_adds),log10((10^links_3_c_mm)+1),cex=0.1,type="l",col=colz[3])
      # abline(h=log10(2-links_3_DD_b_mm)[1],cex=0.1,col=colz[3],lty=3)
      
      lines(log10(1:number_adds),log10((10^links_1_c_mmm)+1),col=colz[1],lty=2)
      lines(log10(1:number_adds),log10((10^links_2_c_mmm)+1),cex=0.1,col=colz[2],type="l",lty=2)
      lines(log10(1:number_adds),log10((10^links_3_c_mmm)+1),cex=0.1,type="l",col=colz[3],lty=2)
      legend("topright",legend=c(expression(a[i]==10^-5.5),expression(a[i]==10^-4.75),expression(a[i]==10^-4),"Deconstructed model","Full model"),
             col=c(colz[3],colz[2],colz[1],"black","black"),lty=c(1,1,1,2))
    }
    dev.off()},silent=T)
}

# Figure 2.7 #
{
  mult<-1
  cex.plot<-2
  load("Pre_Processed/Data_deconstructed.RData")
  pdf("Images/Validation_exploitative.pdf",height=10*mult,width=10*mult)
  #par(mfrow=c(1,2))
  comm_cexx<-1.5*mult
  mtext_cex<-1*mult
  amount_list_2_diff_temp<-read.csv(list.files(pattern="diff",path="/Data/Raw/Validation/")[1])
  plot(as.numeric(strsplit(list.files(pattern="diff",path="/Data/Raw/Validation/")[1],"amount")[[1]][1]),cex.axis=comm_cexx,cex.lab=comm_cexx,
       mean(amount_list_2_diff_temp$d1),xlim=c(0,500),ylim=c(0,1),xlab="Number of consumers",ylab="Accuracy",cex=cex.plot)
  for (j in list.files(pattern="diff",path="/Data/Raw/Validation/")){
    
    amount_list_2_diff_temp<-read.csv(j)
    points(as.numeric(strsplit(j,"amount")[[1]][1]),
           mean(amount_list_2_diff_temp$d1),xlim=c(0,200),ylim=c(0,1),cex=cex.plot)
    
  }
  dev.off()
}

# Figure 2.8 #
{
  load("Pre_Processed/Data_deconstructed.RData")
  mult<-1
  pdf("Images/Validation_phyrric.pdf",height=10*mult,width=10*mult)
  {
    type_Is<-list.files(pattern="phyrric_tests",path="/Data/Raw/Validation/")
    type_IIs<-list.files(pattern="phyrric_2_tests",path="/Data/Raw/Validation/")
    comm_cexx<-1.5*mult
    mtext_cex<-1*mult
    amount_list_1_diff_temp<-read.csv(type_Is[1])
    amount_list_2_diff_temp<-read.csv(type_IIs[1])
    
    alive_ws_dead<-(1-mean(unlist(amount_list_2_diff_temp$d1)))*length(agg_all)
    correct_alive<-mean(unlist(amount_list_2_diff_temp$d1))*length(agg_all)
    dead_ws_alive<-(1-mean(unlist(amount_list_1_diff_temp$d1)))*length(agg_all)
    correct_dead<-mean(unlist(amount_list_1_diff_temp$d1))*length(agg_all)
    dead_alive<-matrix(c(correct_alive, alive_ws_dead, dead_ws_alive, correct_dead),
                       nrow = 2,
                       dimnames = list(Guess = c("Alive", "Dead"),
                                       Truth = c("Alive", "Dead")))
    if(fisher.test(dead_alive, alternative = "greater")$p.value<0.01){
      colz="green"
    }else{
      colz="red"
    }
    
    Sr<-as.numeric(strsplit(type_Is[1],"phyrric_tests")[[1]][1])
    plot(Sr,cex.axis=comm_cexx,cex.lab=comm_cexx,pch=5,
         mean(amount_list_1_diff_temp$d1),xlim=c(0,600),ylim=c(0,1),xlab="Number of resources",ylab="Accuracy")#,col=colz)
    points(Sr,mean(amount_list_2_diff_temp$d1),pch=16) #,col=colz
    
    for (j in 2:length(type_Is)){
      
      amount_list_1_diff_temp<-read.csv(type_Is[j])
      amount_list_2_diff_temp<-read.csv(type_IIs[j])
      
      
      alive_ws_dead<-(1-mean(unlist(amount_list_2_diff_temp$d1)))*length(agg_all)
      correct_alive<-mean(unlist(amount_list_2_diff_temp$d1))*length(agg_all)
      dead_ws_alive<-(1-mean(unlist(amount_list_1_diff_temp$d1)))*length(agg_all)
      correct_dead<-mean(unlist(amount_list_1_diff_temp$d1))*length(agg_all)
      dead_alive<-matrix(c(correct_alive, alive_ws_dead, dead_ws_alive, correct_dead),
                         nrow = 2,
                         dimnames = list(Guess = c("Alive", "Dead"),
                                         Truth = c("Alive", "Dead")))
      if(fisher.test(dead_alive, alternative = "greater")$p.value<0.01){
        colz="green"
      }else{
        colz="red"
      }
      
      
      Sr<-as.numeric(strsplit(type_Is[j],"phyrric_tests")[[1]][1])
      points( Sr, mean(amount_list_1_diff_temp$d1),pch=5) #col=colz,
      points( Sr, mean(amount_list_2_diff_temp$d1),pch=16) #col=colz,
      legend("bottomright",legend=c("Type I","Type II"),
             #paste("Failed to reject",expression(H[0])), paste("Failed",expression(H[0]))),
             pch=c(5,16,NA,NA),lty=c(NA,NA,1,1),
             col=c("black","black","green","red"),bty="n")
      
    }
  }
  dev.off()
}

########################################
            # Chapter Three #
#################################£######

# Figure 3.1 #
{
  mult<-0.85
  try({pdf("../Images/All_deconstructed.pdf",width=12.5*mult,height=12.5*mult)
    
    load("Pre_Processed/Data_full.RData") #O
    comm_cex_lab<-2.2*mult
    cubic<-function(x,c,b,a,d){c+b*x+a*x^2+d*x^3}
    layout.matrix <- matrix(c(1, 2, 3, 4, 5,
                              6, 7, 8, 9, 5), nrow = 5, ncol = 2)
    heights <- c(0.5, 1,1.2,0.7,0.45)
    hScale <- 0.8
    layout(mat = layout.matrix,heights = heights) # Widths of the two columns
    div_f_d<-1.6
    bin<-1000
    {
      # a #
      agg<-rollmean(sort(species_fitness_equilibrium$agg),k=bin)
      comm_width<-5
      par(mar=c(0,comm_width,1,0))
      
      x<-rollmean(sort(log10(species_fitness_equilibrium$agg)),k=bin)
      plot(density(log10(agg),bw=0.1),xaxt="n",main="",lwd=2,col="grey",ylim=c(0,2),xaxs="i",yaxs="i",
           xlim=range(x),cex.lab=comm_cex_lab,cex.axis=comm_cex_lab)
      
      mtext("(a)",side=3,line=-2,adj=0.1,cex=comm_cex_lab)
      abline(v=mean(x),lty="dashed")
      mtext("Full model",side=3,line=-1.2,adj=0.95,cex=comm_cex_lab/div_f_d)
      par(mar=c(0,comm_width,0,0))
      
      # b #
      x<-rollmean(sort(log10(species_fitness_equilibrium$agg)),k=bin)
      xx<-log10(all$agg)
      yy<-log10(all$b_p_Axel_c)
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2 + C * xx^3,start=c(Const=1,A=1,B=1,C=1))
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2+coefs[4]*x^3
      rangee<-1
      
      minn<-min(with(species_fitness_equilibrium,log10(rollmean(fitness[order(agg)],k=bin))))
      refff<-minn
      range_x<-range(x)
      ylim=minn-0.05+c(0,hScale*heights[2])
      with(species_fitness_equilibrium,plot(rollmean(log10(sort(agg)),k=bin),xaxs="i",yaxs="i",
                                            log10(rollmean(fitness[order(agg)],k=bin)),xaxt="n",
                                            ylab=expression(log[10](R(log[10](a[i])))),
                                            type="l",cex=0.1,xlim=range_x,ylim=ylim,yaxt="n",xlab="n",
                                            cex.lab=comm_cex_lab,cex.axis=comm_cex_lab,col="grey",lwd=2))
      mtext("(b)",side=3,line=-2,adj=0.1,cex=comm_cex_lab)
      axis(2,at=c(-0.2,0,0.2),cex.axis=comm_cex_lab)
      legend("bottomright", legend=c("Simulation (parametric)","Simulation (non-parametric)","Mutation bias slope",
                                     expression(log[10](b(a[i])*L(a[i])))),
             pch=c(NA), lty=c(1,1,1,2), col=c("black","grey","blue","red"),cex=comm_cex_lab,bty="n",lwd=c(2,2,2,2))
      pred<-unlist(predict(f_nls_quadratic))
      lines(sort(log10(species_fitness_equilibrium$agg)),log10(pred),lty=1,col="black",lwd=2)
      corrr<-log10(pred)[which.min(abs(sort(log10(species_fitness_equilibrium$agg))-mean(x)))]
      mutation_bias<-log10(a1)
      diff<--mean(x)*mutation_bias/var(log10(species_fitness_equilibrium$agg))/log(10)
      x<-rollmean(sort(log10(species_fitness_equilibrium$agg)),k=bin/4)
      lines(x,-x*mutation_bias/var(log10(species_fitness_equilibrium$agg))/log(10)-diff+corrr/1.5,col="blue",lwd=2)
      x<-rollmean(sort(log10(species_fitness_equilibrium$agg)),k=bin)
      lines(x=rep(mean(x),10),seq(0,10,length.out = 10),lty="dashed")
      lines(y=rep(0,10),x=seq(min(x)-1,mean(x),length.out = 10),lty="dashed")
      
      for_nls_birth_m<-data.frame(with(for_nls_birth,rollmean(birth_rate[order(agg)],k=bin)))
      colnames(for_nls_birth_m)<-"birth_rate"
      for_nls_birth_m$log10agg<-with(for_nls_birth,rollmean(sort(log10(agg)),k=bin))
      b_nls_cubic<-nls(log10(birth_rate) ~ c+b*(log10agg)+a*(log10agg)^2+d*(log10agg)^3,for_nls_birth_m,
                       start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3],d=coef(lm_birth_agg)[3]))
      # b_nls_quadratic<-nls(log10(birth_rate) ~ c+b*(log10agg)+a*(log10agg)^2,for_nls_birth_m,
      #                      start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3]))
      paramz61<-coef(b_nls_cubic)
      prob<-cubic(x,as.numeric(paramz61[1]),as.numeric(paramz61[2]),as.numeric(paramz61[3]),as.numeric(paramz61[4]))
      xx<-with(species_fitness_equilibrium,rollmean(sort(log10(agg)),k=bin))
      yy<-with(species_fitness_equilibrium,rollmean(lifetime[order(agg)],k=bin))
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2 + C * xx^3,start=c(Const=1,A=1,B=1,C=1))  
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2+coefs[4]*x^3
      infered_fitness<-prob+log10(y)
      lines(x,infered_fitness,col="red",lty=2,lwd=2)
      
      # c #
      bin<-1000
      x<-seq(min(log10(agg)),max(log10(agg)),length.out = 1000)
      for_nls_birth_m<-data.frame(with(for_nls_birth,rollmean(birth_rate[order(agg)],k=bin)))
      colnames(for_nls_birth_m)<-"birth_rate"
      for_nls_birth_m$log10agg<-with(for_nls_birth,rollmean(sort(log10(agg)),k=bin))
      b_nls_cubic<-nls(log10(birth_rate) ~ c+b*(log10agg)+a*(log10agg)^2+d*(log10agg)^3,for_nls_birth_m,
                       start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3],d=coef(lm_birth_agg)[3]))
      paramz61<-coef(b_nls_cubic)
      prob<-cubic(x,as.numeric(paramz61[1]),as.numeric(paramz61[2]),as.numeric(paramz61[3]),as.numeric(paramz61[4]))
      
      minn<-min(prob)
      ylim=minn-0.05+c(0,hScale*heights[3])
      
      with(for_nls_birth,plot(rollmean(log10(sort(agg)),k=bin),col="grey",
                              log10(rollmean(birth_rate[order(agg)],k=bin)),
                              type="l",ylab=expression(log[10](b(log[10](a[i])))),xaxt="n",
                              xlim=range_x,cex.lab=comm_cex_lab,lwd=2,xaxs="i",yaxs="i",
                              cex.axis=comm_cex_lab,ylim=ylim))
      lines(x,prob,lwd=2)
      xx<-log10(all$agg)
      yy<-log10(all$b_p_Axel)
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2 + C * xx^3 ,start=c(Const=1,A=1,B=1,C=1)) #+ C * xx^3 ,C=1
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2+coefs[4]*x^3
      lines(x, y,lwd=2.5,col="green")
      
      xx<-log10(all$agg)
      yy<-log10(all$b_p_Axel_c)
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2+ C * xx^3 ,start=c(Const=1,A=1,B=1,C=1))  
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2+coefs[4]*x^3
      lines(x, y,lwd=2.5,col="green",lty="dashed")
      legend("bottomright", legend=c("Prediction","Prediction (corrected)"),
             pch=c(NA), lty=c(1,2), col=c("green","green"),cex=comm_cex_lab,bty="n",lwd=c(2),seg.len = 4)
      mtext("(c)",side=3,line=-2,adj=0.1,cex=comm_cex_lab)
      
      # d #
      par(mar=c(0,comm_width,0,0))
      xx<-with(species_fitness_equilibrium,rollmean(sort(log10(agg)),k=bin))
      yy<-with(species_fitness_equilibrium,rollmean(lifetime[order(agg)],k=bin))
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2+ C * xx^3 ,start=c(Const=1,A=1,B=1,C=1)) 
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2+coefs[4]*x^3
      minn<-min(log10(y))
      ylim=minn+c(0,hScale*heights[4])
      with(species_fitness_equilibrium,plot(rollmean(log10(sort(agg)),k=bin),col="grey",
                                            log10(rollmean(lifetime[order(agg)],k=bin)),
                                            lwd = 2,cex=0.1,xaxs="i",yaxs="i",
                                            type="l",ylab=expression(log[10](L(log[10](a[i])))),main="",yaxt="n",
                                            xlab=expression(log10(a[i])),xlim = range_x, ylim=ylim,cex.lab=comm_cex_lab,
                                            cex.axis=comm_cex_lab))
      lines(x, log10(y), lwd = 2,cex=0.1,xaxs="i",yaxs="i")
      mtext("(d)",side=1,line=-2,adj=0.1,cex=comm_cex_lab)
      axis(2,at=c(2.1,2.3),cex.lab=comm_cex_lab,cex.axis=comm_cex_lab)
    }
    
    plot.new()
    par(mar=c(4,4,2,1))
    mtext(expression(log[10](a[i])),side=1,cex=comm_cex_lab/1.3,line=1.3)
    load("Pre_Processed/Data_deconstructed.RData")
    {
      agg<-rollmean(sort(species_fitness_equilibrium$agg),k=bin)
      par(mar=c(0,4,1,1))
      x<-rollmean(sort(log10(species_fitness_equilibrium$agg)),k=bin)
      plot(density(log10(agg),bw=0.1),xaxt="n",main="",lwd=2,col="grey",ylim=c(0,2),xaxs="i",yaxs="i",
           xlim=range(x),cex.lab=comm_cex_lab,cex.axis=comm_cex_lab,ylab="")
      mtext("(e)",side=3,line=-2,adj=0.1,cex=comm_cex_lab)
      abline(v=mean(x),lty="dashed")
      mtext("Deconstructed model",side=3,line=-1.2,adj=0.95,cex=comm_cex_lab/div_f_d)
      par(mar=c(0,4,0,1))
      
      minn<-min(with(species_fitness_equilibrium,log10(rollmean(fitness[order(agg)],k=bin))))
      refff<-minn
      range_x<-range(x)
      ylim=minn-0.05+c(0,hScale*heights[2])
      with(species_fitness_equilibrium,plot(rollmean(log10(sort(agg)),k=bin),xaxs="i",yaxs="i",
                                            log10(rollmean(fitness[order(agg)],k=bin)),xaxt="n",ylab="",
                                            type="l",cex=0.1,xlim=range_x,ylim=ylim,yaxt="n",xlab="n",
                                            cex.lab=comm_cex_lab,cex.axis=comm_cex_lab,col="grey",lwd=2))
      mtext("(f)",side=3,line=-2,adj=0.1,cex=comm_cex_lab)
      axis(2,at=c(-0.4,-0.2,0,0.2),cex.axis=comm_cex_lab)
      rem_plz<-which(predict(f_nls_quadratic)<0)
      pred<-unlist(predict(f_nls_quadratic)[-rem_plz])
      lines(sort(log10(species_fitness_equilibrium$agg[-rem_plz])),log10(pred),lty=1,col="black",lwd=2)
      
      for_nls_birth_m<-data.frame(with(for_nls_birth,rollmean(birth_rate[order(agg)],k=bin)))
      colnames(for_nls_birth_m)<-"birth_rate"
      for_nls_birth_m$log10agg<-with(for_nls_birth,rollmean(sort(log10(agg)),k=bin))
      b_nls_cubic<-nls(log10(birth_rate) ~ c+b*(log10agg)+a*(log10agg)^2+d*(log10agg)^3,for_nls_birth_m,
                       start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3],d=coef(lm_birth_agg)[3]))
      b_nls_quadratic<-nls(log10(birth_rate) ~ c+b*(log10agg)+a*(log10agg)^2,for_nls_birth_m,
                           start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3]))
      paramz61<-coef(b_nls_quadratic)
      prob<-cubic(x,as.numeric(paramz61[1]),as.numeric(paramz61[2]),as.numeric(paramz61[3]),0)
      xx<-with(species_fitness_equilibrium,rollmean(sort(log10(agg)),k=bin))
      yy<-with(species_fitness_equilibrium,rollmean(lifetime[order(agg)],k=bin))
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2,start=c(Const=1,A=1,B=1)) 
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2
      infered_fitness<-prob+log10(y)
      lines(x,infered_fitness,col="red",lty=2,lwd=2)
      
      corrr<-log10(pred)[which.min(abs(sort(log10(species_fitness_equilibrium$agg))-mean(x)))]
      mutation_bias<-log10(a1)
      diff<--mean(x)*mutation_bias/var(log10(species_fitness_equilibrium$agg))/log(10)
      x<-rollmean(sort(log10(species_fitness_equilibrium$agg)),k=bin/4)
      lines(x,-x*mutation_bias/var(log10(species_fitness_equilibrium$agg))/log(10)-diff+corrr/1.5,col="blue",lwd=2)
      x<-rollmean(sort(log10(species_fitness_equilibrium$agg)),k=bin)
      lines(x=rep(mean(x),10),seq(0,10,length.out = 10),lty="dashed")
      lines(y=rep(0,10),x=seq(min(x)-1,mean(x),length.out = 10),lty="dashed")
      
      bin<-1000
      x<-seq(min(log10(agg)),max(log10(agg)),length.out = 1000)
      for_nls_birth_m<-data.frame(with(for_nls_birth,rollmean(birth_rate[order(agg)],k=bin)))
      colnames(for_nls_birth_m)<-"birth_rate"
      for_nls_birth_m$log10agg<-with(for_nls_birth,rollmean(sort(log10(agg)),k=bin))
      b_nls_cubic<-nls(log10(birth_rate) ~ c+b*(log10agg)+a*(log10agg)^2+d*(log10agg)^3,for_nls_birth_m,
                       start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3],d=coef(lm_birth_agg)[3]))
      b_nls_quadratic<-nls(log10(birth_rate) ~ c+b*(log10agg)+a*(log10agg)^2,for_nls_birth_m,
                           start=list(c=coef(lm_birth_agg)[1],b=coef(lm_birth_agg)[2],a=coef(lm_birth_agg)[3]))
      
      # paramz61<-coef(b_nls_quadratic)
      # prob<-cubic(x,as.numeric(paramz61[1]),as.numeric(paramz61[2]),as.numeric(paramz61[3]),0)
      paramz61<-coef(b_nls_cubic)
      prob<-cubic(x,as.numeric(paramz61[1]),as.numeric(paramz61[2]),as.numeric(paramz61[3]),as.numeric(paramz61[4]))
      minn<-min(prob)
      ylim=minn-0.05+c(0,hScale*heights[3])
      with(for_nls_birth,plot(rollmean(log10(sort(agg)),k=bin),col="grey",
                              log10(rollmean(birth_rate[order(agg)],k=bin)),
                              type="l",ylab="",xaxt="n",
                              xlim=range_x,cex.lab=comm_cex_lab,lwd=2,xaxs="i",yaxs="i",
                              cex.axis=comm_cex_lab,ylim=ylim))
      lines(x,prob,lwd=2)
      
      xx<-log10(all$agg)
      yy<-log10(all$b_p_Axel)
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2  + C * xx^3,start=c(Const=1,A=1,B=1,C=1))  
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2+coefs[4]*x^3
      lines(x, y,lwd=2.5,col="green")
      
      xx<-log10(all$agg)
      yy<-log10(all$b_p_Axel_c)
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2 + C * xx^3 ,start=c(Const=1,A=1,B=1,C=1)) #
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2+coefs[4]*x^3
      lines(x, y,lwd=2.5,col="green",lty="dashed")
      mtext("(g)",side=3,line=-2,adj=0.1,cex=comm_cex_lab)
      
      par(mar=c(0,4,0,1))
      xx<-with(species_fitness_equilibrium,rollmean(sort(log10(agg)),k=bin))
      yy<-with(species_fitness_equilibrium,rollmean(lifetime[order(agg)],k=bin))
      nlmod <- nls(yy ~  Const + A * xx + B * xx^2,start=c(Const=1,A=1,B=1)) # 
      coefs<-coef(nlmod)
      y<-coefs[1]+coefs[2]*x+coefs[3]*x^2#+coefs[4]*x^3
      
      minn<-min(log10(y))
      ylim=minn+c(0,hScale*heights[4])
      with(species_fitness_equilibrium,plot(rollmean(log10(sort(agg)),k=bin),col="grey",
                                            log10(rollmean(lifetime[order(agg)],k=bin)),
                                            lwd = 2,cex=0.1,xaxs="i",yaxs="i",
                                            type="l",ylab="",main="",yaxt="n",
                                            xlab=expression(log10(a[i])),xlim = range_x, ylim=ylim,cex.lab=comm_cex_lab,
                                            cex.axis=comm_cex_lab))
      lines(x, log10(y), lwd = 2,cex=0.1,xaxs="i",yaxs="i")    
      mtext("(h)",side=1,line=-2,adj=0.1,cex=comm_cex_lab)
      axis(2,at=c(1.9,2.1,2.3),cex.lab=comm_cex_lab,cex.axis=comm_cex_lab)
    }
    
    dev.off()},silent=T)
}

# Figure 3.2: Please go to Invasion_Probability.R

# Figure 3.3: Please go to Invasion_Probability.R

# Figure 3.4: Please go to Invasion_Probability.R

# Figure 3.5 #
{
  load("Pre_Processed/Data_deconstructed.RData") 
  species_fitness_equilibrium$l_agg<-log10(species_fitness_equilibrium$agg)
  mean_l10_a<-mean(species_fitness_equilibrium$l_agg);sd_l10_a<-var(species_fitness_equilibrium$l_agg)^0.5
  min_l10_agg<-mean_l10_a-1.96*sd_l10_a;max_l10_agg<-mean_l10_a+1.96*sd_l10_a
  increment_percentage_of_range <- 0.05
  buffer<-5
  agg<-species_fitness_equilibrium$agg
  bins<-seq(min_l10_agg,max_l10_agg,length.out=(1/increment_percentage_of_range+1)) # 
  size<-bins[2]-bins[1]
  
  respp<-unlist(llply(bins,function(x) sum(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),]$extinction == 2)))
  spec<-unlist(llply(bins,function(x) sum(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),]$extinction == 1)))
  comp<-unlist(llply(bins,function(x) sum(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),]$extinction %in% c(61,62,63))))
  over<-unlist(llply(bins,function(x) sum(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),]$extinction == 3)))
  
  # Invader kills residents
  comp_61<-unlist(llply(bins,function(x) with(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),],sum(extinction == 61))))
  # Resident kills resident (after resource invasion)
  comp_62<-unlist(llply(bins,function(x) with(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),],sum(extinction == 62 &
                                                                                                                                                                    lifetime != 0  ))))
  
  comp_62_l0<-unlist(llply(bins,function(x) with(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),],sum(extinction == 62 &
                                                                                                                                                                       lifetime == 0  ))))
  # Resident kills resident (after consumer invasion)
  comp_63<-unlist(llply(bins,function(x) with(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),],sum(extinction == 63 &
                                                                                                                                                                    lifetime != 0  ))))
  comp_63_l0<-unlist(llply(bins,function(x) with(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),],sum(extinction == 63 &
                                                                                                                                                                       lifetime == 0  ))))
  
  comp_62_63<-unlist(llply(bins,function(x) with(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),],
                                                 sum(extinction %in% c(62,63)))))
  
  comp_62_63_l0<-unlist(llply(bins,function(x) with(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),],
                                                    sum(extinction %in% c(62,63) & lifetime == 0  ))))
  
  lifetimes<-unlist(llply(bins,function(x) sum(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),]$lifetime)))+1
  
  all<-unlist(llply(bins,function(x) sum(species_fitness_equilibrium[which(species_fitness_equilibrium$l_agg>x & species_fitness_equilibrium$l_agg<x+size),]$extinction !=0 )))

  mult<-0.5 
  leg.size<-0.5
  sizee<-0.5
  inc_t<-0.2
  pdf("../Images/activation_frequency.pdf",height=10*0.5*mult,width=10*0.5*mult)
  leg.pos<-c(0.8,0.71)
  totall<-spec/lifetimes+respp/lifetimes+over/lifetimes+comp/lifetimes
  uspopage_prop<-data.frame(spec/(lifetimes*totall));colnames(uspopage_prop)<-"Val";uspopage_prop$type<-"(3) Pyrrhic"
  uspopage_prop$type2<-4
  tempp<-data.frame(respp/(lifetimes*totall));colnames(tempp)<-"Val";tempp$type<-"(2) Respiration";tempp$type2<-1
  uspopage_prop<-rbind(uspopage_prop,tempp)
  tempp<-data.frame(over/(lifetimes*totall));colnames(tempp)<-"Val";tempp$type<-"(4) Overexploitation";tempp$type2<-3
  uspopage_prop<-rbind(uspopage_prop,tempp)
  tempp<-data.frame(comp/(lifetimes*totall));colnames(tempp)<-"Val";tempp$type<-"(1) Exploitative";tempp$type2<-2
  uspopage_prop<-rbind(uspopage_prop,tempp)
  uspopage_prop$l_agg<-bins
  
  ggplot(uspopage_prop[order(uspopage_prop$type2),], aes(x=l_agg, y=Val)) + geom_area(colour="black", size=.2, alpha=.4, aes(fill=type)) +
    scale_fill_brewer(palette="Blues",name="Extinction\nMechanism") +
    labs(x = expression(log(a[i])), y = "Frequency of activation (%)" ) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       axis.title.y = element_text(size = rel(1.35*mult), angle = 90),
                       axis.text.y = element_text(size = rel(1.35*mult), angle = 0,color="black"),
                       axis.title.x = element_text(size = rel(1.35*mult), angle = 0),
                       axis.text.x = element_text(size = rel(1.35*mult), angle = 0,color="black"),
                       legend.title = element_text(size=rel(inc_t+leg.size*mult)),
                       legend.text = element_text(size = rel(inc_t+leg.size*mult)),legend.position =leg.pos,
                       legend.key.size = unit(sizee, "cm"),
                       legend.background = element_rect(fill="NA",size=0.5,colour ="NA"))+
                       scale_x_continuous( expand=c(0,0)) + 
    scale_y_continuous( expand=c(0,0),breaks=c(0,0.25,0.5,0.75,1),labels=c("0","25","50","75","100"))
  dev.off()
}

# Figure 3.6 #
{
  load("Pre_Processed/Data_deconstructed.RData") 
  load("Pre_Processed/Numerical_serial_extinction_sub_model/Numerical_serial_extinction_sub_model.RData")
  try({pdf("../Images/components_with_without_serial_scaled_older.pdf",width=7.5,height=7.5)
    points_to_plot<-50
    cut_of<--10
    comm_cex_lab<-1.3
    agg_1<-10^-5.5;agg_2<-10^-4.75;agg_3<-10^-4;
    {
      par(mfrow=c(2,2))
      par(mar=c(0,4.5,4,0))
      bww<-1
      blue_red<-colorRampPalette(c("blue", "red"))
      colz<-blue_red(3)
      
      
      maxi<-max(with(all_gran_table_com,lty="dashed",
                     -density(log10(A2_com),bw=1.5)$x))
      mini<-min(with(all_gran_table_com,lty="dashed",
                     -density(log10(A2_com),bw=1)$x))
      yrange<-range(log10(density(log10(all_gran_table_com$ine_com))$y))
      yrange<-c(yrange[1]*0.25,yrange[2]*0.25)
      all_gran_table_com$l10_agg <- log10(exp(all_gran_table_com$l_agg))
      cut_offf<-0
      cur<-all_gran_table_com[which(all_gran_table_com$agg<agg_1),]
      cur_dens<-data.frame(density(log10(cur$ine_com),bw=bww)$x);colnames(cur_dens)<-"x"
      cur_dens$y<-density(log10(cur$ine_com),bw=bww)$y
      cur_dens[which(cur_dens$x<0),]$y<-exp(cut_of)
      with(cur_dens,plot(x,log10(y),col=colz[1],type="l",ylim=yrange,yaxs="i",xaxs="i",
                         main= "",cex.axis=1.2,cex.lab=1.2,#signif(c(mean(log10(ine_com)),var(log10(ine_com))),2),
                         xlab="",xlim=c(mini,maxi),xaxt="n",ylab=expression(log[10](Density))))
      mtext("(a)",side=3,line=-2,adj=0,cex=comm_cex_lab)
      axis(3,cex.axis=1.2,cex.lab=1.2,at=c(-4,-2,0,2,4,6,8))
      cur<-all_gran_table_com[which(all_gran_table_com$agg>agg_1 & 
                                      all_gran_table_com$agg<agg_2),]
      cur_dens<-data.frame(density(log10(cur$ine_com),bw=bww)$x);colnames(cur_dens)<-"x"
      cur_dens$y<-density(log10(cur$ine_com),bw=bww)$y
      cur_dens[which(cur_dens$x<0),]$y<-exp(cut_of)
      with(cur_dens,lines(x,log10(y),col=colz[2],xlim=c(mini,maxi)))
      
      cur<-all_gran_table_com[which(all_gran_table_com$agg>agg_2 & 
                                      all_gran_table_com$agg<agg_3),]
      cur_dens<-data.frame(density(log10(cur$ine_com),bw=bww)$x);colnames(cur_dens)<-"x"
      cur_dens$y<-density(log10(cur$ine_com),bw=bww)$y
      cur_dens[which(cur_dens$x<0),]$y<-exp(cut_of)
      with(cur_dens,lines(x,log10(y),col=colz[3],xlim=c(mini,maxi)))
      
      # cur<-all_gran_table_com[which(all_gran_table_com$l_agg>(range(plot_x)[1]+cc*3) &
      #                                 all_gran_table_com$l_agg<(range(plot_x)[1]+cc*4)),]
      # cur_dens<-data.frame(density(log10(cur$ine_com),bw=bww)$x);colnames(cur_dens)<-"x"
      # cur_dens$y<-density(log10(cur$ine_com),bw=bww)$y
      # cur_dens[which(cur_dens$x<0),]$y<-exp(cut_of)
      # with(cur_dens,lines(x,log10(y),col=colz[4],xlim=c(mini,maxi)))
      
      mtext(text = expression(log[10](C)-log[10](D)-log[10](A)+log[10](B)),
            side = 3,line = 2,cex=0.75)
      
      par(mar=c(0,0,4,4))
      cut_offf<-0
      cur<-all_gran_table_com[which(all_gran_table_com$agg<agg_1),]
      cur_dens<-density(log10(cur$A1_com-1),bw=bww)
      with(cur_dens,plot(x,log10(y),col=colz[1],type="l",ylim=yrange,yaxs="i",xaxs="i",
                         main= "",cex.axis=1.2,cex.lab=1.2,#signif(c(mean(log10(A1_com-1)),var(log10(A1_com-1))),2),
                         xlab="",xlim=c(mini,maxi),xaxt="n",yaxt="n"))
      mtext("(b)",side=3,line=-2,adj=0,cex=comm_cex_lab)
      axis(3,cex.axis=1.2,cex.lab=1.2)
      cur<-all_gran_table_com[which(all_gran_table_com$agg>agg_1 & 
                                      all_gran_table_com$agg<agg_2),]
      cur_dens<-density(log10(cur$A1_com-1),bw=bww)
      with(cur_dens,lines(x,log10(y),col=colz[2],xlim=c(mini,maxi)))
      
      cur<-all_gran_table_com[which(all_gran_table_com$agg>agg_2 & 
                                      all_gran_table_com$agg<agg_3),]
      cur_dens<-density(log10(cur$A1_com-1),bw=bww)
      with(cur_dens,lines(x,log10(y),col=colz[3],xlim=c(mini,maxi)))
      
      # cur<-all_gran_table_com[which(all_gran_table_com$l_agg>(range(plot_x)[1]+cc*3) & 
      #                                 all_gran_table_com$l_agg<(range(plot_x)[1]+cc*4)),]
      # cur_dens<-density(log10(cur$A1_com-1),bw=bww)
      # with(cur_dens,lines(x,log10(y),col=colz[4],xlim=c(mini,maxi)))
      axis(4,cex.axis=1.2,cex.lab=1.2)
      mtext(text = expression(log[10](C)),
            side = 3,line = 2,cex=1)
      
      par(mar=c(4,4.5,0,0))
      
      cut_offf<-0
      cur<-all_gran_table_com[which(all_gran_table_com$agg<agg_1),]
      cur_dens<-density(-log10(cur$A2_com),bw=bww)
      with(cur_dens,plot(x,log10(y),col=colz[1],type="l",ylim=yrange,yaxs="i",xaxs="i",yaxt="n",xaxt="n",
                         main= "",cex.axis=1.2,cex.lab=1.2,#signif(c(mean(log10(A2_com)),var(log10(A2_com))),2),
                         xlab=expression(-log[10](D)),xlim=c(mini,maxi),ylab=expression(log[10](Density)),yaxt="n"))
      mtext("(c)",side=3,line=-2,adj=0,cex=comm_cex_lab)
      axis(2,cex.axis=1.2,cex.lab=1.2,at=c(-3.25,-2.25,-1.25))
      axis(1,cex.axis=1.2,cex.lab=1.2,at=c(-4,-2,0,2,4,6,8))
      cur<-all_gran_table_com[which(all_gran_table_com$agg>agg_1 & 
                                      all_gran_table_com$agg<agg_2),]
      cur_dens<-density(-log10(cur$A2_com),bw=bww)
      with(cur_dens,lines(x,log10(y),col=colz[2],xlim=c(mini,maxi)))
      
      cur<-all_gran_table_com[which(all_gran_table_com$agg>agg_2 & 
                                      all_gran_table_com$agg<agg_3),]
      cur_dens<-density(-log10(cur$A2_com),bw=bww)
      with(cur_dens,lines(x,log10(y),col=colz[3],xlim=c(mini,maxi)))
      
      # cur<-all_gran_table_com[which(all_gran_table_com$l_agg>(range(plot_x)[1]+cc*3) & 
      #                                 all_gran_table_com$l_agg<(range(plot_x)[1]+cc*4)),]
      # cur_dens<-density(-log10(cur$A2_com),bw=bww)
      # with(cur_dens,lines(x,log10(y),col=colz[4],xlim=c(mini,maxi)))
      
      blue_red<-colorRampPalette(c("blue", "red"))
      # colz<-blue_red(3)
      par(mar=c(4,0,0,4))
      
      maxi<-max(links_1_c_mm[1],links_2_c_mm[1],links_3_c_mm[1])
      mini<-min(log10(2-links_1_DD_b_mm),log10(2-links_2_DD_b_mm),log10(2-links_3_DD_b_mm))
      number_adds<-length(links_1_c_mm)
      
      minnn<-which.min(c(agg_1,agg_2,agg_3))
      maxxx<-which.max(c(agg_1,agg_2,agg_3))
      
      sample_n<-100
      plot(log10(1:number_adds),links_1_c_mm,col=colz[1],cex=0.1,type="l",cex.axis=1.2,cex.lab=1.2,xaxs="i",yaxs="i",
           ylim=c(-1,1.5),xlab="Number of potential serial extinctions",ylab="",yaxt="n",xaxt="n")
      mtext("(d)",side=3,line=-2,adj=0.1,cex=comm_cex_lab)
      # abline(h=log10(1-links_1_DD_b_mm)[1],cex=0.1,col=colz[1],lty=3)
      axis(4,cex.axis=1.2,cex.lab=1.2,at=c(0,1,2,3))
      axis(1,cex.axis=1.2,cex.lab=1.2,at=log10(c(1,2,5,10,50,100,250)),labels=c("Invasion",1,5,10,50,100,250))
      mtext(text = expression(log[10](C)),
            side = 4,line = 2.5,cex=1)
      
      lines(log10(1:number_adds),links_2_c_mm,cex=0.1,col=colz[2],type="l")
      # abline(h=log10(1-links_2_DD_b_mm)[1],cex=0.1,col=colz[2],lty=3)
      
      lines(log10(1:number_adds),links_3_c_mm,cex=0.1,type="l",col=colz[3])
      # abline(h=log10(1-links_3_DD_b_mm)[1],cex=0.1,col=colz[3],lty=3)
      
      lines(log10(1:number_adds),links_1_c_mmm,col=colz[1],lty=2)
      lines(log10(1:number_adds),links_2_c_mmm,cex=0.1,col=colz[2],type="l",lty=2)
      lines(log10(1:number_adds),links_3_c_mmm,cex=0.1,type="l",col=colz[3],lty=2)
      # legend("topright",legend=c("High aggressiveness","Low aggressiveness","Deconstructed model","Full model"),
      #        col=c("red","blue","black","black"),lty=c(1,1,1,2),bty="n")
      legend("topright",
             legend=c(expression(a[i]==10^-5.5),expression(a[i]==10^-4.75),expression(a[i]==10^-4),
                      "Deconstructed model","Full model"),
             col=c(colz[3],colz[2],colz[1],"black","black"),lty=c(1,1,1,1,2),bty="n")
      
    }
    dev.off()},silent=T)
}

# Figure 3.7 #
{
load("Pre_Processed/Data_deconstructed.RData")
load("Pre_Processed/Numerical_serial_extinction_sub_model/Numerical_serial_extinction_sub_model.RData") #O
links_3_AB_D <- 10^(log10(links_3_ine)-log10(links_3_cc-1))
links_3_AB_D_m <- rowMeans(log10(links_3_AB_D))
links_1_AB_D <- 10^(log10(links_1_ine)-log10(links_1_cc-1))
links_1_AB_D_m <- rowMeans(log10(links_1_AB_D))

maxi<-max(with(all_gran_table_com,lty="dashed",
               -density(log10(A2_com),bw=1.5)$x))
mini<-min(with(all_gran_table_com,lty="dashed",
               -density(log10(A2_com),bw=1)$x))

library(RColorBrewer)
mult<-0.6
comm_cex_lab<-1.3
comm_cexx<-3*mult
mtext_cex<-2.3*mult
pdf("../Images/variation.pdf",width=8*1.5,height=5.5*1.5)
print("always run this twice for some reason")
colz_colour_blind<-brewer.pal(10, "Dark2")
points_to_plot<-40
sample_n<-10
randd<-105#110
thik<-2
lab_size<-comm_cexx
point_displace<-1.5
opac_points<-0.3
opac_random<-1
mean_l_type<-1
random_l_type<-5
col_D<-colz_colour_blind[1]
col_C<-colz_colour_blind[2]
col_AB<-colz_colour_blind[3]
col_ine<-colz_colour_blind[4]
col_AB_D<-colz_colour_blind[5]
legend_size<-2
layout.matrix <- matrix(c(1, 2, 3, 4,5,3,6,6,7), nrow = 3, ncol = 3,byrow = T)
layout(mat = layout.matrix,widths = c(2,2,1.5),height=c(3.5,3.5,0.75)) # Widths of the two columns
sample_pp<-unique(round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))
sample_pp_diff<-diff(sample_pp)/4
col_AB_D<-colz_colour_blind[8]
plot_until<- 200
axis_lab<-c(1,10,100)
cex.pp=0.4
{
  layout.matrix <- matrix(c(1, 2, 3, 4,5,3,6,6,7), nrow = 3, ncol = 3,byrow = T)
  layout(mat = layout.matrix,widths = c(2,2,1.5),height=c(3.5,3.5,0.75)) # Widths of the two columns
  divv<-20
  {
    par(mar=c(0,4,4,0))
    number_adds<-dim(links_3_cc2)[1]
    
    whoo<-sample(1:dim(links_3_cc_2)[2],sample_n)
    
    plot(log10(rep(1.2,sample_n)),-log10((links_3_cc_2[1,whoo]-1)/links_3_cc2_2[1,whoo]),xaxs="i",
         col=alpha("white",opac_points),cex=cex.pp, #
         ylim=c(-1,4),xlim=c(0,log10(plot_until)),yaxt="n",
         xaxt="n",ylab="",xlab="",main=expression(log[10](a[i])==-4.25),cex.main=comm_cexx) #xlab=expression(log[10](Number~of~resource~additions))
    axis(2,ylab="",at=c(0,1,2,3,4,5),cex.axis=lab_size,cex.lab=lab_size)
    for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
      points(log10(rep(i,sample_n)),-log10((links_3_cc_2[i,whoo]-1)/links_3_cc2_2[i,whoo]),col=alpha("white",opac_points),cex=cex.pp)
    }
    
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),lty=mean_l_type,
          -links_3_AB_m[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],
          cex=0.1,col=col_C,lwd=thik)
    # lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),
    #       -log10(links_3_AB[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]),
    #       col=alpha(col_C,opac_random),cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          -log10(links_3_AB[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]),
          col=alpha(col_C,opac_random),cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    
    points(log10(rep(1.2,sample_n))*0.66,-log10(links_3_cc2[1,whoo]),col=alpha("white",opac_points),cex=cex.pp,
           ylim=c(-1,5),xlim=c(0,2.5),xlab=expression(log[10](Number~of~resource~additions)))
    for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
      points(log10(rep(i,sample_n)),-log10(links_3_cc2[i,whoo]),col=alpha("white",opac_points),cex=cex.pp)
    }
    
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),lty=mean_l_type,
          -links_3_c2_mm[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],
          cex=0.1,col=col_D,lwd=thik)
    
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          -log10(links_3_cc2[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]),
          col=alpha(col_D,opac_random),cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    
    # abline(h=-log10(median(links_3_cc2)))
    
    
    points(log10(rep(1.2,sample_n))*point_displace,log10(links_3_cc[1,whoo]-1),col=alpha("white",opac_points),cex=cex.pp)
    for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
      points(log10(rep(i,sample_n)),log10(links_3_cc[i,whoo]-1),col=alpha("white",opac_points),cex=cex.pp)
    }
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),col=alpha(col_AB,opac_random),
          log10(links_3_cc[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]-0.5),
          cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    lines(log10(c(1.2,2:number_adds)),links_3_c_mm,cex=0.1,lwd=thik,lty=mean_l_type,col=col_AB)
    # abline(h=log10(median(links_3_cc-1)))
  }
  mtext("(a)",side=1,line=-2,adj=0,cex=comm_cex_lab)
  
  {
    # par(mfrow=c(2,2))
    par(mar=c(0,0,4,4))
    number_adds<-dim(links_1_cc2)[1]
    
    plot(log10(rep(1.2,sample_n)),-log10((links_1_cc_2[1,whoo]-1)/links_1_cc2_2[1,whoo]),col=alpha("white",opac_points),cex=cex.pp,xaxs="i",
         ylim=c(-1,4),xlim=c(0,log10(plot_until)),yaxt="n",xaxt="n",ylab="",xlab="",
         main=expression(log[10](a[i])==-5.25),cex.main=comm_cexx) #xlab=expression(log[10](Number~of~resource~additions))
    # for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
    #   points(log10(rep(i,sample_n)),-log10((links_1_cc_2[i,whoo]-1)/links_1_cc2_2[i,whoo]),col=alpha(col_C,opac_points),cex=cex.pp)
    # }
    # 
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),
          -links_1_AB_m[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],lty=mean_l_type,
          cex=0.1,col=col_C,lwd=thik)
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          -log10(links_1_AB[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]),
          col=alpha(col_C,opac_random),cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    # abline(h=-log10(median((links_1_cc_2[1,]-1)/links_1_cc2_2[1,])))
    
    # points(log10(rep(1.2,sample_n))*0.66,-log10(links_1_cc2[1,whoo]),col=alpha(col_D,opac_points),cex=cex.pp,
    #        ylim=c(-1,5),xlim=c(0,2.5),xlab=expression(log[10](Number~of~resource~additions)))
    # for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
    #   points(log10(rep(i,sample_n)),
    #          -log10(links_1_cc2[i,whoo]),col=alpha(col_D,opac_points),cex=cex.pp)
    # }
    
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),lty=mean_l_type,
          -links_1_c2_mm[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],lwd=thik,
          cex=0.1,col=col_D)
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          -log10(links_1_cc2[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]),
          col=alpha(col_D,opac_random),cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    
    # abline(h=-log10(median(links_1_cc2)))
    
    
    # points(log10(rep(1.2,sample_n))*point_displace,log10(links_1_cc[1,whoo]-1),col=alpha(col_AB,opac_points),cex=cex.pp)
    # for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
    #   points(log10(rep(i,sample_n)),log10(links_1_cc[1,whoo]-1),col=alpha(col_AB,opac_points),cex=cex.pp)
    # }
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          log10(links_1_cc[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]-0.5),
          col=alpha(col_AB,opac_random),cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    lines(log10(c(1.2,2:number_adds)),links_1_c_mm,cex=0.1,col=col_AB,lwd=thik,lty=mean_l_type)
    # abline(h=log10(median(links_1_cc-1)))
  }
  mtext("(b)",side=1,line=-2,adj=0,cex=comm_cex_lab)
  
  plot.new()
  par(mar=c(0,0,0,0))
  legend("center",legend=c("Mean","","Random Consumer","","Sampled Distribution","",expression(log[10](C/D)-log[10](A/B)),
                           "",expression(log[10](C)),"",expression(-log[10](D)),"",expression(log[10](A/B)),"",expression(log[10](AD/B)),
                           "","Extinction threshold"),
         col=c("black",NA,"black",NA,"black",NA,col_ine,NA,col_AB,NA,col_D,NA,col_C,NA,col_AB_D,NA,"black"),
         pch=c(NA,NA,NA,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),lty=c(mean_l_type,NA,random_l_type,NA,NA,NA,1,NA,1,NA,1,NA,1,NA,1,NA,4),bty="n",cex=legend_size)
  

  {
    par(mar=c(4,4,0,0))
    plott<-log10(links_3_ine)[,whoo]
  
    
   
    plot(log10(rep(1.2,sample_n)),plott[1,],col=alpha("white",opac_points),cex=cex.pp,xlim=c(0,log10(plot_until)),ylim=c(-1,5),xaxt="n",xaxs="i",xaxt="n",ylab="",xlab="",
         lwd=thik,cex.axis=lab_size,cex.lab=lab_size)
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),type="l",
         links_3_AB_D_m[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],lty=mean_l_type,
         cex=0.1,col=col_AB_D,lwd=thik)

    # for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
    #   points(log10(rep(i,sample_n)),plott[i,],col=alpha(col_ine,opac_points),cex=cex.pp)
    # }
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),
          links_3_ine_m[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],lty=mean_l_type,
          cex=0.1,col=col_ine,lwd=thik)
    abline(h=0,lty=4)
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          log10(links_3_ine[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]+0.5),
          cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",col=alpha(col_ine,opac_random))
    par(mar=c(4,0,0,4))
    axis(1,cex.axis=lab_size,cex.lab=lab_size,at=log10(axis_lab),labels=axis_lab)
     

    plottt<-data.frame(log10(links_3_AB_D)[,whoo])
    # points(log10(rep(1.2,sample_n)),plottt[1,],col=alpha(col_AB_D,opac_points),cex=cex.pp,xlim=c(0,2.5),ylim=c(-1,5),xaxt="n",xaxs="i",xaxt="n",ylab="",xlab="",
    #      lwd=thik,cex.axis=lab_size,cex.lab=lab_size)
    # for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
    #   points(log10(rep(i,sample_n)),plottt[i,],col=alpha(col_AB_D,opac_points),cex=cex.pp)
    # }
    
    abline(h=0,lty=4)
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          log10(links_3_AB_D[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]+0.5),
          cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",col=alpha(col_AB_D,opac_random))
    par(mar=c(4,0,0,4))
    axis(1,cex.axis=lab_size,cex.lab=lab_size,at=log10(axis_lab),labels=axis_lab)
    mtext("(c)",side=1,line=-2,adj=0,cex=comm_cex_lab) 
    
    plott<-log10(links_1_ine)[,whoo]
    plot(log10(rep(1.2,sample_n)),plott[1,],col=alpha("white",opac_points),cex=cex.pp,xlim=c(0,log10(plot_until)),ylim=c(-1,5),xaxt="n",xaxs="i",xaxt="n",ylab="",xlab="",
         lwd=thik,cex.axis=lab_size,cex.lab=lab_size,yaxt="n")
    # for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
    #   points(log10(rep(i,sample_n)),plott[i,],col=alpha(col_ine,opac_points),cex=cex.pp)
    # }
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),
          links_1_ine_m[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],lty=mean_l_type,
          cex=0.1,col=col_ine,lwd=thik)
    
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),col=alpha(col_ine,opac_random),
          log10(links_1_ine[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]),
          cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",yaxt="n")
    abline(h=0,lty=4)
    axis(1,cex.axis=lab_size,cex.lab=lab_size,at=log10(axis_lab),labels=axis_lab)
    
    plott<-log10(links_1_AB_D)[,whoo]
    # points(log10(rep(1.2,sample_n)),plott[1,],col=alpha("white",opac_points),cex=cex.pp,xlim=c(0,2.5),ylim=c(-1,5),xaxt="n",xaxs="i",xaxt="n",ylab="",xlab="",
    #        lwd=thik,cex.axis=lab_size,cex.lab=lab_size)
    # for (i in c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))){
    #   points(log10(rep(i,sample_n)),plott[i,],col=alpha("white",opac_points),cex=cex.pp)
    # }
    lines(log10(c(1.2,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))),
          links_1_AB_D_m[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot)))],lty=mean_l_type,
          cex=0.1,col=col_AB_D,lwd=thik)
    abline(h=0,lty=4)
    lines(log10(c(1.2,round(seq(2,number_adds/50,length.out=points_to_plot)))),
          log10(links_1_AB_D[c(1,round(10^seq(log10(2),log10(number_adds),length.out=points_to_plot))),randd]+0.5),
          cex=0.1,lty=random_l_type,xlab=expression(log[10](Number~of~resource~additions)),ylab="",col=alpha(col_AB_D,opac_random))
    mtext("(d)",side=1,line=-2,adj=0,cex=comm_cex_lab)
  }
  par(mar=c(0,0,0,0))
  plot.new()
  mtext("Number of resource additions",side=3,line=-1,adj=0.5,cex=mtext_cex)
}
dev.off()
}

########################################
            # Chapter Four  #
########################################

# Figure 4.3 #
{
  pdf("../Images/Meta_com_time_series_no_intra.pdf",height=10*0.5*mult,width=10*mult)
  all_filez<-list.files(pattern="FALSE_meta", path = "Metacommunity/Structured_Topology/",full.names = T) #O
  colz<-colfunc(length(all_filez))#rainbow(dim(all_spec)[2])
  opac<-0.6
  disp_val<-as.numeric(strsplit(all_filez[1],"_")[[1]][1])
  par(mfrow=c(2,1))
  tabel<-read.table(all_filez[1])
  par(mar=c(0,5.5,1,1))
  with(tabel,plot(log10(1:length(occupancy_m)),occupancy_m,ylab="Occupancy %",xaxt="n",
                  type="l",xlab="",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[1],opac),ylim=c(0,100)))
  for (loop in 2:length(all_filez)){
    tabel<-read.table(all_filez[loop])
    with(tabel,lines(log10(1:length(occupancy_m)),occupancy_m,ylab="Occupancy %",
                     type="l",xlab="",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[loop],opac)))
  }
  
  mtext("(a)",side=3,line=-1,adj=0.01,cex=mtext_cex)
  mtext("Dispersal \n Rate",side=3,line=-1,adj=0.99,cex=mtext_cex)
  gradientLegend(valRange=c(0.005,0.1), color = alpha(colz,opac), coords=TRUE, depth=0.1,pos=0.75,pos.num=2,dec=3,
                 side=4,border.col = "black")
  
  par(mar=c(4.5,5.5,1,1))
  tabel<-read.table(all_filez[1])
  yy<-with(tabel,agg_m[which(agg_m!=0)])
  with(tabel,plot(log10(1:length(occupancy_m)),rep(0,length(occupancy_m)),ylab=expression(bar(log[10](a[ij]))),
                  xlab="Time steps",xaxt="n",
                  type="l",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[1],opac),ylim=c(-6,-3)))
  axis(1,labels=c(1,10,100,10000),at=log10(c(1,10,100,10000)),cex.axis=comm_cexx)
  mtext("(b)",side=3,line=-1,adj=0.01,cex=mtext_cex)
  for (loop in 1:length(all_filez)){
    tabel<-read.table(all_filez[loop])
    yy<-with(tabel,agg_m[which(agg_m!=0)])
    lines(log10(1:length(tabel$occupancy_m))[1:length(yy)],yy,ylab="Occupancy %",
          type="l",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[loop],opac))
  }
  legend("topright",legend = "Assembly model \n steady state",lty=c(1),col="black",bty="n",lwd=2)
  abline(h=-5.4)
  dev.off()
}

# Figure 4.4 #
{
library(plotfunctions);library(RColorBrewer)
colfunc <- colorRampPalette(c("blue", "red"))
mult<-1
comm_cexx<-1.5*mult
mtext_cex<-1*mult
pdf("../Images/Meta_com_time_series_2.pdf",height=10*0.5*mult,width=10*mult) #O
all_filez<-list.files(pattern="TRUE_meta",path  = "Metacommunity/Structured_Topology/",full.names = T) #O
colz<-colfunc(length(all_filez))
opac<-0.6
disp_val<-as.numeric(strsplit(all_filez[1],"_")[[1]][1])
par(mfrow=c(2,1))
tabel<-read.table(all_filez[1])
par(mar=c(0,5.5,1,1))
with(tabel,plot(1:length(occupancy_m),occupancy_m,ylab="Occupancy %",
                type="l",xlab="",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[1],opac),ylim=c(0,100)))
for (loop in 2:length(all_filez)){
  tabel<-read.table(all_filez[loop])
  with(tabel,lines(1:length(occupancy_m),occupancy_m,ylab="Occupancy %",
                   type="l",xlab="",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[loop],opac)))
}
mtext("(a)",side=3,line=-1,adj=0.01,cex=mtext_cex)
mtext("Dispersal \n Rate",side=3,line=-1,adj=0.99,cex=mtext_cex)
gradientLegend(valRange=c(0.005,0.1), color = alpha(colz,opac), coords=TRUE, depth=0.1,pos=0.75,pos.num=2,dec=3,
               side=4,border.col = "black")
par(mar=c(4.5,5.5,1,1))
tabel<-read.table(all_filez[1])
yy<-with(tabel,agg_m[which(agg_m!=0)])
plot(1:length(yy),yy,ylab=expression(bar(log[10](a[ij]))),xlim=c(0,dim(tabel)[1]),xlab="Time steps (x100)",
     type="l",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[1],opac),ylim=c(-6,-3))
mtext("(b)",side=3,line=-1,adj=0.01,cex=mtext_cex)
for (loop in 2:length(all_filez)){
  tabel<-read.table(all_filez[loop])
  yy<-with(tabel,agg_m[which(agg_m!=0)])
  lines(1:length(yy),yy,ylab="Occupancy %",
        type="l",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[loop],opac))
}
legend("topright",legend = "Assembly model \n steady state",lty=c(1),col="black",bty="n",lwd=2)
abline(h=-5.4)
dev.off()
}

# Figure 4.5 #
{
pdf("../Images/Meta_com_time_series_random.pdf",height=10*0.5*mult,width=10*mult)
all_filez<-list.files(pattern="TRUE_meta",path  = "Metacommunity/Random_Topology/",full.names = T)
colz<-colfunc(length(all_filez))
opac<-0.6
disp_val<-as.numeric(strsplit(all_filez[1],"_")[[1]][1])
par(mfrow=c(2,1))
tabel<-read.table(all_filez[1])
par(mar=c(0,5.5,1,1))
with(tabel,plot(1:length(occupancy_m),occupancy_m,ylab="Occupancy %",
                type="l",xlab="",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[1],opac),ylim=c(0,100)))
for (loop in 2:length(all_filez)){
  tabel<-read.table(all_filez[loop])
  with(tabel,lines(1:length(occupancy_m),occupancy_m,ylab="Occupancy %",
                   type="l",xlab="",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[loop],opac)))
}

mtext("(a)",side=3,line=-1,adj=0.01,cex=mtext_cex)
mtext("Dispersal \n Rate",side=3,line=-1,adj=0.99,cex=mtext_cex)
gradientLegend(valRange=c(0.005,0.1), color = alpha(colz,opac), coords=TRUE, depth=0.1,pos=0.75,pos.num=2,dec=3,
               side=4,border.col = "black")

par(mar=c(4.5,5.5,1,1))
tabel<-read.table(all_filez[1])
yy<-with(tabel,agg_m[which(agg_m!=0)])
plot(1:length(yy),yy,ylab=expression(bar(log[10](a[ij]))),xlim=c(0,dim(tabel)[1]),xlab="Time steps (x100)",
     type="l",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[1],opac),ylim=c(-6,-3))
mtext("(b)",side=3,line=-1,adj=0.01,cex=mtext_cex)
for (loop in 2:length(all_filez)){
  tabel<-read.table(all_filez[loop])
  yy<-with(tabel,agg_m[which(agg_m!=0)])
  lines(1:length(yy),yy,ylab="Occupancy %",
        type="l",xaxt="n",cex.axis=comm_cexx,cex.lab=comm_cexx,col=alpha(colz[loop],opac))
}
legend("topright",legend = "Assembly model \n steady state",lty=c(1),col="black",bty="n",lwd=2)
abline(h=-5.4)
dev.off()
}


################# End ##################