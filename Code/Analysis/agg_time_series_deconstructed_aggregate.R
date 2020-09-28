all_dirz<-list.dirs(full.names =  F,recursive = F)[-1]
setwd(all_dirz[1])
colz<-rainbow(length(all_dirz))

try(load("all_data.RData"),silent=T)
all_species<-all_species[which(all_species$id!="id"),]
all_species<-all_species[which(all_species$id!=0),]
all_species[]<-lapply(all_species,function(x)as.numeric(as.character(x)))
all_species<-all_species[order(all_species$id),]
all_species_keep<-all_species[which(!all_species$outi<0),]

all_species_aggr<-as.data.frame(all_species_keep$id);colnames(all_species_aggr)<-"id"
all_species_aggr$m_l_a<-all_species_keep$m_l_a
xx<-round(seq(1,length(all_species_aggr$m_l_a),length.out=40))
all_species_aggr<-all_species_aggr[xx,]

for (i in 2:length(all_dirz)){
  setwd(paste("../",all_dirz[i],sep=""))
  try(load("all_data.RData"),silent=T) 
  all_species<-all_species[which(all_species$id!="id"),]
  all_species<-all_species[which(all_species$id!=0),]
  all_species[]<-lapply(all_species,function(x)as.numeric(as.character(x)))
  all_species<-all_species[order(all_species$id),]
  all_species_keep<-all_species[which(!all_species$outi<0),]
  all_species_aggr_temp<-all_species_keep$m_l_a
  xx<-round(seq(1,length(all_species_aggr_temp),length.out=40))
  all_species_aggr_temp<-all_species_aggr_temp[xx]
  all_species_aggr<-cbind(all_species_aggr,all_species_aggr_temp)
}

write.table(all_species_aggr,"Data/Pre_processed/agg_time_series_deconstructed_aggregate.txt")