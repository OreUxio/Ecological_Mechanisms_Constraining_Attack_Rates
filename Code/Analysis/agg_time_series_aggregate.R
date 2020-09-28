all_dirz<-list.dirs(full.names =  F)[-1]

setwd(all_dirz[1])
colz<-rainbow(length(all_dirz))

all_species<-read.table("Data/Full_assembly_model/list_of_aggressivity_values.dat") # Run 
all_species_agg<-data.frame(all_species$V1);colnames(all_species_agg)<-"V1";all_species_agg$V2<-all_species$V2
xx<-round(seq(1,length(all_species_agg$V2),length.out=40))
all_species_agg<-all_species_agg[xx,]

for (i in 2:length(all_dirz)){
  setwd(paste("../",all_dirz[i],sep=""))
  all_species<-read.table("Data/Full_assembly_model/list_of_aggressivity_values.dat")
  all_species_temp<-all_species$V2

  xx<-round(seq(1,length(all_species_temp),length.out=40))
  all_species_temp<-all_species_temp[xx]
  all_species_agg<-cbind(all_species_agg,all_species_temp)
}

write.table(all_species_agg,"Data/Pre_processed/agg_time_series_aggregate.txt")