### Parameters #####
sigma<-5.5
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
number_adds<-500
nreps<-10
rlx_t2<-500
rlx_t<-500
agg_1<-exp(-11.75);agg_2<-exp(-11);agg_3<-exp(-10.25)


setwd("Data/Raw/Serial_extinction_sub_model/")
links_1_cc_2<-as.matrix(read.table(list.files(pattern = "links_1_cc_2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_1_cc_2.txt",full.names = T, recursive = T)[-1]){
  links_1_cc_2<-cbind(links_1_cc_2,read.table(i)*(epsilon*K/resp))
}
links_1_cc_2<-links_1_cc_2[, colSums(links_1_cc_2 != 0) > 0]

links_2_cc_2<-as.matrix(read.table(list.files(pattern = "links_2_cc_2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_2_cc_2.txt",full.names = T, recursive = T)[-1]){
  links_2_cc_2<-cbind(links_2_cc_2,read.table(i)*(epsilon*K/resp))
}
links_2_cc_2<-links_2_cc_2[, colSums(links_2_cc_2 != 0) > 0]

links_3_cc_2<-as.matrix(read.table(list.files(pattern = "links_3_cc_2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_3_cc_2.txt",full.names = T, recursive = T)[-1]){
  links_3_cc_2<-cbind(links_3_cc_2,read.table(i)*(epsilon*K/resp))
}
links_3_cc_2<-links_3_cc_2[, colSums(links_3_cc_2 != 0) > 0]

links_1_cc2_2<-as.matrix(read.table(list.files(pattern = "links_1_cc2_2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_1_cc2_2.txt",full.names = T, recursive = T)[-1]){
  links_1_cc2_2<-cbind(links_1_cc2_2,read.table(i)*(epsilon*K/resp))
}
links_1_cc2_2<- links_1_cc2_2[, colSums( links_1_cc2_2 != 0) > 0]

links_2_cc2_2<-as.matrix(read.table(list.files(pattern = "links_2_cc2_2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_2_cc2_2.txt",full.names = T, recursive = T)[-1]){
  links_2_cc2_2<-cbind(links_2_cc2_2,read.table(i)*(epsilon*K/resp))
}
links_2_cc2_2<- links_2_cc2_2[, colSums(links_2_cc2_2 != 0) > 0]

links_3_cc2_2<-as.matrix(read.table(list.files(pattern = "links_3_cc2_2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_3_cc2_2.txt",full.names = T, recursive = T)[-1]){
  links_3_cc2_2<-cbind(links_3_cc2_2,read.table(i)*(epsilon*K/resp))
}
links_3_cc2_2<- links_3_cc2_2[, colSums(links_3_cc2_2 != 0) > 0]

links_1_cc<-as.matrix(read.table(list.files(pattern = "links_1_cc.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_1_cc.txt",full.names = T, recursive = T)[-1]){
  links_1_cc<-cbind(links_1_cc,read.table(i)*(epsilon*K/resp))
}
links_1_cc<- links_1_cc[, colSums(links_1_cc != 0) > 0]

links_1_cc2<-as.matrix(read.table(list.files(pattern = "links_1_cc2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_1_cc2.txt",full.names = T, recursive = T)[-1]){
  links_1_cc2<-cbind(links_1_cc2,read.table(i)*(epsilon*K/resp))
}
links_1_cc2<- links_1_cc2[, colSums(links_1_cc2 != 0) > 0]


links_2_cc<-as.matrix(read.table(list.files(pattern = "links_2_cc.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_2_cc.txt",full.names = T, recursive = T)[-1]){
  links_2_cc<-cbind(links_2_cc,read.table(i)*(epsilon*K/resp))
}
links_2_cc<- links_2_cc[, colSums(links_2_cc != 0) > 0]


links_2_cc2<-as.matrix(read.table(list.files(pattern = "links_2_cc2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_2_cc2.txt",full.names = T, recursive = T)[-1]){
links_2_cc2<-cbind(links_2_cc2,read.table(i)*(epsilon*K/resp))
}
links_2_cc2<- links_2_cc2[, colSums(links_2_cc2 != 0) > 0]

links_3_cc<-as.matrix(read.table(list.files(pattern = "links_3_cc.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_3_cc.txt",full.names = T, recursive = T)[-1]){
  try({links_3_cc<-cbind(links_3_cc,read.table(i)*(epsilon*K/resp))},silent=T)
}
links_3_cc<- links_3_cc[, colSums(links_3_cc != 0) > 0]

links_3_cc2<-as.matrix(read.table(list.files(pattern = "links_3_cc2.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_3_cc2.txt",full.names = T, recursive = T)[-1]){
  try({links_3_cc2<-cbind(links_3_cc2,read.table(i)*(epsilon*K/resp))},silent=T)
}
links_3_cc2<- links_3_cc2[, colSums(links_3_cc2 != 0) > 0]

links_1_ccc<-as.matrix(read.table(list.files(pattern = "links_1_ccc.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_1_ccc.txt",full.names = T, recursive = T)[-1]){
  links_1_ccc<-cbind(links_1_ccc,read.table(i)*(epsilon*K/resp))
}
links_1_ccc<- links_1_ccc[, colSums(links_1_ccc != 0) > 0]

links_2_ccc<-as.matrix(read.table(list.files(pattern = "links_2_ccc.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_2_ccc.txt",full.names = T, recursive = T)[-1]){
  links_2_ccc<-cbind(links_2_ccc,read.table(i)*(epsilon*K/resp))
}
links_2_ccc<- links_2_ccc[, colSums(links_2_ccc != 0) > 0]

links_3_ccc<-as.matrix(read.table(list.files(pattern = "links_3_ccc.txt",full.names = T, recursive = T)[1])*(epsilon*K/resp))
for (i in list.files(pattern = "links_3_ccc.txt",full.names = T, recursive = T)[-1]){
  links_3_ccc<-cbind(links_3_ccc,read.table(i)*(epsilon*K/resp))
}
links_3_ccc<- links_3_ccc[, colSums(links_3_ccc != 0) > 0]

links_1_DD_bb<-as.matrix(read.table(list.files(pattern = "links_1_DD_bb.txt",full.names = T, recursive = T)[1]))
for (i in list.files(pattern = "links_1_DD_bb.txt",full.names = T, recursive = T)[-1]){
  links_1_DD_bb<-cbind(links_1_DD_bb,read.table(i))
}
links_1_DD_bb<- links_1_DD_bb[, colSums(links_1_DD_bb != 0) > 0]

links_2_DD_bb<-as.matrix(read.table(list.files(pattern = "links_2_DD_bb.txt",full.names = T, recursive = T)[1]))
for (i in list.files(pattern = "links_2_DD_bb.txt",full.names = T, recursive = T)[-1]){
  links_2_DD_bb<-cbind(links_2_DD_bb,read.table(i))
}
links_2_DD_bb<- links_2_DD_bb[, colSums(links_2_DD_bb != 0) > 0]

links_3_DD_bb<-as.matrix(read.table(list.files(pattern = "links_3_DD_bb.txt",full.names = T, recursive = T)[1]))
for (i in list.files(pattern = "links_3_DD_bb.txt",full.names = T, recursive = T)[-1]){
  links_3_DD_bb<-cbind(links_3_DD_bb,read.table(i))
}
links_3_DD_bb<- links_3_DD_bb[, colSums(links_3_DD_bb != 0) > 0]

links_1_ine<-((links_1_cc-1)/links_1_cc2)/((links_1_cc_2-1)/links_1_cc2_2)
links_2_ine<-((links_2_cc-1)/links_2_cc2)/((links_2_cc_2-1)/links_2_cc2_2)
links_3_ine<-((links_3_cc-1)/links_3_cc2)/((links_3_cc_2-1)/links_3_cc2_2)

links_1_AB<-((links_1_cc_2-1)/links_1_cc2_2)
links_2_AB<-((links_2_cc_2-1)/links_2_cc2_2)
links_3_AB<-((links_3_cc_2-1)/links_3_cc2_2)

links_1_AB_m<-rowMeans(log10(links_1_AB),na.rm=T)
links_2_AB_m<-rowMeans(log10(links_2_AB),na.rm=T)
links_3_AB_m<-rowMeans(log10(links_3_AB),na.rm=T)

links_1_ine_m<-rowMeans(log10(links_1_ine),na.rm=T)
links_2_ine_m<-rowMeans(log10(links_2_ine),na.rm=T)
links_3_ine_m<-rowMeans(log10(links_3_ine),na.rm=T)

links_1_c_mmm<-rowMeans(log10(links_1_ccc-1),na.rm=T)
links_2_c_mmm<-rowMeans(log10(links_2_ccc-1),na.rm=T)
links_3_c_mmm<-rowMeans(log10(links_3_ccc-1),na.rm=T)

links_1_c_mm<-rowMeans(log10(links_1_cc-1),na.rm=T)
links_2_c_mm<-rowMeans(log10(links_2_cc-1),na.rm=T)
links_3_c_mm<-rowMeans(log10(links_3_cc-1),na.rm=T)

links_1_c2_mm<-rowMeans(log10(links_1_cc2))
links_2_c2_mm<-rowMeans(log10(links_2_cc2))
links_3_c2_mm<-rowMeans(log10(links_3_cc2))

links_1_DD_b_mm<-rowMeans(links_1_DD_bb)
links_2_DD_b_mm<-rowMeans(links_2_DD_bb)
links_3_DD_b_mm<-rowMeans(links_3_DD_bb)

links_1_cc2_2_m<-rowMeans(log10(links_1_cc2_2))
links_1_cc_2_m<-rowMeans(log10(links_1_cc_2))
links_2_cc2_2_m<-rowMeans(log10(links_2_cc2_2))
links_2_cc_2_m<-rowMeans(log10(links_2_cc_2))
links_3_cc2_2_m<-rowMeans(log10(links_3_cc2_2))
links_3_cc_2_m<-rowMeans(log10(links_3_cc_2))

save.image("../../Pre_Processed/Numerical_serial_extinction_sub_model/Numerical_serial_extinction_sub_model.RData")
