#AddMeta: re_grouped_clonotypes  to OSCC.Rda

ChZhX.T.N_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoChZhX.T.N.txt",header=T,sep = "\t")
HKJ.T.Ca_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoHKJ.T.Ca.txt",header=T,sep = "\t")
HKJ.T.N_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoHKJ.T.N.txt",header = T,sep = "\t")
LWH.T.Ca_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoLWH.T.Ca.txt",header = T,sep = "\t")
LWH.T.N_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoLWH.T.N.txt",header = T,sep = "\t")
LYX.T.Ca_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoLYX.T.Ca.txt",header = T,sep = "\t")
PHD.B.Ca_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoPHD.B.Ca.txt",header = T,sep = "\t")
WRH.G.Ca_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoWRH.G.Ca.txt",header = T,sep = "\t")
WRH.G.N_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoWRH.G.N.txt",header = T,sep = "\t")
XZB.OLK_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoXZB.OLK.txt",header = T,sep = "\t")
ZhChH.T.Ca_clone <- read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/re_grouped_clonotypes.phenoZhChH.T.Ca.txt",header = T,sep = "\t")

barcoder <- function(df, prefix ,sample,  trim="\\-1"){
  df<-df[which(df$alpha.C!="NA" & df$beta.C!="NA" & df$clonotype!="NA" &
              df$alpha.C!="" & df$beta.C!="" & df$clonotype!=""),]
  df$barcode <- gsub(trim, "", df$barcode)
  df$barcode <- paste0(df$barcode, prefix)
  df$sample <- sample
  df
}

ChZhX.T.N_clone <- barcoder(ChZhX.T.N_clone, prefix = "-ChZhX.T.N","ChZhX.T.N")
HKJ.T.Ca_clone <- barcoder(HKJ.T.Ca_clone, prefix = "-HKJ.T.Ca","HKJ.T.Ca")
HKJ.T.N_clone <- barcoder(HKJ.T.N_clone, prefix = "-HKJ.T.N","HKJ.T.N")
LWH.T.Ca_clone <- barcoder(LWH.T.Ca_clone, prefix = "-LWH.T.Ca","LWH.T.Ca")
LWH.T.N_clone <- barcoder(LWH.T.N_clone, prefix = "-LWH.T.N","LWH.T.N")
LYX.T.Ca_clone <- barcoder(LYX.T.Ca_clone, prefix = "-LYX.T.Ca","LYX.T.Ca")
PHD.B.Ca_clone <- barcoder(PHD.B.Ca_clone, prefix = "-PHD.B.Ca","PHD.B.Ca")
WRH.G.Ca_clone <- barcoder(WRH.G.Ca_clone, prefix = "-WRH.G.Ca","WRH.G.Ca")
WRH.G.N_clone <- barcoder(WRH.G.N_clone, prefix = "-WRH.G.N","WRH.G.N")
XZB.OLK_clone <- barcoder(XZB.OLK_clone, prefix = "-XZB.OLK","XZB.OLK")
ZhChH.T.Ca_clone <- barcoder(ZhChH.T.Ca_clone, prefix = "-ZhChH.T.Ca","ZhChH.T.Ca")


OSCC_re_grouped_clonotypes<-rbind.data.frame(ChZhX.T.N_clone,HKJ.T.Ca_clone,HKJ.T.N_clone,LWH.T.Ca_clone,LWH.T.N_clone,
                                             LYX.T.Ca_clone,PHD.B.Ca_clone,WRH.G.Ca_clone,WRH.G.N_clone,XZB.OLK_clone,
                                             ZhChH.T.Ca_clone)
row.names(OSCC_re_grouped_clonotypes)<-OSCC_re_grouped_clonotypes$barcode
OSCC_re_grouped_clonotypes$re_grouped_clonotype<-OSCC_re_grouped_clonotypes$clonotype
OSCC_re_grouped_clonotypes<-OSCC_re_grouped_clonotypes[,-c(9:14)]
OSCC_re_grouped_clonotypes$re_grouped_clonotype<-as.character(OSCC_re_grouped_clonotypes$re_grouped_clonotype)
save(OSCC_re_grouped_clonotypes,file="E:/OSCC/cellranger-4_0_0 result/5/TCR/OSCC_re_grouped_clonotypes.Rda")
write.table(OSCC_re_grouped_clonotypes,file = "E:/OSCC/cellranger-4_0_0 result/5/TCR/OSCC_re_grouped_clonotypes.txt",row.names = F,quote = F,sep = "\t")


clone.state.meta<-read.table("E:/OSCC/cellranger-4_0_0 result/5/TCR/barcode_clonalTypes.txt",header = T,stringsAsFactors = F)
table(clone.state.meta$duplicatedTypes,clone.state.meta$oligoTypes)
row.names(clone.state.meta)<-clone.state.meta$barcode
norCD8T<-AddMetaData(norCD8T,metadata = clone.state.meta)
table(norCD8T@meta.data$duplicatedTypes)
norCD8T@meta.data$duplicatedTypes[which(is.na(norCD8T@meta.data$duplicatedTypes))]<-"NoTCR"
norCD8T@meta.data$duplicatedTypes<-factor(norCD8T@meta.data$duplicatedTypes,levels = c("NoTCR", "uniqueTypes", "duplicatedTypes", "clonalTypes" ))
table(norCD8T@meta.data$oligoTypes)
norCD8T@meta.data$oligoTypes[which(is.na(norCD8T@meta.data$oligoTypes))]<-"NoTCR"
norCD8T@meta.data$oligoTypes<-factor(norCD8T@meta.data$oligoTypes,levels = c("NoTCR", "uniqueTypes", "oligoTypes", "clonalTypes" ))


OSCC_re_grouped_clonotypes$sample<-as.factor(OSCC_re_grouped_clonotypes$sample)
a<-tapply(OSCC_re_grouped_clonotypes$re_grouped_clonotype,OSCC_re_grouped_clonotypes$sample,function(x){
  length(table(x))
})

#?????????????????????
load("E:/MF/scRNA-seq/MF_TCR.meta.Rda")
table(MF_TCR.meta$samples)
MF_TCR.meta$re_grouped_clonotype<-as.character(MF_TCR.meta$re_grouped_clonotype)
head(MF_TCR.meta)
b<-tapply(MF_TCR.meta$re_grouped_clonotype,MF_TCR.meta$samples,function(x){
  length(table(x))
})
       
#????????????????????????clonotypes
#MF17
MF17_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF17" &
                                                             MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF17_re_grouped_clonotypes<-MF17_re_grouped_clonotypes[!duplicated(MF17_re_grouped_clonotypes)]
MF17_re_grouped_clonotypes[1:10]
#MF18
MF18_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF18" &
                                                             MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF18_re_grouped_clonotypes<-MF18_re_grouped_clonotypes[!duplicated(MF18_re_grouped_clonotypes)]
MF18_re_grouped_clonotypes[1:10]
#MF21_1
MF21_1_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF21-1" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF21_1_re_grouped_clonotypes<-MF21_1_re_grouped_clonotypes[!duplicated(MF21_1_re_grouped_clonotypes)]
MF21_1_re_grouped_clonotypes[1:100]
#MF21_2
MF21_2_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF21-2" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF21_2_re_grouped_clonotypes<-MF21_2_re_grouped_clonotypes[!duplicated(MF21_2_re_grouped_clonotypes)]
MF21_2_re_grouped_clonotypes[1:100]
#MF21
MF21_re_grouped_clonotypes<-intersect(MF21_1_re_grouped_clonotypes,MF21_2_re_grouped_clonotypes)
write.csv(MF21_re_grouped_clonotypes,"E:/OSCC/cellranger-4_0_0 result/5/TCR/MF21_intersect_re_grouped_clonotypes.csv")
#MF22
MF22_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF22" &
                                                             MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF22_re_grouped_clonotypes<-MF22_re_grouped_clonotypes[!duplicated(MF22_re_grouped_clonotypes)]
MF22_re_grouped_clonotypes[1:10]
#MF26
MF26_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF26" &
                                                             MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF26_re_grouped_clonotypes<-MF26_re_grouped_clonotypes[!duplicated(MF26_re_grouped_clonotypes)]
MF26_re_grouped_clonotypes[1:10]
#MF27_1
MF27_1_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF27-1" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF27_1_re_grouped_clonotypes<-MF27_1_re_grouped_clonotypes[!duplicated(MF27_1_re_grouped_clonotypes)]
MF27_1_re_grouped_clonotypes[1:10]
#MF27_2
MF27_2_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF27-2" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF27_2_re_grouped_clonotypes<-MF27_2_re_grouped_clonotypes[!duplicated(MF27_2_re_grouped_clonotypes)]
MF27_2_re_grouped_clonotypes[1:10]
#MF27
MF27_re_grouped_clonotypes<-intersect(MF27_1_re_grouped_clonotypes,MF27_2_re_grouped_clonotypes)
write.csv(MF27_re_grouped_clonotypes,"E:/OSCC/cellranger-4_0_0 result/5/TCR/MF27_intersect_re_grouped_clonotypes.csv")

#MF28_1
MF28_1_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF28-1" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF28_1_re_grouped_clonotypes<-MF28_1_re_grouped_clonotypes[!duplicated(MF28_1_re_grouped_clonotypes)]
MF28_1_re_grouped_clonotypes[1:10]
#MF28_2
MF28_2_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF28-2" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF28_2_re_grouped_clonotypes<-MF28_2_re_grouped_clonotypes[!duplicated(MF28_2_re_grouped_clonotypes)]
MF28_2_re_grouped_clonotypes[1:10]
#MF28
MF28_re_grouped_clonotypes<-intersect(MF28_1_re_grouped_clonotypes,MF28_2_re_grouped_clonotypes)
write.csv(MF28_re_grouped_clonotypes,"E:/OSCC/cellranger-4_0_0 result/5/TCR/MF28_intersect_re_grouped_clonotypes.csv")

#MF30_1
MF30_1_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF30-1" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF30_1_re_grouped_clonotypes<-MF30_1_re_grouped_clonotypes[!duplicated(MF30_1_re_grouped_clonotypes)]
MF30_1_re_grouped_clonotypes[1:10]
#MF30_2
MF30_2_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="MF30-2" &
                                                               MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
MF30_2_re_grouped_clonotypes<-MF30_2_re_grouped_clonotypes[!duplicated(MF30_2_re_grouped_clonotypes)]
MF30_2_re_grouped_clonotypes[1:10]
#MF30
MF30_re_grouped_clonotypes<-intersect(MF30_1_re_grouped_clonotypes,MF30_2_re_grouped_clonotypes)
write.csv(MF30_re_grouped_clonotypes,"E:/OSCC/cellranger-4_0_0 result/5/TCR/MF30_intersect_re_grouped_clonotypes.csv")
#pcALCL_1
pcALCL_1_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="pcALCL1" &
                                                                 MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
pcALCL_1_re_grouped_clonotypes<-pcALCL_1_re_grouped_clonotypes[!duplicated(pcALCL_1_re_grouped_clonotypes)]
pcALCL_1_re_grouped_clonotypes[1:10]
#pcALCL_2
pcALCL_2_re_grouped_clonotypes<-as.character(MF_TCR.meta[which(MF_TCR.meta$samples=="pcALCL2" &
                                                                 MF_TCR.meta$re_grouped_clonotype!="NA"),"re_grouped_clonotype"])
pcALCL_2_re_grouped_clonotypes<-pcALCL_2_re_grouped_clonotypes[!duplicated(pcALCL_2_re_grouped_clonotypes)]
pcALCL_2_re_grouped_clonotypes[1:10]

list1 <- list(MF17_CD4 = MF17_re_grouped_clonotypes, MF18_CD4 = MF18_re_grouped_clonotypes, 
              MF21_CD4 = MF21_re_grouped_clonotypes, MF22_CD4 = MF22_re_grouped_clonotypes, 
              MF26_CD4 = MF26_re_grouped_clonotypes, MF27_CD4 = MF27_re_grouped_clonotypes, 
              MF28_CD4 = MF28_re_grouped_clonotypes, MF30_CD4 = MF30_re_grouped_clonotypes,
              pcALCL_1_CD8 = pcALCL_1_re_grouped_clonotypes,pcALCL_2_CD4 = pcALCL_2_re_grouped_clonotypes)
#2????????????????????????
list2<-list()
list2.name<-c()
len<-length(list1)-1
c=0;
for(i in 1:len){
  k=i+1;
  for(j in k:length(list1)){
    c=c+1;
    list2[[c]]<-data.frame(t(intersect(list1[[i]], list1[[j]])))
    name<-c(names(list1[i]),names(list1[j]))
    name<-name[order(name)]
    list2.name<-c(list2.name,paste(name,collapse=","))
  }
}
intersect2 <- rbind.fill(list2)
names(list2)<-list2.name
write.table(intersect2,file="E:/OSCC/cellranger-4_0_0 result/5/TCR/2sample.intersect.re_grouped_clonotypes.txt",sep = "\t",col.names = F,row.names = list2.name,quote = F)

