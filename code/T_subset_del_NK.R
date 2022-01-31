library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)
library(umap)
library(RColorBrewer)
library(cowplot)
library(randomcoloR)
library(corrplot)
library(SingleCellExperiment)
load("E:\\MF\\scRNA-seq\\cell_type_cols.Rda")
cell_type_cols <- c(brewer.pal(9, "Set1"),"#8B008B",
                    "#FF34B3","#CD5C5C","#BC8F8F","#20B2AA","#00F5FF","#FF3030","#ADFF2F","#FFA500","#FF6A6A","#7FFFD4",
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1",
                    "#7CFC00","#000000","#708090")

setwd("E:\\MF\\scRNA-seq\\MF-ALL-20200425\\MF-TCR-20200516\\T_cell-20200522\\T_subset_del_26\\T_subset_del_NK")

load("T_subset_del_NK_celltype.Rda")


pdf(file="Fig S1F.1.pdf", width = 10,height = 3)
FeaturePlot(T_subset_del_NK,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("CD3E", "CD7", "SELPLG"))
dev.off()
pdf(file="Fig S1F.2.pdf", width = 14,height = 6)
FeaturePlot(T_subset_del_NK,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 4,
            features = c("PTPRC", "CD4", "CD8A", "FOXP3", "TBX21", "GATA3", "RORC", "KLRG1"))
dev.off()

pdf(file="Fig 4B.pdf", width = 16,height = 6)
FeaturePlot(T_subset_del_NK,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 5,
            features = c("CD4","CD8B","SELL","CD27","CCR7", "GZMA","GZMB","PRF1","GNLY","IFNG"))
dev.off()

Idents(T_subset_del_NK)<-T_subset_del_NK@meta.data$cytotoxic_Tcm
T_subset_del_NK.del18.27<-subset(T_subset_del_NK,idents=c("Tcm","cytotoxic_CD4","cytotoxic_CD8"))
save(T_subset_del_NK.del18.27,file = "del18.27/T_subset_del_NK.del18.27.Rda")



T_subset_del_NK@meta.data$samples<-as.character(T_subset_del_NK@meta.data$samples)
T_subset_del_NK@meta.data[which(T_subset_del_NK@meta.data$samples=="MF14_v3"),"samples"]<-"MF14"
T_subset_del_NK@meta.data[which(T_subset_del_NK@meta.data$samples=="MF15_v3"),"samples"]<-"MF15"
save(T_subset_del_NK,file="E:/MF/scRNA-seq/MF-ALL/10X/T_subset/T_subset_del_NK.Rda")

Idents(T_subset_del_NK)<-T_subset_del_NK@meta.data$cytotoxic_Tcm
Tcm_T<-subset(T_subset_del_NK,idents=c("Tcm"))
save(Tcm_T,file="Tcm T/Tcm_T.Rda")

T_subset_del_NK@meta.data[3000:3500,]
table(T_subset_del_NK@meta.data$tu.type)
Idents(T_subset_del_NK)<-T_subset_del_NK@meta.data$tu.type

T_subset_del_NK@meta.data$umis[3500:4000]


norCD4T<-subset(T_subset_del_NK,idents=c("CD4T"))
save(norCD4T,file="norCD4T.Rda")
table(norCD4T@meta.data$tu.type)

CD4<-subset(T_subset_del_NK,idents=c("CD4Tu","CD4T","Treg"))
save(CD4,file="CD4.Rda")
table(CD4@meta.data$tu.type)

CD8<-subset(T_subset_del_NK,idents=c("CD8Tu","CD8T"))
table(CD8@meta.data$tu.type)
save(CD8,file="CD8.Rda")

Idents(T_subset_del_NK)<-T_subset_del_NK@meta.data$Type
Tumor_T<-subset(T_subset_del_NK,idents=c("Tumor T"))
save(Tumor_T,file="Tumor_T.Rda")
normal_T<-subset(T_subset_del_NK,idents=c("Normal T"))
save(normal_T,file="normal_T.Rda")

Idents(T_subset_del_NK)<-T_subset_del_NK@meta.data$tu.type
Treg<-subset(T_subset_del_NK,idents=c("Treg"))
save(Treg,file="Treg.Rda")  
  
load("D:/HSM/MF/MF-ALL-20200425/MF-TCR-20200516/T_cell-20200522/T_subset_del_26/T_subset.clone.state.meta.Rda")
T_subset_del_NK<-AddMetaData(T_subset_del_NK,metadata = T_subset.clone.state.meta)
table(T_subset_del_NK@meta.data$clone.state.duplicated)
meta<-T_subset_del_NK@meta.data
meta<-subset(meta,meta$clone.state.duplicated!="NoTCR")
table(meta$clone.state.duplicated)
table(meta$clone.state.oligo)

pdf(file="bar_clone.state_split.pdf", width = 8,height = 4)
ggplot(T_subset_del_NK@meta.data, aes(x=seurat_clusters,fill=clone.state.duplicated))+ 
  xlab(NULL) + ylab(NULL) + 
  labs(title = "")+
  geom_bar(position = "fill")+
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=rev(c("#007371","#45D04F","#FFEA00","grey80")))
ggplot(meta, aes(x=seurat_clusters,fill=clone.state.duplicated))+ 
  xlab(NULL) + ylab(NULL) + 
  labs(title = "")+
  geom_bar(position = "fill")+
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=rev(c("#007371","#45D04F","#FFEA00")))

ggplot(T_subset_del_NK@meta.data, aes(x=seurat_clusters,fill=clone.state.oligo))+ 
  xlab(NULL) + ylab(NULL) + 
  labs(title = "")+
  geom_bar(position = "fill")+
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=rev(c("#007371","#45D04F","#FFEA00","grey80")))
ggplot(meta, aes(x=seurat_clusters,fill=clone.state.oligo))+ 
  xlab(NULL) + ylab(NULL) + 
  labs(title = "")+
  geom_bar(position = "fill")+
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=rev(c("#007371","#45D04F","#FFEA00")))

dev.off()

head(T_subset_del_NK@meta.data)
table(T_subset_del_NK@meta.data$seurat_clusters,T_subset_del_NK@meta.data$samples)

factor_level<-names(table(T_subset_del_NK@meta.data$cellType))[c(-13,-17)]
label<-c("0.tu.MF28.exh","1.CD4mT","2.CD8eT","3.tu.TCRloss.MF28","4.tu.CD8.p1",
         "5.tu.MF17","6.CD4Treg","7.tu.MF21","8.tu.MF26","9.tu.MF30-2.CD79A",
         "10.tu.TCRloss.Mix.MS4A1","11.tu.TCRloss.Mix","12.tu.CD8.MF14",
         "13.tu.MF15","14.tu.CD8.p1","15.mT.CD79A","16.tu.TCRloss.p1","17.tu.p2",
         "18.nT.MF14/15","19.tu.MF30-1","20.tu.MF22","21.tu.MF30-2.CD79A",
         "22.tu.TCRloss.MF17","23.tu.MF28.exh"
)
T_subset_del_NK@meta.data$cellType<-factor(T_subset_del_NK@meta.data$cellType,levels = factor_level,labels = label)

T_subset_del_NK@meta.data$seurat_clusters<-factor(T_subset_del_NK@meta.data$seurat_clusters,levels = c(0:11,13:15,17:25),labels = 0:23)

current.cluster.ids <- c(0:23)
new.cluster.ids <- c("Tumor T", "CD4 T", "CD8 T", "Tumor T", "Tumor T","Tumor T","Treg",
                     "Tumor T", "Tumor T", "Tumor T", "Tumor T", "Tumor T",  "Tumor T", 
                     "Tumor T", "Tumor T",  "DN T", "Tumor T", "Tumor T", "CD8 T",
                     "Tumor T","Tumor T","Tumor T","Tumor T","Tumor T")

T_subset_del_NK@meta.data$Type <- plyr::mapvalues(x = T_subset_del_NK@meta.data$seurat_clusters,
                                           from = current.cluster.ids,
                                           to = new.cluster.ids)
count_10X<-data.frame(table(T_subset_del_NK@meta.data$Type,T_subset_del_NK@meta.data$samples))

load("E:\\MF\\scRNA-seq\\MF_plate\\R-workspace\\MF_plates_cell_type.Rda")
current.cluster.ids <- c(0:8)
new.cluster.ids <- c("Tumor T", "CD4 T", "CD8 T", "Treg", "Tumor T","Fibroblast","B cell",
                     "Tumor T", "CD8 T")
MF_plates@meta.data$Type <- plyr::mapvalues(x = MF_plates@meta.data$seurat_clusters,
                                            from = current.cluster.ids,
                                            to = new.cluster.ids)
count_plates<-data.frame(table(MF_plates@meta.data$Type,MF_plates@meta.data$samples)[c("Tumor T","CD4 T", "CD8 T", "Treg"),])

count<-rbind(count_10X,count_plates)
colnames(count)<-c("Type","Samples","CellNumber")

sum<-tapply(count$CellNumber, count$Samples, sum)
sum<-rep(sum,table(count$Samples))
count$sum<-sum
count$proportion<-count$CellNumber/count$sum
write.csv(count,file = "Type.Samples.count.csv",row.names = F)
count$Samples<-factor(count$Samples,levels =c("MF18","MF27-1","MF27-2", "pcALCL2","MF30-1", "MF7","MF22","MF4",
                                              "MF21-1","MF6", "MF30-2", "MF28-2","MF21-2","MF15_v3","MF28-1",
                                              "MF14_v3","MF17", "MF26","pcALCL1"))
pdf(file="Type.Samples.count.pdf", width = 10,height = 6)
ggplot(count, aes(x=Samples,y=CellNumber,fill=Type))+ 
  xlab(NULL) + ylab(NULL) + 
  theme(legend.position = "right") +
  geom_bar(stat = "identity",position = "fill")+
  scale_fill_manual(values=cell_type_cols) + theme_bw() +
  theme(axis.text.x =  element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 12), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

current.cluster.ids <- c(0:23)
new.cluster.ids <- c("CD4Tu","CD4T","CD8T","CD4Tu","CD8Tu",
                     "CD4Tu","Treg","CD4Tu","CD4Tu","CD4Tu",
                     "CD4Tu","CD4Tu","CD8Tu","CD4Tu","CD8Tu","DNT","CD8Tu","CD4Tu",
                     "CD8T","CD4Tu","CD4Tu","CD4Tu","CD4Tu","CD4Tu")

T_subset_del_NK@meta.data$tu.type <- plyr::mapvalues(x = T_subset_del_NK@meta.data$seurat_clusters,
                                                  from = current.cluster.ids,
                                                  to = new.cluster.ids)

table(T_subset_del_NK@meta.data$tu.type)

current.cluster.ids <- names(table(T_subset_del_NK@meta.data$samples))
new.cluster.ids <- c("cytotoxic_CD8","cytotoxic_CD4","cytotoxic_CD4","no_tumor","Tcm",
                     "Tcm","cytotoxic_CD4","cytotoxic_CD4","no_tumor","no_tumor",
                     "Tcm","Tcm","Tcm","Tcm","cytotoxic_CD8","cytotoxic_CD4")

T_subset_del_NK@meta.data$cytotoxic_Tcm <- plyr::mapvalues(x = T_subset_del_NK@meta.data$samples,
                                                   from = current.cluster.ids,
                                                   to = new.cluster.ids)

table(T_subset_del_NK@meta.data$cytotoxic_Tcm)

T_subset_del_NK@meta.data$cytotoxic_Tcm_2<-factor(T_subset_del_NK@meta.data$cytotoxic_Tcm,
                                                       levels = c("cytotoxic_CD8","cytotoxic_CD4","Tcm"),
                                                       labels = c("cytotoxic","cytotoxic","Tcm"))


save(T_subset_del_NK,file="T_subset_del_NK_celltype.Rda")

pdf(file="Fraction_by_tu.type_cytotoxic_Tcm.pdf", width = 8,height = 6)
ggplot(T_subset_del_NK@meta.data, aes(cytotoxic_Tcm,fill=tu.type))+ 
  xlab(NULL) + ylab(NULL) + 
  theme(legend.position = "right") +
  geom_bar(position = "fill")+
  scale_fill_manual(values=cell_type_cols) + theme_bw() +
  theme(axis.text.x =  element_text( vjust = 0.5, hjust = 0.5, size = 12), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(T_subset_del_NK@meta.data, aes(tu.type,fill=cytotoxic_Tcm))+ 
  xlab(NULL) + ylab(NULL) + 
  theme(legend.position = "right") +
  geom_bar(position = "fill")+
  scale_fill_manual(values=cell_type_cols) + theme_bw() +
  theme(axis.text.x =  element_text( vjust = 0.5, hjust = 0.5, size = 12), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

cytotoxic_samples<-c("MF14_v3","MF15_v3","MF17","MF22","MF26","pcALCL1","pcALCL2")
Tcm_samples<-c("MF21-1","MF21-2","MF28-1","MF28-2","MF30-1","MF30-2")

Idents(T_subset_del_NK)<-T_subset_del_NK@meta.data$tu.type

CD4T<-subset(T_subset_del_NK,idents=c("CD4Tu","CD4T","Treg"))
save(CD4T,file="CD4T.Rda")

Idents(CD4T)<-CD4T@meta.data$samples
table(Idents(CD4T))

CD4_cytotoxic<-subset(CD4T,idents=cytotoxic_samples)
Idents(CD4_cytotoxic)<-CD4_cytotoxic@meta.data$tu.type
tumor_vs_normal.markers<-FindMarkers(CD4_cytotoxic,ident.1 = c("CD4Tu"), ident.2 = c("CD4T","Treg"))
x <-tumor_vs_normal.markers[tumor_vs_normal.markers[,5] < 0.01 & tumor_vs_normal.markers[,2] > 0.1 | tumor_vs_normal.markers[,2] < -0.1,]
write.csv(x,file="CD4_cytotoxic_tumor_VS._CD4normal.markers.csv")

CD4_Tcm<-subset(CD4T,idents=Tcm_samples)
Idents(CD4_Tcm)<-CD4_Tcm@meta.data$tu.type
tumor_vs_normal.markers<-FindMarkers(CD4_Tcm,ident.1 = c("CD4Tu"), ident.2 = c("CD4T","Treg"))
x <-tumor_vs_normal.markers[tumor_vs_normal.markers[,5] < 0.01 & tumor_vs_normal.markers[,2] > 0.1 | tumor_vs_normal.markers[,2] < -0.1,]
write.csv(x,file="CD4_Tcm_tumor_VS._CD4normal.markers.csv")


CD8T<-subset(T_subset_del_NK,idents=c("CD8T","CD8Tu"))
save(CD8T,file="CD8T.Rda")

Idents(CD8T)<-CD8T@meta.data$samples
table(Idents(CD8T))

CD8_cytotoxic<-subset(CD8T,idents=cytotoxic_samples)
Idents(CD8_cytotoxic)<-CD8_cytotoxic@meta.data$tu.type
tumor_vs_normal.markers<-FindMarkers(CD8_cytotoxic,ident.1 = c("CD8Tu"), ident.2 = c("CD8T"))
x <-tumor_vs_normal.markers[tumor_vs_normal.markers[,5] < 0.01 & tumor_vs_normal.markers[,2] > 0.1 | tumor_vs_normal.markers[,2] < -0.1,]
write.csv(x,file="CD8_cytotoxic_tumor_VS._CD8normal.markers.csv")

CD8_Tcm<-subset(CD8T,idents=Tcm_samples)
table(CD8_Tcm@meta.data$samples)
Idents(CD8_Tcm)<-CD8_Tcm@meta.data$tu.type
table(Idents(CD8_Tcm))
tumor_vs_normal.markers<-FindMarkers(CD8_Tcm,ident.1 = c("CD8Tu"), ident.2 = c("CD8T"))
x <-tumor_vs_normal.markers[tumor_vs_normal.markers[,5] < 0.01 & tumor_vs_normal.markers[,2] > 0.1 | tumor_vs_normal.markers[,2] < -0.1,]
write.csv(x,file="CD8_Tcm_tumor_VS._CD8normal.markers.csv")

T_subset_del_NK<-subset(T_subset_del_NK,idents=c("CD8T"))
save(T_subset_del_NK,file="T_subset_del_NK/T_subset_del_NK.Rda")

load("D:\\HSM\\MF\\MF-ALL-20200425\\MF-TCR-20200516\\T_cell-20200522\\T_subset_markers_pos_neg.Rdata")
table(T_subset.markers$cluster)
head(T_subset.markers)
T_subset.markers<- T_subset.markers[T_subset.markers$cluster!=12|
                                      T_subset.markers$cluster!=16|
                                      T_subset.markers$cluster!=26,]
T_subset.markers$cluster<-factor(T_subset.markers$cluster,levels = c(0:11,13:15,17:25),labels = 0:23)
x <- T_subset.markers[T_subset.markers$avg_logFC > 0.5,] %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
specific_genes <- unique(x$gene)

pdf(file = "Dotplot_T_subset_features_genes.pdf", width = 20, height = 6)
# Dotplot_cols <- RColorBrewer::brewer.pal(4,"Set2")

DotPlot(object = T_subset_del_NK, cols = c("#87CEFF", "#FF7F00"),
        features = rev(specific_genes), col.max = 2, col.min = -1.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5))

dev.off()

Idents(T_subset_del_NK) <- "seurat_clusters"
pdf(file="UMAP_number.pdf", width = 8, height = 6)
DimPlot(object = T_subset_del_NK, 
        reduction = "tsne",
        label = TRUE, 
        label.size = 6,
        pt.size = 0.1,
        cols= cell_type_cols)
DimPlot(object = T_subset_del_NK, 
        label = TRUE, 
        label.size = 6,
        pt.size = 0.1,
        cols= cell_type_cols)
dev.off()

pdf(file="UMAP_by_cellType.pdf", width = 12, height = 7)
DimPlot(object = T_subset_del_NK,
        #reduction = "tsne",
        label = T,
        label.size = 3,
        pt.size = 0.1,
        group.by = "cellType",
        cols=cell_type_cols)
DimPlot(object = T_subset_del_NK,
        #reduction = "tsne",
        label = F,
        label.size = 3,
        pt.size = 0.1,
        group.by = "cellType",
        cols=cell_type_cols)
dev.off()

meta<-T_subset_del_NK@meta.data[which(T_subset_del_NK@meta.data$samples=="MF14_v3"|
                                        T_subset_del_NK@meta.data$samples=="MF15_v3"|
                                        T_subset_del_NK@meta.data$samples=="MF17"|
                                        T_subset_del_NK@meta.data$samples=="MF21-1"|
                                        T_subset_del_NK@meta.data$samples=="MF21-2"|
                                        T_subset_del_NK@meta.data$samples=="MF22"|
                                        T_subset_del_NK@meta.data$samples=="MF26"|
                                        T_subset_del_NK@meta.data$samples=="MF28-1"|
                                        T_subset_del_NK@meta.data$samples=="MF28-2"|
                                        T_subset_del_NK@meta.data$samples=="MF30-1"|
                                        T_subset_del_NK@meta.data$samples=="MF30-2"|
                                        T_subset_del_NK@meta.data$samples=="pcALCL1"|
                                        T_subset_del_NK@meta.data$samples=="pcALCL2"),]
meta$samples<-factor(meta$samples,levels = c("MF14_v3","MF15_v3","MF17","MF22","MF26","pcALCL1","pcALCL2","MF21-1","MF21-2","MF28-1","MF28-2","MF30-1","MF30-2"))

pdf(file="bar_percent_Type_samples_tu_nor.pdf", width = 8, height = 6)
ggplot(meta,aes(x=samples,fill=Type))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())
dev.off()

pdf(file="bar_percent_Type_samples.pdf", width = 8, height = 6)
ggplot(T_subset_del_NK@meta.data,aes(x=Type,fill=samples))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())

ggplot(T_subset_del_NK@meta.data,aes(x=samples,fill=Type))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())
dev.off()

table(T_subset_del_NK@meta.data$Subtype)
pdf(file="UMAP_group.by_Subtype.pdf", width = 8, height = 6)
DimPlot(object = T_subset_del_NK, 
        reduction = "umap",
        label = TRUE, 
        label.size = 6,
        group.by = "Subtype",
        pt.size = 0.1,
        cols= cell_type_cols)
DimPlot(object = T_subset_del_NK, 
        reduction = "umap",
        label = TRUE, 
        label.size = 4,
        group.by = "Subtype_2",
        pt.size = 0.1,
        cols= cell_type_cols[5:10])
dev.off()

pdf(file="UMAP_group.by_Type.pdf", width = 8, height = 6)
DimPlot(object = T_subset_del_NK, 
        reduction = "umap",
        label = TRUE, 
        label.size = 6,
        group.by = "Type",
        pt.size = 0.1,
        cols= cell_type_cols)
dev.off()

pdf(file="UMAP_group.by_sample.pdf", width = 8, height = 6)
DimPlot(T_subset_del_NK, 
        label = F, 
        label.size = 4,
        pt.size = 0.1,
        group.by = "samples",
        order = rev(c("MF14_v3","MF15_v3", "MF17", "MF18", "MF21-1", "MF21-2","MF22", "MF26", 
                      "MF27-1", "MF27-2", "MF28-1", "MF28-2","MF30-1", "MF30-2", "pcALCL1", "pcALCL2")),
        cols= cell_type_cols)
DimPlot(T_subset_del_NK, 
        label = T, 
        label.size = 4,
        pt.size = 0.1,
        group.by = "samples",
        order = rev(c("MF14_v3","MF15_v3", "MF17", "MF18", "MF21-1", "MF21-2","MF22", "MF26", 
                      "MF27-1", "MF27-2", "MF28-1", "MF28-2","MF30-1", "MF30-2", "pcALCL1", "pcALCL2")),
        cols= cell_type_cols)
dev.off()

pdf(file="Umap_clonotype-all.pdf", width = 7, height = 6)
DimPlot(object = T_subset_del_NK,
        reduction="umap",
        label = F, 
        pt.size = 0.1,
        group.by = "raw_clonotype_id")+NoLegend()
dev.off()

pdf(file="Fraction_by_is_cell_cluster.pdf", width = 8,height = 6)
ggplot(T_subset_del_NK@meta.data, aes(seurat_clusters,fill=is_cell))+ 
  xlab(NULL) + ylab(NULL) + 
  theme(legend.position = "right") +
  geom_bar(position = "fill")+
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title=element_text(hjust=0.5),legend.position = "right")+
  labs(fill="TRC information",title="Fraction of TCR in T cells") +  #改变图例的标题
  scale_fill_manual(values=rev(c("#529DC6", "grey80")),labels=c("TCR loss", "With TCR")) #设置图例的标签

ggplot(T_subset_del_NK@meta.data, aes(is_cell,fill=seurat_clusters))+ 
  xlab(NULL) + ylab(NULL) + 
  theme(legend.position = "right") +
  geom_bar(position = "fill") +
  scale_fill_manual(values=cell_type_cols) + theme_bw() +
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

meta<-T_subset_del_NK@meta.data[,c("samples","is_cell")]
meta<-subset(meta,meta$samples!="MF14_v3" & meta$samples!="MF15_v3")

pdf(file="Fraction_by_is_cell_samples.pdf", width = 8,height = 6)
ggplot(meta, aes(samples,fill=is_cell))+ 
  xlab(NULL) + ylab(NULL) +
  geom_bar(position = "fill")+
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title=element_text(hjust=0.5), legend.position = "right")+
  labs(fill="TRC information",title="Fraction of TCR in T cells") +  #改变图例的标题
  scale_fill_manual(values=rev(c("#529DC6", "grey80")),labels=c("TCR loss", "With TCR")) #设置图例的标签
dev.off()

pdf(file="Fraction_by_cellType_samples.pdf", width = 10,height = 8)
ggplot(T_subset_del_NK@meta.data, aes(cellType,fill=samples))+ 
  xlab(NULL) + ylab(NULL) + 
  theme(legend.position = "right") +
  geom_bar(position = "fill")+
  scale_fill_manual(values=cell_type_cols) + theme_bw() +
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), axis.text.y = element_text(size = 12), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + coord_flip()
dev.off()


pdf(file="vlnplot_features_count.pdf", width = 20, height = 4)
VlnPlot(T_subset_del_NK, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        pt.size = 0,
        cols = cell_type_cols)

VlnPlot(T_subset_del_NK, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "samples",
        pt.size = 0,
        cols = cell_type_cols)
dev.off()

pdf(file="FeaturePlot_features_count.pdf", width = 22, height = 6)
FeaturePlot(T_subset_del_NK, 
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
            ncol = 3, 
            pt.size = 0.1,
            cols = c("#377EB8","#F781BF"))
dev.off()


genes<-c("PTPRC","CD3E","CD4","CD8A","CD8B","CD7","TOX","FOXP3","TCF7","CCR7", "IL7R","IFNG",
         "GZMB","NKG7","GNLY","EPCAM","PDCD1","CTLA4","TIGIT","CD79A","MS4A1")

pdf(file = "Feature_specific_markers.pdf", width = 18, height = 12)
FeaturePlot(T_subset_del_NK, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=5,
            features =genes)
dev.off()

pdf(file="./Vlnplot_specific_markers.pdf", width =23, height = 14)
VlnPlot(object = T_subset_del_NK, 
        cols = cell_type_cols,
        ncol = 4,
        pt.size = 0,
        features = genes)  + NoLegend()

dev.off()

pdf(file = "Feature_HLA_markers.pdf", width = 10, height = 3)
FeaturePlot(T_subset_del_NK, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=3,
            features =c("HLA-DPB1","HLA-E","MIF"))
dev.off()

pdf(file="./Vlnplot_HLA_markers.pdf", width =10, height = 3)
VlnPlot(object = T_subset_del_NK, 
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0,
        group.by = "tu.type",
        features = c("HLA-DPB1","HLA-E","MIF"))  + NoLegend()
dev.off()

p1 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              #group.by = "RNA_snn_res.0.6",
              features = genes[1]) + 
  theme_bw() + 
  coord_flip() +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank() , 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  NoLegend() + panel_border(color = "black")

p2 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[2]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.y = element_blank() , 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


p3 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[3]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.y = element_blank() , 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p4 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[4]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p5 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[5]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p6 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[6]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


p7 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[7]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p8 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[8]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p9 <- VlnPlot(T_subset_del_NK, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[9]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p10 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[10]) +
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p11 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[11]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),  
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  NoLegend() + panel_border(color = "black")


p12 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[12]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p13 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[13]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p14 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[14]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p15 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[15]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p16 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[16]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p17 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[17]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p18 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[18]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p19 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[19]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p20 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[20]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p21 <- VlnPlot(T_subset_del_NK, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[21]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

pdf(file = "plot_grid_cluster_marker.pdf", width = 22, height =6)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, ncol = 21)
dev.off()

Naive <- list(c('CD62L','CD27','LEF1','IL7R','CCR7','CD28','SELL','LRRN3','NPM1'))
Effect <- list(c("CCL4","CCL5","CCL4L2","GZMK","GZMH","NKG7","GNLY","TNFRSF4","MX1","PRF1","GZMB"))
Memory <- list(c('CD62L','CD27','LEF1','IL7R','CCR7','CD28','SELL'))
Exhausted <- list(c('LAG3', 'TIGIT', 'PDCD1', 'HAVCR2','EOMES','ENTPD1','CD38','CD101', 'CTLA4','CD244','BTLA','CD160','TIM3'))#co-inhibitory receptors
Cytotoxic <- list(c('NKG7', 'PRF1', 'GZMA', 'GZMB', 'GZMK', 'IFNG', 'CCL4', 'CST7',"GZMH","GZMM"))
Apoptosis <- list(c('RIPK1','LMBR1L','BCL10','KDELR1','FASLG','CLC','RIPK1','RIPK1','LMBR1L','LMBR1L','DNAJA3','DNAJA3','AKT1','FAS','IL2RA','SIVA1'))#GO_BP annotations: T cell apoptotic process & activation-induced cell death of T cells
Stem <- list(c('CD34','CD38','TRAF2','TNIK','KIT','BMI1','CD44','CD47','CD96','GLI1','GLI2','IL3RA','CLEC12A','HAVCR2','MSI2','LY6E'))


levels(T_subset_del_NK)
Idents(T_subset_del_NK) <- "cellType"
# T_subset_del_NK <- subset(T_subset_del_NK, idents = c("MF14_v3", "pcALCL"), invert = TRUE)
# Idents(T_subset_del_NK) <- "Type"
# T_subset_del_NK <- subset(T_subset_del_NK, idents = c("T_normal"), invert = TRUE)
head(T_subset_del_NK@meta.data)

T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Effect, name = "Effect_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Memory, name = "Memory_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Exhausted, name = "Exhausted_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Cytotoxic, name = "Cytotoxic_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Apoptosis, name = "Apoptosis_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Stem, name = "Stem_score")


Th1 <- list(c("KLRD1","CXCR3","CCR5","TBX21","STAT1","STAT4","IL2","IFNG","TNF"))
Th2 <- list(c("CCR3","CCR4","CXCR4","GATA3","STAT6","IL4","IL5","IL6","IL13"))
Th9 <- list(c("IRF4","GATA3","CCR6","STAT6","SPI1","IL9","IL10"))
Th17 <- list(c("KLRB1","CCR4","CCR6","RORC","STAT3","IL17A","IL17F","IL21","IL22"))
Th22 <- list(c("CCR4","CCR6","CCR10","AHR","BNC2","FOXO4","STAT3","IL13","IL22","TNF"))
Tfh <- list(c("CD69","CXCR5","CCR7","STAT3","BCL6","IL6","IL10","IL12","IL21"))

T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Th1, name = "Th1_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Th2, name = "Th2_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Th9, name = "Th9_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Th17, name = "Th17_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Th22, name = "Th22_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = Tfh, name = "Tfh_score")

head(T_subset_del_NK@meta.data)

myColors1 <- c(brewer.pal(9, "Set1"),"#8B008B",
               "#FF34B3","#CD5C5C","#BC8F8F","#20B2AA","#00F5FF","#ADFF2F","#FFA500","#FF6A6A","#7FFFD4",
               "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#FF3030",
               "#7CFC00","#000000","#708090")
load("D:\\HSM\\MF\\cell_type_cols.Rda")


pdf("./T_subset_del_NK_score_plots_cellType.pdf", width = 8, height = 4.5)
ggplot(T_subset_del_NK@meta.data, aes(x=Cyto1_score1,y=Exha1_score1,color = cellType))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=myColors1)+
  #geom_smooth(method=lm)+
  theme_bw()
dev.off()


pdf("./T_subset_del_NK_score_plots_Subtype.pdf", width = 6, height = 4.5)
ggplot(T_subset_del_NK@meta.data, aes(x=Cyto1_score1,y=Exha1_score1,color = Tumor_type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=brewer.pal(5,"Set2"))+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=Cyto1_score1,y=Exha1_score1,color = Subtype))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none", legend.title = element_text(size = 12, colour = "blue", face = "bold" ))+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=Cyto1_score1,y=Exha1_score1,color = Subtype_2))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none", legend.title = element_text(size = 12, colour = "blue", face = "bold" ))+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
dev.off()

pdf("./T_subset_del_NK_score_boxplots.pdf", width = 12, height = 6)
ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Effect_score1,fill=cellType,))+
  scale_fill_manual(values=cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-1.5,1.5))+ theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Memory_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-1.0,1.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Exhausted_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,0.5))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Cytotoxic_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-1.3,2.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Apoptosis_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,0.5))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Stem_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,0.5))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=S.Score,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,1.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=G2M.Score,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,1.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Th1_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,1.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Th2_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,1.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Th9_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,0.5))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Th17_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,0.5))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Th22_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,1.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=seurat_clusters,y=Tfh_score1,fill=cellType),)+
  scale_fill_manual(values= cell_type_cols)+
  geom_boxplot(aes(group=seurat_clusters),outlier.colour = NA, width=0.3,color="black")+
  #geom_violin(width = 0.6, color="black")+
  #geom_point(position=position_jitter(width = 0.1,height = 0),alpha=0.1,shape=20,color="c(brewer.pal(9, "Set1")")+
  coord_cartesian(ylim = c(-0.5,1.0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

dev.off()



pdf("./T_subset_del_NK_score_plots_line_facet.pdf", width = 10, height = 8)
ggplot(T_subset_del_NK@meta.data, aes(x=Cyto1_score1,y=Exha1_score1)) +
  stat_density2d(aes(fill=..level..), geom="polygon") +
  facet_wrap(. ~ samples, ncol = 4) + 
  #scale_fill_gradient(low = "white", high = terrain.colors(100))+
  scale_fill_gradient(low = "white", high = rev(heat.colors(100)))+
  #scale_fill_gradientn(colours = terrain.colors(100))+
  #scale_fill_gradient(limits=c(0,2))+
  #xlim(-1.5,3)+
  #ylim(-0.6,2.3)+
  geom_vline(xintercept = 0.35,color="black",lty=2,lwd=0.5)+
  geom_hline(yintercept = -0.16 ,color="black",lty=2,lwd=0.5)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=Cyto1_score1,y=Exha1_score1)) +
  stat_density2d(aes(fill=..level..), geom="polygon") +
  #scale_fill_gradient(low = "white", high = terrain.colors(100))+
  scale_fill_gradient(low = "#F5F5F5", high = "#B22222")+
  #scale_fill_gradientn(colours = terrain.colors(100))+
  #scale_fill_gradient(limits=c(0,2))+
  #xlim(-1.5,3)+
  #ylim(-0.6,2.3)+
  geom_vline(xintercept = -0,color="black",lty=2,lwd=0.5)+
  geom_hline(yintercept = -0.15 ,color="black",lty=2,lwd=0.5)+
  theme_bw()
dev.off()

head(T_subset_del_NK@meta.data)
table(T_subset_del_NK@meta.data$tu.type)
table(T_subset_del_NK@meta.data$Type)
Idents(T_subset_del_NK)<-T_subset_del_NK@meta.data$Type
library(clusterProfiler)
library(org.Hs.eg.db)

tumor_vs_normal.markers<-FindMarkers(T_subset_del_NK,ident.1 = c("Tumor T"), ident.2 = c("Normal T"),logfc.threshold = 0.1)
save(tumor_vs_normal.markers,file="tumor_vs_normal.markers.Rda")
x <-tumor_vs_normal.markers[tumor_vs_normal.markers[,5] < 0.01 & tumor_vs_normal.markers[,2] > 0.1 | tumor_vs_normal.markers[,2] < -0.1,]
write.csv(x,file="tumor VS normal.markers.csv")

pdf(file="EnhancedVolcano_tumor VS normal.pdf", width = 10, height = 8)
EnhancedVolcano::EnhancedVolcano(markers,lab = rownames(markers),
                                 x = "avg_logFC",
                                 y = "p_val_adj",
                                 titleLabSize = 12,
                                 subtitle = '',
                                 #subtitleLabSize = 10,
                                 captionLabSize = 14,
                                 pCutoff=10e-6,
                                 FCcutoff= 0.5,
                                 transcriptPointSize = 1.0,
                                 transcriptLabSize = 4,
                                 title = "tumor VS normal",
                                 colAlpha = 1,
                                 #legend = c("NS","Av_log2FC","Adj_p","Av_log2FC & Adj_p"),
                                 legend = c("NS","Log2 FC","P","P & Log2 FC"),
                                 xlim = c(-1.5, 1.5),
                                 legendPosition = "top",
                                 legendLabSize = 10,
                                 legendIconSize = 3.0,
                                 xlab = bquote("Average"~Log[2]~"FC"),
                                 ylab = bquote(~-Log[10]~adjusted~italic(P))) +
  theme_bw()
dev.off()
#R-下载某一条通路的所有基因名字（GO）
#GO:0050870: positive regulation of T cell activation 
#GO:0002286: T cell activation involved in immune response
#GO:0050862: positive regulation of T cell receptor signaling pathway
#GO:0050868: negative regulation of T cell activation
#GO:0042098: T cell proliferation
#GO:0070234: positive regulation of T cell apoptotic process
library(org.Hs.eg.db)
library(clusterProfiler)
getGOgenes<-function(GO_id){
  e<-mget(GO_id,org.Hs.egGO2ALLEGS)
  a<-e[[1]]
  length(a)
  b<-a[!duplicated(a)]
  length(b)
  eg <- bitr(b, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
  head(eg)
  genelist<-eg$SYMBOL
  return(genelist)
}
T.activation<-getGOgenes("GO:0002286")
write.csv(T.activation,file = "GO_0002286.T cell activation involved in immune response.csv")
T.apoptotic<-getGOgenes("GO:0070234")

intersect(T.activation,T.proliferation)
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = list(T.activation), name = "T.activation.gene_score")
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = list(T.apoptotic), name = "T.apoptotic.gene_score")

load("tumor_vs_normal.markers.Rda")
tumor.high.gene<-tumor_vs_normal.markers[tumor_vs_normal.markers[,5] < 0.01 & tumor_vs_normal.markers[,2] > 0.25,]
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = list(rownames(tumor.high.gene)), name = "tumor.high.gene_score")

normal.high.gene<-tumor_vs_normal.markers[tumor_vs_normal.markers[,5] < 0.01 & tumor_vs_normal.markers[,2] < -0.25,]
T_subset_del_NK <- AddModuleScore(object = T_subset_del_NK, features = list(rownames(normal.high.gene)), name = "normal.high.gene_score")

pdf("./normal.tumor.high.gene_score_plots.pdf", width = 6, height = 4.5)
ggplot(T_subset_del_NK@meta.data, aes(x=normal.high.gene_score1,y=T.activation.gene_score1,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=tumor.high.gene_score1,y=T.activation.gene_score1,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=tumor.high.gene_score1,y=normal.high.gene_score1,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=normal.high.gene_score1,y=G2M.Score,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=tumor.high.gene_score1,y=G2M.Score,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=T.apoptotic.gene_score1,y=G2M.Score,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data, aes(x=T.activation.gene_score1,y=G2M.Score,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
dev.off()
pdf("./T.activation_G2M_score_plots.pdf", width = 6, height = 4.5)
ggplot(T_subset_del_NK@meta.data, aes(x=G2M.Score,y=T.activation.gene_score1,color = Type))+
  geom_point(cex=1, alpha=1)+
  #scale_color_manual()+
  theme(legend.position = "none")+
  scale_color_manual(values=cell_type_cols)+
  #geom_smooth(method=lm)+
  theme_bw()
ggplot(T_subset_del_NK@meta.data,aes(x=Type,y=T.activation.gene_score1,fill=Type))+
  scale_fill_manual(values=c(cell_type_cols))+
  geom_boxplot(aes(group=Type),outlier.colour = NA, width=0.3,color="black")+
  coord_cartesian(ylim = c(-0.1,0.1))+ theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10), 
        legend.title = element_text(size = 12, colour = "blue", face = "bold" ))

ggplot(T_subset_del_NK@meta.data,aes(x=Type,y=G2M.Score,fill=Type))+
  scale_fill_manual(values=c(cell_type_cols))+
  geom_boxplot(aes(group=Type),outlier.colour = NA, width=0.3,color="black")+
  coord_cartesian(ylim = c(-0.15,0.1))+ theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10), 
        legend.title = element_text(size = 12, colour = "blue", face = "bold" ))
dev.off()