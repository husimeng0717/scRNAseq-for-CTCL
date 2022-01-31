library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)

load("E:\\MF\\scRNA-seq\\MF-ALL-20200425\\MF.Rda")
load("D:\\HSM\\MF\\cell_type_cols.Rda")
setwd("D:\\HSM\\MF\\MF-ALL-20200425\\MF-TCR-20200516")
head(MF@meta.data)
MF_TCR@meta.data$raw_clonotype_id[MF_TCR@meta.data$raw_clonotype_id=="None"]<-NA
table(MF_TCR@meta.data$samples)
table(MF_TCR@meta.data$cell_type)
table(MF_TCR@meta.data$cellType2)
table(MF_TCR@meta.data$raw_clonotype_id)

Idents(MF_TCR)<-MF_TCR@meta.data$seurat_clusters
Idents(MF_TCR)


p1<-FeaturePlot(MF_TCR,
                 reduction="umap",
                 cols = c("grey80","red"),
                 features = c("PTPRC"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p2<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("CD3E"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p3<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("ITK"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p4<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("CD7"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p5<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("CD4"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p6<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("CD8A"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p7<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("CD8B"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p8<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("CD79A"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p9<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("MS4A1"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p10<-FeaturePlot(MF_TCR,
                reduction="umap",
                cols = c("grey80","red"),
                features = c("NCAM1"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p11<-FeaturePlot(MF_TCR,
                 reduction="umap",
                 cols = c("grey80","red"),
                 features = c("FCGR3A"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p12<-FeaturePlot(MF_TCR,
                 reduction="umap",
                 cols = c("grey80","red"),
                 features = c("COL1A1"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p13<-FeaturePlot(MF_TCR,
                 reduction="umap",
                 cols = c("grey80","red"),
                 features = c("KRT14"))+
  NoLegend()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p14<-FeaturePlot(MF_TCR,
                 reduction="umap",
                 cols = c("grey80","red"),
                 features = c("EPCAM"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p15<-FeaturePlot(MF_TCR,
                 reduction="umap",
                 cols = c("grey80","red"),
                 features = c("PMEL"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
p16<-FeaturePlot(MF_TCR,
                 reduction="umap",
                 cols = c("grey80","red"),
                 features = c("CD14"))+
  NoLegend()+
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

p<-list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)
pdf(file = "FeaturePlot_cellType_specialGenes.pdf", width = 12, height =10)
CombinePlots(p,ncol = 4)
dev.off()

Idents(MF_TCR)<-MF_TCR$cellType
markers.to.plot <- c("PTPRC","CD3E","CD4","CD40LG","CD8A","CD8B","CD79A","MS4A1","CD14","LYZ","ACTA2" ,"DCN","ENG","VWF",
                     "GNLY","NKG7","MZB1","IGLL5","KRT14","KRT6A","IRF7","LILRA4","TPSAB1","CPA3","PMEL","TYRP1")
pdf(file="DotPlot.pdf", width = 10, height = 4)
DotPlot(MF_TCR, features =rev(markers.to.plot),cols =c("#4B60E1","#FF555E"))+ RotatedAxis()+ #x轴刻度标签倾斜45度
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
      axis.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()

pdf(file="Umap_split.by_samples.pdf", width = 12, height = 10)
DimPlot(MF_TCR,
        reduction="umap",
        label = F,
        split.by = "samples",
        cols = cell_type_cols)+NoLegend()
dev.off()

pdf(file="Umap_clonotype-all.pdf", width = 8, height = 6)
DimPlot(object = MF_TCR,
        reduction="umap",
        label = F, 
        label.size = 8,
        pt.size = 0.5,
        #cols= c("#E41A1C" ,"#377EB8" ,"#4DAF4A", "#984EA3" ,"#FF7F00","#FFFF33", "#A65628", "#F781BF", "#008B8B" ,"#8B008B", "#FF34B3"),
        group.by = "raw_clonotype_id")+NoLegend()
dev.off()

head(MF_TCR@meta.data)
table(MF_TCR@meta.data$clonotype)
MF_TCR@meta.data$clonotype<-MF_TCR@meta.data$raw_clonotype_id
MF_TCR@meta.data[which(MF_TCR@meta.data$raw_clonotype_id!="clonotype1"| 
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype2"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype3"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype4"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype5"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype6"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype7"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype8"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype9"|
                         MF_TCR@meta.data$raw_clonotype_id!="clonotype10"|
                         MF_TCR@meta.data$raw_clonotype_id!="NA"),'clonotype']<-"other clonotype"

MF_TCR@meta.data[which(MF_TCR@meta.data$raw_clonotype_id=="clonotype1"| 
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype2"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype3"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype4"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype5"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype6"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype7"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype8"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype9"|
                         MF_TCR@meta.data$raw_clonotype_id=="clonotype10"),'clonotype']<-MF_TCR@meta.data[which(MF_TCR@meta.data$raw_clonotype_id=="clonotype1"| 
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype2"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype3"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype4"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype5"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype6"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype7"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype8"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype9"|
                                                                                                                  MF_TCR@meta.data$raw_clonotype_id=="clonotype10"),'raw_clonotype_id']
factor_level<-c(paste("clonotype",1:10,sep=""),"other clonotype")
MF_TCR@meta.data$clonotype<-factor(MF_TCR@meta.data$clonotype,levels = factor_level)

pdf(file="Umap_clonotype.pdf", width = 8, height = 6)
DimPlot(object = MF_TCR,
        reduction="umap",
        label = F, 
        label.size = 8,
        pt.size = 0.5,
        cols= c("#E41A1C" ,"#377EB8" ,"#4DAF4A", "#984EA3" ,"#FF7F00","#FFFF33", "#A65628", "#F781BF", "#008B8B" ,"#8B008B","#6AB5E4", "#FF34B3"),
        group.by = "clonotype")
dev.off()

pdf(file="bar_percent_cellType_samples.pdf", width = 8, height = 6)
ggplot(MF_TCR@meta.data,aes(x=cellType,fill=samples))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())

ggplot(MF_TCR@meta.data,aes(x=samples,fill=cellType))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())

ggplot(MF_TCR@meta.data,aes(x=cellType2,fill=samples))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())

ggplot(MF_TCR@meta.data,aes(x=samples,fill=cellType2))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())
dev.off()

MF_TCR.meta<-MF_TCR@meta.data
save(MF_TCR.meta,file="MF_TCR.meta.Rda")
MF_TCR.meta.subset<-MF_TCR.meta[which(MF_TCR.meta$cellType2=="CD4T"| 
                             MF_TCR.meta$cellType2=="CD8T"|
                             MF_TCR.meta$cellType2=="DNT"|
                             MF_TCR.meta$cellType2=="tuT"),]
MF_TCR.meta.subset$cellType2<-factor(MF_TCR.meta.subset$cellType2,levels = c("CD4T","CD8T","DNT","tuT"))
table(MF_TCR.meta.subset$cellType2)  
MF_TCR.meta.subset$samples<-factor(MF_TCR.meta.subset$samples,levels = c("MF14_v3","MF15_v3", "MF17", "MF18", "MF21-1", "MF21-2","MF22", "MF26", 
                                                                         "MF27-1", "MF27-2", "MF28-1", "MF28-2","MF30-1", "MF30-2", "pcALCL1", "pcALCL2"))

pdf(file="bar_percent_cellType_Tumor_type.pdf", width = 8, height = 4)
ggplot(MF_TCR.meta.subset,aes(x=samples,fill=cellType2))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols[5:8])+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())
dev.off()

pdf(file="bar_percent_cellType_Tumor_type.pdf", width = 8, height = 4)
ggplot(MF_TCR@meta.data,aes(x=cellType,fill=Tumor_type))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=c("#98F251","#FFDF00"))+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())
dev.off()

pdf(file="bar_percent_seurat_clusters_samples.pdf", width = 10, height = 6)
ggplot(MF_TCR@meta.data,aes(x=seurat_clusters,fill=samples))+
  geom_bar(position = 'fill',width=0.8)+
  labs(x='',y='')+
  scale_fill_manual(values=cell_type_cols)+
  theme(text = element_text(size=16))+ 
  RotatedAxis()+
  theme(panel.background = element_blank())
dev.off()

pdf(file = "FeaturePlot_Feature_nCount_percent.mt.pdf", width = 14, height = 4)
FeaturePlot(MF_TCR,
            reduction="umap",
            ncol = 3,
            cols = c("#377EB8","#F781BF"),
            features = c("nFeature_RNA","nCount_RNA","percent.mt"))
dev.off()

pdf(file = "VlnPlot_Feature_nCount_percent.mt.pdf", width = 8, height = 8)
VlnPlot(MF_TCR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1,
        group.by = "seurat_clusters",pt.size = 0,cols = cell_type_cols)
dev.off()

#immune cell special genes
pdf(file = "FeaturePlot_immune.pdf", width = 12, height =10)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            features = c("PTPRC","CD2","CD7", "CD3D","CD3E","CD3G","CD4","CD8A","CD8B"))
dev.off()
pdf(file="VlnPlot_immune.pdf", width = 12, height = 10)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features = c("PTPRC","CD2","CD7", "CD3D","CD3E","CD3G","CD4","CD8A","CD8B"))+NoLegend()
dev.off()

#γδT
pdf(file = "FeaturePlot_γδT.pdf", width = 10, height =3)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("TRGC1","TRGC2","TRDC"))
dev.off()
pdf(file="VlnPlot_γδT.pdf", width = 10, height = 3)
VlnPlot(object = MF_TCR, 
        reduction="umap",
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features = c("TRGC1","TRGC2","TRDC"))
dev.off()

#immuno_repress
pdf(file = "FeaturePlot_immuno_repress.pdf", width = 10, height = 6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("PDCD1","CTLA4","EPCAM","FOXP3","CXCL13","CD274"))
dev.off()
pdf(file="VlnPlot_immuno_repress.pdf", width = 10, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        reduction="umap",
        pt.size = 0,
        ncol = 3,
        features = c("PDCD1","CTLA4","EPCAM","FOXP3","CXCL13","CD274"))
dev.off()

#naive
pdf(file = "FeaturePlot_naive.pdf", width = 9, height = 6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("IL7R","CCR7","LEF1","CD27","CD28","SELL"))
dev.off()
pdf(file = "VlnPlot_Naive.pdf", width = 9, height = 6)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("IL7R","CCR7","LEF1","CD27","CD28","SELL"))
dev.off()

#Activated
pdf(file = "FeaturePlot_Activated.pdf", width = 8, height = 4)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("IFNG","FOS","JUN"))
dev.off()
pdf(file = "VlnPlot_Activated.pdf", width = 8, height = 4)
VlnPlot(MF_TCR,
        reduction="umap",
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("IFNG","FOS","JUN"))
dev.off()

#Effector
pdf(file = "FeaturePlot_Effector.pdf", width = 12, height = 8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("CCL4","CCL5","CCL4L2","GZMK","GZMH","NKG7","GNLY","TNFRSF4","MX1","PRF1","GZMB"))
dev.off()
pdf(file = "VlnPlot_effector.pdf", width = 12, height = 8)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("CCL4","CCL5","CCL4L2","GZMK","GZMH","NKG7","GNLY","TNFRSF4","MX1","PRF1","GZMB"))
dev.off()

#Exhausted
pdf(file = "FeaturePlot_Exhausted.pdf", width = 12, height = 8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("PDCD1","CTLA4","LAG3","EOMES","HAVCR2","ENTPD1","CD38","CD101","TIGIT"))
dev.off()
pdf(file = "VlnPlot_exhausted.pdf", width = 12, height = 8)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("PDCD1","CTLA4","LAG3","EOMES","HAVCR2","ENTPD1","CD38","CD101","TIGIT","LAYN"))
dev.off()

#Treg
pdf(file = "FeaturePlot_Treg.pdf", width = 8, height = 6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("CTLA4","IL2RA","FOXP3","TNFRSF18","IKZF2","RUNX1"))
dev.off()
pdf(file = "VlnPlot_Treg.pdf", width = 8, height = 6)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("CTLA4","IL2RA","FOXP3","TNFRSF18","IKZF2","RUNX1"))
dev.off()

#CD8+ cytotoxic T
pdf(file = "FeaturePlot_CD8+ cytotoxic.pdf", width = 8, height = 6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("CD8A","LAMP1","FASLG","GZMB","IFNG","PRF1"))
dev.off()
pdf(file = "VlnPlot_CD8+ cytotoxic.pdf", width = 8, height = 6)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("CD8A","LAMP1","FASLG","GZMB","IFNG","PRF1"))
dev.off()

#Th1
pdf(file = "FeaturePlot_Th1.pdf", width = 12, height = 8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("KLRD1","CXCR3","CCR5","TBX21","STAT1","STAT4","IL2","IFNG","TNF"))
dev.off()
pdf(file = "VlnPlot_Th1.pdf", width = 12, height = 8)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("KLRD1","CXCR3","CCR5","TBX21","STAT1","STAT4","IL2","IFNG","TNF"))
dev.off()

#Th2
pdf(file = "FeaturePlot_Th2.pdf", width = 12, height = 8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("CCR3","CCR4","CXCR4","GATA3","STAT6","IL4","IL5","IL6","IL13"))
dev.off()
pdf(file = "VlnPlot_Th2.pdf", width = 12, height = 8)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("CCR3","CCR4","CXCR4","GATA3","STAT6","IL4","IL5","IL6","IL13"))
dev.off()

#Th17
pdf(file = "FeaturePlot_Th17.pdf", width = 12, height = 8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("KLRB1","CCR4","CCR6","RORC","STAT3","IL17A","IL17F","IL21","IL22"))
dev.off()
pdf(file = "VlnPlot_Th17.pdf", width = 12, height = 8)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("KLRB1","CCR4","CCR6","RORC","STAT3","IL17A","IL17F","IL21","IL22"))
dev.off()

#Th9
pdf(file = "FeaturePlot_Th9.pdf", width = 12, height = 8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("IRF4","GATA3","CCR6","STAT6","SPI1","IL9","IL10"))
dev.off()
pdf(file = "VlnPlot_Th9.pdf", width = 12, height = 8)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features = c("IRF4","GATA3","CCR6","STAT6","SPI1","IL9","IL10"))
dev.off()

#Th22
pdf(file = "FeaturePlot_Th22.pdf", width = 12, height = 12)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("CCR4","CCR6","CCR10","AHR","BNC2","FOXO4","STAT3","IL13","IL22","TNF"))
dev.off()
pdf(file = "VlnPlot_Th22.pdf", width = 12, height = 12)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features =c("CCR4","CCR6","CCR10","AHR","BNC2","FOXO4","STAT3","IL13","IL22","TNF"))
dev.off()

#Tfh
pdf(file = "FeaturePlot_Tfh.pdf", width = 12, height = 8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("CD69","CXCR5","CCR7","STAT3","BCL6","IL6","IL10","IL12","IL21"))
dev.off()
pdf(file = "VlnPlot_Tfh.pdf", width = 12, height = 8)
VlnPlot(MF_TCR,
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0, 
        features =c("CD69","CXCR5","CCR7","STAT3","BCL6","IL6","IL10","IL12","IL21"))
dev.off()

#Fibroblasts cell special genes
pdf(file = "FeaturePlot_Fibroblasts.pdf", width = 10, height =8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("FAP","CXCL12", "PDPN", "COL1A2","DCN", "COL3A1", "COL6A1"))
dev.off()
pdf(file="VlnPlot_Fibroblasts.pdf", width = 10, height = 8)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features = c("FAP","CXCL12", "PDPN", "COL1A2","DCN", "COL3A1", "COL6A1"))
dev.off()

#Myeloid cells
pdf(file = "FeaturePlot_Myeloid.pdf", width = 8, height =4)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("S100A8","S100A9"))
dev.off()
pdf(file="VlnPlot_Myeloid.pdf", width = 8, height = 4)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features =c("S100A8","S100A9"))
dev.off()

#MoMø special genes
pdf(file = "FeaturePlot_MoMac.pdf", width = 10, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("LYZ","CD14", "MS4A7", "CD68","FCGR2A", "CSF1R"))
dev.off()
pdf(file="VlnPlot_MoMac.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features =c("LYZ","CD14", "MS4A7", "CD68","FCGR2A", "CSF1R"))
dev.off()

#Dendritic cells special genes
pdf(file = "FeaturePlot_DC.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("CD40", "CD80","CD83", "CCR7"))
dev.off()
pdf(file="VlnPlot_DC.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("CD40", "CD80","CD83", "CCR7"))
dev.off()

#cDC1s
pdf(file = "FeaturePlot_cDC1s.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("XCR1", "CADM1"))
dev.off()
pdf(file="VlnPlot_cDC1s.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("XCR1", "CADM1"))
dev.off()

#cDC2s
pdf(file = "FeaturePlot_cDC2s.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("XCR1", "CADM1"))
dev.off()
pdf(file="VlnPlot_cDC2s.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("CD1A", "CD172A"))
dev.off()

#Plasmacytoid dendritic cells (pDCs) special genes
pdf(file = "FeaturePlot_pDCs.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("LILRA4","IRF7", "TLR7","IRF4"))
dev.off()
pdf(file="VlnPlot_pDCs.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("LILRA4","IRF7", "TLR7","IRF4"))
dev.off()

#Mast cells special genes
pdf(file = "FeaturePlot_Mast_cells.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("CMA1", "MS4A2","TPSAB1", "TPSB2","FCER1A"))
dev.off()
pdf(file="VlnPlot_Mast_cells.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("CMA1", "MS4A2","TPSAB1", "TPSB2","FCER1A"))
dev.off()

#Myocytes special genes
pdf(file = "FeaturePlot_Myocytes.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("ACTA1", "ACTN2","MYL2", "MYH2"))
dev.off()
pdf(file="VlnPlot_Myocytes.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("ACTA1", "ACTN2","MYL2", "MYH2"))
dev.off()

#Endothelial cells special genes
pdf(file = "FeaturePlot_Endothelial.pdf", width = 10, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("PECAM1","VWF", "ENG","IFI27","HLA-DRB1","NPDC1"))
dev.off()
pdf(file="VlnPlot_Endothelial_cells.pdf", width = 10, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features = c("PECAM1","VWF", "ENG","IFI27","HLA-DRB1","NPDC1"))
dev.off()


#B cells special genes
pdf(file = "FeaturePlot_B_cells.pdf", width = 6, height =4)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("CD19","CD79A", "MS4A1","CD79B"))
dev.off()
pdf(file="VlnPlot_B_cells.pdf", width = 6, height = 4)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("CD19","CD79A", "MS4A1","CD79B"))
dev.off()

#Plasma cells special genes
pdf(file = "FeaturePlot_Plasma_cells.pdf", width = 6, height =4)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("IGKC","IGLC2", "IGHG1"))
dev.off()
pdf(file="VlnPlot_Plasma_cells.pdf", width = 6, height = 4)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features = c("IGKC","IGLC2", "IGHG1"))
dev.off()

#Proliferating
pdf(file = "FeaturePlot_Proliferating.pdf", width = 6, height =2)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("TOP2A", "MKI67"))
dev.off()
pdf(file="VlnPlot_Proliferating.pdf", width = 6, height = 2)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("TOP2A", "MKI67"))
dev.off()

#NK cells special genes
pdf(file = "FeaturePlot_NK.pdf", width = 12, height =8)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("NKG7","GNLY","NCR1","NCAM1"))
dev.off()
pdf(file="VlnPlot_NK.pdf", width = 12, height = 8)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("NKG7","GNLY","NCR1","NCAM1"))
dev.off()

#NKT cells special genes
pdf(file = "FeaturePlot_NKT.pdf", width = 4, height =3)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            features = c("KLRB1"))
dev.off()
pdf(file="VlnPlot_NKT.pdf", width = 4, height = 3)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        features = c("KLRB1"))+NoLegend()
dev.off()

#Epithelial cells special genes
pdf(file = "FeaturePlot_Epithelial.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("EPCAM","JUP","KRT5","KRT14"))
dev.off()
pdf(file="VlnPlot_Epithelial.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("EPCAM","JUP","KRT5","KRT14","EGFR"))
dev.off()

#Neutrophils 
pdf(file = "FeaturePlot_Neutrophils.pdf", width = 8, height =6)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 3,
            features = c("G0S2","FCGR3B","CSF3R","CXCL8","S100A8","S100A9"))
dev.off()
pdf(file="VlnPlot_Neutrophils.pdf", width = 8, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features = c("G0S2","FCGR3B","CSF3R","CXCL8","S100A8","S100A9"))+NoLegend()
dev.off()

#Melanocyte
pdf(file = "FeaturePlot_Melanocyte.pdf", width = 7, height =3)
FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("PMEL","TYRP1"))
dev.off()
pdf(file="VlnPlot_Melanocyte.pdf", width = 7, height = 3)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 3,
        features = c("PMEL","TYRP1"))+NoLegend()
dev.off()

FeaturePlot(MF_TCR,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("IL10"))
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        pt.size = 0,
        ncol = 2,
        features = c("CLEC12B","SFRP1","SFRP2","GZMK"))+NoLegend()

#T: 0,2,3,5,6,7,9,11,13,14,17,18,19,22,24,28,29
#  Normal T: 3,6,7,13,17,19,22,29
#   CD8+: 3,7,17,19,22
#   CD4+: 6
#     Treg: 13
#   CD- CD8-: 29
#  Tumor T: 0,2,5,9,11,14,18,24,28,
# B: 1,10,20,21,30,31,
# plasma B: 16
# Mac: 4,26
# NKT: 15
# pDC: 25  Myofibroblasts: 8,32  Endothelial: 12  Epidermis: 23  Melanocytes/Mast cells: 27

Idents(MF_TCR)<-MF_TCR@meta.data$seurat_clusters
cellType1<-c("0.T cells","1.B cells","2.T cells","3.T cells","4.Macrophages","5.T cells","6.T cells","7.T cells",
            "8.Myofibroblasts","9.T cells","10.B cells","11.T cells","12.Endothelial cells","13.T cells","14.T cells",
            "15.NK","16.plasma B","17.T cells","18.T cells","19.T cells","20.B cells","21.B cells",
            "22.T","23.Epidermis","24.T cells","25.pDC","26.Macrophages","27.Mast/Melanocytes","28.T cells",
            "29.T cells","30.B cells","31.B cells","32.Myofibroblasts")
MF_TCR$cellType1<-factor(MF_TCR@meta.data$seurat_clusters,levels =0:32,labels =cellType1)
table(MF_TCR@meta.data$cellType1,MF_TCR@meta.data$seurat_clusters)
table(MF_TCR@meta.data$cellType1)
pdf(file="Umap_cellType1.pdf", width = 11, height = 5)
DimPlot(MF_TCR,
        label = T,
        reduction="umap",
        cols = cell_type_cols,
        label.size= 4,
        group.by = "cellType1")
dev.off()

cellType2<-c("tuT","B","tuT","CD8T","Mac","tuT","CD4/8T","CD8T",
             "MyoFb","tuT","B","tuT","Endo","Treg","tuT",
             "NK","pB","CD8T","tuT","CD8T","B","B",
             "CD8T","Epi","tuT","pDC","Mac","M/M","tuT",
             "DNT","B","B","MyoFb")
MF_TCR$cellType2<-factor(MF_TCR@meta.data$seurat_clusters,levels =0:32,labels =cellType2)
table(MF_TCR@meta.data$cellType2,MF_TCR@meta.data$seurat_clusters)
table(MF_TCR@meta.data$cellType2)
pdf(file="Umap_cellType2.pdf", width = 7, height = 5)
DimPlot(MF_TCR,
        label = T,
        reduction="umap",
        cols = cell_type_cols,
        label.size= 4,
        group.by = "cellType2")
dev.off()

Idents(MF_TCR)<-MF_TCR@meta.data$cellType2
MF_TCR.cellType2.markers <- FindAllMarkers(MF_TCR, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
save(MF_TCR, file = "MF_TCR_cellType2_markers_pos_neg.Rdata")
load("MF_TCR_cellType2_markers_pos_neg.Rdata")
x <- MMF_TCR_cellType2 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(x,file="MF_TCR_cellType2_top50_markers_pos_neg.csv")

x <- T_subset.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv(x,file="MF_TCR_cellType2_top100_markers_pos_neg.csv")

x <- T_subset.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
write.csv(x,file="MF_TCR_cellType2_top200_markers_pos_neg.csv")

x <- T_subset.markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_logFC)
write.csv(x,file="MF_TCR_cellType2_top500_markers_pos_neg.csv")


#Normal T: 3,6,7,13,17,19,22,29
cellType3<-c("tuT","B","tuT","norT","Mac","tuT","norT","norT",
             "MyoFb","tuT","B","tuT","Endo","norT","tuT",
             "NK","pB","norT","tuT","norT","B","B",
             "norT","Epi","tuT","pDC","Mac","M/M","tuT",
             "norT","B","B","MyoFb")
MF_TCR$cellType3<-factor(MF_TCR@meta.data$seurat_clusters,levels =0:32,labels =cellType3)
table(MF_TCR@meta.data$cellType3,MF_TCR@meta.data$seurat_clusters)
table(MF_TCR@meta.data$cellType3)
pdf(file="Umap_cellType3.pdf", width = 7, height = 5)
DimPlot(MF_TCR,
        label = T,
        reduction="umap",
        cols = cell_type_cols,
        label.size= 4,
        group.by = "cellType3")
dev.off()


cellType<-c("T cells","B cells","T cells","T cells","Macrophages","T cells","T cells","T cells",
            "Myofibroblasts","T cells","B cells","T cells","Endothelial cells","T cells","T cells",
            "NK","plasma B","T cells","T cells","T cells","B cells","B cells",
            "T cells","Epidermis","T cells","pDC","Macrophages","Mast/Melanocytes","T cells",
            "T cells","B cells","B cells","Myofibroblasts")
MF_TCR$cellType<-factor(MF_TCR@meta.data$seurat_clusters,levels =0:32,labels =cellType)
table(MF_TCR@meta.data$cellType,MF_TCR@meta.data$seurat_clusters)
table(MF_TCR@meta.data$cellType)
pdf(file="Umap_cellType.pdf", width = 7.5, height = 5)
DimPlot(MF_TCR,
        label = T,
        reduction="umap",
        cols = cell_type_cols,
        label.size= 4,
        group.by = "cellType")
dev.off()

save(MF_TCR,file="MF_TCR.Rda")


p1 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              #group.by = "RNA_snn_res.0.6",
              features = c("PTPRC")) + 
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

p2 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("CD3E")) + 
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


p3 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("CD4")) + 
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

p4 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("CD8A")) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p5 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("CD8B")) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p6 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("FOXP3")) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


p7 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("ITK")) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p8 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("CD7")) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p9 <- VlnPlot(MF_TCR, 
              cols = cell_type_cols,
              pt.size = 0,
              features = c("CD2")) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p10 <- VlnPlot(MF_TCR, 
               cols = cell_type_cols,
               pt.size = 0,
               features = c("TOX")) +
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p11 <- VlnPlot(MF_TCR, 
               cols = cell_type_cols,
               pt.size = 0,
               features = c("KLRD1")) + 
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


p12 <- VlnPlot(MF_TCR, 
               cols = cell_type_cols,
               pt.size = 0,
               features = c("FCGR3A")) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


pdf(file = "plot_grid_cluster_marker.pdf", width = 18, height =6)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol = 12)
dev.off()

Idents(MF_TCR)<-MF_TCR@meta.data$seurat_clusters
c6<-subset(MF_TCR,idents = 6)
table(MF_TCR@meta.data$cellType2)
DimPlot(c6,
        label = T,
        reduction="umap",
        cols = cell_type_cols,
        label.size= 4,
        group.by = "cellType2")
c6 <- FindNeighbors(object = c6, dims = 1:30)
c6 <- FindClusters(object = c6, resolution = 0.2)
c6 <- RunTSNE(object = c6, dims = 1:30)
c6 <- RunUMAP(c6, dims = 1:30)

save(c6,file="c6.Rda")

DimPlot(c6,
        label = T,
        reduction="umap",
        cols = cell_type_cols,
        label.size= 4)
FeaturePlot(c6,
            reduction="umap",
            cols = c("grey80","red"),
            ncol = 2,
            features = c("CD4","CD8A"))
head(c6@meta.data)
table(c6@meta.data$RNA_snn_res.0.2)
cellType2<-c("CD4T","CD4T","CD8T")
c6$cellType2<-factor(c6@meta.data$seurat_clusters,levels =0:2,labels =cellType2)
table(c6$cellType2)

table(MF_TCR@meta.data$cellType2)
MF_TCR_meta<- MF_TCR@meta.data 
MF_TCR_meta$cellType2<-as.character(MF_TCR_meta$cellType2)
c6_meta<-c6@meta.data
c6_meta$cellType2<-as.character(c6_meta$cellType2)
MF_TCR_meta[which(MF_TCR_meta$seurat_clusters==6),'cellType2']<-c6_meta$cellType2
MF_TCR_meta$cellType2<-as.factor(MF_TCR_meta$cellType2)
head(MF_TCR_meta)
table(MF_TCR_meta$cellType2)

MF_TCR@meta.data<-MF_TCR_meta
pdf(file="Umap_cellType2.pdf", width = 7, height = 5)
DimPlot(MF_TCR,
        label = T,
        reduction="umap",
        cols = cell_type_cols,
        label.size= 4,
        group.by = "cellType2")
dev.off()

head(MF_TCR@meta.data)

meta<-data.frame(Cell=row.names(MF_TCR@meta.data),cell_type=MF_TCR@meta.data$cellType2)
write.csv(meta,file="cellphonedb/meta.csv",row.names = F)
count<-MF_TCR@assays$RNA@data
row.names(count)<-row.names(MF_TCR)
count[1:6,1:10]
write.csv(count,file="cellphonedb/count.csv")


meta<-data.frame(Cell=row.names(MF_TCR@meta.data),cell_type=MF_TCR@meta.data$cellType3)
write.csv(meta,file="cellphonedb/meta.csv",row.names = F)
###########cellphonedb 
db.mean<-read.table("D:\\HSM\\MF\\MF-ALL-20200425\\MF-TCR-20200516\\cellphonedb\\MT_project\\means.txt",sep = "\t",header = T)
db.pvalue<-read.table("D:\\HSM\\MF\\MF-ALL-20200425\\MF-TCR-20200516\\cellphonedb\\MT_project\\pvalues.txt",sep = "\t",header = T)
db.deconvoluted<-read.table("D:\\HSM\\MF\\MF-ALL-20200425\\MF-TCR-20200516\\cellphonedb\\MT_project\\deconvoluted.txt",sep = "\t",header = T)
db.significant_means<-read.table("D:\\HSM\\MF\\MF-ALL-20200425\\MF-TCR-20200516\\cellphonedb\\MT_project\\significant_means.txt",sep = "\t",header = T)

#dot_plot: 
# cellphonedb plot dot_plot 
# --means-path=./out/MT_project/means.txt 
# --pvalues-path=./out/MT_project/pvalues.txt

Idents(MF_TCR) <- "seurat_clusters"

M1<-c('CD80', 'CD86', 'CD68', 'IL1R1', 'TLR2', 'TLR4', 'TNF', 'IL1A', 'IL1B', 'IL6', 'IL12A', 'IL12B',
      'IL23A', 'IL23R', 'PTGS2', 'IL27', 'CXCL9', 'CXCL10', 'CXCL11', 'IL15')

pdf(file = "Feature_M1.pdf", width = 18, height = 12)
FeaturePlot(MF_TCR, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=5,
            features =M1)
dev.off()

pdf(file="./Vlnplot_M1.pdf", width =30, height = 12)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        ncol = 5,
        pt.size = 0,
        features = M1)  + NoLegend()

dev.off()



M2<-c('MRC1', 'TLR1', 'TLR8', 'CCL17', 'CCL24', 'IL10', 'TGFB1', 'CCL13', 'CD209', 'ALOX5')

pdf(file = "Feature_M2.pdf", width = 18, height = 6)
FeaturePlot(MF_TCR, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=5,
            features =M2)
dev.off()

pdf(file="./Vlnplot_M2.pdf", width =30, height = 6)
VlnPlot(object = MF_TCR, 
        cols = cell_type_cols,
        ncol = 5,
        pt.size = 0,
        features = M2)  + NoLegend()

dev.off()

