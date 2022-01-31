library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(pheatmap)
setwd("E:/MF/scRNA-seq/MF-ALL-20200425/MF-TCR-20200516/other.cell.non.T")
load("other.cell.non.T.cellType.Rda")
save(other.cell.non.T,file = "other.cell.non.T.cellType.Rda")
load("E:/MF/scRNA-seq/MF-ALL/10X/MF_cellType.Rda")
cell_type_cols <- c(brewer.pal(9, "Set1"),
                    "#FF34B3","#CD5C5C","#BC8F8F","#20B2AA","#00F5FF","#FF3030","#ADFF2F","#FFA500","#FF6A6A","#7FFFD4",
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1",
                    "#7CFC00","#000000","#708090")
samples<-MF_TCR@meta.data[match(row.names(other.cell.non.T@meta.data),row.names(MF_TCR@meta.data)),c("orig.ident","samples")]
other.cell.non.T<-AddMetaData(other.cell.non.T,samples)

current.cluster.ids <- names(table(other.cell.non.T@meta.data$samples))
new.cluster.ids <- c("TCYEM","TCYEM","TCYEM","no_tumor","TCM",
                     "TCM","TCYEM","TCYEM","no_tumor","no_tumor",
                     "TCM","TCM","TCM","TCM","TCYEM","TCYEM")

other.cell.non.T@meta.data$cytotoxic_Tcm <- plyr::mapvalues(x = other.cell.non.T@meta.data$samples,
                                                           from = current.cluster.ids,
                                                           to = new.cluster.ids)
Idents(other.cell.non.T)<-"cytotoxic_Tcm"
other.cell.non.T.del.no.tumor<-subset(other.cell.non.T,idents=c("TCYEM","TCM"))
save(other.cell.non.T.del.no.tumor,file = "del.no.tumor/other.cell.non.T.del.no.tumor.Rda")

pdf(file="UMAP_cytotoxic_Tcm.pdf", width = 8, height = 6)
DimPlot(other.cell.non.T, 
        label = F, 
        label.size = 6,
        pt.size = 0.1,
        group.by = "cytotoxic_Tcm",
        cols= cell_type_cols)
dev.off()

head(MF_TCR@meta.data)

current.cluster.ids <- names(table(MF_TCR@meta.data$cellType))
new.cluster.ids <- c("Tumor T","Normal T","B cells","Macrophages","Myofibroblasts",
                     "Endothelial cells","plasma cells","ILC1","NK","Epidermis","pDC","Mast/Melanocytes")

MF_TCR@meta.data$cellType.new <- plyr::mapvalues(x = MF_TCR@meta.data$cellType,
                                                  from = current.cluster.ids,
                                                  to = new.cluster.ids)
MF_TCR@meta.data$cellType.new <- factor(MF_TCR@meta.data$cellType.new,
                                         levels = c("Tumor T","Normal T","B cells","Macrophages","Myofibroblasts",
                                                    "Endothelial cells","plasma cells","ILC1","NK","Epidermis","pDC","Mast/Melanocytes"))
table(MF_TCR@meta.data$cellType.new)

Idents(MF_TCR)<-"cellType.new"
table(Idents(MF_TCR))
other.cell.non.T<-subset(MF_TCR,idents = c("B cells","Macrophages","Myofibroblasts",
                                           "Endothelial cells","plasma cells","ILC1",
                                           "NK","Epidermis","pDC","Mast/Melanocytes"))
table(other.cell.non.T@meta.data$cellType.new)

#re-group
other.cell.non.T <- RunPCA(object = other.cell.non.T,  
                        features = VariableFeatures(object = other.cell.non.T), 
                        npcs = 50, 
                        verbose = FALSE)
ElbowPlot(other.cell.non.T,ndims = 50)
# uMAP/t-SNE and Clustering
other.cell.non.T <- FindNeighbors(other.cell.non.T, reduction = "pca", dims = 1:30)
other.cell.non.T <- FindClusters(other.cell.non.T, resolution = 0.2)
other.cell.non.T <- RunUMAP(other.cell.non.T, reduction = "pca", dims = 1:30)
other.cell.non.T <- RunTSNE(other.cell.non.T, reduction = "pca", dims = 1:30)
pdf(file="UMAP_number.pdf", width = 8, height = 6)
DimPlot(object = other.cell.non.T, 
        label = TRUE, 
        label.size = 6,
        pt.size = 0.1,
        cols= cell_type_cols)
dev.off()

#B: 0,1,3,16;Mac:2,7;DC:12;MyoFib:4;Fib:9;Endo:5;Plasma:8;ILC1/NK:6,14;Epi:10;pDC:11;Mast:13;Melanocyte:15
current.cluster.ids <- names(table(other.cell.non.T@meta.data$seurat_clusters))
new.cluster.ids <- c("B cells","B cells","Macrophages","B cells","Myofibroblasts","Endothelial cells","ILC1/NK",
                     "Macrophages","Plasma cells","Fibroblasts","Epithelial cells","pDC","DC","Mast cells",
                     "ILC1/NK","Melanocyte","B cells")
other.cell.non.T@meta.data$cellType.new <- plyr::mapvalues(x = other.cell.non.T@meta.data$seurat_clusters,
                                                  from = current.cluster.ids,
                                                  to = new.cluster.ids)
table(other.cell.non.T@meta.data$cellType.new)


NK.ILC1<-subset(other.cell.non.T,idents = c(6,14))
table(Idents(NK.ILC1))
NK.ILC1 <- RunPCA(object = NK.ILC1,  
                           features = VariableFeatures(object = other.cell.non.T), 
                           npcs = 50, 
                           verbose = FALSE)
ElbowPlot(NK.ILC1,ndims = 50)
# uMAP/t-SNE and Clustering
NK.ILC1 <- FindNeighbors(NK.ILC1, reduction = "pca", dims = 1:30)
NK.ILC1 <- FindClusters(NK.ILC1, resolution = 0.2)
NK.ILC1 <- RunUMAP(NK.ILC1, reduction = "pca", dims = 1:30)
NK.ILC1 <- RunTSNE(NK.ILC1, reduction = "pca", dims = 1:30)
pdf(file="NK.ILC1.UMAP_number.pdf", width = 4, height = 3)
DimPlot(object = NK.ILC1, 
        label = TRUE, 
        label.size = 6,
        pt.size = 0.1,
        cols= cell_type_cols)
dev.off()
pdf(file = "NK.ILC1.Feature_specific_markers.pdf", width = 10, height = 8)
FeaturePlot(NK.ILC1, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=3,
            features =c("NKG7","GNLY","KLRD1","NCAM1","FCGR3A","IL7R"))
dev.off()
#0:ILC1;  1,2:NK
current.cluster.ids <- names(table(NK.ILC1@meta.data$seurat_clusters))
new.cluster.ids <- c("ILC1","NK","NK")
NK.ILC1@meta.data$cellType.new <- plyr::mapvalues(x = NK.ILC1@meta.data$seurat_clusters,
                                                 from = current.cluster.ids,
                                                 to = new.cluster.ids)
table(NK.ILC1@meta.data$cellType.new)

NK.ILC1.meta.data<-NK.ILC1@meta.data[,c("samples","cellType.new")]
NK.ILC1.meta.data$cellType.new<-as.character(NK.ILC1.meta.data$cellType.new)
table(NK.ILC1.meta.data$cellType.new)


other.cell.non.T.meta.data<-other.cell.non.T@meta.data[,c("samples","cellType.new")]
other.cell.non.T.meta.data$cellType.new<-as.character(other.cell.non.T.meta.data$cellType.new)
table(other.cell.non.T.meta.data$cellType.new)

other.cell.non.T.meta.data[match(row.names(NK.ILC1.meta.data),row.names(other.cell.non.T.meta.data)),"cellType.new"]<-NK.ILC1.meta.data$cellType.new
table(other.cell.non.T.meta.data$cellType.new)
other.cell.non.T<-AddMetaData(other.cell.non.T,other.cell.non.T.meta.data)
table(other.cell.non.T@meta.data$cellType.new)
other.cell.non.T@meta.data$cellType.new <- factor(other.cell.non.T@meta.data$cellType.new,
                                        levels = c("B cells","Macrophages","Myofibroblasts","Endothelial cells","Plasma cells",
                                                   "Fibroblasts","ILC1","NK","Epithelial cells","pDC","DC","Mast cells","Melanocyte"))
save(other.cell.non.T,file = "other.cell.non.T.cellType.Rda")

pdf(file="UMAP_number.cellType.pdf", width = 8, height = 6)
DimPlot(object = other.cell.non.T, 
        label = TRUE, 
        label.size = 4,
        group.by = "cellType.new",
        pt.size = 0.1,
        cols= cell_type_cols)
dev.off()

table(other.cell.non.T@meta.data$cellType.new,other.cell.non.T@meta.data$seurat_clusters)

Idents(other.cell.non.T)<-"cellType.new"
other.cell.non.T.markers<-FindAllMarkers(other.cell.non.T, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
save(other.cell.non.T.markers, file = "other.cell.non.T_markers_pos_neg.cellType.new.Rdata")

top10 <- other.cell.non.T.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf(file="DoHeatmap.cellType.new.pdf", width = 10, height = 8)
DoHeatmap(other.cell.non.T, features = top10$gene,group.by="cellType.new", label = F,group.colors =cell_type_cols) +
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                       mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")
dev.off()

x <- other.cell.non.T.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(x,file="other.cell.non.T_top50_markers_pos_neg.cellType.new.csv")

x <- other.cell.non.T.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv(x,file="other.cell.non.T_top100_markers_pos_neg.cellType.new.csv")

x <- other.cell.non.T.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
write.csv(x,file="other.cell.non.T_top200_markers_pos_neg.cellType.new.csv")

x <- other.cell.non.T.markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_logFC)
write.csv(x,file="other.cell.non.T_top500_markers_pos_neg.cellType.new.csv")


genes<-c("PTPRC","CD19","CD79A","MS4A1","CD40","CD14","LAMP3","CD68","LYZ","NKG7","GNLY",
         "KLRD1","FCGR3A","IL7R","MZB1","IGLL5","IRF7","LILRA4","CMA1", "MS4A2",
         "VWF","ENG","EPCAM","JUP","KRT5","KRT14","PMEL","TYRP1","FAP","ACTA2")

pdf(file = "Feature_specific_markers.pdf", width = 18, height = 18)
FeaturePlot(other.cell.non.T, 
            #reduction = "tsne",
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol=5,
            features =genes)
dev.off()

pdf(file="./Vlnplot_specific_markers.pdf", width =23, height = 18)
VlnPlot(object = other.cell.non.T, 
        cols = cell_type_cols,
        ncol = 5,
        pt.size = 0,
        features = genes)  + NoLegend()

dev.off()


p1 <- VlnPlot(other.cell.non.T, 
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

p2 <- VlnPlot(other.cell.non.T, 
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


p3 <- VlnPlot(other.cell.non.T, 
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

p4 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[4]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p5 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[5]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p6 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[6]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


p7 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[7]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p8 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[8]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p9 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[9]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p10 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[10]) +
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p11 <- VlnPlot(other.cell.non.T, 
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


p12 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[12]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p13 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[13]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p14 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[14]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p15 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[15]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p16 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[16]) + 
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

p17 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[17]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p18 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[18]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p19 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[19]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p20 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[20]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


p21 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[21]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p22 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[22]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p23 <- VlnPlot(other.cell.non.T, 
              cols = cell_type_cols,
              pt.size = 0,
              features = genes[23]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p24 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[24]) +
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p25 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[25]) + 
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


p26 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[26]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p27 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[27]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p28 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[28]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p29 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[29]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")

p30 <- VlnPlot(other.cell.non.T, 
               cols = cell_type_cols,
               pt.size = 0,
               features = genes[30]) + 
  theme_bw() + 
  coord_flip() + 
  theme(axis.line.x = element_blank(),          axis.ticks.x = element_blank(), axis.text.x = element_blank(),          axis.title.x = element_blank(),         axis.line.y = element_blank(),          axis.ticks.y = element_blank(),axis.text.y = element_blank() ,          axis.title.y = element_blank(),         plot.title = element_text(hjust = 0.5))+ 
  NoLegend() + panel_border(color = "black")


pdf(file = "plot_grid_cluster_marker.pdf", width = 16, height =6)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,
          p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30, ncol = 15)
dev.off()

