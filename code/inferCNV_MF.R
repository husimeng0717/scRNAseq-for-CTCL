library(NGCHM)
library(infercnv)
library(Seurat)
library(tidyverse)

# data = Read10X(data.dir = "10x_data_dir/") 
# seurat_obj = CreateSeuratObject(raw.data=data, min.cells=3, min.genes=200)
# singleCell.counts.matrix = as.matrix(seurat_obj@raw.data[,seurat_obj@cell.names])

# use more palatable column names (cell identifiers)            
# cell.names <- sapply(seq_along(colnames(singleCell.counts.matrix)), function(i) paste0("cell_", i), USE.NAMES = F)      
# colnames(singleCell.counts.matrix) = cell.names    

# write the output table 
write.table(round(singleCell.counts.matrix, digits=3), file='singleCell.counts.matrix.txt', quote=F, sep="\t")   


load("C:/Users/aking/Desktop/MF-workspace/MF_all_integrated/res0.8/T_cell_subset/20191127-Tumor_vs_normal/Tcells_subset_TCR_20191127.Rda")
DefaultAssay(Tcells_subset_TCR) <- "RNA"
singleCell.counts.matrix <- as.matrix(Tcells_subset_TCR@assays$RNA[,])
write.table(singleCell.counts.matrix, file = "singleCell.counts.matrix.txt", sep = "\t")
cellAnnotations <- Tcells_subset_TCR@meta.data[,c(2,11)]
names(cellAnnotations)[1] <- "Single_cell"
head(cellAnnotations)
cellAnnotations$Single_cell <- row.names(cellAnnotations)
head(cellAnnotations)
write.table(cellAnnotations, file = "cellAnnotations.txt",sep = "\t", row.names = F, col.names = F)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="singleCell.counts.matrix.txt",
                                    annotations_file="cellAnnotations.txt",
                                    delim="\t",
                                    gene_order_file="C:/Users/aking/Desktop/MF-workspace/MF_all_integrated/res0.8/T_cell_subset/INFERcnv/gene_position.txt",
                                    ref_group_names=c("5","8","9"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff= 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups= T,   # cluster
                             denoise=T,
                             HMM=T)
###########################################################################################################
###########################################################################################################
