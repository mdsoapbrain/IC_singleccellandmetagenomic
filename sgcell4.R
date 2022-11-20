rm(list=ls())
# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("/Users/danny/sgcell/GSE160936_RAW_simple")


# load libraries
#install.packages('Seurat')

library(dplyr)
library(patchwork)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SeuratObject)
library(cowplot)

# get data location
dirs <- list.dirs(path = '/Users/danny/sgcell/GSE160936_RAW_simple',
                  recursive = F, full.names = F)
for(x in dirs){
  name <-  x
  
  cts <- ReadMtx(mtx = paste0('/Users/danny/sgcell/GSE160936_RAW_simple/',x,'/matrix.mtx.gz'),
                 features = paste0('/Users/danny/sgcell/GSE160936_RAW_simple/',x,'/features.tsv.gz'),
                 cells = paste0('/Users/danny/sgcell/GSE160936_RAW_simple/',x,'/barcodes.tsv.gz'))
  
# create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}
ls()
# merge datasets
merged_seurat <- merge(S44_SSC_ND, y = c(S45_SSC_AD, S46_EC_ND, S47_EC_AD),
                       add.cell.ids = ls()[4:7],
                       project = 'ADND')
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('sample','Tissue','Disease','Barcode'), 
                                    sep = '_')



# 1. QC -------
View(merged_seurat@meta.data)
# % MT reads
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
View(merged_seurat@meta.data)

VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)

# split the object by dataset
brain.list <- SplitObject(merged_seurat, split.by = "Disease")

# perform standard preprocessing on each object

for (i in 1:length(brain.list)) {
  brain.list[[i]] <- NormalizeData(brain.list[[i]])
  brain.list[[i]] <- FindVariableFeatures(
    brain.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )}
# find anchors
brain.anchors <- FindIntegrationAnchors(object.list = brain.list, dims = 1:20)
# integrate combine data
brain.integrated <- IntegrateData(anchorset = brain.anchors, dims = 1:20)

#
DefaultAssay(brain.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
brain.integrated <- ScaleData(brain.integrated)

brain.integrated <- RunPCA(brain.integrated, npcs = 30)
brain.integrated <- RunUMAP(brain.integrated, reduction = "pca", dims = 1:30)
brain.integrated <- FindNeighbors(brain.integrated, reduction = "pca", dims = 1:20)
brain.integrated <- FindClusters(brain.integrated, resolution = 0.5)

p1 <- DimPlot(brain.integrated, reduction = "umap", group.by = "Disease")
p2 <- DimPlot(brain.integrated, reduction = "umap", label = TRUE, 
              repel = TRUE) 
plot_grid(p1, p2)
DimPlot(brain.integrated, reduction = "umap", split.by = "Disease")

DefaultAssay(brain.integrated) <- 'RNA'

# findConserved markers -------------
markers_cluster0 <- FindConservedMarkers(brain.integrated,
                                         ident.1 = 0,
                                         grouping.var = 'Disease')

head(markers_cluster0)
brain.integrated <- RenameIdents(brain.integrated, `0` = 'Microglia-1')
DimPlot(brain.integrated, reduction = 'umap', label = T)
{
brain.integrated <- RenameIdents(brain.integrated, `2` = 'Microglia-2')
brain.integrated <- RenameIdents(brain.integrated, `3` = 'Microglia-3')


markers_cluster1 <- FindConservedMarkers(brain.integrated,
                                         ident.1 = 1,
                                         grouping.var = 'Disease')

head(markers_cluster1)
brain.integrated <- RenameIdents(brain.integrated, `1` = 'Astrocytes-1')
DimPlot(brain.integrated, reduction = 'umap', label = T)

brain.integrated <- RenameIdents(brain.integrated, `4` = 'Astrocytes-2')
brain.integrated <- RenameIdents(brain.integrated, `5` = 'Astrocytes-3')
DimPlot(brain.integrated, reduction = 'umap', label = T)
######
markers_cluster6 <- FindConservedMarkers(brain.integrated,
                                         ident.1 = 6,
                                         grouping.var = 'Disease')
head(markers_cluster6)
brain.integrated <- RenameIdents(brain.integrated, `6` = 'Oligodendrocytes')

markers_cluster7 <- FindConservedMarkers(brain.integrated,
                                         ident.1 = 7,
                                         grouping.var = 'Disease')
head(markers_cluster7)
brain.integrated <- RenameIdents(brain.integrated, `Oligodendrocytes` = 'Oligodendrocytes-1')
brain.integrated <- RenameIdents(brain.integrated, `7` = 'Oligodendrocytes-2')
DimPlot(brain.integrated, reduction = 'umap', label = T)


markers_cluster9 <- FindConservedMarkers(brain.integrated,
                                         ident.1 = 9,
                                         grouping.var = 'Disease')
head(markers_cluster9)
brain.integrated <- RenameIdents(brain.integrated, `9` = 'Endothelial cells')
DimPlot(brain.integrated, reduction = 'umap', label = T)
######
markers_cluster10 <- FindConservedMarkers(brain.integrated,
                                         ident.1 = 10,
                                         grouping.var = 'Disease')
head(markers_cluster10)

brain.integrated <- RenameIdents(brain.integrated, `10` = 'Fibroblasts')
DimPlot(brain.integrated, reduction = 'umap', label = T)
#######
markers_cluster11 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 11,
                                          grouping.var = 'Disease')
head(markers_cluster11)

brain.integrated <- RenameIdents(brain.integrated, `11` = 'Erythroid cells')
DimPlot(brain.integrated, reduction = 'umap', label = T)

#######
markers_cluster12 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 12,
                                          grouping.var = 'Disease')
head(markers_cluster12)

brain.integrated <- RenameIdents(brain.integrated, `12` = 'NK-cells')
DimPlot(brain.integrated, reduction = 'umap', label = T)

#######
markers_cluster13 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 13,
                                          grouping.var = 'Disease')
head(markers_cluster13)

brain.integrated <- RenameIdents(brain.integrated, `13` = 'Oligodendrocytes-3')
DimPlot(brain.integrated, reduction = 'umap', label = T)

#######
markers_cluster14 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 14,
                                          grouping.var = 'Disease')
head(markers_cluster14)

brain.integrated <- RenameIdents(brain.integrated, `14` = 'Smooth muscle cells')
DimPlot(brain.integrated, reduction = 'umap', label = T)


#######
markers_cluster15 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 15,
                                          grouping.var = 'Disease')
head(markers_cluster15)

brain.integrated <- RenameIdents(brain.integrated, `15` = 'Neurons')
DimPlot(brain.integrated, reduction = 'umap', label = T)

#######
markers_cluster16 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 16,
                                          grouping.var = 'Disease')
head(markers_cluster16)

brain.integrated <- RenameIdents(brain.integrated, `OPCs` = 'OPCs-1')
DimPlot(brain.integrated, reduction = 'umap', label = T)

#######
markers_cluster17 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 17,
                                          grouping.var = 'Disease')
head(markers_cluster17)

brain.integrated <- RenameIdents(brain.integrated, `\tSmooth muscle cells-2` = 'Smooth muscle cells-2')
DimPlot(brain.integrated, reduction = 'umap', label = T)

#######
markers_cluster8 <- FindConservedMarkers(brain.integrated,
                                          ident.1 = 8,
                                          grouping.var = 'Disease')
head(markers_cluster8)

brain.integrated <- RenameIdents(brain.integrated, `8` = 'OPCs-2')
DimPlot(brain.integrated, reduction = 'umap', label = T)

}

# cells already have annotations provided in the metadata
View(brain.integrated@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
#Idents(brain.integrated) <- brain.integrated@meta.data$seurat_annotations
#Idents(brain.integrated)

#DimPlot(brain.integrated, reduction = 'umap', label = TRUE)

brain.integrated$celltype.Disease <- paste(Idents(brain.integrated), brain.integrated$Disease, sep = "_")
brain.integrated$celltype <- Idents(brain.integrated)
Idents(brain.integrated) <- brain.integrated$celltype.Disease
Microglia01.response <- FindMarkers(brain.integrated, ident.1 = "Microglia-1_ND", ident.2 = "Microglia-1_AD", verbose = FALSE)
head(Microglia01.response, n = 5)
{
DimPlot(brain.integrated, reduction = 'umap', label = TRUE)

Microglia02.response <- FindMarkers(brain.integrated, ident.1 
                                    = "Microglia-2_ND", ident.2 = "Microglia-2_AD", verbose = FALSE)
head(Microglia02.response, n = 5)

Microglia03.response <- FindMarkers(brain.integrated, ident.1 
                                    = "Microglia-3_ND", ident.2 = "Microglia-3_AD", verbose = FALSE)
head(Microglia03.response, n = 5)

Astrocytes01.response <- FindMarkers(brain.integrated, ident.1 
                                    = "Astrocytes-1_ND", ident.2 = "Astrocytes-1_AD", verbose = FALSE)
head(Astrocytes01.response, n = 5)

Astrocytes02.response <- FindMarkers(brain.integrated, ident.1 
                                     = "Astrocytes-2_ND", ident.2 = "Astrocytes-2_AD", verbose = FALSE)
head(Astrocytes02.response, n = 5)

Astrocytes03.response <- FindMarkers(brain.integrated, ident.1 
                                     = "Astrocytes-3_ND", ident.2 = "Astrocytes-3_AD", verbose = FALSE)
head(Astrocytes03.response, n = 5)
}
##plot
plots <- VlnPlot(brain.integrated, features = c("XIST","UTY", "SYNDIG1"), split.by = "Disease", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
{
plots <- VlnPlot(brain.integrated, features = c("PTPRE","HSPA1A", "USP9Y"), split.by = "Disease", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

plots <- VlnPlot(brain.integrated, features = c("BCL6","NLGN4Y", "APOE"), split.by = "Disease", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

plots <- VlnPlot(brain.integrated, features = c("NEAT1","GAPDH", "GFAP"), split.by = "Disease", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

plots <- VlnPlot(brain.integrated, features = c("AHCYL1"), split.by = "Disease", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
}
##Feature plot
Idents(brain.integrated) <- factor(Idents(brain.integrated), levels = c("Microglia-1", "Microglia-2", 
                                                                      "Microglia-3", "Astrocytes-1", "Astrocytes-2", "Astrocytes-3"))
markers.to.plot <- c("XIST", "UTY", "SYNDIG1", "PTPRE", "HSPA1A", "USP9Y", "BCL6", "NLGN4Y", "APOE", 
                     "NEAT1", "GAPDH", "GFAP", "AHCYL1")
DotPlot(brain.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Disease") + RotatedAxis()







