#setting the directory
setwd("C:/Users/S236282/Desktop/Itaconate data/Young Aged Data")

set.seed(1234)

#loading the libraries
library(Seurat)
library("hdf5r")
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggpubr)
library(rstatix)
library(dittoSeq)
library(wesanderson)
library(celldex)
library(SingleR)

#reading the rds file
young_aged_data= readRDS("Young_Aged.rds")
View(young_aged_data)

View(young_aged_data@meta.data)

# Extract gene names 
gene_names_sct <- rownames(GetAssayData(young_data, assay = "SCT", layer = "data"))

# View the first few gene names
head(gene_names_sct)
write.csv(gene_names_sct, file = "gene_names.csv")

Idents(young_aged_data) <- "TissueID"

head(young_aged_data@meta.data$location)
table(young_aged_data@meta.data$location)

#Checking whether the data is preprocessed and labelled
str(young_aged_data)
head(young_aged_data)
names(young_aged_data@assays)

young_aged_data@assays$Spatial@counts  # Raw counts for Spatial assay
young_aged_data@assays$Spatial@data    # Normalized data for Spatial assay

#checking the metadata
if ("meta.data" %in% slotNames(young_aged_data)) {
  head(young_aged_data@meta.data)
} else {
  cat("No meta.data slot found in the object.")
}

#Creating UMAP and Dimplots
VariableFeatures(young_aged_data)
ScaleData(young_aged_data)

#find no. of variance
ElbowPlot(young_aged_data)

# Use 9 PCs for clustering and UMAP
young_aged_data <- FindNeighbors(young_aged_data, reduction = "pca", dims = 1:9)
young_aged_data <- FindClusters(young_aged_data,resolution =c(0.1,0.3,0.5,0.7,1), verbose = FALSE)
young_aged_data <- RunUMAP(young_aged_data, dims = 1:9, min.dist = 0.5)


p1 <- DimPlot(young_aged_data, reduction = "umap", label = T, pt.size = 1.5, label.size = 6) +
  ggtitle("Cluster Identification: Young and Aged")
p1

SpatialDimPlot(young_aged_data, images = "slice1.5", label = TRUE)

#featires to plot
genes = c("Bnip3","Bnip3l")

#Dotplot of the gene of interest
DotPlot(young_aged_data, features = genes, assay = "SCT", group.by = "Groups", col.min = 0) 

#DotPlot(object = young_aged_data, features = genes,
        #cols = c("blue3", "orange"), 
        #group.by = "TissueID", 
        #scale = T,  
        #split.by = 'age') 

####################################### Subsetting the data #############################################

young_data = subset(young_aged_data, subset = age == "young")
aged_data = subset(young_aged_data, subset = age == "aged")

###################################### ONLY YOUNG DATA ##########################################

#Creating UMAP and Dimplots
VariableFeatures(young_data)
ScaleData(young_data)

#find no. of variance
ElbowPlot(young_data)

Idents(young_data) <- "TissueID"
Idents(young_data)

# Use 9 PCs for clustering and UMAP
young_data <- FindNeighbors(young_data, reduction = "pca", dims = 1:9)
young_data <- FindClusters(young_data,resolution = 0.5, verbose = FALSE)
young_data <- RunUMAP(young_data, dims = 1:9, min.dist = 0.5)


DimPlot(young_data, reduction = "umap", group.by = "TissueID", label = T, pt.size = 2.0, label.size = 6) +
  ggtitle("Cluster Identification: Young")

SpatialDimPlot(young_data, images = "slice1.3", label = F, group.by = "TissueID", pt.size.factor = 4.0 )


##################################### ONLY AGED DATA #####################################################

#Creating UMAP and Dimplots
VariableFeatures(aged_data)
ScaleData(aged_data)

#find no. of variance
ElbowPlot(aged_data)

Idents(aged_data) <- "TissueID"
Idents(aged_data)

# Use 9 PCs for clustering and UMAP
aged_data <- FindNeighbors(aged_data, reduction = "pca", dims = 1:9)
aged_data <- FindClusters(aged_data,resolution = 0.5, verbose = FALSE)
aged_data <- RunUMAP(aged_data, dims = 1:9, min.dist = 0.5)


DimPlot(aged_data, reduction = "umap", group.by = "TissueID", label = T, pt.size = 2.0, label.size = 6) +
  ggtitle("Cluster Identification: Aged")

SpatialDimPlot(aged_data, images = "slice1.6", label = TRUE, group.by = "TissueID")

######################################## MODULE ENRICHMENT #############################################

mito = c("Mff","Ambra1","Atg13","Atg14","Becn1","Cdc37","Clec16a","Gba","Hdac6", "Htra2","Huwe1","Map1lc3b","Optn","Phb2","Pink1",
         "Retreg1","Rnf41","Slc25a4","Slc25a5","Smurf1","Sqstm1","Taz","Tigar","Tomm7","Trp53","Ulk1","Usp30","Vdac1","Vps13c","Vps13d")

fiss = c("Gdap1","Opa3","Triap1","Bak1","Bax","Dnm1l","Fis1","Sumo2","Sumo3","Sumo1","Dnm1")

fus = c("Mfn1","Mfn2","Opa1","Mtch2","Mtch1","Mief2","Mief1")

#module list
mito_list = list(Mitophagy = mito, Fission = fiss, Fusion = fus)

# Run module analysis
module_all <- AddModuleScore(young_aged_data, features = mito_list, assay = "SCT", name = "cluster")

head(module_all@meta.data)

DotPlot(module_all, features = c("cluster1", "cluster2", "cluster3"), group.by = "Groups", col.min = 0, scale.min = 0) +
  theme(axis.text.x = element_text(angle = 30)) +
  ggtitle("Module Enrichment: Young and Aged") +
  scale_x_discrete(labels = c("Mitophagy", "Fission", "Fusion"))


SpatialFeaturePlot(module_young, features = unlist(mito_list), images = "slice1.3")


