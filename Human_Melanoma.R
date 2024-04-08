#Loading Libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(dplyr)

path <- "C:/Users/jcust/OneDrive/Documents/RPractice"
list.files(path)

print("Time to run a Single Cell RNA seq Analysis!!")

#import the data (Human Melanoma)
HM_Data.data <- Read10X_h5("10k_Human_DTC_Melanoma_5p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")

#create seurat object
HM_data <- CreateSeuratObject(counts = HM_Data.data, project = "HU", min.cells = 3, min.features = 200)
print(HM_data)

#quality control
HM_data[["percent.mt"]] <- PercentageFeatureSet(HM_data, pattern = "^MT-")
#visualize HM data up until this point
#Visualize QC metrics as a violin plot
VlnPlot(HM_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(HM_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HM_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filtering of the data
HM_data <- subset(HM_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalizing the data the values in the normalized line are default
HM_data <- NormalizeData(HM_data, normalization.method = "LogNormalize", scale.factor = 10000)


#Identification of highly variable features
HM_data <- FindVariableFeatures(HM_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(HM_data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(HM_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#scaling of the data
all.genes <- rownames(HM_data)
HM_data <- ScaleData(HM_data, features = all.genes)


#PCA
HM_data <- RunPCA(HM_data, features = VariableFeatures(object = HM_data))
print(HM_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HM_data, dims = 1:2, reduction = "pca")
DimPlot(HM_data, reduction = "pca") + NoLegend()
DimHeatmap(HM_data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(HM_data, dims = 1:17, cells = 500, balanced = TRUE)


#determine dimensionality of the data
ElbowPlot(HM_data)


#clustering
HM_data <- FindNeighbors(HM_data, dims = 1:20)
HM_data <- FindClusters(HM_data, resolution = 0.5)


#run UMAP
HM_data <- RunUMAP(HM_data, dims = 1:20)
DimPlot(HM_data, reduction = "umap")

#save
saveRDS(HM_data, file = 'Human_Melanoma.RDS')


#markers to aid with cluster identification
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
HM_data.markers <- FindAllMarkers(HM_data), only.pos = TRUE)
HM_data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
