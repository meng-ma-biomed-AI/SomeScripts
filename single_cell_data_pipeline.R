# scripts to analyze single cell RNA-Seq data 
# data: 10X Human Genomics

library(Seurat)
library(tidyverse)

# load the NSCLC dataset 
nsclc.sparse.m <- Read10X_h5(filename = 'data/20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m)
cts <- nsclc.sparse.m$'Gene Expression'


# initialize the Seurat object with the row (non-normalized data)
nsclc.seurat.obj <- CreateSeuratObject(counts=cts, project="NSCLC", min.cells = 3, min.features=200)
str(nsclc.sparse.obj)
nsclc.sparse.obj
# 32978 features across 71880 samples within 1 assay 

# 1 QC- filter low quality cells 

View(nsclc.seurat.obj@meta.data)

# % MT reads 
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method='lm')

# 2. Filtering 
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# 3 normalize data 
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor =  10000)
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)

# 4 Identify highly variable features 
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# iswnridfy rhw 10 most highly variable genes
top50 <- head(VariableFeatures(nsclc.seurat.obj), 50)

# plot variable features with and without labels 
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top50, repel = TRUE)

# 5 scaling 
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
str(nsclc.seurat.obj)

# 6 perform linear dimensionality reduction 
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results 
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5 )
DimHeatmap(nsclc.seurat.obj, dims=1, cells=500, balanced = TRUE)

# determine dimensionality of the data 
ElbowPlot(nsclc.seurat.obj)


# 7 clustering 
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj,dims = 1:15)

# understanding resolution 
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.3", label = TRUE)

# setting identity of clusters 
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.3"
Idents(nsclc.seurat.obj)

# non-linear dimensionality redduction 
# if you haven't install UMAP, you can do so via reticulate py_install(packges = 
# 'umap-learn')
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set 'label = TRUE' or use the LabelClusters function to help label individual clusters
DimPlot(nsclc.seurat.obj, reduction = 'umap')




