# Remove gene from gene matrix prior to clustering

# Load library

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

# load each of the seurat objects of filtered and processed eye data
# we load in the seurat objects where R7 and R8 clusters are annotated as Rh3, Rh4 and Rh5, Rh6 clusters. 
# dorsal third R7s were clustered with Rh4 and Dorsal rim area were clustered as Rh3
# a1dr = adult 1 day eye data set
# a3dr = adult 3 day eye data set
# a7dr = adult 7 day eye data set

# We will extract the SoupX processed counts matrices from each of these data set and then remove Rh3, Rh4, Rh5 and/or Rh6

a1dr.data <- as.matrix(GetAssayData(a1dr, slot = "counts", assay = "RNA"))
a1drh34.data <- a1dr.data[-(which(rownames(a1dr.data) %in% c("Rh3","Rh4"))),]
a1drh56.data <- a1dr.data[-(which(rownames(a1dr.data) %in% c("Rh5","Rh6"))),]
a1dra.data <- a1dr.data[-(which(rownames(a1dr.data) %in% c("Rh3","Rh4","Rh5","Rh6"))),]
a1drninae.data <- a1dr.data[-(which(rownames(a1dr.data) %in% c("ninaE"))),]
a1drall.data <- a1dr.data[-(which(rownames(a1dr.data) %in% c("ninaE","Rh3","Rh4","Rh5","Rh6"))),]

a3dr.data <- as.matrix(GetAssayData(a3dr, slot = "counts", assay = "RNA"))
a3drh34.data <- a3dr.data[-(which(rownames(a3dr.data) %in% c("Rh3","Rh4"))),]
a3drh56.data <- a3dr.data[-(which(rownames(a3dr.data) %in% c("Rh5","Rh6"))),]
a3dra.data <- a3dr.data[-(which(rownames(a3dr.data) %in% c("Rh3","Rh4","Rh5","Rh6"))),]
a3drninae.data <- a3dr.data[-(which(rownames(a1dr.data) %in% c("ninaE"))),]
a3drall.data <- a3dr.data[-(which(rownames(a1dr.data) %in% c("ninaE","Rh3","Rh4","Rh5","Rh6"))),]

a7dr.data <- as.matrix(GetAssayData(a7dr, slot = "counts", assay = "RNA"))
a7drh34.data <- a7dr.data[-(which(rownames(a7dr.data) %in% c("Rh3","Rh4"))),]
a7drh56.data <- a7dr.data[-(which(rownames(a7dr.data) %in% c("Rh5","Rh6"))),]
a7dra.data <- a7dr.data[-(which(rownames(a7dr.data) %in% c("Rh3","Rh4","Rh5","Rh6"))),]
a7drninae.data <- a7dr.data[-(which(rownames(a1dr.data) %in% c("ninaE"))),]
a7drall.data <- a7dr.data[-(which(rownames(a1dr.data) %in% c("ninaE","Rh3","Rh4","Rh5","Rh6"))),]

# We run Seurat scTransform workflow on these count matrices

# create seurat object
a1dr34 <- CreateSeuratObject(counts = a1drh34.data, project = "a1d_remove_Rh34")
a1dr56 <- CreateSeuratObject(counts = a1drh56.data, project = "a1d_remove_Rh56")
a1dra <- CreateSeuratObject(counts = a1dra.data, project = "a1d_remove_Rh")
a1drninae <- CreateSeuratObject(counts = a1drninae.data, project = "a1d_remove_ninaE")
a1drall <- CreateSeuratObject(counts = a1drall.data, project = "a1d_remove_all_Rh")

a3dr34 <- CreateSeuratObject(counts = a3drh34.data, project = "a3d_remove_Rh34")
a3dr56 <- CreateSeuratObject(counts = a3drh56.data, project = "a3d_remove_Rh56")
a3dra <- CreateSeuratObject(counts = a3dra.data, project = "a3d_remove_Rh")
a3drninae <- CreateSeuratObject(counts = a3drninae.data, project = "a3d_remove_ninaE")
a3drall <- CreateSeuratObject(counts = a3drall.data, project = "a3d_remove_all_Rh")

a7dr34 <- CreateSeuratObject(counts = a7drh34.data, project = "a7d_remove_Rh34")
a7dr56 <- CreateSeuratObject(counts = a7drh56.data, project = "a7d_remove_Rh56")
a7dra <- CreateSeuratObject(counts = a7dra.data, project = "a7d_remove_Rh")
a7drninae <- CreateSeuratObject(counts = a7drninae.data, project = "a7d_remove_ninaE")
a7drall <- CreateSeuratObject(counts = a7drall.data, project = "a7d_remove_all_Rh")

# Add in percent mitochondrial genes metadata
a1dr34 <- PercentageFeatureSet(a1dr34, pattern = "mt:", col.name = "percent.mt")
a1dr56 <- PercentageFeatureSet(a1dr56, pattern = "mt:", col.name = "percent.mt")
a1dra <- PercentageFeatureSet(a1dra, pattern = "mt:", col.name = "percent.mt")
a1drninae <- PercentageFeatureSet(a1drninae, pattern = "mt:", col.name = "percent.mt")
a1drall <- PercentageFeatureSet(a1drall, pattern = "mt:", col.name = "percent.mt")

a3dr34 <- PercentageFeatureSet(a3dr34, pattern = "mt:", col.name = "percent.mt")
a3dr56 <- PercentageFeatureSet(a3dr56, pattern = "mt:", col.name = "percent.mt")
a3dra <- PercentageFeatureSet(a3dra, pattern = "mt:", col.name = "percent.mt")
a3drninae <- PercentageFeatureSet(a3drninae, pattern = "mt:", col.name = "percent.mt")
a3drall <- PercentageFeatureSet(a3drall, pattern = "mt:", col.name = "percent.mt")

a7dr34 <- PercentageFeatureSet(a7dr34, pattern = "mt:", col.name = "percent.mt")
a7dr56 <- PercentageFeatureSet(a7dr56, pattern = "mt:", col.name = "percent.mt")
a7dra <- PercentageFeatureSet(a7dra, pattern = "mt:", col.name = "percent.mt")
a7drninae <- PercentageFeatureSet(a7drninae, pattern = "mt:", col.name = "percent.mt")
a7drall <- PercentageFeatureSet(a7drall, pattern = "mt:", col.name = "percent.mt")

# SCTransform
a1dr34 <- SCTransform(a1dr34, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a1dr56 <- SCTransform(a1dr56, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a1dra <- SCTransform(a1dra, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a1drninae <- SCTransform(a1drninae, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a1drall <- SCTransform(a1drall, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)

a3dr34 <- SCTransform(a3dr34, variable.features.n = 2800, vars.to.regress = "percent.mt", verbose = TRUE)
a3dr56 <- SCTransform(a3dr56, variable.features.n = 2800, vars.to.regress = "percent.mt", verbose = TRUE)
a3dra <- SCTransform(a3dra, variable.features.n = 2800, vars.to.regress = "percent.mt", verbose = TRUE)
a3drninae <- SCTransform(a3drninae, variable.features.n = 2800, vars.to.regress = "percent.mt", verbose = TRUE)
a3drall <- SCTransform(a3drall, variable.features.n = 2800, vars.to.regress = "percent.mt", verbose = TRUE)

a7dr34 <- SCTransform(a7dr34, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)
a7dr56 <- SCTransform(a7dr56, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)
a7dra <- SCTransform(a7dra, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)
a7drninae <- SCTransform(a7drninae, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)
a7drall <- SCTransform(a7drall, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)

# clustering with 50 dimensions
a1dr34 <- RunPCA(a1dr34, verbose = TRUE)
a1dr34 <- RunUMAP(a1dr34, dims = 1:50, verbose = TRUE)
a1dr34 <- FindNeighbors(a1dr34, dims = 1:50, verbose = TRUE)
a1dr34 <- FindClusters(a1dr34, verbose = TRUE)

a1dr56 <- RunPCA(a1dr56, verbose = TRUE)
a1dr56 <- RunUMAP(a1dr56, dims = 1:50, verbose = TRUE)
a1dr56 <- FindNeighbors(a1dr56, dims = 1:50, verbose = TRUE)
a1dr56 <- FindClusters(a1dr56, verbose = TRUE)

a1dra <- RunPCA(a1dra, verbose = TRUE)
a1dra <- RunUMAP(a1dra, dims = 1:50, verbose = TRUE)
a1dra <- FindNeighbors(a1dra, dims = 1:50, verbose = TRUE)
a1dra <- FindClusters(a1dra, verbose = TRUE)

a1drninae <- RunPCA(a1drninae, verbose = TRUE)
a1drninae <- RunUMAP(a1drninae, dims = 1:50, verbose = TRUE)
a1drninae <- FindNeighbors(a1drninae, dims = 1:50, verbose = TRUE)
a1drninae <- FindClusters(a1drninae, verbose = TRUE)

a1drall <- RunPCA(a1drall, verbose = TRUE)
a1drall <- RunUMAP(a1drall, dims = 1:50, verbose = TRUE)
a1drall <- FindNeighbors(a1drall, dims = 1:50, verbose = TRUE)
a1drall <- FindClusters(a1drall, verbose = TRUE)

#
a3dr34 <- RunPCA(a3dr34, verbose = TRUE)
a3dr34 <- RunUMAP(a3dr34, dims = 1:50, verbose = TRUE)
a3dr34 <- FindNeighbors(a3dr34, dims = 1:50, verbose = TRUE)
a3dr34 <- FindClusters(a3dr34, verbose = TRUE)

a3dr56 <- RunPCA(a3dr56, verbose = TRUE)
a3dr56 <- RunUMAP(a3dr56, dims = 1:50, verbose = TRUE)
a3dr56 <- FindNeighbors(a3dr56, dims = 1:50, verbose = TRUE)
a3dr56 <- FindClusters(a3dr56, verbose = TRUE)

a3dra <- RunPCA(a3dra, verbose = TRUE)
a3dra <- RunUMAP(a3dra, dims = 1:50, verbose = TRUE)
a3dra <- FindNeighbors(a3dra, dims = 1:50, verbose = TRUE)
a3dra <- FindClusters(a3dra, verbose = TRUE)

a3drninae <- RunPCA(a3drninae, verbose = TRUE)
a3drninae <- RunUMAP(a3drninae, dims = 1:50, verbose = TRUE)
a3drninae <- FindNeighbors(a3drninae, dims = 1:50, verbose = TRUE)
a3drninae <- FindClusters(a3drninae, verbose = TRUE)

a3drall <- RunPCA(a3drall, verbose = TRUE)
a3drall <- RunUMAP(a3drall, dims = 1:50, verbose = TRUE)
a3drall <- FindNeighbors(a3drall, dims = 1:50, verbose = TRUE)
a3drall <- FindClusters(a3drall, verbose = TRUE)

#
a7dr34 <- RunPCA(a7dr34, verbose = TRUE)
a7dr34 <- RunUMAP(a7dr34, dims = 1:50, verbose = TRUE)
a7dr34 <- FindNeighbors(a7dr34, dims = 1:50, verbose = TRUE)
a7dr34 <- FindClusters(a7dr34, verbose = TRUE)

a7dr56 <- RunPCA(a7dr56, verbose = TRUE)
a7dr56 <- RunUMAP(a7dr56, dims = 1:50, verbose = TRUE)
a7dr56 <- FindNeighbors(a7dr56, dims = 1:50, verbose = TRUE)
a7dr56 <- FindClusters(a7dr56, verbose = TRUE)

a7dra <- RunPCA(a7dra, verbose = TRUE)
a7dra <- RunUMAP(a7dra, dims = 1:50, verbose = TRUE)
a7dra <- FindNeighbors(a7dra, dims = 1:50, verbose = TRUE)
a7dra <- FindClusters(a7dra, verbose = TRUE)

a7drninae <- RunPCA(a7drninae, verbose = TRUE)
a7drninae <- RunUMAP(a7drninae, dims = 1:50, verbose = TRUE)
a7drninae <- FindNeighbors(a7drninae, dims = 1:50, verbose = TRUE)
a7drninae <- FindClusters(a7drninae, verbose = TRUE)

a7drall <- RunPCA(a7drall, verbose = TRUE)
a7drall <- RunUMAP(a7drall, dims = 1:50, verbose = TRUE)
a7drall <- FindNeighbors(a7drall, dims = 1:50, verbose = TRUE)
a7drall <- FindClusters(a7drall, verbose = TRUE)

# we can plot the cluster plots at this point but the cluster labels will be seurat auto generated numbers
DimPlot(a7dr34, label = TRUE)

# transfer labels from unaltered data sets (a1dr, a3dr, a7dr)
a1d_ident <- Idents(a1dr)
a1dr34[["old.ident"]] <- a1d_ident
Idents(a1dr34)<- "old.ident"
a1dr56[["old.ident"]] <- a1d_ident
Idents(a1dr56)<- "old.ident"
a1dra[["old.ident"]] <- a1d_ident
Idents(a1dra)<- "old.ident"
a1drninae[["old.ident"]] <- a1d_ident
Idents(a1drninae)<- "old.ident"
a1drall[["old.ident"]] <- a1d_ident
Idents(a1drall)<- "old.ident"

a3d_ident <- Idents(a3dr)
a3dr34[["old.ident"]] <- a3d_ident
Idents(a3dr34)<- "old.ident"
a3dr56[["old.ident"]] <- a3d_ident
Idents(a3dr56)<- "old.ident"
a3dra[["old.ident"]] <- a3d_ident
Idents(a3dra)<- "old.ident"
a3drninae[["old.ident"]] <- a3d_ident
Idents(a3drninae)<- "old.ident"
a3drall[["old.ident"]] <- a3d_ident
Idents(a3drall)<- "old.ident"

a7d_ident <- Idents(a7dr)
a7dr34[["old.ident"]] <- a7d_ident
Idents(a7dr34)<- "old.ident"
a7dr56[["old.ident"]] <- a7d_ident
Idents(a7dr56)<- "old.ident"
a7dra[["old.ident"]] <- a7d_ident
Idents(a7dra)<- "old.ident"
a7drninae[["old.ident"]] <- a7d_ident
Idents(a7drninae)<- "old.ident"
a7drall[["old.ident"]] <- a7d_ident
Idents(a7drall)<- "old.ident"

# We can generate clusterplots and featureplots at this point
# Remember that Rh3/Rh4/Rh5/Rh6 were removed from the objects and we could not generate FeaturePlots for the removed genes

#saveRDS
saveRDS(a1dr34, file = "a1d_Rh34_remove.rds")
saveRDS(a1dr56, file = "a1d_Rh56_remove.rds")
saveRDS(a1dra, file = "a1d_RhAll_remove.rds")
saveRDS(a1drninae, file = "a1d_ninaE_remove.rds")
saveRDS(a1drall, file = "a1d_ninaE_RhAll_remove.rds")

saveRDS(a3dr34, file = "a3d_Rh34_remove.rds")
saveRDS(a3dr56, file = "a3d_Rh56_remove.rds")
saveRDS(a3dra, file = "a3d_RhAll_remove.rds")
saveRDS(a3drninae, file = "a3d_ninaE_remove.rds")
saveRDS(a3drall, file = "a3d_ninaE_RhAll_remove.rds")

saveRDS(a7dr34, file = "a7d_Rh34_remove.rds")
saveRDS(a7dr56, file = "a7d_Rh56_remove.rds")
saveRDS(a7dra, file = "a7d_RhAll_remove.rds")
saveRDS(a7drninae, file = "a7d_ninaE_remove.rds")
saveRDS(a7drall, file = "a7d_ninaE_RhAll_remove.rds")