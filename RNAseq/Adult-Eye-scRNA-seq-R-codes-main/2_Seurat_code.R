#Seurat SCTransform for all data sets

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

# We had already ran SoupX to correct for ambient RNA in each data set
# SoupX output was saved as: A1D_SoupX_out, A3D_SoupX_out, A7D_SoupX_out, A1F_SoupX_out
# for adult male 1 day, 3 days, 7 days old adult eyes and 1 day old adult female eyes, respectively.

# Create seurat objects for each data set with cut offs, and run Seurat SCTransform


### adult male 1 day
a1dall <- CreateSeuratObject(counts = A1D_SoupX_out, project = "A1D_male", min.cells = 3, min.features = 200)
a1dall <- PercentageFeatureSet(a1dall, pattern = "mt:", col.name = "percent.mt")

# Generate some QC plots
plot1 <- FeatureScatter(a1dall, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(a1dall, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Chose the following nFeature and percent mito cutoffs after viewing QC. 
# Then run the rest of Seurat
a1dall <- subset(a1dall, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 20)
a1dall <- SCTransform(a1dall, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a1dall <- RunPCA(a1dall, verbose = TRUE)
a1dall <- RunUMAP(a1dall, dims = 1:50, verbose = TRUE)
a1dall <- FindNeighbors(a1dall, dims = 1:50, verbose = TRUE)
a1dall <- FindClusters(a1dall, verbose = TRUE)
DimPlot(a1dall, label = TRUE)
# To remove brain / non-eye tissue, we used feature plots of fne, moody and repo
# Note that moody is also seen in pigment cell cluster (Pdh, w positive). the pigment cell cluster is left alone
FeaturePlot(a1dall, features = c("fne","moody","repo"), order = T)

# Subset such that we only get eye cells. Get raw counts data and recluster
# fne, moody and/or repo positive clusters are removed
# These are the cluster numbers assigned by Seurat during the time when the analysis was ran
# When Seurat was ran again the cluster numbers may change, but the criteria of selecting against fne, moody and/or repo still holds true
a1dsubset <- subset(a1dall, idents = c(31,19,22,18,20,16,26,15,10,24,27,30,28,23,19), invert = TRUE)

a1dsubset.data <- as.matrix(GetAssayData(a1dsubset, slot = "counts", assay = "RNA"))

a1d <- CreateSeuratObject(counts = a1dsubset.data, project = "A1D_male")
a1d <- PercentageFeatureSet(a1d, pattern = "mt:", col.name = "percent.mt")
a1d <- SCTransform(a1d, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a1d <- RunPCA(a1d, verbose = TRUE)
a1d <- RunUMAP(a1d, dims = 1:50, verbose = TRUE)
a1d <- FindNeighbors(a1d, dims = 1:50, verbose = TRUE)
a1d <- FindClusters(a1d, verbose = TRUE)
DimPlot(a1d, label = TRUE)

# The first subset (a1dsubset) may not be perfect as some fne/moody/repo positive cells were not removed
# These cells form their own clusters so we can subset them out
# after subsetting again, RunPCA, RunUMAP, FindNeighbors and FindClusters were repeated with the same conditions

a1d2 <- subset(a1d, idents = c(19), invert = TRUE)
a1d2 <- RunPCA(a1d2, verbose = TRUE)
a1d2 <- RunUMAP(a1d2, dims = 1:50, verbose = TRUE)
a1d2 <- FindNeighbors(a1d2, dims = 1:50, verbose = TRUE)
a1d2 <- FindClusters(a1d2, verbose = TRUE)
DimPlot(a1d2, label = TRUE)

# we then annotate the clusters with known cell type markers:
# ninaE for R1-6; Rh3, Rh4 for R7; Rh5, Rh6 for R8; Pdh, w for pigment cells;ct, Crys for cone cells
# Further in vivo analyses show that wrapper is a marker for primary pigment cells and santa-maria is a marker for secondary and tertiary pigment cells
# Therefore we renamed the initial pigment cell clusters into primary pigment cells and secondary and tertiary pigment cells
# Assign cluster identity:
# a1d_ct has clusters annotated as cell types
# a1d_rh has clusters annotated with simplified Rh types (dorsal third R7 are grouped with Rh4, dorsal rim are grouped with Rh3)
# a1d_rhc has clusters annotated with Rh types, dorsal rim and dorsal third R7
a1d_ct <- a1d2
a1d_rh <- a1d2
a1d_rhc <- a1d2
r7.name <- WhichCells(object = a1d_ct, idents = c(16,10,5,8,12,20))
r8.name <- WhichCells(object = a1d_ct, idents = c(18,6))
r16.name <- WhichCells(object = a1d_ct, idents = c(2,3,5,7,13,14))
cone.name <- WhichCells(object = a1d_ct,idents = c(21))
pigm.name <- WhichCells(object = a1d_ct, idents = c(0,1,4,11,19,9))
pigm1.name <- WhichCells(object = a1d_ct, idents = c(17))
a1d_ct <- SetIdent(object = a1d_ct, cells = r16.name, value = "R1-6")
a1d_ct <- SetIdent(object = a1d_ct, cells = r7.name, value = "R7")
a1d_ct <- SetIdent(object = a1d_ct, cells = r8.name, value = "R8")
a1d_ct <- SetIdent(object = a1d_ct, cells = cone.name, value = "Cone")
a1d_ct <- SetIdent(object = a1d_ct, cells = pigm1.name, value = "1_Pigment")
a1d_ct <- SetIdent(object = a1d_ct, cells = pigm.name, value = "2,3_Pigments")

rh3.name <- WhichCells(object = a1d_rh, idents = c(8,12,20))
rh4.name <- WhichCells(object = a1d_rh, idents = c(16,10,15))
rh5.name <- WhichCells(object = a1d_rh, idents = c(18))
rh6.name <- WhichCells(object = a1d_rh, idents = c(6))
r16.name <- WhichCells(object = a1d_rh, idents = c(2,3,5,7,13,14))
cone.name <- WhichCells(object = a1d_rh,idents = c(21))
pigm.name <- WhichCells(object = a1d_rh, idents = c(0,1,4,11,19,9))
pigm1.name <- WhichCells(object = a1d_rh, idents = c(17))
a1d_rh <- SetIdent(object = a1d_rh, cells = r16.name, value = "R1-6")
a1d_rh <- SetIdent(object = a1d_rh, cells = rh3.name, value = "Rh3")
a1d_rh <- SetIdent(object = a1d_rh, cells = rh4.name, value = "Rh4")
a1d_rh <- SetIdent(object = a1d_rh, cells = rh5.name, value = "Rh5")
a1d_rh <- SetIdent(object = a1d_rh, cells = rh6.name, value = "Rh6")
a1d_rh <- SetIdent(object = a1d_rh, cells = cone.name, value = "Cone")
a1d_rh <- SetIdent(object = a1d_rh, cells = pigm1.name, value = "1_Pigment")
a1d_rh <- SetIdent(object = a1d_rh, cells = pigm.name, value = "2,3_Pigments")

D3.name <- WhichCells(object = a1d_rhc, idents = c(16))
DRA.name <- WhichCells(object = a1d_rhc, idents = c(20))
rh3.name <- WhichCells(object = a1d_rhc, idents = c(8,12))
rh4.name <- WhichCells(object = a1d_rhc, idents = c(10,15))
rh5.name <- WhichCells(object = a1d_rhc, idents = c(18))
rh6.name <- WhichCells(object = a1d_rhc, idents = c(6))
r16.name <- WhichCells(object = a1d_rhc, idents = c(2,3,5,7,13,14))
cone.name <- WhichCells(object = a1d_rhc,idents = c(21))
pigm.name <- WhichCells(object = a1d_rhc, idents = c(0,1,4,11,19,9))
pigm1.name <- WhichCells(object = a1d_rhc, idents = c(17))
a1d_rhc <- SetIdent(object = a1d_rhc, cells = D3.name, value = "Dorsal_3rd_R7")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = DRA.name, value = "DRA")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = r16.name, value = "R1-6")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = rh3.name, value = "Rh3")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = rh4.name, value = "Rh4")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = rh5.name, value = "Rh5")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = rh6.name, value = "Rh6")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = cone.name, value = "Cone")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = pigm1.name, value = "1_Pigment")
a1d_rhc <- SetIdent(object = a1d_rhc, cells = pigm.name, value = "2,3_Pigments")

saveRDS(a1d_ct , file = "Adult_1_day_celltype.rds")
saveRDS(a1d_rh, file = "Adult_1_day_Rh.rds")
saveRDS(a1d_rhc, file = "Adult_1_day_Rh_complete.rds")

### adult male 3 days
# The procedure is the same as 1 day males
a3dall <- CreateSeuratObject(counts = A3D_SoupX_out, project = "A1D_male", min.cells = 3, min.features = 200)
a3dall <- PercentageFeatureSet(a3dall, pattern = "mt:", col.name = "percent.mt")

#QC plots
plot1 <- FeatureScatter(a3dall, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(a3dall, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Chose the following nFeature and percent mito cutoffs after viewing QC. 
# Then run the rest of Seurat
a3dall <- subset(a3dall, subset = nFeature_RNA > 200 & nFeature_RNA < 2800 & percent.mt < 25)
a3dall <- SCTransform(a3dall, variable.features.n = 2800, vars.to.regress = "percent.mt", verbose = TRUE)
a3dall <- RunPCA(a3dall, verbose = TRUE)
a3dall <- RunUMAP(a3dall, dims = 1:50, verbose = TRUE)
a3dall <- FindNeighbors(a3dall, dims = 1:50, verbose = TRUE)
a3dall <- FindClusters(a3dall, verbose = TRUE)
DimPlot(a3dall, label = TRUE)
# To remove brain / non-eye tissue, we used feature plots of fne, moody and repo
# Note that moody is also seen in pigment cell cluster (Pdh, w positive). the pigment cell cluster is left alone
FeaturePlot(a3dall, features = c("fne","moody","repo"), order = T)

# Subset such that we only get eye cells. Get raw counts data and recluster
# fne, moody and/or repo positive clusters are removed
# These are the cluster numbers assigned by Seurat during the time when the analysis was ran
# When Seurat was ran again the cluster numbers may change, but the criteria of selecting against fne, moody and/or repo still holds true
a3dsubset <- subset(a3dall, idents = c(17,18,21,23,15,24,19,22,16), invert = TRUE)

a3dsubset.data <- as.matrix(GetAssayData(a3dsubset, slot = "counts", assay = "RNA"))

# Rerun sctransform on the subset raw count data
# FindCluster resolution was increased to 5 to more accurately assign cluster identity
a3d <- CreateSeuratObject(counts = a3dsubset.data, project = "A3D_male")
a3d <- PercentageFeatureSet(a3d, pattern = "mt:", col.name = "percent.mt")
a3d <- SCTransform(a3d, variable.features.n = 2800, vars.to.regress = "percent.mt", verbose = TRUE)
a3d <- RunPCA(a3d, verbose = TRUE)
a3d <- RunUMAP(a3d, dims = 1:50, verbose = TRUE)
a3d <- FindNeighbors(a3d, dims = 1:50, verbose = TRUE)
a3d <- FindClusters(a3d, verbose = TRUE, resolution = 5)
DimPlot(a3d, label = TRUE)
# The subset step above removed all the non eye tissue. No need to subset again. 
# Cluster identities were assigned similarly to 1day data
a3d_ct <- a3d
a3d_rh <- a3d

r7.name <- WhichCells(object = a3d_ct, idents = c(20,32,45,52,14,22,49,23,42,40,9,46,50,51))
r8.name <- WhichCells(object = a3d_ct, idents = c(43,27,38,35,25,29))
r16.name <- WhichCells(object = a3d_ct, idents = c(47,18,15,34,37,26,17,28,39,19,16,30,12,3,53,24,41,13,21,5,10,7,4,31,33))
cone.name <- WhichCells(object = a3d_ct,idents = c(44))
pigm.name <- WhichCells(object = a3d_ct, idents = c(36,0,1,8,2,11,6))
pigm1.name <- WhichCells(object = a3d_ct, idents = c(48))

a3d_ct <- SetIdent(object = a3d_ct, cells = r16.name, value = "R1-6")
a3d_ct <- SetIdent(object = a3d_ct, cells = r7.name, value = "R7")
a3d_ct <- SetIdent(object = a3d_ct, cells = r8.name, value = "R8")
a3d_ct <- SetIdent(object = a3d_ct, cells = cone.name, value = "Cone")
a3d_ct <- SetIdent(object = a3d_ct, cells = pigm1.name, value = "1_Pigment")
a3d_ct <- SetIdent(object = a3d_ct, cells = pigm.name, value = "2,3_Pigments")

rh3.name <- WhichCells(object = a3d_rh, idents = c(20,32,45,52,14,22,49))
rh4.name <- WhichCells(object = a3d_rh, idents = c(23,42,40,9,46,50,51))
rh5.name <- WhichCells(object = a3d_rh, idents = c(43))
rh6.name <- WhichCells(object = a3d_rh, idents = c(27,38,35,25,29))
r16.name <- WhichCells(object = a3d_rh, idents = c(47,18,15,34,37,26,17,28,39,19,16,30,12,3,53,24,41,13,21,5,10,7,4,31,33))
cone.name <- WhichCells(object = a3d_rh,idents = c(44))
pigm.name <- WhichCells(object = a3d_rh, idents = c(36,0,1,8,2,11,6))
pigm1.name <- WhichCells(object = a3d_rh, idents = c(48))

a3d_rh <- SetIdent(object = a3d_rh, cells = r16.name, value = "R1-6")
a3d_rh <- SetIdent(object = a3d_rh, cells = rh3.name, value = "Rh3")
a3d_rh <- SetIdent(object = a3d_rh, cells = rh4.name, value = "Rh4")
a3d_rh <- SetIdent(object = a3d_rh, cells = rh5.name, value = "Rh5")
a3d_rh <- SetIdent(object = a3d_rh, cells = rh6.name, value = "Rh6")
a3d_rh <- SetIdent(object = a3d_rh, cells = cone.name, value = "Cone")
a3d_rh <- SetIdent(object = a3d_rh, cells = pigm1.name, value = "1_Pigment")
a3d_rh <- SetIdent(object = a3d_rh, cells = pigm.name, value = "2,3_Pigments")

D3.name <- WhichCells(object = a3d_rhc, idents = c(50,51))
DRA.name <- WhichCells(object = a3d_rhc, idents = c(52))
rh3.name <- WhichCells(object = a3d_rhc, idents = c(20,32,45,14,22,49))
rh4.name <- WhichCells(object = a3d_rhc, idents = c(23,42,40,9,46))
rh5.name <- WhichCells(object = a3d_rhc, idents = c(43))
rh6.name <- WhichCells(object = a3d_rhc, idents = c(27,38,35,25,29))
r16.name <- WhichCells(object = a3d_rhc, idents = c(47,18,15,34,37,26,17,28,39,19,16,30,12,3,53,24,41,13,21,5,10,7,4,31,33))
cone.name <- WhichCells(object = a3d_rhc,idents = c(44))
pigm.name <- WhichCells(object = a3d_rhc, idents = c(36,0,1,8,2,11,6))
pigm1.name <- WhichCells(object = a3d_rhc, idents = c(48))

a3d_rhc <- SetIdent(object = a3d_rhc, cells = D3.name, value = "Dorsal_3rd_R7")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = DRA.name, value = "DRA")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = r16.name, value = "R1-6")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = rh3.name, value = "Rh3")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = rh4.name, value = "Rh4")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = rh5.name, value = "Rh5")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = rh6.name, value = "Rh6")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = cone.name, value = "Cone")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = pigm1.name, value = "1_Pigment")
a3d_rhc <- SetIdent(object = a3d_rhc, cells = pigm.name, value = "2,3_Pigments")

saveRDS(a3d_ct, file = "Adult_3day_celltype.rds")
saveRDS(a3d_rh, file = "Adult_3day_rh.rds")
saveRDS(a3d_rhc, file = "Adult_3day_rh_complete.rds")

### Adult 7 days
# same procedure as 1day and 3days eyes
a7dall <- CreateSeuratObject(counts = A7D_SoupX_out, project = "A7D_male", min.cells = 3, min.features = 200)
a7dall <- PercentageFeatureSet(a7dall, pattern = "mt:", col.name = "percent.mt")

plot1 <- FeatureScatter(a7dall, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(a7dall, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Chose the following nFeature and percent mito cutoffs after viewing QC. 
# Then run the rest of Seurat
a7dall <- subset(a7dall, subset = nFeature_RNA > 200 & nFeature_RNA < 2700 & percent.mt < 18)
a7dall <- SCTransform(a7dall, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)
a7dall <- RunPCA(a7dall, verbose = TRUE)
a7dall <- RunUMAP(a7dall, dims = 1:50, verbose = TRUE)
a7dall <- FindNeighbors(a7dall, dims = 1:50, verbose = TRUE)
a7dall <- FindClusters(a7dall, verbose = TRUE)
DimPlot(a7dall, label = TRUE)
# To remove brain / non-eye tissue, we used feature plots of fne, moody and repo
# Note that moody is also seen in pigment cell cluster (Pdh, w positive). the pigment cell cluster is left alone
FeaturePlot(a7dall, features = c("fne","moody","repo"), order = T)

# Subset such that we only get eye cells. Get raw counts data and recluster
# fne, moody and/or repo positive clusters are removed
# These are the cluster numbers assigned by Seurat during the time when the analysis was ran
# When Seurat was ran again the cluster numbers may change, but the criteria of selecting against fne, moody and/or repo still holds true
a7dsubset <- subset(a7dall, idents = c(21,18,22,11,13,15,14,20,12), invert = TRUE)

a7dsubset.data <- as.matrix(GetAssayData(a7dsubset, slot = "counts", assay = "RNA"))

# Rerun sctransform on the subset raw count data
a7d <- CreateSeuratObject(counts = a7dsubset.data, project = "A7D_male")
a7d <- PercentageFeatureSet(a7d, pattern = "mt:", col.name = "percent.mt")
a7d <- SCTransform(a7d, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)
a7d <- RunPCA(a7d, verbose = TRUE)
a7d <- RunUMAP(a7d, dims = 1:50, verbose = TRUE)
a7d <- FindNeighbors(a7d, dims = 1:50, verbose = TRUE)
a7d <- FindClusters(a7d, verbose = TRUE)
DimPlot(a7d, label = TRUE)
# The subset step above removed all the non eye tissue. No need to subset again. 
# Cluster identities were assigned similarly to 1day data
a7d_ct <- a7d
a7d_rh <- a7d
a7d_rhc <- a7d

r7.name <- WhichCells(object = a7d_ct, idents = c(11,12,10,8,15))
r8.name <- WhichCells(object = a7d_ct, idents = c(13,3))
r16.name <- WhichCells(object = a7d_ct, idents = c(0,1,6,4,5))
cone.name <- WhichCells(object = a7d_ct,idents = c(16))
pigm.name <- WhichCells(object = a7d_ct, idents = c(2,9,7))
pigm1.name <- WhichCells(object = a7d_ct, idents = c(14))
a7d_ct <- SetIdent(object = a7d_ct, cells = r16.name, value = "R1-6")
a7d_ct <- SetIdent(object = a7d_ct, cells = r7.name, value = "R7")
a7d_ct <- SetIdent(object = a7d_ct, cells = r8.name, value = "R8")
a7d_ct <- SetIdent(object = a7d_ct, cells = cone.name, value = "Cone")
a7d_ct <- SetIdent(object = a7d_ct, cells = pigm1.name, value = "1_Pigment")
a7d_ct <- SetIdent(object = a7d_ct, cells = pigm.name, value = "2,3_Pigments")

rh3.name <- WhichCells(object = a7d_rh, idents = c(8,10,15))
rh4.name <- WhichCells(object = a7d_rh, idents = c(11,12))
rh5.name <- WhichCells(object = a7d_rh, idents = c(13))
rh6.name <- WhichCells(object = a7d_rh, idents = c(3))
r16.name <- WhichCells(object = a7d_rh, idents = c(0,1,6,4,5))
cone.name <- WhichCells(object = a7d_rh,idents = c(16))
pigm.name <- WhichCells(object = a7d_rh, idents = c(2,9,7))
pigm1.name <- WhichCells(object = a7d_rh, idents = c(14))
a7d_rh <- SetIdent(object = a7d_rh, cells = r16.name, value = "R1-6")
a7d_rh <- SetIdent(object = a7d_rh, cells = rh3.name, value = "Rh3")
a7d_rh <- SetIdent(object = a7d_rh, cells = rh4.name, value = "Rh4")
a7d_rh <- SetIdent(object = a7d_rh, cells = rh5.name, value = "Rh5")
a7d_rh <- SetIdent(object = a7d_rh, cells = rh6.name, value = "Rh6")
a7d_rh <- SetIdent(object = a7d_rh, cells = cone.name, value = "Cone")
a7d_rh <- SetIdent(object = a7d_rh, cells = pigm1.name, value = "1_Pigment")
a7d_rh <- SetIdent(object = a7d_rh, cells = pigm.name, value = "2,3_Pigments")

a7d_rhc <- FindClusters(a7d_rhc, verbose = TRUE, resolution = 2) #need to increase resolution to isolate dorsal 3rd R7s
D3.name <- WhichCells(object = a7d_rhc, idents = c(25))
DRA.name <- WhichCells(object = a7d_rhc, idents = c(24))
rh3.name <- WhichCells(object = a7d_rhc, idents = c(9,15,18))
rh4.name <- WhichCells(object = a7d_rhc, idents = c(10,12))
rh5.name <- WhichCells(object = a7d_rhc, idents = c(22))
rh6.name <- WhichCells(object = a7d_rhc, idents = c(11,14,19))
r16.name <- WhichCells(object = a7d_rhc, idents = c(1,2,4,7,13,16,20,3,5,17))
cone.name <- WhichCells(object = a7d_rhc,idents = c(27))
pigm.name <- WhichCells(object = a7d_rhc, idents = c(0,6,26,8,21))
pigm1.name <- WhichCells(object = a7d_rhc, idents = c(23))
a7d_rhc <- SetIdent(object = a7d_rhc, cells = D3.name, value = "Dorsal_3rd_R7")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = DRA.name, value = "DRA")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = r16.name, value = "R1-6")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = rh3.name, value = "Rh3")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = rh4.name, value = "Rh4")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = rh5.name, value = "Rh5")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = rh6.name, value = "Rh6")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = cone.name, value = "Cone")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = pigm1.name, value = "1_Pigment")
a7d_rhc <- SetIdent(object = a7d_rhc, cells = pigm.name, value = "2,3_Pigments")
saveRDS(a7d_ct, file = "Adult_7day_celltype.rds")
saveRDS(a7d_rh, file = "Adult_7day_rh.rds")
saveRDS(a7d_rhc, file = "Adult_7day_rh_complete.rds")


### Adult 1 day Female
# same procedure as above
a1fall <- CreateSeuratObject(counts = A1D_F_SoupX_out, project = "A1D_female", min.cells = 3, min.features = 200)
a1fall <- PercentageFeatureSet(a1fall, pattern = "mt:", col.name = "percent.mt")

plot1 <- FeatureScatter(a1fall, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(a1fall, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Chose the following nFeature and percent mito cutoffs after viewing QC. 
# Then run the rest of Seurat
a1fall <- subset(a1fall, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)
a1fall <- SCTransform(a1fall, variable.features.n = 4000, vars.to.regress = "percent.mt", verbose = TRUE)
a1fall <- RunPCA(a1fall, verbose = TRUE)
a1fall <- RunUMAP(a1fall, dims = 1:50, verbose = TRUE)
a1fall <- FindNeighbors(a1fall, dims = 1:50, verbose = TRUE)
a1fall <- FindClusters(a1fall, verbose = TRUE)
DimPlot(a1fall, label = TRUE)
# To remove brain / non-eye tissue, we used feature plots of fne, moody and repo
# Note that moody is also seen in pigment cell cluster (Pdh, w positive). the pigment cell cluster is left alone
FeaturePlot(a1fall, features = c("fne","moody","repo"), order = T)

# Subset such that we only get eye cells. Get raw counts data and recluster
# fne, moody and/or repo positive clusters are removed
# These are the cluster numbers assigned by Seurat during the time when the analysis was ran
# When Seurat was ran again the cluster numbers may change, but the criteria of selecting against fne, moody and/or repo still holds true
a1fsubset <- subset(a1fall, idents = c(22,24,21,15,26,20,19,23,13,18), invert = TRUE)

a1fsubset.data <- as.matrix(GetAssayData(a1fsubset, slot = "counts", assay = "RNA"))

# Rerun sctransform on the subset raw count data
# cluster resolution of 8 was used
a1f <- CreateSeuratObject(counts = a1fsubset.data, project = "A1D_female")
a1f <- PercentageFeatureSet(a1f, pattern = "mt:", col.name = "percent.mt")
a1f <- SCTransform(a1f, variable.features.n = 2700, vars.to.regress = "percent.mt", verbose = TRUE)
a1f <- RunPCA(a1f, verbose = TRUE)
a1f <- RunUMAP(a1f, dims = 1:50, verbose = TRUE)
a1f <- FindNeighbors(a1f, dims = 1:50, verbose = TRUE)
a1f <- FindClusters(a1f, verbose = TRUE, resolution = 8)
DimPlot(a1f, label = TRUE)
# The subset step above removed all the non eye tissue. 
# There are 1 small cluster of Rh3 positive cells (29) and 1 small cluster of Rh6 positive cells (54) near/within R1-6
# These cells do not express pros nor CG2082 but express ninaE strongly
# We suspect these are droplets with multiple cells or cells with ambient RNA contamination that were not removed by SoupX.
# therefore we remove them and reran PCA clustering

a1f2 <- subset(a1f, idents = c(29,54), invert = TRUE)
a1f2 <- RunPCA(a1f2, verbose = TRUE)
a1f2 <- RunUMAP(a1f2, dims = 1:50, verbose = TRUE)
a1f2 <- FindNeighbors(a1f2, dims = 1:50, verbose = TRUE)
a1f2 <- FindClusters(a1f2, verbose = TRUE, resolution = 5)
DimPlot(a1f2, label = TRUE)

# Cluster identities were assigned similarly to 1day data
# cone cells and primary pigment cells are absent in this data set and were not included in cluster annotation
a1f_ct <- a1f2
a1f_rh <- a1f2

r7.name <- WhichCells(object = a1f_ct, idents = c(43,27,34,16,36,13,6,9,44))
r8.name <- WhichCells(object = a1f_ct, idents = c(21,39,26,37))
r16.name <- WhichCells(object = a1f_ct, idents = c(28,24,46,19,2,5,42,23,11,35,1,7,17,41,32,40,15,4,8,14,20,18,3,29))
pigm.name <- WhichCells(object = a1f_ct, idents = c(38,22,25,12,10,33,45,31,0,30))
a1f_ct <- SetIdent(object = a1f_ct, cells = r16.name, value = "R1-6")
a1f_ct <- SetIdent(object = a1f_ct, cells = r7.name, value = "R7")
a1f_ct <- SetIdent(object = a1f_ct, cells = r8.name, value = "R8")
a1f_ct <- SetIdent(object = a1f_ct, cells = pigm.name, value = "2,3_Pigments")

rh3.name <- WhichCells(object = a1f_rh, idents = c(36,13,6,9,44))
rh4.name <- WhichCells(object = a1f_rh, idents = c(43,27,34,16))
rh5.name <- WhichCells(object = a1f_rh, idents = c(37))
rh6.name <- WhichCells(object = a1f_rh, idents = c(21,39,26))
r16.name <- WhichCells(object = a1f_rh, idents = c(28,24,46,19,2,5,42,23,11,35,1,7,17,41,32,40,15,4,8,14,20,18,3,29))
pigm.name <- WhichCells(object = a1f_rh, idents = c(38,22,25,12,10,33,45,31,0,30))
a1f_rh <- SetIdent(object = a1f_rh, cells = r16.name, value = "R1-6")
a1f_rh <- SetIdent(object = a1f_rh, cells = rh3.name, value = "Rh3")
a1f_rh <- SetIdent(object = a1f_rh, cells = rh4.name, value = "Rh4")
a1f_rh <- SetIdent(object = a1f_rh, cells = rh5.name, value = "Rh5")
a1f_rh <- SetIdent(object = a1f_rh, cells = rh6.name, value = "Rh6")
a1f_rh <- SetIdent(object = a1f_rh, cells = pigm.name, value = "2,3_Pigments")

saveRDS(a1f_ct, file = "Adult_1day_Female_celltype.rds")
saveRDS(a1f_rh, file = "Adult_1day_Female_rh.rds")

# Markers for each cell type were called using FindAllMarkers function

a1d.markers <- FindAllMarkers(a1d_ct, only.pos = TRUE, min.pct = 0.25, logfc.threshold=0.25)
a3d.markers <- FindAllMarkers(a3d_ct, only.pos = TRUE, min.pct = 0.25, logfc.threshold=0.25)
a7d.markers <- FindAllMarkers(a7d_ct, only.pos = TRUE, min.pct = 0.25, logfc.threshold=0.25)
a1f.markers <- FindAllMarkers(a1f_ct, only.pos = TRUE, min.pct = 0.25, logfc.threshold=0.25)

