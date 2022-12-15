# Harmony Integration of Adult 1 Day Male and Females

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(SeuratData)
library(SeuratWrappers)
library(harmony)

#load in Seurat objects of clustered adult 1 day male and female
# a1d is the adult 1 day male data set
# a1f is the adult 1 day female data set

# Get the counts from each and create a fresh seurat object with them
a1d.rawdata <- as.matrix(GetAssayData(a1d, slot = "counts", assay = "RNA"))
a1f.rawdata <- as.matrix(GetAssayData(a1f, slot = "counts", assay = "RNA"))

a1dm <- CreateSeuratObject(counts = a1d.rawdata, project = "A1d_male")
a1dfm <- CreateSeuratObject(counts = a1f.rawdata, project = "A1d_female")
# we save the original cluster annotation to each of these seurat objects as "old.ident"
a1dm[["old.ident"]] <- Idents(object = a1d)
a1fm[["old.ident"]] <- Idents(object = a1f)

# merge and add %mito, then SCT and RunPCA
a1dh <- merge(a1dm, y = a1dfm, add.cell.ids = c("M", "F"), project = "a1d_merge")
a1dh <- PercentageFeatureSet(a1dh, pattern = "mt:", col.name = "percent.mt")
a1dh <- SCTransform(a1dh, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a1dh <- RunPCA(a1dh, verbose = TRUE)

# Run Harmony. we will use SCT as the assay use
# resolution of 2 is used for FindCluster step
a1dh <- RunHarmony(a1dh, group.by.vars = "orig.ident", assay.use = "SCT")

a1dh <- RunUMAP(a1dh, reduction = "harmony", dims = 1:50)
a1dh <- FindNeighbors(a1dh, reduction = "harmony", dims = 1:50)
a1dh <- FindClusters(a1dh, resolution = 2)
DimPlot(a1dh, group.by = "orig.ident")
DimPlot(a1dh, label = T)

# We assigned cluster identities based on marker genes similarly to how the rest of the data sets were assigned.
a1dh_ct <- a1dh

r7.name <- WhichCells(object = a1dh_ct, idents = c(22,33,26,10,32,17,5,20))
r8.name <- WhichCells(object = a1dh_ct, idents = c(27,25,15,18))
r16.name <- WhichCells(object = a1dh_ct, idents = c(3,4,6,8,19,12,13,7,21,24,16,9,14,31))
cone.name <- WhichCells(object = a1dh_ct,idents = c(35))
pigm.name <- WhichCells(object = a1dh_ct, idents = c(36,34,29,23,28,11,0,1,2,28))
pigm1.name <- WhichCells(object = a1dh_ct, idents = c(30))
a1dh_ct <- SetIdent(object = a1dh_ct, cells = r16.name, value = "R1-6")
a1dh_ct <- SetIdent(object = a1dh_ct, cells = r7.name, value = "R7")
a1dh_ct <- SetIdent(object = a1dh_ct, cells = r8.name, value = "R8")
a1dh_ct <- SetIdent(object = a1dh_ct, cells = cone.name, value = "Cone")
a1dh_ct <- SetIdent(object = a1dh_ct, cells = pigm1.name, value = "1_Pigment")
a1dh_ct <- SetIdent(object = a1dh_ct, cells = pigm.name, value = "2,3_Pigments")

a1dh_ct[["cell_type"]] <- Idents(object=a1dh_ct)

saveRDS(a1dh_ct, file = "A1D_F_M_Harmony_Cell_Type.rds")

# To find male or female specific genes in each cell type cluster
MvFa1dR7.markers <- FindMarkers(a1dh_ct, ident.1 = "A1d_male", ident.2 = "A1d_female",group.by = "orig.ident",subset.ident = "R7",min.pct = 0.25)
MvFa1dR8.markers <- FindMarkers(a1dh_ct, ident.1 = "A1d_male", ident.2 = "A1d_female",group.by = "orig.ident",subset.ident = "R8",min.pct = 0.25)
MvFa1dR16.markers <- FindMarkers(a1dh_ct, ident.1 = "A1d_male", ident.2 = "A1d_female",group.by = "orig.ident",subset.ident = "R1-6",min.pct = 0.25)
MvFa1dpigm.markers <- FindMarkers(a1dh_ct, ident.1 = "A1d_male", ident.2 = "A1d_female",group.by = "orig.ident",subset.ident = "2,3_Pigments",min.pct = 0.25)

# to plot ClusterPlot with newly assigned cluster identity
DimPlot(a1dh_ct, group.by = "cell_type")
# to plot ClusterPlot with cluster identity from male alone and female alone data sets
DimPlot(a1dh_ct, group.by = "old.ident")