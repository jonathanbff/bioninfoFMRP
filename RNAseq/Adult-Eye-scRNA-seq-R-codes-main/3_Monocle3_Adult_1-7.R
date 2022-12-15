#Monocle3 for soupX adult eye day 1 - day 7

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(sctransform)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

#load in Seurat objects of processed adult 1-day, 3-day and 7-day old eye data
a1d <- readRDS("~/A1D.rds")
a3d <- readRDS("~/A3D.rds")
a7d <- readRDS("~/A7D.rds")

#Get raw count data from each Seurat object
a1d.rawdata <- as.matrix(GetAssayData(a1d, slot = "counts", assay = "RNA"))
a3d.rawdata <- as.matrix(GetAssayData(a3d, slot = "counts", assay = "RNA"))
a7d.rawdata <- as.matrix(GetAssayData(a7d, slot = "counts", assay = "RNA"))

#Create new seurat objects from raw count data
a1dm <- CreateSeuratObject(counts = a1d.rawdata, project = "a1d")
a3dm <- CreateSeuratObject(counts = a3d.rawdata, project = "a3d")
a7dm <- CreateSeuratObject(counts = a7d.rawdata, project = "a7d")

#transform each seurat object into cell_data_set objects for Monocle3 to process
cds1 <- as.cell_data_set(a1dm)
cds3 <- as.cell_data_set(a3dm)
cds7 <- as.cell_data_set(a7dm)

#as.cell_data_set did not bring the gene name information from the seurat object
#We need to estimate size factors and add back gene names into CDS
cds1 <- estimate_size_factors(cds1)
cds3 <- estimate_size_factors(cds3)
cds7 <- estimate_size_factors(cds7)

cds1@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(a1d[["RNA"]])
cds3@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(a3d[["RNA"]])
cds7@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(a7d[["RNA"]])

# Can now process cds as shown in Monocle3 tutorial
big_cds <- combine_cds(list(cds1, cds3,cds7))
adcds <- preprocess_cds(big_cds, num_dim = 100)
adcds <- reduce_dimension(adcds)
adcds <- cluster_cells(adcds, resolution=1e-5)
plot_cells(adcds)

#to check for any batch effects
plot_cells(adcds, color_cells_by="orig.ident", label_cell_groups=FALSE)

#there are significant batch effects in our merged CDS, therefore we need to perform batch effects correction.
#renamed new cds object as adcds2
adcds2 <- align_cds(adcds, num_dim = 100, alignment_group = "orig.ident")
adcds2 <- reduce_dimension(adcds2)
plot_cells(adcds2, color_cells_by="orig.ident", label_cell_groups=FALSE)
#Batch effects corrected

#cluster cell with higher resolution for cell cluster identity annotation
adcds2 <- cluster_cells(adcds2, resolution=1e-3)
#assign cell cluster identity based on marker genes
colData(adcds2)$assigned_cell_type <- as.character(clusters(adcds2))
colData(adcds2)$assigned_cell_type = dplyr::recode(colData(adcds2)$assigned_cell_type,
                                                   "31"="R1-6",
                                                   "8"="R1-6",
                                                   "16"="R1-6",
                                                   "11"="R1-6",
                                                   "32"="R1-6",
                                                   "24"="R1-6",
                                                   "25"="R1-6",
                                                   "17"="R1-6",
                                                   "27"="R1-6",
                                                   "19"="R1-6",
                                                   "4"="R1-6",
                                                   "14"="R1-6",
                                                   "10"="R1-6",
                                                   "13"="R1-6",
                                                   "3"="R1-6",
                                                   "9"="R1-6",
                                                   "5"="R1-6",
                                                   "22"="R7",
                                                   "23"="R7",
                                                   "15"="R7",
                                                   "30"="R7",
                                                   "34"="R7",
                                                   "6"="R8",
                                                   "29"="R8",
                                                   "35"="R8",
                                                   "1"="2,3_pigments",
                                                   "21"="2,3_pigments",
                                                   "18"="2,3_pigments",
                                                   "28"="2,3_pigments",
                                                   "2"="2,3_pigments",
                                                   "12"="2,3_pigments",
                                                   "7"="2,3_pigments",
                                                   "20"="2,3_pigments",
                                                   "26"="2,3_pigments",
                                                   "37"="2,3_pigments",
                                                   "33"="1_pigments",
                                                   "36"="Cone")

#To assign cell clusters with colors that correspond with other figures
cell_type_color <- c("R1-6" = "#f263df",
                     "R7" = "#639dfc",
                     "R8" = "#1bbec3",
                     "Cone" = "#17b93f",
                     "1_pigments" = "#b69e1c",
                     "2,3_pigments" = "#f6756f")
plot_cells(adcds2, color_cells_by="assigned_cell_type", label_cell_groups = F)+scale_color_manual(values = cell_type_color)


#trajectory analyses
#learngraph, tells monocle to learn the trajectory
adcds2 <- learn_graph(adcds2)

#order cells based on pseudotime.  first pick the root node (early in development)
adcds2 <- order_cells(adcds2)

#plot by pseudotime
plot_cells(adcds2, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
