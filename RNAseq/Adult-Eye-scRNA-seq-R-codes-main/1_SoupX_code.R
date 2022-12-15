# SoupX processing for adult data sets

# Run Seurat on unfiltered data sets without any cut offs to obtain cluster information

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

a7d.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
a7d <- CreateSeuratObject(counts = a7d.data, project = "a7d_male")
a7d <- PercentageFeatureSet(a7d, pattern = "mt:", col.name = "percent.mt")
a7d <- SCTransform(a7d, variable.features.n = 3500, vars.to.regress = "percent.mt", verbose = TRUE)
a7d <- RunPCA(a7d, verbose = TRUE)
a7d <- RunUMAP(a7d, dims = 1:50, verbose = TRUE)
a7d <- FindNeighbors(a7d, dims = 1:50, verbose = TRUE)
a7d <- FindClusters(a7d, verbose = TRUE)

#repeated for 1day, 3days and 1day female data sets, to make a1d, a3d and a1f, respectively

#SoupX
library(SoupX)
#load table of counts (toc, for the cells) and table of droplets (tod, for raw droplets)
toc = Seurat::Read10X(data.dir = "~/filtered_feature_bc_matrix/")
tod = Seurat::Read10X(data.dir = "~/raw_feature_bc_matrix/")
#load them into soup
sc = SoupChannel(tod, toc)

# set cluster information for the soupchannel
# the setnames are the seurat cluster # generated above
# colnames to bring the names for each cell associated with the cluster info
sc1 = sc
sc1 = setClusters(sc1, setNames(a1dn$seurat_clusters, colnames(a1dn)))

# Run SoupX with default parameters, auto estimate the contamination
sc1 = autoEstCont(sc1)
# output the adjusted matrix of counts
a7d_soupx_out = adjustCounts(sc1)

saveRDS(a7d_soupx_out, file = "A7D_SoupX_out.rds")

# This was repeated for the remaining data sets to generate the following rds files
# A1D_SoupX_out.rds
# A3D_SoupX_out.rds
# A1F_SoupX_out.rds

# These are loaded into R as objects and can be directly analyzed by Seurat