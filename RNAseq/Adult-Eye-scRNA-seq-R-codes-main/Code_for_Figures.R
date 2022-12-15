# R code for Figure generation

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

# Load in processed seurat objects for 1 day, 3 days and 7 days old males where clusters are annotated with major cell types (e.g. R7, R8), as a1d, a3d and a7d, respectively
# Load in processed seurat objects for 1 day old females, as a1f
# Load in processed seurat objects for 1 day male and female Harmony integrated data, as a1dh

# load in processed seurat objects for 1 day, 3 days and 7 days old males where clusters are annotated with Rh subtypes instead of R7 or R8 as a1dr, a3dr and a7dr, respectively. 
# load in processed seurat objects for Rh removed 1 day, 3 days and 7 days old male eye datasets as the following:
# 1 day old male: a1dr34, a1dr56 and a1dra, as 1 day old Rh3,4 removed, Rh5,6 removed and Rh3-6 removed datasets respectively
# 3 day old male: a3dr34, a3dr56 and a3dra, as 1 day old Rh3,4 removed, Rh5,6 removed and Rh3-6 removed datasets respectively
# 7 day old male: a7dr34, a7dr56 and a7dra, as 1 day old Rh3,4 removed, Rh5,6 removed and Rh3-6 removed datasets respectively

# Figure 1 F - H
DimPlot(a1d, label = F, pt.size = 2) + NoAxes() + NoLegend() #1F
DimPlot(a3d, label = F, pt.size = 2) + NoAxes() + NoLegend() #1G
DimPlot(a7d, label = F, pt.size = 2) + NoAxes() + NoLegend() #1H

# Figure 2 A
DimPlot(a7d, label = F, pt.size = 2, cols = c("grey","grey","grey","red","grey","grey")) + NoAxes() + NoLegend()
# Figure 2 B-E
FeaturePlot(a7d, features = "Rh5", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 B
FeaturePlot(a7d, features = "Rh6", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 C
FeaturePlot(a7d, features = "sens", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 D
FeaturePlot(a7d, features = "CG2082", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 E
# Figure 2 I
DimPlot(a7d, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","red","grey")) + NoAxes() + NoLegend()
# Figure 2 J-M
FeaturePlot(a7d, features = "Rh3", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 J
FeaturePlot(a7d, features = "Rh4", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 K
FeaturePlot(a7d, features = "pros", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 L
FeaturePlot(a7d, features = "igl", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #2 M

# Figure 3 A
DimPlot(a7d, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","grey","red")) + NoAxes() + NoLegend()
# Figure 3 B-D
FeaturePlot(a7d, features = "ninaE", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #3 B
FeaturePlot(a7d, features = "oc", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #3 C
FeaturePlot(a7d, features = "Zasp66", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #3 D
# Figure 3 F
VlnPlot(a7d, features = "ninaE") + NoLegend()

# Figure 4 A
DimPlot(a7d, label = F, pt.size = 3, cols = c("grey","grey","red","grey","grey","grey")) + NoAxes() + NoLegend()
# Figure 4 B-D
FeaturePlot(a7d, features = "ct", pt.size = 3, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #4 B
FeaturePlot(a7d, features = "Crys", pt.size = 3, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #4 C
FeaturePlot(a7d, features = "CG5597", pt.size = 3, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #4 D

# Figure 5 A
DimPlot(a7d, label = F, pt.size = 2, cols = c("red","red","grey","grey","grey","grey")) + NoAxes() + NoLegend()
# Figure 5 B-E
FeaturePlot(a7d, features = "w", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #5 B
FeaturePlot(a7d, features = "Pdh", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #5 C
FeaturePlot(a7d, features = "santa-maria", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #5 D
FeaturePlot(a7d, features = "wrapper", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #5 E
# Figure 5 F
VlnPlot(a7d, features = "Pdh") + NoLegend()

# Figure 6 A
DimPlot(a1f, label = F, pt.size = 2) + NoAxes() + NoLegend()
# Figure 6 B-D
DimPlot(a1dh, label = F, pt.size = 2) + NoAxes() + NoLegend() #6 B
DimPlot(a1dh, label = F, pt.size = 2, split.by = "orig.ident") + NoAxes() + NoLegend() #6 C, D
# Figure 6 E-G
FeaturePlot(a1dh, features = "CG6999", pt.size = 2, order = T, split.by = "orig.ident") & NoAxes() & NoLegend() & theme(plot.title= element_blank(), axis.title.y.right = element_blank()) #6E, F
FeaturePlot(a1dh, features = "t", pt.size = 2, order = T, split.by = "orig.ident") & NoAxes() & NoLegend() & theme(plot.title= element_blank(), axis.title.y.right = element_blank()) #6G, H
# Figure 6 I, J
a1dh_subset_vlnplot <- subset(a1dh, idents = c("1_Pigment","Cone"), invert = TRUE) # remove primary pigment and cone cells
VlnPlot(a1dh_subset_vlnplot, features = "CG6999", split.by = "orig.ident", split.plot = T) & theme(plot.title = element_text(face = "italic")) & NoLegend() #6 I
VlnPlot(a1dh_subset_vlnplot, features = "t", split.by = "orig.ident", split.plot = T) & theme(plot.title = element_text(face = "italic")) & NoLegend() #6 J

# Figure 7
DimPlot(a7dr, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7A
DimPlot(a7dr34, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7B
DimPlot(a7dr56, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7C
DimPlot(a7dra, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7D

DimPlot(a3dr, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7E
DimPlot(a3dr34, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7F
DimPlot(a3dr56, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7G
DimPlot(a3dra, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7H

DimPlot(a1dr, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7I
DimPlot(a1dr34, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7J
DimPlot(a1dr56, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7K
DimPlot(a1dra, pt.size = 2, label = F, cols = c("grey","grey","grey","red","blue","green3","orchid","grey")) & NoAxes() & NoLegend() #7L