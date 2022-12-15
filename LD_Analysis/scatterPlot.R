if(!require("ggplot2")) {
   install.packages("ggplot2", dependencies=TRUE)
   if(!require("ggplot2")) stop ("ggplot2 not found")
}

if(!require("ggrepel")) {
   install.packages("ggrepel", dependencies=TRUE)
   if(!require("ggrepel")) stop ("ggrepel not found")
}

if(!require("tidyr")) {
   install.packages("tidyr", dependencies=TRUE)
   if(!require("tidyr")) stop ("tidyr not found")
}

if(!require("dplyr")) {
   install.packages("dplyr", dependencies=TRUE)
   if(!require("dplyr")) stop ("dplyr not found")
}

if(!require("scales")) {
   install.packages("scales", dependencies=TRUE)
   if(!require("scales")) stop ("scales not found")
}

##########################################################
## CODE FOR GENERATING VOLCANO PLOTS                    ##
## ADAPTED FROM JULIA -- SCRIPT AUTOMATICALLY GENERATED ##
## LAST MODIFIED 3/20/2019 - KBRIN BIOINFORMATICS CORE  ##
##########################################################

Sys.setenv("DISPLAY"=":0")                     ## NEED SO X11 WORKS PROPERLY
df <- read.csv("LOD_scatter.txt", header=TRUE, sep="\t")
FN <- "LOD_scatterPlot.png"                                # output file name
title <-"LOD vs. Distance"     # plot title 

df$sig <- "LOD <= 3"
df$sig[df$LOD > 3] <- "LOD > 3"
df$sig[df$LOD > 200] <- "LOD > 200"

maxCoord <- 30000
minCoord <- 0
minPVal <- min(df$p_value)
minCoord
maxCoord
maxYCoord <- max(df$LOD)
df$gene <- rownames(df)
input <- df
topGenes <- head(input, 20)                                                    # top 20 genes
topGenes

png(FN, units="in", width=10, height=5, res=300, type="cairo")

############################
## CREATE LOD scatterplot ##
############################
p = ggplot(input, aes(DISTANCE, LOD)) +
  geom_point(aes(col=sig)) +                                                         # add points colored by significance
  scale_color_manual(values=c("black", "pink", "red")) + 
  labs(title=title, x ="Distance", y = "LOD")  +
  scale_x_continuous(limits = c(minCoord, maxCoord),breaks=pretty_breaks(n=11)) +                  # limits printing every other value
  scale_y_continuous(limits = c(0, maxYCoord),breaks=pretty_breaks(n=11)) +                 # limits printing every value

  geom_text_repel(data=topGenes, aes(label=gene), size = 3) +                 # adding text for the top 20 genes
  theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))
p
dev.off()

############################
## CREATE R^2 Scatterplot ##
############################
df <- read.csv("r2_scatter.txt", header=TRUE, sep="\t")
FN <- "r2_scatterPlot.png"                                # output file name
title <-"r^2 vs. Distance"     # plot title

df$sig <- "r^2 <= 0.8"
df$sig[df$r2 > 0.8] <- "r^2 > 0.8"
df$sig[df$r2 > 0.9] <- "r^2 > 0.9"

maxCoord <- 30000
minCoord <- 0
minCoord
maxCoord
maxYCoord <- max(df$r2)
df$gene <- rownames(df)
input <- df
topGenes <- head(input, 20)                                                    # top 20 genes
topGenes

png(FN, units="in", width=10, height=5, res=300, type="cairo")

############################
## CREATE LOD scatterplot ##
############################
p = ggplot(input, aes(DISTANCE, r2)) +
  geom_point(aes(col=sig)) +                                                         # add points colored by significance
  scale_color_manual(values=c("black", "pink", "red")) +
  labs(title=title, x ="Distance", y = "r^2")  +
  scale_x_continuous(limits = c(minCoord, maxCoord),breaks=pretty_breaks(n=11)) +                  # limits printing every other value
  scale_y_continuous(limits = c(0, maxYCoord),breaks=pretty_breaks(n=11)) +                 # limits printing every value

  geom_text_repel(data=topGenes, aes(label=gene), size = 3) +                 # adding text for the top 20 genes
  theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))
p
dev.off()


#####################################
## END OF R CODE FOR VOLCANO PLOTS ##
#####################################
sessionInfo()
