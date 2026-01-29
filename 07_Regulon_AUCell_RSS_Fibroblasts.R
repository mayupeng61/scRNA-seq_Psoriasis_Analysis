############################################################
## Regulon activity analysis in fibroblasts
## AUCell integration, visualization and RSS analysis
############################################################

## ===============================
## Load required libraries
## ===============================
library(data.table)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggsci)
library(reshape2)
library(SCENIC)
library(BiocParallel)

## ===============================
## Load fibroblast Seurat object
## ===============================
Fibroblasts <- readRDS("Fibroblasts.rds")
data<-as.data.frame(Fibroblasts@assays$RNA$data)
write.csv(data, file = "data_fibroblasts.csv", row.names = TRUE)
## ===============================
## Load AUCell score matrix
## ===============================
aucell_scores <- read.csv(
  "aucell_scores_t.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE
)

# Set regulon names as rownames
rownames(aucell_scores) <- aucell_scores$Regulon
aucell_scores <- aucell_scores[, -1]

## ===============================
## Add AUCell scores as a new assay
## ===============================
Fibroblasts[["AUC"]] <- CreateAssayObject(
  counts = aucell_scores
)
DefaultAssay(Fibroblasts) <- "AUC"

## ===============================
## Define color palette
## ===============================
my.colors <- colorRampPalette(
  c("lightblue", "darkred")
)(2)

## ===============================
## Define feature lists
## ===============================
Features <- rownames(aucell_scores)

TF_list <- c(
  "Atf3(+)","Atf4(+)","Bach1(+)","Bcl3(+)","Bhlhe40(+)",
  "Cebpb(+)","Cebpd(+)","Egr1(+)","Ets2(+)",
  "Fos(+)","Fosb(+)","Fosl1(+)","Fosl2(+)",
  "Foxp1(+)","Irf1(+)",
  "Jun(+)","Junb(+)","Jund(+)",
  "Klf4(+)","Mafk(+)","Mitf(+)",
  "Nfe2l1(+)","Nfe2l2(+)","Nfkb1(+)",
  "Nr1d1(+)","Sfpq(+)",
  "Stat3(+)","Tbx15(+)","Tcf4(+)","Twist1(+)"
)

## ===============================
## DotPlot: all regulons across samples
## ===============================
DotPlot(
  Fibroblasts,
  features = Features,
  group.by = "sample",
  cols = my.colors
) + RotatedAxis()

## ===============================
## DotPlot: selected TF regulons across samples
## ===============================
DotPlot(
  Fibroblasts,
  features = TF_list,
  group.by = "sample",
  cols = my.colors
) + RotatedAxis()

## ===============================
## Regulon Specificity Score (RSS) analysis
## ===============================
anno_col <- Fibroblasts@meta.data
head(anno_col)

rss <- calcRSS(
  AUC = aucell_scores,
  cellAnnotation = anno_col$cell_type
)

## ===============================
## RSS visualization for Tnc-Fibroblasts
## ===============================
plotRSS_oneSet(
  rss,
  "Tnc-Fibroblasts",
  n = 10
)

## ===============================
## Global RSS heatmap
## ===============================
rssPlot <- plotRSS(
  rss,
  col.low = "#FBE4D8",
  col.mid = "#522B5B",
  col.high = "#100019",
  verbose = FALSE
)

rssPlot

