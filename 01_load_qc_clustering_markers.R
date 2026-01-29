############################################
## Single-cell RNA-seq analysis main script
## Load data, QC, clustering and markers
############################################

## ========== Load packages ==========
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ComplexHeatmap)
library(RColorBrewer)

## ========== Load raw 10X data ==========
nl <- Read10X(data.dir = "NORMAL")
nl <- CreateSeuratObject(
  counts = nl,
  project = "nl",
  min.cells = 3,
  min.features = 200
)

nl


imq <- Read10X(data.dir = "IMQ")
imq <- CreateSeuratObject(
  counts = imq,
  project = "imq",
  min.cells = 3,
  min.features = 200
)

imq


ce <- Read10X(data.dir = "CE")
ce <- CreateSeuratObject(
  counts = ce,
  project = "ce",
  min.cells = 3,
  min.features = 200
)

ce


## ========== Merge datasets ==========
obj <- merge(
  nl,
  y = list(imq, ce),
  add.cell.ids = c("nl", "imq", "ce"),
  project = "obj"
)


## ========== Quality control ==========
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^Rp[ls]")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[ab]")

VlnPlot(
  obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.rp",
    "percent.hb"
  ),
  ncol = 3,
  pt.size = 0,
  group.by = "orig.ident"
)

obj <- subset(
  obj,
  percent.mt < 10 &
    nFeature_RNA < 7500 &
    nCount_RNA < 100000 &
    percent.rp < 40 &
    percent.hb < 1
)


## ========== Normalization & integration ==========
obj <- NormalizeData(obj) %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.)) %>%
  IntegrateLayers(HarmonyIntegration) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  RunTSNE(reduction = "harmony", dims = 1:15)


## ========== Visualization ==========
DimPlot(obj, reduction = "umap", group.by = "Cell_Type")
ElbowPlot(obj)
UMAPPlot(obj) + TSNEPlot(obj)
DimPlot(obj, split.by = "orig.ident")



## ========== Marker gene identification ==========
obj <- JoinLayers(obj)

obj.markers <- FindAllMarkers(
  obj,
  only.pos = TRUE
)

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


## ========== Save object ==========
saveRDS(obj.markers, file = "imq.markers.rds")
saveRDS(obj, file = "imq.rds")