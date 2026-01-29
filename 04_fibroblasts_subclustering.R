################################################
## 1. Fibroblast subset extraction
################################################
obj <- readRDS("imq.rds")
# Subset fibroblasts
Fibroblasts <- subset(obj, idents = "Fibroblasts")

################################################
## 2. Re-normalization and feature selection
################################################

# Data normalization
Fibroblasts <- NormalizeData(Fibroblasts)
Fibroblasts <- FindVariableFeatures(
  Fibroblasts,
  selection.method = "vst",
  nfeatures = 2000
)

################################################
## 3. PCA and UMAP dimensional reduction
################################################

# PCA
Fibroblasts <- ScaleData(Fibroblasts)
Fibroblasts <- RunPCA(
  Fibroblasts,
  features = VariableFeatures(object = Fibroblasts)
)
ElbowPlot(Fibroblasts)

# UMAP
Fibroblasts <- RunUMAP(Fibroblasts, dims = 1:15)

################################################
## 4. Clustering analysis
################################################

Fibroblasts <- FindNeighbors(Fibroblasts, dims = 1:15)
Fibroblasts <- FindClusters(Fibroblasts, resolution = 0.1)
Fibroblasts <- JoinLayers(Fibroblasts)

################################################
## 5. Marker gene identification
################################################

Fibroblasts.markers <- FindAllMarkers(
  Fibroblasts,
  only.pos = TRUE
)

Fibroblasts.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


################################################
## 6. Marker expression validation (VlnPlot)
################################################

VlnPlot(Fibroblasts, features = c("Thy1","Tgfbi"))   # 0
VlnPlot(Fibroblasts, features = c("Tnc"))            # 1
VlnPlot(Fibroblasts, features = c("Cygb"))           # 2
VlnPlot(Fibroblasts, features = c("Col11a1","Acan")) # 3
VlnPlot(Fibroblasts, features = c("Cxcl13"))         # 4
VlnPlot(Fibroblasts, features = c("Rspo3"))          # 5
VlnPlot(Fibroblasts, features = c("Mki67"))          # 6
VlnPlot(Fibroblasts, features = c("Pax7"))           # 7

################################################
## 7. Fibroblast subcluster annotation
################################################

# Define new cluster identities
celltype <- c(
  "Thy1-Fibroblasts",
  "Tnc-Fibroblasts",
  "Cygb-Fibroblasts",
  "Col11a1-Fibroblasts",
  "Cxcl13-Fibroblasts",
  "Col23a1-Fibroblasts",
  "Mki67-Fibroblasts",
  "Pax7-Fibroblasts"
)

# Map cluster IDs to cell types
names(celltype) <- as.character(0:7)

# Update cell identity labels
Fibroblasts <- RenameIdents(Fibroblasts, celltype)
cell_type <- Idents(Fibroblasts)
Fibroblasts <- AddMetaData(
  Fibroblasts,
  metadata = cell_type,
  col.name = "cell_type"
)

table(Fibroblasts$cell_type)

################################################
## 8. Add sample metadata
################################################

Fibroblasts <- Fibroblasts %>%
  AddMetaData(
    metadata = ifelse(
      Fibroblasts[["orig.ident"]] == "ce",
      "Celastrol",
      ifelse(
        Fibroblasts[["orig.ident"]] == "imq",
        "Psoriasis",
        "Normal"
      )
    ),
    col.name = "sample"
  )

Fibroblasts$sample <- factor(
  Fibroblasts$sample,
  levels = c("Normal", "Psoriasis", "Celastrol")
)

################################################
## 9. Save fibroblast object
################################################

saveRDS(Fibroblasts, file = "Fibroblasts.rds")
