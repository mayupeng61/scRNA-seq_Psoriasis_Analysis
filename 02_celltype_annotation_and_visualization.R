############################################
## Cell type annotation and sample metadata
############################################

## ========== Load required packages ==========
library(Seurat)
library(dplyr)
library(scRNAtoolVis)
obj <- readRDS("imq.rds")

## ========== Define cell type annotations ==========
# Define cluster-to-cell-type mapping
CellType <- c(
  "Fibroblasts", "Fibroblasts", "Fibroblasts",
  "Keratinocytes", "Keratinocytes",
  "Smooth Muscle Cells",
  "Keratinocytes",
  "Myeloid Cells",
  "Keratinocytes",
  "Fibroblasts",
  "Keratinocytes",
  "Keratinocytes",
  "T Cells", "T Cells",
  "Endothelial Cells",
  "Fibroblasts",
  "Schwann Cells",
  "Myeloid Cells",
  "Fibroblasts",
  "Myeloid Cells",
  "T Cells",
  "Smooth Muscle Cells"
)

# Map original cluster IDs (0–21) to cell types
names(CellType) <- as.character(0:21)


## ========== Rename cluster identities ==========
# Seurat object is assumed to be named "obj"
obj <- RenameIdents(obj, CellType)


## ========== Add cell type metadata ==========
Cell_Type <- Idents(obj)

obj <- AddMetaData(
  obj,
  metadata = Cell_Type,
  col.name = "Cell_Type"
)

# Check cell counts per cell type
table(obj$Cell_Type)


## ========== Add sample group metadata ==========
# Map orig.ident to biological groups
obj <- obj %>%
  AddMetaData(
    metadata = ifelse(
      obj[["orig.ident"]] == "ce",
      "Celastrol",
      ifelse(
        obj[["orig.ident"]] == "imq",
        "Psoriasis",
        "Normal"
      )
    ),
    col.name = "sample"
  )


## ========== UMAP visualization by cell type ==========
clusterCornerAxes(
  object = obj,
  reduction = "umap",
  clusterCol = "Cell_Type",
  noSplit = TRUE,
  cellLabel = TRUE,
  cellLabelSize = 5
)


## ========== Marker gene expression validation ==========
features <- c(
  # Myeloid cells
  "Fcer1g", "Tyrobp", "Il1b", "Cd68", "Aif1",
  # Keratinocytes
  "Krt5", "Krt14", "Krt1", "Krt10",
  # T cells
  "Cd3d", "Cd3g", "Cd8a",
  # Smooth muscle cells
  "Acta2", "Myh11", "Cnn1",
  # Schwann cells
  "Mpz", "Pmp22", "S100b",
  # Fibroblasts
  "Col1a1", "Col1a2",
  # Endothelial cells
  "Pecam1"
)

jjDotPlot(
  object = obj,
  gene = features,
  id = "Cell_Type"
)


############################################
## Cell type proportion analysis
############################################

## ========== Calculate cell type ratios ==========
Cellratio <- prop.table(
  table(obj$Cell_Type, obj$sample),
  margin = 2
)

Cellratio <- as.data.frame(Cellratio)

# Set sample order
Cellratio$Var2 <- factor(
  Cellratio$Var2,
  levels = c("Normal", "Psoriasis", "Celastrol")
)


## ========== Define color palette ==========
allcolour <- c(
  "#20B2AA", "#FFA500", "#9370DB", "#98FB98", "#F08080",
  "#1E90FF", "#7CFC00", "#FFFF00", "#808000", "#FF00FF",
  "#FA8072", "#7B68EE", "#9400D3", "#800080", "#A0522D",
  "#D2B48C", "#D2691E", "#87CEEB", "#40E0D0", "#5F9EA0",
  "#FF1493", "#0000CD", "#008B8B", "#FFE4B5", "#8A2BE2",
  "#228B22", "#E9967A", "#4682B4", "#32CD32", "#F0E68C",
  "#FFFFE0", "#EE82EE", "#FF6347", "#6A5ACD", "#9932CC",
  "#8B008B", "#8B4513", "#DEB887"
)


## ========== Visualization ==========
library(ggplot2)

ggplot(Cellratio) +
  geom_bar(
    aes(x = Var2, y = Freq, fill = Var1),
    stat = "identity",
    width = 0.7,
    linewidth = 0.5,
    colour = "#222222"
  ) +
  theme_classic() +
  labs(
    x = "Sample",
    y = "Ratio"
  ) +
  scale_fill_manual(values = allcolour) +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5,
      linetype = "solid"
    )
  )

## ========== Save object ==========
saveRDS(obj.markers, file = "imq.markers.rds")
saveRDS(obj, file = "imq.rds")