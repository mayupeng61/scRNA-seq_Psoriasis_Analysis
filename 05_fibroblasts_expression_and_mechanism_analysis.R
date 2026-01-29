################################################
## 1. Load required R packages (duplicates are intentionally retained)
################################################

library(ggplot2)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(scRNAtoolVis)

obj <- readRDS("imq.rds")
Fibroblasts <- readRDS("Fibroblasts.rds")

################################################
## 2. Subsetting fibroblasts and myeloid cells
################################################

Tnc_Fibroblasts <- subset(
  Fibroblasts,
  idents = "Tnc-Fibroblasts"
)

Myeloid_Cells <- subset(
  obj,
  idents = "Myeloid Cells"
)

################################################
## 3. Defining sample factor levels
################################################

obj$sample <- factor(
  obj$sample,
  levels = c("Normal", "Psoriasis", "Celastrol")
)
Fibroblasts$sample <- factor(
  Fibroblasts$sample,
  levels = c("Normal", "Psoriasis", "Celastrol")
)
Myeloid_Cells$sample <- factor(
  Myeloid_Cells$sample,
  levels = c("Normal", "Psoriasis", "Celastrol")
)

################################################
## 4. Tnc expression and statistical analysis in fibroblasts
################################################

vln <- VlnPlot(
  Fibroblasts,
  features = "Tnc",
  group.by = "sample",
  pt.size = 0
) +
  stat_compare_means(
    comparisons = list(
      c("Normal", "Psoriasis"),
      c("Celastrol", "Psoriasis")
    ),
    label = "p.signif",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#F1ECE1",
      "Psoriasis" = "#b34644",
      "Celastrol" = "#518d6b"
    )
  )

# Normal vs Psoriasis
Tnc_normal_psoriasis <- FindMarkers(
  Fibroblasts,
  ident.1 = "Normal",
  ident.2 = "Psoriasis",
  features = "Tnc",
  group.by = "sample",
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)


# Celastrol vs Psoriasis
Tnc_celastrol_psoriasis <- FindMarkers(
  Fibroblasts,
  ident.1 = "Celastrol",
  ident.2 = "Psoriasis",
  features = "Tnc",
  group.by = "sample",
  test.use = "wilcox"
)

################################################
## 5. Tnf expression and statistical analysis in myeloid cells
################################################

Myeloid_Cells <- subset(obj, idents = "Myeloid Cells")

vln <- VlnPlot(
  Myeloid_Cells,
  features = "Tnf",
  group.by = "sample",
  pt.size = 0
) +
  stat_compare_means(
    comparisons = list(
      c("Normal", "Psoriasis"),
      c("Celastrol", "Psoriasis")
    ),
    label = "p.signif",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#F1ECE1",
      "Psoriasis" = "#b34644",
      "Celastrol" = "#518d6b"
    )
  )

# Normal vs Psoriasis
Tnf_normal_psoriasis <- FindMarkers(
  Myeloid_Cells,
  ident.1 = "Normal",
  ident.2 = "Psoriasis",
  features = "Tnf",
  group.by = "sample",
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)


# Celastrol vs Psoriasis
Tnf_celastrol_psoriasis <- FindMarkers(
  Myeloid_Cells,
  ident.1 = "Celastrol",
  ident.2 = "Psoriasis",
  features = "Tnf",
  group.by = "sample",
  test.use = "wilcox"
)

################################################
## 6. Tnfrsf1a expression and statistical analysis in fibroblasts
################################################

vln <- VlnPlot(
  Fibroblasts,
  features = "Tnfrsf1a",
  group.by = "sample",
  pt.size = 0
) +
  stat_compare_means(
    comparisons = list(
      c("Normal", "Psoriasis"),
      c("Celastrol", "Psoriasis")
    ),
    label = "p.signif",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#F1ECE1",
      "Psoriasis" = "#b34644",
      "Celastrol" = "#518d6b"
    )
  )

# Normal vs Psoriasis
Tnfrsf1a_normal_psoriasis <- FindMarkers(
  Fibroblasts,
  ident.1 = "Normal",
  ident.2 = "Psoriasis",
  features = "Tnfrsf1a",
  group.by = "sample",
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)

# Celastrol vs Psoriasis
Tnfrsf1a_celastrol_psoriasis <- FindMarkers(
  Fibroblasts,
  ident.1 = "Celastrol",
  ident.2 = "Psoriasis",
  features = "Tnfrsf1a",
  group.by = "sample",
  test.use = "wilcox"
)

################################################
## 7. featureCornerAxes visualization
################################################

featureCornerAxes(
  object = obj,
  reduction = "umap",
  groupFacet = "sample",
  relLength = 0.5,
  relDist = 0.2,
  features = c("Tnc")
)

################################################
## 8. Identification of all markers (preparation for volcano plot)
################################################


Fibroblasts.markers1 <- FindAllMarkers(
  Fibroblasts,
  only.pos = F
)

Fibroblasts.markers1 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> TOP10

################################################
## 9. UMAP corner label visualization of fibroblast subclusters
################################################

clusterCornerAxes(
  object = Fibroblasts,
  reduction = "umap",
  noSplit = TRUE,
  clusterCol = "cell_type"
)

clusterCornerAxes(
  object = Fibroblasts,
  reduction = "umap",
  clusterCol = "cell_type",
  noSplit = FALSE,
  groupFacet = "sample",
  aspect.ratio = 1,
  relLength = 0.5,
  cornerTextSize = 2,
  themebg = "bwCorner"
)




################################################
## 10. Volcano plot analysis
################################################

jjVolcano(
  diffData = Fibroblasts.markers1,
  tile.col = corrplot::COL2("RdYlBu", 15)[4:12],
  size = 3.5,
  fontface = "italic",
  polar = TRUE
)

################################################
## 11. Gene set definition
################################################

annoGene <- c(
  "Ccl2", "Cxcl12",
  "Il6", "Saa3"
)

immuneReceptorGenes <- c(
  "Ccr2",
  "Cxcr4",
  "Il6ra",
  "Il6st",
  "Tlr2",
  "Scarb1",
  "Tlr4"
)

immuneReceptorGenes2 <- c(
  "Ccr2",
  "Cxcr4",
  "Il6ra",
  "Tlr2",
  "Scarb1"
)

################################################
## 12. DotPlot visualization across different cell types
################################################

jjDotPlot(
  object = Fibroblasts,
  gene = annoGene,
  id = "cell_type"
)

jjDotPlot(
  object = obj,
  gene = immuneReceptorGenes,
  id = "Cell_Type"
)

jjDotPlot(
  object = Fibroblasts,
  gene = "Tlr4",
  id = "sample"
) + ggtitle("Fibroblasts")

jjDotPlot(
  object = Myeloid_Cells,
  gene = immuneReceptorGenes2,
  id = "sample"
) + ggtitle("Myeloid Cells")

jjDotPlot(
  object = Tnc_Fibroblasts,
  gene = annoGene,
  id = "sample",
  ytree = TRUE
) + ggtitle("Tnc-Fibroblasts")

jjDotPlot(
  object = Fibroblasts,
  gene = "Tnfrsf1a",
  id = "celltype",
  ytree = TRUE
) + ggtitle("Fibroblasts")

################################################
## 13. Cell proportion analysis
################################################

table(
  Idents(Fibroblasts),
  Fibroblasts$cell_type
)

Cellratio <- prop.table(
  table(Fibroblasts$cell_type, Fibroblasts$sample),
  margin = 2
)

print(Cellratio)

Cellratio <- as.data.frame(Cellratio)

Cellratio$Var2 <- factor(
  Cellratio$Var2,
  levels = c("Normal", "Psoriasis", "Celastrol")
)

allcolour <- c(
  "#20B2AA","#FFA500","#9370DB","#98FB98","#F08080",
  "#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF",
  "#FA8072","#7B68EE","#9400D3","#800080","#A0522D",
  "#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
  "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2",
  "#228B22","#E9967A","#4682B4","#32CD32","#F0E68C",
  "#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC",
  "#8B008B","#8B4513","#DEB887"
)

ggplot(Cellratio) +
  geom_bar(
    aes(x = Var2, y = Freq, fill = Var1),
    stat = "identity",
    width = 0.7,
    linewidth = 0.5,
    colour = "#222222"
  ) +
  theme_classic() +
  labs(x = "Sample", y = "Ratio") +
  scale_fill_manual(values = allcolour) +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5,
      linetype = "solid"
    )
  )
