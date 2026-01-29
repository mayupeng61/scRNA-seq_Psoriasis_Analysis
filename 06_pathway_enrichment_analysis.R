## =====================================================
## Load required packages
## =====================================================

## Base / data handling
library(methods)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(magrittr)
library(data.table)

## Visualization
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(ggsci)
library(aplot)
library(ggfun)
library(ggplotify)
library(ggridges)
library(gghalves)

## Seurat
library(Seurat)
library(SeuratObject)

## Differential & pathway analysis
library(limma)

library(GSVA)
library(fgsea)
library(msigdbr)

## irGSEA core dependencies
library(GSEABase)
library(AUCell)
library(UCell)
library(singscore)
library(decoupleR)
library(Nebulosa)
library(RobustRankAggreg)

## Parallel
library(doParallel)
library(doRNG)

## irGSEA
library(irGSEA)

## ======================================
## fgsea: All cells (Celastrol vs Psoriasis)
## ======================================
obj <- readRDS("imq.rds")

Idents(obj) <- "sample"

C_P <- FindMarkers(
  obj,
  ident.1 = "Celastrol",
  ident.2 = "Psoriasis",
  test.use = "wilcox",
  min.pct = 0.1,
  verbose = FALSE
)

C_P$ID <- rownames(C_P)

ranks <- C_P$avg_log2FC
names(ranks) <- rownames(C_P)

fgseaRes <- fgsea(
  pathways = genesets,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

print(head(fgseaResTidy))

draw_result <- fgseaResTidy %>%
  mutate(label = case_when(
    NES >  1 & padj < 0.05 ~ "up",
    NES < -1 & padj < 0.05 ~ "down",
    TRUE ~ "no"
  )) %>%
  arrange(desc(NES))

draw_result$label <- factor(draw_result$label, levels = c("down", "no", "up"))

ggplot(draw_result, aes(reorder(pathway, NES), NES, fill = label)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(
    values = c("down" = "#008020", "no" = "gray", "up" = "#08519C")
  ) +
  coord_flip() +
  labs(
    x = "Pathways",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw()
## ======================================
## fgsea: All cells (Psoriasis vs Normal)
## ======================================

P_N <- FindMarkers(
  obj,
  ident.1 = "Psoriasis",
  ident.2 = "Normal",
  test.use = "wilcox",
  min.pct = 0.1,
  verbose = FALSE
)

ranks2 <- P_N$avg_log2FC
names(ranks2) <- rownames(P_N)

fgseaRes <- fgsea(
  pathways = genesets,
  stats = ranks2,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
draw_result <- fgseaResTidy %>%
  mutate(label = case_when(
    NES >  1 & padj < 0.05 ~ "up",
    NES < -1 & padj < 0.05 ~ "down",
    TRUE ~ "no"
  )) %>%
  arrange(desc(NES))

draw_result$label <- factor(draw_result$label, levels = c("down", "no", "up"))

ggplot(draw_result, aes(reorder(pathway, NES), NES, fill = label)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(
    values = c("down" = "#008020", "no" = "gray", "up" = "#08519C")
  ) +
  coord_flip() +
  labs(
    x = "Pathways",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw()
## ======================================
## fgsea: Fibroblasts (Celastrol vs Psoriasis)
## ======================================
Fibroblasts <- readRDS("Fibroblasts.rds")

Idents(Fibroblasts) <- "sample"

C_P <- FindMarkers(
  Fibroblasts,
  ident.1 = "Celastrol",
  ident.2 = "Psoriasis",
  test.use = "wilcox",
  min.pct = 0.1,
  verbose = FALSE
)

C_P$ID <- rownames(C_P)

ranks <- C_P$avg_log2FC
names(ranks) <- rownames(C_P)

fgseaRes <- fgsea(
  pathways = genesets,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

print(head(fgseaResTidy))
draw_result <- fgseaResTidy %>%
  mutate(label = case_when(
    NES >  1 & padj < 0.05 ~ "up",
    NES < -1 & padj < 0.05 ~ "down",
    TRUE ~ "no"
  )) %>%
  arrange(desc(NES))

draw_result$label <- factor(draw_result$label, levels = c("down", "no", "up"))

ggplot(draw_result, aes(reorder(pathway, NES), NES, fill = label)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(
    values = c("down" = "#008020", "no" = "gray", "up" = "#08519C")
  ) +
  coord_flip() +
  labs(
    x = "Pathways",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw()

## ======================================
## fgsea: Fibroblasts (Psoriasis vs Normal)
## ======================================

P_N <- FindMarkers(
  Fibroblasts,
  ident.1 = "Psoriasis",
  ident.2 = "Normal",
  test.use = "wilcox",
  min.pct = 0.1,
  verbose = FALSE
)

ranks2 <- P_N$avg_log2FC
names(ranks2) <- rownames(P_N)

fgseaRes <- fgsea(
  pathways = genesets,
  stats = ranks2,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
draw_result <- fgseaResTidy %>%
  mutate(label = case_when(
    NES >  1 & padj < 0.05 ~ "up",
    NES < -1 & padj < 0.05 ~ "down",
    TRUE ~ "no"
  )) %>%
  arrange(desc(NES))

draw_result$label <- factor(draw_result$label, levels = c("down", "no", "up"))

ggplot(draw_result, aes(reorder(pathway, NES), NES, fill = label)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(
    values = c("down" = "#008020", "no" = "gray", "up" = "#08519C")
  ) +
  coord_flip() +
  labs(
    x = "Pathways",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw()

## ======================================
## fgsea: Tnc-Fibroblasts (Celastrol vs Psoriasis)
## ======================================
Tnc_Fibroblasts <- subset(
  Fibroblasts,
  idents = "Tnc-Fibroblasts"
)

Idents(Tnc_Fibroblasts) <- "sample"

C_P <- FindMarkers(
  Tnc_Fibroblasts,
  ident.1 = "Celastrol",
  ident.2 = "Psoriasis",
  test.use = "wilcox",
  min.pct = 0.1,
  verbose = FALSE
)

C_P$ID <- rownames(C_P)

ranks <- C_P$avg_log2FC
names(ranks) <- rownames(C_P)

fgseaRes <- fgsea(
  pathways = genesets,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

print(head(fgseaResTidy))
draw_result <- fgseaResTidy %>%
  mutate(label = case_when(
    NES >  1 & padj < 0.05 ~ "up",
    NES < -1 & padj < 0.05 ~ "down",
    TRUE ~ "no"
  )) %>%
  arrange(desc(NES))

draw_result$label <- factor(draw_result$label, levels = c("down", "no", "up"))

ggplot(draw_result, aes(reorder(pathway, NES), NES, fill = label)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(
    values = c("down" = "#008020", "no" = "gray", "up" = "#08519C")
  ) +
  coord_flip() +
  labs(
    x = "Pathways",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw()

## ======================================
## fgsea: Tnc-Fibroblasts (Psoriasis vs Normal)
## ======================================

P_N <- FindMarkers(
  Tnc_Fibroblasts,
  ident.1 = "Psoriasis",
  ident.2 = "Normal",
  test.use = "wilcox",
  min.pct = 0.1,
  verbose = FALSE
)

ranks2 <- P_N$avg_log2FC
names(ranks2) <- rownames(P_N)

fgseaRes <- fgsea(
  pathways = genesets,
  stats = ranks2,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

draw_result <- fgseaResTidy %>%
  mutate(label = case_when(
    NES >  1 & padj < 0.05 ~ "up",
    NES < -1 & padj < 0.05 ~ "down",
    TRUE ~ "no"
  )) %>%
  arrange(desc(NES))

draw_result$label <- factor(draw_result$label, levels = c("down", "no", "up"))

ggplot(draw_result, aes(reorder(pathway, NES), NES, fill = label)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(
    values = c("down" = "#008020", "no" = "gray", "up" = "#08519C")
  ) +
  coord_flip() +
  labs(
    x = "Pathways",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw()

## ======================================
## irGSEA score: pathway activity scoring
## ======================================

Fibroblasts <- irGSEA.score(
  object = Fibroblasts,
  assay = "RNA",
  slot = "data",
  seeds = 123,
  ncores = 4,
  min.cells = 3,
  min.feature = 0,
  species = "Mus musculus",
  category = "H",
  custom = FALSE,
  geneset = NULL,
  msigdb = TRUE,
  subcategory = NULL,
  geneid = "symbol",
  method = c(
    "AUCell",
    "UCell",
    "singscore",
    "ssgsea",
    "JASMINE",
    "viper"
  ),
  aucell.MaxRank = NULL,
  ucell.MaxRank = NULL,
  kcdf = "Gaussian"
)

## ======================================
## irGSEA integrate: subtype-level comparison
## ======================================

result.dge <- irGSEA.integrate(
  object = Fibroblasts,
  group.by = "cell_type",
  metadata = NULL,
  col.name = NULL,
  method = c(
    "AUCell",
    "UCell",
    "singscore",
    "ssgsea",
    "JASMINE",
    "viper"
  )
)

## ======================================
## irGSEA heatmap (RRA integration)
## ======================================

irGSEA.heatmap.plot <- irGSEA.heatmap(
  object = result.dge,
  method = "RRA",
  top = 50,
  show.geneset = NULL
)

irGSEA.heatmap.plot

## ======================================
## irGSEA density scatterplot (UMAP)
## ======================================

scatterplot <- irGSEA.density.scatterplot(
  object = Fibroblasts,
  method = "UCell",
  show.geneset = "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  reduction = "umap"
)

scatterplot

scatterplot <- irGSEA.density.scatterplot(
  object = Fibroblasts,
  method = "UCell",
  show.geneset = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  reduction = "umap"
)

scatterplot

## ======================================
## irGSEA half violin plot
## ======================================

halfvlnplot <- irGSEA.halfvlnplot(
  object = Fibroblasts,
  method = "UCell",
  show.geneset = "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

halfvlnplot

halfvlnplot <- irGSEA.halfvlnplot(
  object = Fibroblasts,
  method = "UCell",
  show.geneset = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

halfvlnplot
