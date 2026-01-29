############################################################
## CellChat multi-sample analysis (Normal / IMQ / Celastrol)
## Species: Mouse
## Purpose:
##  - Infer cell–cell communication networks
##  - Compare interaction number and strength across conditions
##  - Focus on TNF and TENASCIN signaling pathways
############################################################

## =========================
## 1. Load required packages
## =========================
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(patchwork)
library(gridExtra)

obj <- readRDS("imq.rds")

## =========================
## 2. Split Seurat object by condition
## =========================
Idents(obj) <- 'sample'

NL  <- subset(obj, idents = "Normal")
CE  <- subset(obj, idents = "Celastrol")
IMQ <- subset(obj, idents = "Psoriasis")

## =========================================================
## 3. CellChat analysis: Normal (NL)
## =========================================================
data.input1 <- GetAssayData(object = NL, layer = 'data', assay = "RNA")

identity1 <- data.frame(
  group = NL$Cell_Type,
  row.names = names(NL$Cell_Type)
)

unique(identity1$group)

cellchat_NL <- createCellChat(object = data.input1)
summary(cellchat_NL)

cellchat_NL <- addMeta(cellchat_NL, meta = identity1, meta.name = "labels")
cellchat_NL <- setIdent(cellchat_NL, ident.use = "labels")
levels(cellchat_NL@idents)

cellchat_NL@DB <- CellChatDB.mouse
cellchat_NL <- subsetData(cellchat_NL)
cellchat_NL <- identifyOverExpressedGenes(cellchat_NL)
cellchat_NL <- identifyOverExpressedInteractions(cellchat_NL)
cellchat_NL <- computeCommunProb(cellchat_NL, raw.use = TRUE, population.size = TRUE)
cellchat_NL <- computeCommunProbPathway(cellchat_NL)
cellchat_NL <- aggregateNet(cellchat_NL)
cellchat_NL <- netAnalysis_computeCentrality(cellchat_NL, slot.name = "netP")

saveRDS(cellchat_NL, file = "cellchat_NL.rds")
cellchat_NL <- readRDS("cellchat_NL.rds")

## =========================================================
## 4. CellChat analysis: Celastrol (CE)
## =========================================================
data.input2 <- GetAssayData(object = CE, layer = 'data', assay = "RNA")

identity2 <- data.frame(
  group = CE$Cell_Type,
  row.names = names(CE$Cell_Type)
)

unique(identity2$group)

cellchat_CE <- createCellChat(object = data.input2)
summary(cellchat_CE)

cellchat_CE <- addMeta(cellchat_CE, meta = identity2, meta.name = "labels")
cellchat_CE <- setIdent(cellchat_CE, ident.use = "labels")
levels(cellchat_CE@idents)

cellchat_CE@DB <- CellChatDB.mouse
cellchat_CE <- subsetData(cellchat_CE)
cellchat_CE <- identifyOverExpressedGenes(cellchat_CE)
cellchat_CE <- identifyOverExpressedInteractions(cellchat_CE)
cellchat_CE <- computeCommunProb(cellchat_CE, raw.use = TRUE, population.size = TRUE)
cellchat_CE <- computeCommunProbPathway(cellchat_CE)
cellchat_CE <- aggregateNet(cellchat_CE)
cellchat_CE <- netAnalysis_computeCentrality(cellchat_CE, slot.name = "netP")

saveRDS(cellchat_CE, file = "cellchat_CE.rds")
cellchat_CE <- readRDS("cellchat_CE.rds")

## =========================================================
## 5. CellChat analysis: Psoriasis (IMQ)
## =========================================================
data.input3 <- GetAssayData(object = IMQ, layer = 'data', assay = "RNA")

identity3 <- data.frame(
  group = IMQ$Cell_Type,
  row.names = names(IMQ$Cell_Type)
)

unique(identity3$group)

cellchat_IMQ <- createCellChat(object = data.input3)
summary(cellchat_IMQ)

cellchat_IMQ <- addMeta(cellchat_IMQ, meta = identity3, meta.name = "labels")
cellchat_IMQ <- setIdent(cellchat_IMQ, ident.use = "labels")
levels(cellchat_IMQ@idents)

cellchat_IMQ@DB <- CellChatDB.mouse
cellchat_IMQ <- subsetData(cellchat_IMQ)
cellchat_IMQ <- identifyOverExpressedGenes(cellchat_IMQ)
cellchat_IMQ <- identifyOverExpressedInteractions(cellchat_IMQ)
cellchat_IMQ <- computeCommunProb(cellchat_IMQ, raw.use = TRUE, population.size = TRUE)
cellchat_IMQ <- computeCommunProbPathway(cellchat_IMQ)
cellchat_IMQ <- aggregateNet(cellchat_IMQ)
cellchat_IMQ <- netAnalysis_computeCentrality(cellchat_IMQ, slot.name = "netP")

saveRDS(cellchat_IMQ, file = "cellchat_IMQ.rds")
cellchat_IMQ <- readRDS("cellchat_IMQ.rds")

## =========================================================
## 6. Overview comparison of interaction numbers
## =========================================================
object.list <- list(
  NL  = cellchat_NL,
  IMQ = cellchat_IMQ,
  CE  = cellchat_CE
)

weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))

par(mfrow = c(1,3), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2],
    edge.width.max = 12,
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

## =========================================================
## 7. Merge CellChat objects for differential analysis
## =========================================================
cellchat_1   <- mergeCellChat(list(cellchat_NL, cellchat_IMQ),
                              add.names = c("NL", "IMQ"), cell.prefix = FALSE)

cellchat_2   <- mergeCellChat(list(cellchat_IMQ, cellchat_CE),
                              add.names = c("IMQ", "CE"), cell.prefix = FALSE)

cellchat_all <- mergeCellChat(list(cellchat_NL, cellchat_IMQ, cellchat_CE),
                              add.names = c("NL", "IMQ", "CE"), cell.prefix = FALSE)

netVisual_diffInteraction(cellchat_1, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_1, weight.scale = TRUE, measure = "weight")
netVisual_diffInteraction(cellchat_2, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_2, weight.scale = TRUE, measure = "weight")

## =========================================================
## 8. TNF signaling analysis
## =========================================================
pathways.show <- c("TNF")

weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathways.show)

par(mfrow = c(1,3), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

pairLR.TNF <- extractEnrichedLR(object.list$CE,
                                signaling = pathways.show,
                                geneLR.return = FALSE)

LR.show <- pairLR.TNF[1, ]

plotGeneExpression(cellchat_all, signaling = "TNF", type = "violin")

cellchat_all@meta$datasets <- factor(
  cellchat_all@meta$datasets,
  levels = c("NL", "IMQ", "CE")
)

plotGeneExpression(cellchat_all,
                   signaling = "TNF",
                   split.by = "datasets",
                   type = "violin")

## =========================================================
## 9. TENASCIN signaling analysis
## =========================================================
pathways.show <- c("TENASCIN")

weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathways.show)

par(mfrow = c(1,3), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

pairLR.TENASCIN <- extractEnrichedLR(object.list$CE,
                                     signaling = pathways.show,
                                     geneLR.return = FALSE)

## =========================================================
## 10. Bubble plot for TENASCIN ligand–receptor pairs
## =========================================================
p <- netVisual_bubble(
  cellchat_all,
  sources.use = c(1),
  targets.use = c(1,2,3,4,5,6,7),
  comparison = c(1,2,3),
  pairLR.use = pairLR.TENASCIN,
  angle.x = 45
)

p
