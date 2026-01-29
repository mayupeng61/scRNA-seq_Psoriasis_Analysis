# ==============================================================================
# Script Name: Sensitivity Analysis and Ligand-Receptor Interaction (CellChat)
# Description: This script performs stratified downsampling on scRNA-seq data 
#              to control for cell proportion shifts, followed by comparative 
#              CellChat analysis across three conditions: Normal, Psoriasis (IMQ), 
#              and Celastrol-treated.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(gridExtra)

# ------------------------------------------------------------------------------
# 2. Data Loading and Preprocessing
# ------------------------------------------------------------------------------
# Load the Seurat object
obj <- readRDS("imq.rds")

# Ensure the active identity class is set to Cell_Type
Idents(obj) <- "Cell_Type" 

# ==============================================================================
# 3. Sensitivity Analysis: Stratified Downsampling
# ==============================================================================
# Goal: To eliminate potential bias caused by shifting cell proportions, we 
# downsample each cell type to the minimum number observed across any group.

# A. Calculate Downsampling Thresholds (Minimum across groups)
# ------------------------------------------------------------------------------
# Generate a count matrix and identify the minimum cell count for each cell type
# across all groups (e.g., if Fibroblasts are 3500 in NL, 3600 in IMQ, 2812 in CE,
# the threshold becomes 2812).
counts_matrix <- table(obj$Cell_Type, obj$sample)
min_per_type <- apply(counts_matrix, 1, min)

message("--- Downsampling thresholds per cell type ---")
print(min_per_type)

# B. Execute Stratified Downsampling Loop
# ------------------------------------------------------------------------------
downsampled_cells <- list()

# Iterate through each cell type
for (ct in names(min_per_type)) {
  threshold <- min_per_type[ct]
  
  # Iterate through each sample group (Normal, Psoriasis, Celastrol)
  for (grp in unique(obj$sample)) {
    
    # Identify target cell IDs for the specific group and cell type
    cells_target <- WhichCells(obj, expression = (sample == grp & Cell_Type == ct))
    
    # Downsampling logic:
    # If the cell count exceeds the threshold, random sampling is performed.
    # Otherwise, all cells are retained.
    if (length(cells_target) > threshold) {
      set.seed(123) # Set random seed for reproducibility
      selected <- sample(cells_target, size = threshold)
    } else {
      selected <- cells_target
    }
    
    downsampled_cells <- c(downsampled_cells, selected)
  }
}

downsampled_cells <- unlist(downsampled_cells)

# C. Generate the Balanced Seurat Object
# ------------------------------------------------------------------------------
obj.balanced <- subset(obj, cells = downsampled_cells)

# D. Verify Results
# ------------------------------------------------------------------------------
# The cell counts for each cell type across different groups should now be identical.
message("--- Downsampling complete. Verifying balanced cell counts ---")
print(table(obj.balanced$Cell_Type, obj.balanced$sample))


# ==============================================================================
# 4. Split Object by Condition
# ==============================================================================
Idents(obj.balanced) <- 'sample'

NL  <- subset(obj.balanced, idents = "Normal")
CE  <- subset(obj.balanced, idents = "Celastrol")
IMQ <- subset(obj.balanced, idents = "Psoriasis")


# ==============================================================================
# 5. CellChat Analysis: Normal (NL)
# ==============================================================================
data.input1 <- GetAssayData(object = NL, layer = 'data', assay = "RNA")

identity1 <- data.frame(
  group = NL$Cell_Type,
  row.names = names(NL$Cell_Type)
)

# Initialize CellChat object
cellchat_NL <- createCellChat(object = data.input1)
summary(cellchat_NL)

# Add meta data and set identities
cellchat_NL <- addMeta(cellchat_NL, meta = identity1, meta.name = "labels")
cellchat_NL <- setIdent(cellchat_NL, ident.use = "labels")

# Set database and run standard processing
cellchat_NL@DB <- CellChatDB.mouse
cellchat_NL <- subsetData(cellchat_NL)
cellchat_NL <- identifyOverExpressedGenes(cellchat_NL)
cellchat_NL <- identifyOverExpressedInteractions(cellchat_NL)

# Compute communication probability
# Note: population.size = TRUE considers the effect of cell abundance, 
# though downsampling has equalized this to some extent.
cellchat_NL <- computeCommunProb(cellchat_NL, raw.use = TRUE, population.size = TRUE)
cellchat_NL <- computeCommunProbPathway(cellchat_NL)
cellchat_NL <- aggregateNet(cellchat_NL)
cellchat_NL <- netAnalysis_computeCentrality(cellchat_NL, slot.name = "netP")


# ==============================================================================
# 6. CellChat Analysis: Celastrol (CE)
# ==============================================================================
data.input2 <- GetAssayData(object = CE, layer = 'data', assay = "RNA")

identity2 <- data.frame(
  group = CE$Cell_Type,
  row.names = names(CE$Cell_Type)
)

# Initialize CellChat object
cellchat_CE <- createCellChat(object = data.input2)
summary(cellchat_CE)

# Add meta data and set identities
cellchat_CE <- addMeta(cellchat_CE, meta = identity2, meta.name = "labels")
cellchat_CE <- setIdent(cellchat_CE, ident.use = "labels")

# Set database and run standard processing
cellchat_CE@DB <- CellChatDB.mouse
cellchat_CE <- subsetData(cellchat_CE)
cellchat_CE <- identifyOverExpressedGenes(cellchat_CE)
cellchat_CE <- identifyOverExpressedInteractions(cellchat_CE)

# Compute communication probability
cellchat_CE <- computeCommunProb(cellchat_CE, raw.use = TRUE, population.size = TRUE)
cellchat_CE <- computeCommunProbPathway(cellchat_CE)
cellchat_CE <- aggregateNet(cellchat_CE)
cellchat_CE <- netAnalysis_computeCentrality(cellchat_CE, slot.name = "netP")


# ==============================================================================
# 7. CellChat Analysis: Psoriasis (IMQ)
# ==============================================================================
data.input3 <- GetAssayData(object = IMQ, layer = 'data', assay = "RNA")

identity3 <- data.frame(
  group = IMQ$Cell_Type,
  row.names = names(IMQ$Cell_Type)
)

# Initialize CellChat object
cellchat_IMQ <- createCellChat(object = data.input3)
summary(cellchat_IMQ)

# Add meta data and set identities
cellchat_IMQ <- addMeta(cellchat_IMQ, meta = identity3, meta.name = "labels")
cellchat_IMQ <- setIdent(cellchat_IMQ, ident.use = "labels")

# Set database and run standard processing
cellchat_IMQ@DB <- CellChatDB.mouse
cellchat_IMQ <- subsetData(cellchat_IMQ)
cellchat_IMQ <- identifyOverExpressedGenes(cellchat_IMQ)
cellchat_IMQ <- identifyOverExpressedInteractions(cellchat_IMQ)

# Compute communication probability
cellchat_IMQ <- computeCommunProb(cellchat_IMQ, raw.use = TRUE, population.size = TRUE)
cellchat_IMQ <- computeCommunProbPathway(cellchat_IMQ)
cellchat_IMQ <- aggregateNet(cellchat_IMQ)
cellchat_IMQ <- netAnalysis_computeCentrality(cellchat_IMQ, slot.name = "netP")


# ==============================================================================
# 8. Comparative Analysis: Interaction Numbers
# ==============================================================================
object.list <- list(
  NL  = cellchat_NL,
  IMQ = cellchat_IMQ,
  CE  = cellchat_CE
)

# Calculate maximum weight for consistent scaling
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


# ==============================================================================
# 9. Differential Interaction Analysis
# ==============================================================================
# Merge objects for pairwise and group comparisons
cellchat_1   <- mergeCellChat(list(cellchat_NL, cellchat_IMQ),
                              add.names = c("NL", "IMQ"), cell.prefix = FALSE)

cellchat_2   <- mergeCellChat(list(cellchat_IMQ, cellchat_CE),
                              add.names = c("IMQ", "CE"), cell.prefix = FALSE)

cellchat_all <- mergeCellChat(list(cellchat_NL, cellchat_IMQ, cellchat_CE),
                              add.names = c("NL", "IMQ", "CE"), cell.prefix = FALSE)

# Visualize differential interactions (Strength and Count)
par(mfrow = c(1,4), xpd = TRUE)
netVisual_diffInteraction(cellchat_1, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_1, weight.scale = TRUE, measure = "weight")
netVisual_diffInteraction(cellchat_2, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_2, weight.scale = TRUE, measure = "weight")


# ==============================================================================
# 10. Pathway Specific Analysis: TNF Signaling
# ==============================================================================
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


# ==============================================================================
# 11. Pathway Specific Analysis: TENASCIN Signaling
# ==============================================================================
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


# ==============================================================================
# 12. Visualization: Bubble Plot for TENASCIN Ligand-Receptor Pairs
# ==============================================================================
# Note: Ensure pairLR.TENASCIN is defined, or use signaling="TENASCIN"
p <- netVisual_bubble(
  cellchat_all,
  sources.use = c(1),                 # Define source cluster IDs
  targets.use = c(1,2,3,4,5,6,7),     # Define target cluster IDs
  comparison = c(1,2,3),
  signaling = c("TENASCIN"),          # Plotting all pairs in TENASCIN pathway
  angle.x = 45
)

print(p)