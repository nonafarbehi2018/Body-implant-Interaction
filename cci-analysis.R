library(ggplot2)
library(Seurat)
library(SeuratData)
library(CellChat)
library(patchwork)
library(dplyr)
library(readr)
library(tiff)
load("~/Desktop/Implants-2025-Feb/new_new_annotation/Visium_data_with_spatial_new.RData")
setwd("/media/bml/USER_DATA/Nona/implant_data_new/new_data_celltypes/new_new_annotation/")

############### Rod W2 ##################
# Create scRNA-seq cellchat object
RodW2_seu <- JoinLayers(RodW2_seu)

cellchat_rod_w2_new <- createCellChat(RodW2_seu, group.by = "cell_type_updated")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat_rod_w2_new@DB <- CellChatDB.use
cellchat_rod_w2_new@meta$samples <- "rod_w2"

# Identify over expressed genes and interactions
cellchat_rod_w2_new <- subsetData(cellchat_rod_w2_new)
cellchat_rod_w2_new <- identifyOverExpressedGenes(cellchat_rod_w2_new)
cellchat_rod_w2_new <- identifyOverExpressedInteractions(cellchat_rod_w2_new)
cellchat_rod_w2_new <- smoothData(cellchat_rod_w2_new, adj = PPI.mouse)

# Compute the communication probability and infer cellular communication network 
future::plan("multisession", workers = 4) 
min_percentage <- 0.1
rod_wk2.new.groupsize <- as.numeric(table(cellchat_rod_w2_new@idents))
uwprod_wk2 <- floor(min_percentage / 100.0 * sum(rod_wk2.new.groupsize))
cellchat_rod_w2_new <- computeCommunProb(cellchat_rod_w2_new, raw.use = FALSE, population.size = FALSE)
cellchat_rod_w2_new <- filterCommunication(cellchat_rod_w2_new, min.cells = uwprod_wk2)
cellchat_rod_w2_new <- computeCommunProbPathway(cellchat_rod_w2_new)
cellchat_rod_w2_new <- aggregateNet(cellchat_rod_w2_new)
cellchat_rod_w2_new.net <- subsetCommunication(cellchat_rod_w2_new)
cellchat_rod_w2_new_pathway <- cellchat_rod_w2_new@netP$pathways

library(svglite)

# Plot overall cci
svglite(filename = "rod_w2_recluster_all_genes.svg", width=12.5, height=8.5)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_circle(cellchat_rod_w2_new@net$count, vertex.weight = rod_wk2.new.groupsize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = cell_colors_RodW2)
gg2 <- netVisual_circle(cellchat_rod_w2_new@net$weight, vertex.weight = rod_wk2.new.groupsize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = cell_colors_RodW2)
dev.off()

# River plot
library(ggalluvial)
library(NMF)
selectK(cellchat_rod_w2_new, pattern = "outgoing")
nPatterns = 7
cellchat_rod_w2_new <- identifyCommunicationPatterns(cellchat_rod_w2_new, pattern = "outgoing", k = nPatterns)
gg4 <- netAnalysis_river(cellchat_rod_w2_new, pattern = "outgoing", color.use = cell_colors_RodW2, color.use.pattern = colors_pattern, cutoff = 0.3)
ggsave(filename= "rodw2_river_7patterns1.svg", plot=gg4, width = 5.5, height = 3.5, units = 'in', dpi = 300)

cell_colors_RodW2 <- c(
  "Keratinocytes1"              = "#00CC00",
  "Activated fibroblasts1"      = "#800080",
  "Activated fibroblasts2"        = "#BF40BF",
  "Macrophages1"                = "#FF3366",
  "Macrophages2"                = "#FFCC00",
  "Immune cells"                = "salmon",
  "Adipocytes1"                 = "#3399FF",
  "Adipocytes2"                 = "#66FFFF",
  "Muscle1"                     = "#0033CC",
  "Hair Follicle keratinocytes" = "#CC6600",
  "Activated fibroblasts & macrophages"="#FF9000"
)


colors_pattern <- c(
  "Pattern 1" = "darksalmon",
  "Pattern 2" = "yellowgreen",
  "Pattern 3" = "green3",
  "Pattern 4" = "deepskyblue",
  "Pattern 5" = "hotpink",
  "Pattern 6" = "gold",
  "Pattern 7" = "mediumpurple"
)

cell_colors_RodW2_list <- sort(names(cell_colors_RodW2))
cell_colors_RodW2 <- cell_colors_RodW2[cell_colors_RodW2_list]


# plot specific pathways 
pathways.show <- c("RANKL","IGF", "PDGF", "PLAU", "GRN", "IGFBP")  # "RANKL", "IGF", "PDGF", "PLAU","PDGF", "YAP"
for (i in 1:length(pathways.show)) {
  print(pathways.show[i])
  # Circle plot
  svglite(filename = paste0(pathways.show[i], "_circle_rod_w2_recluster.svg"), width=8, height=6)
  par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
  netVisual_aggregate(cellchat_rod_w2_new, signaling = pathways.show[i], layout = "circle", color.use = cell_colors_RodW2)
  dev.off()
  
  gg <- plotGeneExpression(cellchat_rod_w2_new, signaling = pathways.show[i], enriched.only = TRUE, type = "violin", color.use = cell_colors_RodW2)
  gg
  ggsave(filename=paste0(pathways.show[i], "_gene_expression_rod_w2_recluster.svg"), plot=gg, width = 5.5, height = 8.5, units = 'in', dpi = 300)
}

# plot signaling roles
cellchat_rod_w2_new <- netAnalysis_computeCentrality(cellchat_rod_w2_new, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_network(cellchat_rod_w2_new, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_rod_w2_new, color.use = cell_colors_RodW2)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat_rod_w2_new, signaling = c("IGF", "PDGF", "PLAU"))
#> Signaling role analysis on the cell-cell communication network from user's input
library(ggpubr)
gg3 <- ggarrange(gg1,gg2, ncol = 2, nrow = 1)
ggsave(filename="rod_w2_signaling_roles.svg", plot=gg1, width = 5.5, height = 5.5, units = 'in', dpi = 300)


####################################RodW4######################################################
cellchat_rod_w4_new <- createCellChat(RodW4_seu, group.by = "cell_type_updated")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

cellchat_rod_w4_new@DB <- CellChatDB.use
cellchat_rod_w4_new@meta$samples <- "rod_w4"

cellchat_rod_w4_new <- subsetData(cellchat_rod_w4_new)
cellchat_rod_w4_new <- identifyOverExpressedGenes(cellchat_rod_w4_new)
cellchat_rod_w4_new <- identifyOverExpressedInteractions(cellchat_rod_w4_new)
cellchat_rod_w4_new <- smoothData(cellchat_rod_w4_new, adj = PPI.mouse)

# Compute the communication probability and infer cellular communication network # 
future::plan("multisession", workers = 4) 

min_percentage <- 0.1
rod_w4.new.groupsize <- as.numeric(table(cellchat_rod_w4_new@idents))
uwpsphere_rod_w4 <- floor(min_percentage / 100.0 * sum(rod_w4.new.groupsize))
cellchat_rod_w4_new <- computeCommunProb(cellchat_rod_w4_new, raw.use = FALSE, population.size = FALSE)
cellchat_rod_w4_new <- filterCommunication(cellchat_rod_w4_new, min.cells = uwpsphere_rod_w4)
cellchat_rod_w4_new <- computeCommunProbPathway(cellchat_rod_w4_new)
cellchat_rod_w4_new <- aggregateNet(cellchat_rod_w4_new)
cellchat_rod_w4_new.net <- subsetCommunication(cellchat_rod_w4_new)
cellchat_rod_w4_new_pathway <- cellchat_rod_w4_new@netP$pathways

library(svglite)

svglite(filename = "rod_w4_recluster_all_genes.svg", width=12.5, height=8.5)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_circle(cellchat_rod_w4_new@net$count, vertex.weight = rod_w4.new.groupsize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = cell_colors_RodW4)
gg2 <- netVisual_circle(cellchat_rod_w4_new@net$weight, vertex.weight = rod_w4.new.groupsize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = cell_colors_RodW4)
dev.off()

library(ggalluvial)
library(NMF)
selectK(cellchat_rod_w4_new, pattern = "outgoing")
nPatterns = 8
cellchat_rod_w4_new <- identifyCommunicationPatterns(cellchat_rod_w4_new, pattern = "outgoing", k = nPatterns)
gg4 <- netAnalysis_river(cellchat_rod_w4_new, pattern = "outgoing", color.use = cell_colors_RodW4, color.use.pattern = colors_pattern, cutoff = 0.3)
ggsave(filename= "RodW4_river_8patterns1.svg", plot=gg4, width = 8, height = 6.5, units = 'in', dpi = 300)

# Example: Defining colors for cell groups
celltype_colors_RodW4 <- c(
  "Keratinocytes1"                = "#00CC00",
  "Keratinocytes2"                = "#FFCCFF",
  "Hair Follicle keratinocytes"   = "#CC6600",
  "Adipocytes1"                   = "#3399FF",
  "Adipocytes2"                   = "#66FFFF",
  "Muscle1"                       = "#0033CC",
  "Activated fibroblasts1"        = "#800080",
  "Macrophages1"                  = "#FF3366",
  "Macrophages2"                  = "#FFCC00",
  "Fibroblasts"                   = "#FF0000"
)

colors_pattern <- c(
  "Pattern 1" = "darksalmon",
  "Pattern 2" = "yellowgreen",
  "Pattern 3" = "green3",
  "Pattern 4" = "deepskyblue",
  "Pattern 5" = "hotpink",
  "Pattern 6" = "gold",
  "Pattern 7" = "mediumpurple",
  "Pattern 8" = "lightblue"
)

cell_colors_RodW4_list <- order(names(celltype_colors_RodW4))
cell_colors_RodW4 <- celltype_colors_RodW4[cell_colors_RodW4_list]





pathways.show <- c("IGF", "PDGF", "PLAU","PDGF", "TGFb")  # "RANKL", "IGF", "PDGF", "PLAU","PDGF", "YAP"
for (i in 1:length(pathways.show)) {
  print(pathways.show[i])
  # Circle plot
  svglite(filename = paste0(pathways.show[i], "_circle_rod_w4_recluster.svg"), width=8, height=6)
  par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
  netVisual_aggregate(cellchat_rod_w4_new, signaling = pathways.show[i], layout = "circle", color.use = cell_colors_RodW4)
  dev.off()
  
  gg <- plotGeneExpression(cellchat_rod_w4_new, signaling = pathways.show[i], enriched.only = TRUE, type = "violin", color.use = cell_colors_RodW4)
  gg
  ggsave(filename=paste0(pathways.show[i], "_gene_expression_rod_w4_recluster.svg"), plot=gg, width = 5.5, height = 8.5, units = 'in', dpi = 300)
}

cellchat_rod_w4_new <- netAnalysis_computeCentrality(cellchat_rod_w4_new, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_rod_w4_new, color.use = cell_colors_RodW4)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
#> Signaling role analysis on the cell-cell communication network from user's input
library(ggpubr)
ggsave(filename="rod_w4_signaling_roles.svg", plot=gg1, width = 5.5, height = 5.5, units = 'in', dpi = 300)

#########################################################SphereW2############################################################
cellchat_sphere_w2_new <- createCellChat(SphereW2_seu, group.by = "cell_type_updated")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

cellchat_sphere_w2_new@DB <- CellChatDB.use
cellchat_sphere_w2_new@meta$samples <- "sphere_w2"

# Preprocessing the expression data for cell-cell communication analysis #
cellchat_sphere_w2_new <- subsetData(cellchat_sphere_w2_new)
cellchat_sphere_w2_new <- identifyOverExpressedGenes(cellchat_sphere_w2_new)
cellchat_sphere_w2_new <- identifyOverExpressedInteractions(cellchat_sphere_w2_new)
cellchat_sphere_w2_new <- smoothData(cellchat_sphere_w2_new, adj = PPI.mouse)

# Compute the communication probability and infer cellular communication network # 
future::plan("multisession", workers = 4) 

min_percentage <- 0.1
sphere_wk2.new.groupsize <- as.numeric(table(cellchat_sphere_w2_new@idents))
uwprod_wk2 <- floor(min_percentage / 100.0 * sum(sphere_wk2.new.groupsize))
cellchat_sphere_w2_new <- computeCommunProb(cellchat_sphere_w2_new, raw.use = FALSE, population.size = FALSE)
cellchat_sphere_w2_new <- filterCommunication(cellchat_sphere_w2_new, min.cells = uwprod_wk2)
cellchat_sphere_w2_new <- computeCommunProbPathway(cellchat_sphere_w2_new)
cellchat_sphere_w2_new <- aggregateNet(cellchat_sphere_w2_new)
cellchat_sphere_w2_new.net <- subsetCommunication(cellchat_sphere_w2_new)
cellchat_sphere_w2_new_pathway <- cellchat_sphere_w2_new@netP$pathways

library(svglite)

svglite(filename = "sphere_w2_recluster_all_genes.svg", width=12.5, height=8.5)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_circle(cellchat_sphere_w2_new@net$count, vertex.weight = sphere_wk2.new.groupsize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = cell_colors_SphereW2)
gg2 <- netVisual_circle(cellchat_sphere_w2_new@net$weight, vertex.weight = sphere_wk2.new.groupsize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = cell_colors_SphereW2)
dev.off()

library(ggalluvial)
library(NMF)
selectK(cellchat_sphere_w2_new, pattern = "outgoing")
nPatterns = 7
cellchat_sphere_w2_new <- identifyCommunicationPatterns(cellchat_sphere_w2_new, pattern = "outgoing", k = nPatterns)
gg4 <- netAnalysis_river(cellchat_sphere_w2_new, pattern = "outgoing", color.use = cell_colors_SphereW2, color.use.pattern = colors_pattern, cutoff=0.3)
ggsave(filename= "sphw2_river_7patterns.svg", plot=gg4, width = 5.5, height = 3.5, units = 'in', dpi = 300)

colors_pattern <- c(
  "Pattern 1" = "darksalmon",
  "Pattern 2" = "yellowgreen",
  "Pattern 3" = "green3",
  "Pattern 4" = "deepskyblue",
  "Pattern 5" = "hotpink",
  "Pattern 6" = "gold",
  "Pattern 7" = "mediumpurple"
)

cell_colors_SphereW2 <- c(
  "Keratinocytes1"                = "#00CC00",
  "Keratinocytes2"                = "#FFCCFF",
  "Adipocytes1"                   = "#3399FF",
  "Muscle1"                       = "#0033CC",
  "Macrophages1"                  = "#FF3366",
  "Macrophages2"                  = "#FFCC00",
  "Activated fibroblasts1"        = "#800080",
  "Activated fibroblasts2"        = "#BF40BF",
  "Activated fibroblasts3"        = "#D98CD9"
)


cell_colors_SphereW2_list <- sort(names(cell_colors_SphereW2))
cell_colors_SphereW2 <- cell_colors_SphereW2[cell_colors_SphereW2_list]

pathways.show <- c("RANKL", "IGF", "PDGF", "PLAU","PDGF","TGFb", "IGFBP") 
for (i in 1:length(pathways.show)) {
  print(pathways.show[i])
  # Circle plot
  svglite(filename = paste0(pathways.show[i], "_circle_sphere_w2_recluster.svg"), width=8, height=6)
  par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
  netVisual_aggregate(cellchat_sphere_w2_new, signaling = pathways.show[i], layout = "circle", color.use = cell_colors_SphereW2)
  dev.off()
  
  gg <- plotGeneExpression(cellchat_sphere_w2_new, signaling = pathways.show[i], enriched.only = TRUE, type = "violin", color.use = cell_colors_SphereW2)
  
  ggsave(filename=paste0(pathways.show[i], "_gene_expression_sphere_w2_recluster.svg"), plot=gg, width = 5.5, height = 8.5, units = 'in', dpi = 300)
}

cellchat_sphere_w2_new <- netAnalysis_computeCentrality(cellchat_sphere_w2_new, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_sphere_w2_new, color.use = cell_colors_SphereW2)

library(ggpubr)
ggsave(filename="sphere_w2_signaling_roles.svg", plot=gg1, width = 5.5, height = 5.5, units = 'in', dpi = 300)
save(cellchat_sphere_w2_new, cellchat_rod_w2_new,　file = "new_annotations_cellchat_analysis.rda")

