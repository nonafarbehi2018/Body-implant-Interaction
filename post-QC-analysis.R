# ============================================================
# Body-implant interaction spatial transcriptomics analysis
#
# Author: Nona Farbehi
# Repository: Body-implant-Interaction
# Purpose:
#   1) Import AnnData objects
#   2) Harmonise cell type labels
#   3) Convert to Seurat objects
#   4) Generate publication-quality UMAP and marker plots
#   5) Run pairwise DE analyses across matched clusters
#   6) Run CellChat communication analysis
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(anndata)
  library(reticulate)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(CellChat)
})

# ---------------------------
# 1) User configuration
# ---------------------------
paths <- list(
  rod_w2_h5ad    = "data/rod_w2_updated_cell_cluster_name.h5ad",
  rod_w4_h5ad    = "data/rod_w4_updated_cell_cluster_name.h5ad",
  sphere_w2_h5ad = "data/sphere_w2_updated_cell_cluster_name.h5ad",
  rdata_main     = "data/Visium_skin_all3.RData",
  rdata_cellchat = "data/Visium_skin_all3_cellChat.RData",
  output_dir     = "output"
)

dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$output_dir, "rds"), showWarnings = FALSE, recursive = TRUE)

analysis_params <- list(
  assay = "RNA",
  meta_col_celltype = "cell_type_updated",
  n_variable_features = 2000,
  npcs = 30,
  min_cells_de = 10,
  de_min_pct = 0.1,
  de_logfc_threshold = 0,
  de_label_logfc = 1,
  cellchat_min_cells = 5,
  cellchat_db_search = "Secreted Signaling"
)

# ---------------------------
# 2) Helper functions
# ---------------------------

assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
}

save_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(paths$output_dir, "figures", paste0(filename, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in"
  )
  ggsave(
    filename = file.path(paths$output_dir, "figures", paste0(filename, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
}

save_base_plot <- function(filename, plot_expr, width = 10, height = 6, res = 300) {
  png(
    filename = file.path(paths$output_dir, "figures", paste0(filename, ".png")),
    width = width,
    height = height,
    units = "in",
    res = res
  )
  plot_expr()
  dev.off()

  pdf(
    file = file.path(paths$output_dir, "figures", paste0(filename, ".pdf")),
    width = width,
    height = height
  )
  plot_expr()
  dev.off()
}

ann_to_seurat <- function(ad_obj, dataset_name, assay = "RNA") {
  X <- py_to_r(ad_obj$X)
  genes <- py_to_r(ad_obj$var_names)
  cells <- py_to_r(ad_obj$obs_names)

  if (inherits(X, "dgRMatrix")) {
    X <- as(X, "dgCMatrix")
  } else if (!inherits(X, "dgCMatrix")) {
    X <- Matrix(as.matrix(X), sparse = TRUE)
  }

  counts <- Matrix::t(X)
  rownames(counts) <- genes
  colnames(counts) <- cells

  meta <- py_to_r(ad_obj$obs) |> as.data.frame()
  rownames(meta) <- cells
  meta$dataset <- dataset_name

  CreateSeuratObject(counts = counts, meta.data = meta, assay = assay)
}

recode_celltypes <- function(seu, mapping, meta_col = analysis_params$meta_col_celltype) {
  stopifnot(meta_col %in% colnames(seu@meta.data))
  seu[[meta_col]][, 1] <- dplyr::recode(seu[[meta_col]][, 1], !!!mapping)
  seu
}

normalize_and_embed <- function(seu, assay = analysis_params$assay, nfeatures = analysis_params$n_variable_features, npcs = analysis_params$npcs) {
  DefaultAssay(seu) <- assay

  if (length(VariableFeatures(seu)) == 0) {
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  }

  if (!"pca" %in% names(seu@reductions)) {
    seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = npcs, verbose = FALSE)
  }

  if (!"umap" %in% names(seu@reductions)) {
    seu <- RunUMAP(seu, dims = seq_len(npcs), reduction = "pca", verbose = FALSE)
  }

  seu
}

make_umap_plot <- function(seu, title, colors, group.by = analysis_params$meta_col_celltype) {
  DimPlot(
    seu,
    reduction = "umap",
    group.by = group.by,
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    scale_color_manual(values = colors, drop = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    ) +
    ggtitle(title)
}

make_marker_dotplot <- function(seu, marker_list, split.by = "dataset", group.by = analysis_params$meta_col_celltype, title = NULL) {
  features <- unique(unlist(marker_list))
  features <- features[features %in% rownames(seu)]

  DotPlot(
    seu,
    features = features,
    group.by = group.by,
    split.by = split.by,
    cols = c("grey90", "steelblue", "firebrick")
  ) +
    RotatedAxis() +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(title %||% "Marker expression across datasets")
}

run_cluster_de <- function(obj, dataset1, dataset2, cluster_id, cluster_col = "cluster") {
  sub_obj <- subset(obj, subset = dataset %in% c(dataset1, dataset2) & .data[[cluster_col]] == cluster_id)

  if (ncol(sub_obj) == 0) return(NULL)
  n_per_dataset <- table(sub_obj$dataset)
  if (!all(c(dataset1, dataset2) %in% names(n_per_dataset))) return(NULL)
  if (any(n_per_dataset[c(dataset1, dataset2)] < analysis_params$min_cells_de)) return(NULL)

  DefaultAssay(sub_obj) <- analysis_params$assay
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  Idents(sub_obj) <- sub_obj$dataset

  de_res <- FindMarkers(
    sub_obj,
    ident.1 = dataset1,
    ident.2 = dataset2,
    logfc.threshold = analysis_params$de_logfc_threshold,
    min.pct = analysis_params$de_min_pct,
    slot = "data",
    verbose = FALSE
  )

  if (nrow(de_res) == 0) return(NULL)

  de_res |>
    rownames_to_column("gene") |>
    mutate(
      cluster = as.character(cluster_id),
      comparison = paste(dataset1, "vs", dataset2),
      neglog10_padj = -log10(p_val_adj + 1e-300),
      highlight = abs(avg_log2FC) >= analysis_params$de_label_logfc & p_val_adj < 0.05
    )
}

make_volcano_plot <- function(de_tbl, title) {
  ggplot(de_tbl, aes(x = avg_log2FC, y = neglog10_padj)) +
    geom_point(aes(color = highlight), alpha = 0.7, size = 1.5) +
    geom_vline(xintercept = c(-analysis_params$de_label_logfc, analysis_params$de_label_logfc), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    ggrepel::geom_text_repel(
      data = subset(de_tbl, highlight),
      aes(label = gene),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
    theme_bw(base_size = 11) +
    labs(
      title = title,
      x = "Average log2 fold-change",
      y = expression(-log[10](adjusted~p)),
      color = NULL
    )
}

run_cellchat_pipeline <- function(seu, group.by = analysis_params$meta_col_celltype, species_db = CellChatDB.mouse) {
  cellchat_db <- subsetDB(species_db, search = analysis_params$cellchat_db_search)
  cc <- createCellChat(object = seu, group.by = group.by)
  cc@DB <- cellchat_db
  cc@idents <- droplevels(cc@idents)

  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = analysis_params$cellchat_min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------
# 3) Colour palette and markers
# ---------------------------
celltype_colors_master <- c(
  "Keratinocytes1" = "#00CC00",
  "Keratinocytes2" = "#FFCCFF",
  "Keratinocytes3" = "#CC33CC",
  "Hair Follicle keratinocytes" = "#CC6600",
  "Adipocytes1" = "#3399FF",
  "Adipocytes2" = "#66FFFF",
  "Muscle1" = "#0033CC",
  "Activated fibroblasts1" = "#800080",
  "Activated fibroblasts2" = "#BF40BF",
  "Activated fibroblasts & macrophages" = "#FF9000",
  "Macrophages1" = "#FF3366",
  "Macrophages2" = "#FFCC00",
  "Immune cells" = "salmon",
  "Fibroblasts1" = "#FF0000"
)

marker_list <- list(
  `Activated fibroblasts2` = c("Fn1", "Lum", "Acta2", "Tagln", "Postn", "Thbs1"),
  `Activated fibroblasts & macrophages` = c("Fn1", "Spp1", "Lgals3", "Il1b", "Tnf", "Timp1"),
  Adipocytes1 = c("Pparg", "Adipoq", "Fabp4", "Plin1", "Lpl", "Cebpa"),
  Adipocytes2 = c("Pparg", "Adipoq", "Fabp4", "Lipe"),
  `Hair Follicle keratinocytes` = c("Krt14", "Krt5", "Krt15", "Krt17", "Sox9", "Lhx2"),
  Keratinocytes1 = c("Krt14", "Krt5", "Krt1", "Krt10"),
  Keratinocytes2 = c("Krt1", "Krt10", "Sprr1a", "Sprr2a", "Lor"),
  `Immune cells` = c("Ptprc", "Cd3e", "Cd3d", "Lyz2", "Ms4a7", "Nkg7"),
  Macrophages1 = c("Adgre1", "Cd68", "Aif1", "C1qa", "C1qb", "Lyz2"),
  Macrophages2 = c("Il1b", "Tnf", "Nos2", "Ccl2", "Spp1", "Lgals3"),
  Muscle1 = c("Acta2", "Tagln", "Myh11", "Cnn1", "Myh6", "Myh7")
)

# ---------------------------
# 4) Load data
# ---------------------------
message("Loading input files...")
purrr::walk(paths[c("rod_w2_h5ad", "rod_w4_h5ad", "sphere_w2_h5ad")], assert_file_exists)

RodW2 <- read_h5ad(paths$rod_w2_h5ad)
RodW4 <- read_h5ad(paths$rod_w4_h5ad)
SphereW2 <- read_h5ad(paths$sphere_w2_h5ad)

# Optional legacy objects
if (file.exists(paths$rdata_main)) load(paths$rdata_main)
if (file.exists(paths$rdata_cellchat)) load(paths$rdata_cellchat)

# ---------------------------
# 5) Convert to Seurat
# ---------------------------
RodW2_seu <- ann_to_seurat(RodW2, "RodW2")
RodW4_seu <- ann_to_seurat(RodW4, "RodW4")
SphereW2_seu <- ann_to_seurat(SphereW2, "SphereW2")

# ---------------------------
# 6) Harmonise labels
# ---------------------------
RodW2_mapping <- c(
  "Inflammatory cells" = "Activated fibroblasts & macrophages",
  "Activated fibroblasts" = "Activated fibroblasts2",
  "Activated fibroblasts2" = "Activated fibroblasts1"
)

SphereW2_mapping <- c(
  "Inflammatory cells" = "Activated fibroblasts & macrophages",
  "Activated fibroblasts & macrophages" = "Macrophages1",
  "Macrophages1" = "Activated fibroblasts2",
  "Keratinocytes3" = "Activated fibroblasts & macrophages"
)

RodW2_seu <- recode_celltypes(RodW2_seu, RodW2_mapping)
SphereW2_seu <- recode_celltypes(SphereW2_seu, SphereW2_mapping)

# ---------------------------
# 7) Normalisation + UMAP
# ---------------------------
RodW2_seu <- normalize_and_embed(RodW2_seu)
RodW4_seu <- normalize_and_embed(RodW4_seu)
SphereW2_seu <- normalize_and_embed(SphereW2_seu)

Idents(RodW2_seu) <- RodW2_seu[[analysis_params$meta_col_celltype]][, 1]
Idents(RodW4_seu) <- RodW4_seu[[analysis_params$meta_col_celltype]][, 1]
Idents(SphereW2_seu) <- SphereW2_seu[[analysis_params$meta_col_celltype]][, 1]

rodw2_colors <- celltype_colors_master[unique(RodW2_seu[[analysis_params$meta_col_celltype]][, 1])]
rodw4_colors <- celltype_colors_master[unique(RodW4_seu[[analysis_params$meta_col_celltype]][, 1])]
spherew2_colors <- celltype_colors_master[unique(SphereW2_seu[[analysis_params$meta_col_celltype]][, 1])]

p_umap_RodW2 <- make_umap_plot(RodW2_seu, "RodW2 – UMAP by cell type", rodw2_colors)
p_umap_RodW4 <- make_umap_plot(RodW4_seu, "RodW4 – UMAP by cell type", rodw4_colors)
p_umap_SphereW2 <- make_umap_plot(SphereW2_seu, "SphereW2 – UMAP by cell type", spherew2_colors)

save_plot(p_umap_RodW2, "RodW2_UMAP_celltypes", width = 7, height = 6)
save_plot(p_umap_RodW4, "RodW4_UMAP_celltypes", width = 7, height = 6)
save_plot(p_umap_SphereW2, "SphereW2_UMAP_celltypes", width = 7, height = 6)

# ---------------------------
# 8) Merge and marker plots
# ---------------------------
skin_combined <- merge(
  x = RodW2_seu,
  y = list(RodW4_seu, SphereW2_seu),
  add.cell.ids = c("RodW2", "RodW4", "SphereW2"),
  project = "Body_implant_interaction"
)
DefaultAssay(skin_combined) <- analysis_params$assay
Idents(skin_combined) <- skin_combined[[analysis_params$meta_col_celltype]][, 1]

p_markers <- make_marker_dotplot(
  skin_combined,
  marker_list = marker_list,
  title = "Skin cell type marker expression across RodW2, RodW4 and SphereW2"
)
save_plot(p_markers, "marker_dotplot_by_dataset", width = 13, height = 7)

vln_genes <- c("Acta2", "Adipoq", "Krt14", "Cd68")
p_vln <- VlnPlot(
  skin_combined,
  features = vln_genes,
  group.by = analysis_params$meta_col_celltype,
  split.by = "dataset",
  pt.size = 0
)
save_plot(p_vln, "marker_violinplots_by_dataset", width = 12, height = 7)

# ---------------------------
# 9) Differential expression
# ---------------------------
combined <- skin_combined
combined$cluster <- if ("seurat_clusters" %in% colnames(combined@meta.data)) {
  as.factor(combined$seurat_clusters)
} else {
  as.factor(combined[[analysis_params$meta_col_celltype]][, 1])
}

comparison_pairs <- list(
  c("RodW2", "RodW4"),
  c("RodW2", "SphereW2"),
  c("RodW4", "SphereW2")
)

all_de <- list()
for (pair in comparison_pairs) {
  d1 <- pair[[1]]
  d2 <- pair[[2]]
  shared_clusters <- intersect(
    unique(as.character(combined$cluster[combined$dataset == d1])),
    unique(as.character(combined$cluster[combined$dataset == d2]))
  )

  for (cl in shared_clusters) {
    de_tbl <- run_cluster_de(combined, d1, d2, cl)
    if (is.null(de_tbl)) next

    tag <- paste0(d1, "_vs_", d2, "_cluster_", cl)
    all_de[[tag]] <- de_tbl
    write.csv(de_tbl, file.path(paths$output_dir, "tables", paste0("DE_", tag, ".csv")), row.names = FALSE)

    p_de <- make_volcano_plot(de_tbl, paste0("DE: ", d1, " vs ", d2, " – cluster ", cl))
    save_plot(p_de, paste0("volcano_", tag), width = 8, height = 6)
  }
}

if (length(all_de) > 0) {
  bind_rows(all_de) |>
    write.csv(file.path(paths$output_dir, "tables", "DE_all_comparisons_combined.csv"), row.names = FALSE)
}

# ---------------------------
# 10) CellChat
# ---------------------------
DefaultAssay(RodW2_seu) <- analysis_params$assay
DefaultAssay(SphereW2_seu) <- analysis_params$assay

cellchat_RodW2 <- run_cellchat_pipeline(RodW2_seu)
cellchat_SphereW2 <- run_cellchat_pipeline(SphereW2_seu)

cellchat_merged <- mergeCellChat(
  list(RodW2 = cellchat_RodW2, SphereW2 = cellchat_SphereW2),
  add.names = c("RodW2", "SphereW2")
)

p_count <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
p_weight <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
save_plot(p_count, "cellchat_compare_interactions_count", width = 6, height = 5)
save_plot(p_weight, "cellchat_compare_interactions_weight", width = 6, height = 5)

save_base_plot(
  filename = "cellchat_circle_interaction_count",
  width = 12,
  height = 6,
  plot_expr = function() {
    par(mfrow = c(1, 2), xpd = TRUE)
    netVisual_circle(
      cellchat_RodW2@net$count,
      vertex.weight = as.numeric(table(cellchat_RodW2@idents)),
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = "RodW2 – Number of interactions"
    )
    netVisual_circle(
      cellchat_SphereW2@net$count,
      vertex.weight = as.numeric(table(cellchat_SphereW2@idents)),
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = "SphereW2 – Number of interactions"
    )
  }
)

saveRDS(RodW2_seu, file.path(paths$output_dir, "rds", "RodW2_seu_processed.rds"))
saveRDS(RodW4_seu, file.path(paths$output_dir, "rds", "RodW4_seu_processed.rds"))
saveRDS(SphereW2_seu, file.path(paths$output_dir, "rds", "SphereW2_seu_processed.rds"))
saveRDS(cellchat_RodW2, file.path(paths$output_dir, "rds", "cellchat_RodW2.rds"))
saveRDS(cellchat_SphereW2, file.path(paths$output_dir, "rds", "cellchat_SphereW2.rds"))

message("Analysis complete. Outputs saved to: ", normalizePath(paths$output_dir))
