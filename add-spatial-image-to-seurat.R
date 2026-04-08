library(ggplot2)
library(Seurat)
library(SeuratData)
library(CellChat)
library(patchwork)
library(dplyr)
library(readr)
library(tiff)
setwd("~/Desktop/Implants-2025-Feb/new_new_annotation")
load("~/Desktop/Implants-2025-Feb/new_new_annotation/Visium_skin_all3.RData")

############RodW2############ 
# read png image 
Idents(RodW2_seu) <- RodW2_seu$cell_type_updated

rod_w2_old <- read_rds("rod_w2.rds")

RodW2_seu_spatial <- CreateSeuratObject(GetAssayData(RodW2_seu,assay = "RNA",layer = "data"), project = "SeuratProject", assay = "Spatial",
                    meta.data = RodW2_seu[[]])



RodW2_seu_spatial@images$slice1 <- rod_w2_old@images$slice1

RodW2_seu_spatialplot <- SpatialPlot(RodW2_seu_spatial, group.by = "cell_type_updated", pt.size.factor = 4, cols = celltype_colors_master)
ggsave("RodW2_seu_spatialplot.pdf", RodW2_seu_spatialplot, width = 8, height = 6, dpi = 300)

###################RodW4#####################
Idents(RodW4_seu) <- RodW4_seu$cell_type_updated

rod_w4_old <- read_rds("rod_w4.rds")

RodW4_seu_spatial <- CreateSeuratObject(GetAssayData(RodW4_seu,assay = "RNA",layer = "data"), project = "SeuratProject", assay = "Spatial",
                                        meta.data = RodW4_seu[[]])

cell_colors_RodW4 <- c(
  "Macrophages1"= "#FF3366","Keratinocytes1"= "#00CC00","Activated fibroblasts2"="#FF9000",
  "Macrophages2"= "#FFCC00","Adipocytes1"="#3399FF","Muscle1"="#0033CC",
  "Adipocytes2"= "#66FFFF","Activated fibroblasts"="#800080","Immune cells"= "salmon",
  "Hair Follicle keratinocytes"= "#CC6600","Keratinocytes2"= "#FFCCFF"
)

RodW4_seu_spatial@images$slice1 <- rod_w4_old@images$slice1

RodW4_seu_spatialplot <- SpatialPlot(RodW4_seu_spatial, group.by = "cell_type_updated", pt.size.factor = 3, cols = celltype_colors_master)
ggsave("RodW4_seu_spatialplot.pdf", RodW4_seu_spatialplot, width = 8, height = 6, dpi = 300)

###################SphereW2#####################
Idents(SphereW2_seu) <- SphereW2_seu$cell_type_updated

sphere_w2_old <- read_rds("sphere_w2.rds")

SphereW2_seu_spatial <- CreateSeuratObject(GetAssayData(SphereW2_seu,assay = "RNA",layer = "data"), project = "SeuratProject", assay = "Spatial",
                                        meta.data = SphereW2_seu[[]])

cell_colors_SphereW2 <- c(
  "Keratinocytes1" = "#00CC00","Adipocytes1"="#3399FF","Activated fibroblasts2"="#FF9000",
  "Muscle1"="#0033CC","Macrophages1"= "#FF3366","Activated fibroblasts"= "#800080",
  "Keratinocytes2"="#FFCCFF","Adipocytes2"= "#66FFFF","Macrophages2"= "#FFCC00",
  "Hair Follicle keratinocytes"= "#CC6600"
)

SphereW2_seu_spatial@images$slice1 <- sphere_w2_old@images$slice1

SphereW2_seu_spatialplot <- SpatialPlot(SphereW2_seu_spatial, group.by = "cell_type_updated", pt.size.factor = 6, cols = celltype_colors_master)
ggsave("SphereW2_seu_spatialplot.pdf", SphereW2_seu_spatialplot, width = 8, height = 6, dpi = 300)

######### update colors and labels #####
library(anndata)
load("Visium_data_with_spatial.RData")
#SphereW2_seu_spatialplot <- SpatialPlot(SphereW2_seu_spatial, group.by = "cell_type_updated", pt.size.factor = 6, cols = celltype_colors_master)
RodW2_seu_spatial$cell_type_updated

cell_colors_SphereW2 <- c(
  "Keratinocytes1"                = "#00CC00",
  "Keratinocytes2"                = "#FFCCFF",
  "Adipocytes1"                   = "#3399FF",
  "Adipocytes2"                   = "#66FFFF",
  "Muscle1"                       = "#0033CC",
  "Macrophages1"                  = "#FF3366",
  "Macrophages2"                  = "#FFCC00",
  "Activated fibroblasts1"        = "#800080",
  "Activated fibroblasts2"        = "#BF40BF",
  "Activated fibroblasts3"        = "#D98CD9",   # NEW SUGGESTED COLOR
  "Hair Follicle keratinocytes"   = "#CC6600"
)


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

cell_colors_RodW4 <- c("Macrophages1"= "#FF3366",
                       "Muscle1"="#0033CC",
                       "Adipocytes1"="#3399FF",
                       "Macrophages2"= "#FFCC00",
                       "Keratinocytes2"="#FFCCFF",
                       "Adipocytes2"= "#66FFFF",
                       "Keratinocytes1"= "#00CC00",
                       "Fibroblasts1"= "#FF0000",
                       "Hair Follicle keratinocytes"="#CC6600",
                       "Activated fibroblasts"="#800080")


levels(RodW2_seu_spatial$cell_type_updated)[levels(RodW2_seu_spatial$cell_type_updated) =="Keratinocytes2"] <- "Activated fibroblasts2"
levels(SphereW2_seu_spatial$cell_type_updated)[levels(RodW2_seu_spatial$cell_type_updated)=="Activated fibroblasts & macrophages"] <- "Activated fibroblasts3"

levels(RodW2_seu$cell_type_updated)[levels(RodW2_seu$cell_type_updated) =="Keratinocytes2"] <- "Activated fibroblasts2"
levels(SphereW2_seu$cell_type_updated)[levels(SphereW2_seu$cell_type_updated)=="Activated fibroblasts & macrophages"] <- "Activated fibroblasts3"

# check the latest label 
RodW2_seu_spatialplot <- SpatialPlot(RodW2_seu_spatial, group.by = "cell_type_updated", pt.size.factor = 4, cols = cell_colors_RodW2)
SphereW2_seu_spatialplot <- SpatialPlot(SphereW2_seu_spatial, group.by = "cell_type_updated", pt.size.factor = 6, cols = cell_colors_SphereW2)

save.image(file = "Visium_data_with_spatial_new.RData")
