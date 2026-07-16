rm(list = ls())
closeAllConnections()
gc()
library(anndata)
library(Seurat)
library(CellChat)
library(dplyr)
library(tibble)
library(ggplot2)
library(future.apply)
library(progressr)
library(Matrix)
library(SparseArray)

setwd("/home/shared/SpatialCrosstalk")
source("./codes/CCCinference.R")
source("./codes/Visualization.R")
options(future.globals.maxSize= 60*1024^3)

# data:https://www.nature.com/articles/s41588-025-02193-3
data.dir <- "/home/shared/SpatialCrosstalk/datasets/VisiumHD/"
deconv <- read.csv("/home/shared/SpatialCrosstalk/datasets/VisiumHD/DeconvolutionResults_P1CRC.csv")
rownames(deconv) <- deconv$barcode

SeuratObj <- Load10X_Spatial(
  data.dir = data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE
)
meta_all <- SeuratObj@meta.data
meta_all$celltype_supervised <- deconv[rownames(meta_all), "DeconvolutionLabel2"]
table(meta_all$celltype_supervised)

# Select the cells
# Load the coordinates
coords <- SeuratObj@images[[1]]@boundaries$centroids@coords %>% as.data.frame()
rownames(coords) <- colnames(SeuratObj)
target_cell <- "s_008um_00556_00532-1"
target_cell %in% colnames(SeuratObj)

thres <- 800
xc <- coords[target_cell, ] %>% as.numeric()
xc_u <- xc[1]+thres; xc_l <- xc[1]-thres
yc_u <- xc[2]+thres; yc_l <- xc[2]-thres

coords_used <- coords[coords$x <= xc_u & coords$x >= xc_l, ]
coords_used <- coords_used[coords_used$y <= yc_u & coords_used$y >= yc_l, ]
cell_used <- rownames(coords_used)
meta <- coords_used
meta$annotation <- deconv[cell_used, "DeconvolutionLabel2"]
table(meta$annotation)

# Determine the receivers
thres <- 750
xc <- coords[target_cell, ] %>% as.numeric()
xc_u <- xc[1]+thres; xc_l <- xc[1]-thres
yc_u <- xc[2]+thres; yc_l <- xc[2]-thres

coords_used <- coords[coords$x <= xc_u & coords$x >= xc_l, ]
coords_used <- coords_used[coords_used$y <= yc_u & coords_used$y >= yc_l, ]
receiver <- rownames(coords_used)
table(deconv[receiver, "DeconvolutionLabel2"])

# SpaXTalk workflow
exp_data <- SeuratObj@assays$Spatial$counts
exp_data <- exp_data[, cell_used]
exp_data <- NormalizeData(exp_data)

saveRDS(exp_data, file = "./datasets/VisiumHD/expression_used.RDS")
write.table(meta, file = "./datasets/VisiumHD/meta_used.csv", quote = F,sep = ",")

LRDB <- readRDS("./databases/LR/CellChat_human.RDS") %>% as.data.frame()
RTDB <- readRDS("./databases/RT/KEGG_PPR_human.RDS") %>% as.data.frame()
RTDB <- RTDB[RTDB$From %in% LRDB$To, ]
RTDB <- RTDB %>% 
    group_by(To) %>%
    arrange(PPR) %>%
    slice_head(n = 10) %>%
    ungroup()

coords_new <- meta %>% as.data.frame()
coords_new$x <- (coords_new$x - min(coords_new$x))/30
coords_new$y <- (coords_new$y - min(coords_new$y))/30
mat <- dist(coords_new) %>% as.matrix()
# diag(mat) <- Inf
# min(mat)
sort(mat[21,])[1:15]

meta$x <- coords_new$x
meta$y <- coords_new$y

neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receiver,dist_threshold = 3,decay_type = "exp",normalize = T)
neighbor_list[[1]]
CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[6]]) %>% sum()
GRNs <- InferGRN(exps = exp_data, neighbor_list = neighbor_list, exp_thres = 0.05,PPR_top = 99,useC = F, RTDB = RTDB, nsim = 0,lambda = 1, nclusters = 64, use_delta_y = T)
dim(GRNs$beta)

allresults <- list(
    # SVG = SVGs,
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/VHD/VHD.RDS")
