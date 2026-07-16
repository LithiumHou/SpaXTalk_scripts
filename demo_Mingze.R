rm(list = ls())
closeAllConnections()
gc()
library(Seurat)
library(CellChat)
library(dplyr)
library(tibble)
library(ggplot2)
library(future.apply)
library(progressr)
library(Matrix)
library(SparseArray)

handlers("default")
handlers("cli")
handlers(global = TRUE)
options(future.globals.maxSize= 60*1024^3)

# Run this in the command line: PKG_LIBS='-llapack -lblas' R CMD SHLIB fit_spatial_local_effect.c
dyn.load("/home/shared/SpatialCrosstalk/codes/fit_spatial_local_effect.so")
fit_spatial_local_effect_batch <- function(data_mat, pair_i, pair_j, neighbors,
                                           lambda = 1, progress = TRUE, update_every = 50) {
  .Call("fit_spatial_local_effect_batch",
        data_mat, as.integer(pair_i), as.integer(pair_j),
        neighbors, as.numeric(lambda),
        as.logical(progress), as.integer(update_every))
}
# dyn.load("/home/shared/SpatialCrosstalk/codes/pairwise_moran_batch.so")
# pairwise_moran_precompute_cols <- function(neighbors, colnames_data) {
#  .Call("pairwise_moran_precompute_cols", neighbors, as.character(colnames_data))
# }

# pairwise_moran_batch <- function(data_mat, pair_i, pair_j, neighbors, precomp,
#                                  nsim, progress = TRUE, update_every = 50) {
#   .Call("pairwise_moran_batch",
#         data_mat, as.integer(pair_i), as.integer(pair_j),
#         neighbors, precomp, as.integer(nsim),
#         as.logical(progress), as.integer(update_every))
# }
setwd("/home/shared/SpatialCrosstalk")
source("./codes/CCCinference.R")
source("./codes/Visualization.R")
options(future.globals.maxSize= 60*1024^3)
LRDB <- readRDS("./databases/LR/CellChat_human.RDS") %>% as.data.frame()
RTDB <- readRDS("./databases/RT/KEGG_PPR_human.RDS") %>% as.data.frame()
RTDB <- RTDB[RTDB$From %in% LRDB$To, ]


SeuratObj <- readRDS("/home/shared/SpatialCrosstalk/datasets/Mingze2025/Mingze_object.RDS")
exp_data_all <- SeuratObj@assays$RNA$data
meta_all <- SeuratObj@meta.data
meta_all$annotation <- meta_all$celltype

# Select one complete response and one disease progression (PD)
# PD sample
CPID <- 31808
meta <- meta_all[meta_all$CPID == CPID, ]
exp_data <- exp_data_all[, rownames(meta)]
exp_data[1:5,1:5]
meta <- meta[,c("x","y","annotation")]
table(meta$annotation)

coords_new <- meta %>% as.data.frame()
coords_new$x <- (coords_new$x)*0.02
coords_new$y <- (coords_new$y)*0.02
mat <- dist(coords_new) %>% as.matrix()
diag(mat) <- Inf
sort(mat[1,])[1:15]

meta$x <- coords_new$x
meta$y <- coords_new$y

neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = NULL,dist_threshold = 5,decay_type = "exp",normalize = 1,nclusters = 1)
neighbor_list[[1]]
CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[3]]) %>% sum()
GRNs <- InferGRN(exp_data, neighbor_list = neighbor_list,  exp_thres = 0.05, PPR_top = 15, useC = F, RTDB = RTDB, nsim = 0,lambda = 1, nclusters = 64, use_delta_y = T)
allresults <- list(
    # SVG = SVGs,
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/Mingze2025/Mingze_PD_new.RDS")

# CR sample
CPID <- 31810
meta <- meta_all[meta_all$CPID == CPID, ]
exp_data <- exp_data_all[, rownames(meta)]
exp_data[1:5,1:5]
meta <- meta[,c("x","y","annotation")]
table(meta$annotation)

coords_new <- meta %>% as.data.frame()
coords_new$x <- (coords_new$x)*0.02
coords_new$y <- (coords_new$y)*0.02
mat <- dist(coords_new) %>% as.matrix()
diag(mat) <- Inf
sort(mat[1,])[1:15]

meta$x <- coords_new$x
meta$y <- coords_new$y

neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = NULL,dist_threshold = 5,decay_type = "exp",normalize = 1,nclusters = 1)
neighbor_list[[10]]
CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[3]]) %>% sum()
GRNs <- InferGRN(exp_data, neighbor_list = neighbor_list,  exp_thres = 0.05,PPR_top = 15, useC = F, RTDB = RTDB, nsim = 0,lambda = 1, nclusters = 64, use_delta_y = T)
GRNs$beta[1:5,1:5]
allresults <- list(
    # SVG = SVGs,
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/Mingze2025/Mingze_CR_new.RDS")
