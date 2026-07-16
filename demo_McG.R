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


handlers("default")
handlers("cli")
handlers(global = TRUE)

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

SeuratObj <- readRDS("./datasets/McGregor2025/McGregor_object.RDS")
clusters <- table(SeuratObj$CellAnnotation.Level0) %>% as.data.frame()
cluster_n <- clusters[clusters$Freq > 100, "Var1"] %>% as.character()
SeuratObj$annotation <- SeuratObj$CellAnnotation.Level0

# SeuratObj <- subset(SeuratObj, annotation %in% cluster_n)
receiver_cluster = "Endothelium"
Prepare_data(SeuratObj, receiver_cluster = receiver_cluster)
exp_data[1:15,1:15]
LRDB <- readRDS("./databases/LR/CellChat_human.RDS") %>% as.data.frame()
RTDB <- readRDS("./databases/RT/KEGG_PPR_human.RDS") %>% as.data.frame()
RTDB <- RTDB[RTDB$From %in% LRDB$To, ]

# mat <- dist(meta[,1:2]) %>% as.matrix()
# diag(mat) <- Inf
coords_new <- meta %>% as.data.frame()
coords_new$x <- (coords_new$x - min(coords_new$x))/10
coords_new$y <- (coords_new$y - min(coords_new$y))/10
# mat <- dist(coords_new) %>% as.matrix()
# diag(mat) <- Inf
# min(mat)
# sort(mat[1,])[1:5]

meta$x <- coords_new$x
meta$y <- coords_new$y

# SVGs <- InferSVG(exp_data, locs = meta[,1:2], receiver = receiver)
neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receiver,dist_threshold = 2.5,decay_type = "exp",normalize = T)
neighbor_list[[1]]

CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[1]]) %>% sum()
gc()
GRNs <- InferGRN(exp_data, neighbor_list = neighbor_list,  exp_thres = 0.05, PPR_top = 99,useC = F, RTDB = RTDB, nsim = 0,lambda = 1, nclusters = 32, use_delta_y = T)
GRNs$beta[1:5,1:10]
allresults <- list(
    # SVG = SVGs,
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/McGregor2025/McGregor_result_new.RDS")
