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
LRDB <- readRDS("./databases/LR/CellChat_mouse.RDS") %>% as.data.frame()
RTDB <- readRDS("./databases/RT/KEGG_PPR_mouse.RDS") %>% as.data.frame()
RTDB <- RTDB[RTDB$From %in% LRDB$To, ]

# E9.5 data
adata <- read_h5ad("./datasets/Chen2022/E9.5_MOSTA.h5ad")
exp_data <- t(adata$X)
meta <- adata$obs %>% as.data.frame()
coords <- adata$obsm[["spatial"]] %>% as.data.frame()
meta <- data.frame(x = coords$V1,
                   y = coords$V2,
                   annotation = meta$annotation %>% as.character(),
                   stringsAsFactors = FALSE)
rownames(meta) <- colnames(exp_data)
table(meta$annotation)

receiver <- rownames(meta[meta$annotation %in% c("Neural crest"), ])

neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receiver,dist_threshold = 3,decay_type = "exp",normalize = T)
neighbor_list[[1]]

CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[1]]) %>% sum()
# neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receiver,dist_threshold = 2,decay_type = "exp",normalize = T)
GRNs <- InferGRN(exps = exp_data, neighbor_list = neighbor_list, PPR_top = 15,exp_thres = 0.05,useC = F, RTDB = RTDB, nsim = 0,lambda = 1, nclusters = 64, use_delta_y = T)
betas <- GRNs$beta
# GRNs$moranI[1:5,1:5]
allresults <- list(
    # SVG = SVGs,
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/Chen2022/E95_new.RDS")

# SCENIC analysis
SeuratObj <- CreateSeuratObject(exp_data, meta.data = meta)
SeuratObj <- subset(SeuratObj, annotation == "Neural crest")
work_dir <- "./codes/scenic/"
tf_list <- unique(RTDB$To) %>% sort()

# E11.5 data
gc()
adata <- read_h5ad("./datasets/Chen2022/E11.5_MOSTA.h5ad")
exp_data <- t(adata$X)
meta <- adata$obs %>% as.data.frame()
coords <- adata$obsm[["spatial"]] %>% as.data.frame()
meta <- data.frame(x = coords$V1,
                   y = coords$V2,
                   annotation = meta$annotation %>% as.character(),
                   stringsAsFactors = FALSE)
rownames(meta) <- colnames(exp_data)
table(meta$annotation)
receivers <- rownames(meta[meta$annotation %in% c("Heart"), ])

mat <- dist(meta[,1:2]) %>% as.matrix()
neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receivers,dist_threshold = 2,decay_type = "insq",normalize = F,nclusters = 1)
neighbor_list[[1]]
CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[3]]>0) %>% sum()
GRNs <- InferGRN(exp_data, neighbor_list = neighbor_list, RTDB = RTDB, nsim = 201)
GRNs$moranI[1:5,1:5]

allresults <- list(
    # SVG = SVGs,
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/Chen2022/E115_new.RDS")

# E13.5 data
adata <- read_h5ad("./datasets/Chen2022/E13.5_MOSTA.h5ad")
exp_data <- t(adata$X)
meta <- adata$obs %>% as.data.frame()
coords <- adata$obsm[["spatial"]] %>% as.data.frame()
meta <- data.frame(x = coords$V1,
                   y = coords$V2,
                   annotation = meta$annotation %>% as.character(),
                   stringsAsFactors = FALSE)
rownames(meta) <- colnames(exp_data)
table(meta$annotation)
receivers <- rownames(meta[meta$annotation %in% c("Heart"), ])

mat <- dist(meta[,1:2]) %>% as.matrix()
neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receivers,dist_threshold = 2,decay_type = "insq",normalize = F,nclusters = 1)
neighbor_list[[1]]
CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[3]]>0) %>% sum()
GRNs <- InferGRN(exp_data, neighbor_list = neighbor_list, RTDB = RTDB, nsim = 201)
GRNs$moranI[1:5,1:5]

allresults <- list(
    # SVG = SVGs,
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/Chen2022/E135_new.RDS")
