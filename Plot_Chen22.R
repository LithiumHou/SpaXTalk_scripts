rm(list = ls())
closeAllConnections()
gc()
library(Seurat)
library(CellChat)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)
library(ggnewscale)
library(ggforce)
library(ggpubr)
library(scatterpie)
library(viridis)
library(Matrix)
library(ggraph)
library(paletteer)

setwd("/home/shared/SpatialCrosstalk")
source("./codes/Process_results.R")
source("./codes/Visualization.R")

adata <- anndata::read_h5ad("./datasets/Chen2022/E9.5_MOSTA.h5ad")
meta <- adata$obs %>% as.data.frame()
coords <- adata$obsm[["spatial"]] %>% as.data.frame()
meta <- data.frame(x = coords$V1,
                   y = coords$V2,
                   annotation = meta$annotation %>% as.character(),
                   row.names = rownames(meta), 
                   stringsAsFactors = FALSE)
exp_data <- t(adata$X)
colnames(exp_data) <- rownames(meta)
rownames(exp_data) <- rownames(adata$var)
rm(adata, coords)
target_type <- "Neural crest"

filen <- "./results/Chen2022/E95_new.RDS"
allresults <- readRDS(filen)
CCC_origin <- allresults$CCC
CCC_mat <- Aggr_CCC_results(CCC_origin, flag = "Sender") # sender*LR
CCC_df <- tidyfst::mat_df(CCC_mat)
GRN_filtered <- Process_GRN(allresults$GRN, beta_thres = 1e-5)
XT_results <- Aggr_XT_results(CCC_df, GRN_filtered, threshold = 1e-8)

# 0. Plot the data
table(meta$annotation)
ggplot() + 
  geom_point(data = meta, aes(x = x, y = y, color = annotation),size = 3) +
  # scale_color_brewer()+
  theme_classic() +
  theme(text = element_text(size = 24))+
  guides(color = guide_legend(override.aes = list(size = 12))) +
  labs(color = "Cell type")
ggsave(filename = "./results/Chen2022/supp/Data.PDF", width = 12, height = 12)

# 1. Dominant senders and receivers
colnames(CCC_mat) <- gsub("="," - ", colnames(CCC_mat))
Plot_DominantLR(results = CCC_mat)
ggsave(filename = "./results/Chen2022/supp/DominantLR.PDF", width = 8, height = 14, device = cairo_pdf)

Ls <- sub(" - .*", "", colnames(CCC_mat))
Rs <- sub(".* - ", "", colnames(CCC_mat))
df <- data.frame(Ligand = Ls, Receptor = Rs, Weight = colSums(CCC_mat))
pdf("./results/Chen2022/supp/DominantLR2.PDF", width = 8,height = 8)
SigXTalkR::PlotCCI_CirclePlot(df, topk = 15)
dev.off()

results <- Aggr_CCC_results(CCC_origin, flag = "LR") # sender*LR
max_sender <- apply(results, 2, which.max)
sender_types <- row.names(results)[max_sender]
sender_types <- meta[sender_types,'annotation']
names(sender_types) <- colnames(results)
all_senders <- unique(sender_types) %>% sort()
df <- meta[meta$annotation == target_type, ]
df$sender <- 0
df[names(sender_types), "sender"] <- sender_types
df$sender <- factor(df$sender, levels = c(all_senders))
colors <- c("#004774", "#FF8000FF", "#979736", "#FFFF33FF", 
"#f4e3ea", "#33FF00FF", "#A6EDFFFF", "#1AB2FFFF", "#CCBFFFFF")
ggplot() +
    geom_point(
        data = df,
        aes(x = x, y = y, color = sender_types[rownames(df)]),
        alpha = 0.9,
        size = 4
    ) +
    scale_color_manual(values = colors) +
    coord_fixed() + 
    theme_bw() +
    theme(text = element_text(size = 24)) +
    labs(color = "Dominant sender cell type") +
    guides(color = guide_legend(override.aes = list(size = 12)))

filen <- paste0("./results/Chen2022/figure0130/CCC_sender.PDF")
ggsave(filename = filen, width = 16, height = 16)

aggr_XT <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_results, sum)
XT_toplot <- aggr_XT %>% arrange(desc(Weight)) %>% slice_head(n=15)
uf <- unique(XT_toplot$Ligand) %>% sort()
cols <- paletteer::paletteer_c("grDevices::Set 3", length(uf))
col_map <- setNames(cols, uf)
feature_colors <- col_map[XT_toplot$Ligand]
# Plot_XTAlluvial(aggr_XT[,c(1,2,3,4)], topk = 10)

filen <- paste0("./results/Chen2022/supp/All_","Allu",".PDF")
pdf(filen, width = 4.5, height = 4)
alluvial::alluvial(XT_toplot[,1:3], gap.width = 0.02,xw = 0.2,cw = 0.2,
                  freq = XT_toplot$Weight, border = NA, alpha = 0.75,
                  col = feature_colors,
                  axis_labels = c("Ligand", "Receptor", "Target"))
dev.off()

data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = NULL)
PlotSpa_SRVector(meta = meta, data = data,threshold = .8)
filen <- paste0("./results/Chen2022/supp/Vector_ALL.PDF")
ggsave(filename = filen, device = cairo_pdf, width = 14, height = 14)

# 2. vector field
features <- c("Mif=Cd44","Mdk=Ncl","Mdk=Itgb1")
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[1])
PlotSpa_SRVector(meta = meta, data = data,threshold = .75)
filen <- paste0("./results/Chen2022/supp/Vector_",features[1],".PDF")
ggsave(filename = filen, width = 14, height = 14,device = cairo_pdf)

data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[2])
PlotSpa_SRVector(meta = meta, data = data,threshold = .75)
filen <- paste0("./results/Chen2022/supp/Vector_",features[2],".PDF")
ggsave(filename = filen, width = 14, height = 14,device = cairo_pdf)

data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[3])
PlotSpa_SRVector(meta = meta, data = data,threshold = .75)
filen <- paste0("./results/Chen2022/supp/Vector_",features[3],".PDF")
ggsave(filename = filen, device = cairo_pdf, width = 14, height = 14)



# 3. subgroup
cluster_results <- Group_signals(results = XT_results, meta = meta, features = NULL, groups = NULL, cluster_k = 3)
PlotSpa_Scatterpie(cluster_results, meta, topk = 5, pie_disp = "target")
filen <- paste0("./results/Chen2022/figure0130/Subgroup_","target",".PDF")
ggsave(filename = filen, width = 14, height = 14)

# Group 1
g1_cells <- cluster_results$Group1$Cells
XT_filtered <- filter(XT_results, Cell %in% g1_cells)
aggr_XT <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_filtered, sum) %>%
    group_by(Ligand) %>%
    mutate(Spe = Weight/sum(Weight)) %>%
    ungroup() %>%
    group_by(SSC) %>%
    mutate(Fid = Weight/sum(Weight)) %>%
    ungroup() 
XT_toplot <- aggr_XT %>% arrange(desc(Weight)) %>% slice_head(n=10)
uf <- unique(XT_toplot$Ligand) %>% sort()
cols <- paletteer::paletteer_c("grDevices::Set 3", length(uf))
col_map <- setNames(cols, uf)
feature_colors <- col_map[XT_toplot$Ligand]
# Plot_XTAlluvial(aggr_XT[,c(1,2,3,4)], topk = 10)

filen <- paste0("./results/Chen2022/figure0130/Subgroup1_","Allu",".PDF")
pdf(filen, width = 4.5, height = 4)
alluvial::alluvial(XT_toplot[,1:3], gap.width = 0.02,xw = 0.2,cw = 0.2,
                  freq = XT_toplot$Weight, border = NA, alpha = 0.75,
                  col = feature_colors,
                  axis_labels = c("Ligand", "Receptor", "Target"))
dev.off()

# Group 2
g2_cells <- cluster_results$Group2$Cells
XT_filtered <- filter(XT_results, Cell %in% g2_cells)
aggr_XT <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_filtered, sum) %>%
    group_by(Ligand) %>%
    mutate(Spe = Weight/sum(Weight)) %>%
    ungroup() %>%
    group_by(SSC) %>%
    mutate(Fid = Weight/sum(Weight)) %>%
    ungroup() 
# Plot_XTAlluvial(aggr_XT[,c(1,2,3,4)], topk = 10)
XT_toplot <- aggr_XT %>% arrange(desc(Weight)) %>% slice_head(n=10)
uf <- unique(XT_toplot$Ligand) %>% sort()
cols <- paletteer::paletteer_c("grDevices::Set 3", length(uf))
col_map <- setNames(cols, uf)
feature_colors <- col_map[XT_toplot$Ligand]
# Plot_XTAlluvial(aggr_XT[,c(1,2,3,4)], topk = 10)

filen <- paste0("./results/Chen2022/figure0130/Subgroup2_","Allu",".PDF")
pdf(filen, width = 4.5, height = 4)
alluvial::alluvial(XT_toplot[,1:3], gap.width = 0.02,xw = 0.2,cw = 0.2,
                  freq = XT_toplot$Weight, border = NA, alpha = 0.75,
                  col = feature_colors,
                  axis_labels = c("Ligand", "Receptor", "Target"))
dev.off()

# Group 3
g3_cells <- cluster_results$Group3$Cells
XT_filtered <- filter(XT_results, Cell %in% g3_cells)
aggr_XT <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_filtered, sum) %>%
    group_by(Ligand) %>%
    mutate(Spe = Weight/sum(Weight)) %>%
    ungroup() %>%
    group_by(SSC) %>%
    mutate(Fid = Weight/sum(Weight)) %>%
    ungroup() 
# Plot_XTAlluvial(aggr_XT[,c(1,2,3,4)], topk = 10)
XT_toplot <- aggr_XT %>% arrange(desc(Weight)) %>% slice_head(n=10)
uf <- unique(XT_toplot$Ligand) %>% sort()
cols <- paletteer::paletteer_c("grDevices::Set 3", length(uf))
col_map <- setNames(cols, uf)
feature_colors <- col_map[XT_toplot$Ligand]
# Plot_XTAlluvial(aggr_XT[,c(1,2,3,4)], topk = 10)

filen <- paste0("./results/Chen2022/figure0130/Subgroup3_","Allu",".PDF")
pdf(filen, width = 4.5, height = 4)
alluvial::alluvial(XT_toplot[,1:3], gap.width = 0.02,xw = 0.2,cw = 0.2,
                  freq = XT_toplot$Weight, border = NA, alpha = 0.75,
                  col = feature_colors,
                  axis_labels = c("Ligand", "Receptor", "Target"))
dev.off()

# Specificity
features <- c("Akt1","Myc","Mapk3")
XT_spe_features <- XT_results %>% 
    group_by(Cell, Ligand) %>% 
    mutate(Spe = Weight/sum(Weight)) %>%
    ungroup() %>% filter(SSC %in% features)
aggr_spe <- aggregate(Spe ~ Cell + SSC, data = XT_spe_features, mean)

for (fe in features){
    agg <- filter(aggr_spe, SSC == fe)
    agg$Cell <- as.character(agg$Cell)
    df <- meta
    df$value <- 0
    df[agg$Cell,"value"] <- agg$Spe
    df <- df[order(df$value, decreasing = T), ]
    ggplot() +
        geom_point(
            data = df[df$value > 0,],
            aes(x = x, y = y, color = value),
            alpha = 0.9,
            size = 4
        ) +
        scale_color_gradientn(
            colors = c("yellow","blue"),
            name = "Specificity"
        ) +        
        ggnewscale::new_scale("color") +   
        coord_fixed() + 
        geom_point(
            data = df[df$value == 0,],
            aes(x = x, y = y),
            color = "lightgray",
            alpha = 0.4,
            size = 2
        ) +
        # Theme
        theme_bw() +
        theme(text = element_text(size = 24)) +
        guides(color = guide_colorbar(order = 1))
    filen <- paste0("./results/Chen2022/figure0130/Specificity_",fe,".PDF")
    ggsave(filename = filen, width = 14, height = 14)
}

##### Center vs. border
cell_dist <- Dist_to_border(meta)
cell_dist <- cell_dist[names(CCC_origin)]

CCC_LR <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = NULL)
contrib <- Sender_contribution(CCC_LR, meta)
cont_self <- data.frame(dist = cell_dist, cont = contrib[names(cell_dist),"Neural crest"])
unique(cont_self$dist)
group <- cell_dist[rownames(contrib)]
names(group) <- rownames(contrib)
distance_threshold <- c(1.01,1.51,2.01,3.01)
group <- cut(group, breaks = c(-Inf, distance_threshold, Inf), labels = c("0-1", "1-1.5", "1.5-2", "2-3", ">3")) %>% as.character()
names(group) <- rownames(contrib)
cont_results <- lapply(unique(group), function(g){
    cells <- names(group)[group == g]
    contrib_sub <- contrib[cells,] %>% tidyfst::mat_df()
    cont_sum <- aggregate(value~col, data = contrib_sub, mean)
    cont_sum$group <- g
    colnames(cont_sum) <- c("Sender","contribution","group")
    cont_sum
}) 
cont_results <- do.call(rbind, cont_results)
cont_results$group <- factor(cont_results$group, levels = c("0-1", "1-1.5", "1.5-2", "2-3", ">3"))
table(cont_results$Sender)
cont_results$Sender <- factor(cont_results$Sender, levels = all_senders)
ggplot(cont_results, aes(x = group, y = contribution, fill = Sender)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    theme(legend.position = "none") +
    labs(x = "Distance to border", y = "Average contribution", fill = "Sender cell type") +
    theme(axis.text.x = element_text(size = 24))
filen <- paste0("./results/Chen2022/supp/Contribution_by_distance.PDF")
ggsave(filename = filen, width = 6, height = 14)

# downstream
ligands <- c("Mif","Fn1","Col1a2","Col2a1","Col4a1")
d1_cells <- names(group[group == "0-1"]) %>% sort()
XT_used <- filter(XT_results,Cell %in% d1_cells) %>% filter(Ligand %in% ligands)
aggr_XT1 <- aggregate(Weight ~ Ligand + SSC, data = XT_used, sum) %>%
  group_by(Ligand) %>% 
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()

d2_cells <- names(group[group == "1.5-2"]) %>% sort()
XT_used <- filter(XT_results,Cell %in% d2_cells) %>% filter(Ligand %in% ligands)
aggr_XT2 <- aggregate(Weight ~ Ligand + SSC, data = XT_used, sum) %>%
  group_by(Ligand) %>% 
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()

d3_cells <- names(group[group == ">3"]) %>% sort()
XT_used <- filter(XT_results,Cell %in% d3_cells) %>% filter(Ligand %in% ligands)
aggr_XT3 <- aggregate(Weight ~ Ligand + SSC, data = XT_used, sum) %>%
  group_by(Ligand) %>% 
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()

g1 <- aggr_XT1 %>% group_by(SSC) %>% summarise(total = sum(Spe)) %>% arrange(desc(total)) %>% slice_head(n = 10) %>% pull(SSC)
g2 <- aggr_XT2 %>% group_by(SSC) %>% summarise(total = sum(Spe)) %>% arrange(desc(total)) %>% slice_head(n = 10) %>% pull(SSC)
g3 <- aggr_XT3 %>% group_by(SSC) %>% summarise(total = sum(Spe)) %>% arrange(desc(total)) %>% slice_head(n = 10) %>% pull(SSC)
genes <- c(g1,g2,g3) %>% unique() %>% sort() 

aggr_XT1 %>% filter(SSC %in% genes) %>% 
  ggplot() +
  geom_point(aes(x = SSC, y = Ligand, color = Spe),size = 8) +
  scale_color_gradientn(colors = c("lightgray","#7bdeff","blue"),limits = c(0,0.12),name = "Specificity",oob = scales::squish) +
  theme_bw() +
  theme(text = element_text(size = 24))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "right") +
  labs(x = "TF", y = "Ligand")
filen <- paste0("./results/Chen2022/supp/Spe_Close.PDF")
ggsave(filename = filen, width = 15, height = 4)

aggr_XT2 %>% filter(SSC %in% genes) %>% 
  ggplot() +
  geom_point(aes(x = SSC, y = Ligand, color = Spe),size = 8) +
  scale_color_gradientn(colors = c("lightgray","#7bdeff","blue"),limits = c(0,0.12),name = "Specificity",oob = scales::squish) +
  theme_bw() +
  theme(text = element_text(size = 24))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "right") +
  labs(x = "TF", y = "Ligand")
filen <- paste0("./results/Chen2022/supp/Spe_medium.PDF")
ggsave(filename = filen, width = 15, height = 4)

aggr_XT3 %>% filter(SSC %in% genes) %>% 
  ggplot() +
  geom_point(aes(x = SSC, y = Ligand, color = Spe),size = 8) +
  scale_color_gradientn(colors = c("lightgray","#7bdeff","blue"),limits = c(0,0.12),name = "Specificity",oob = scales::squish) +
  theme_bw() +
  theme(text = element_text(size = 24))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "right") +
  labs(x = "TF", y = "Ligand")
filen <- paste0("./results/Chen2022/supp/Spe_Far.PDF")
ggsave(filename = filen, width = 15, height = 4)
aggregate(Spe~SSC,data = aggr_XT3,sum)

# sample study
color_use <- c("aquamarine","bisque","pink")
names(color_use) <- c("Ligand","Receptor","Target")
XT_filtered <- XT_results[XT_results$Cell == "185_151", ] %>% arrange(desc(Weight)) %>% slice_head(n = 6)
Plot_MultiLayer_Singlecell(XT_filtered, start_trim = 0.3, end_trim = 0.3, arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,rect_height = 0.3, label_size = 6,layer_fill_colors = color_use)
ggsave(filename = "./results/Chen2022/figure0130/multilayer1.PDF", width = 3600, height = 2400, units = "px")


XT_filtered <- XT_results[XT_results$Cell == "185_153", ] %>% arrange(desc(Weight)) %>% slice_head(n = 6)
Plot_MultiLayer_Singlecell(XT_filtered, start_trim = 0.3, end_trim = 0.3, arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,rect_height = 0.3, label_size = 6,layer_fill_colors = color_use)
ggsave(filename = "./results/Chen2022/figure0130/multilayer2.PDF", width = 3600, height = 2400, units = "px")

####### SCENIC PART ########
# run scenic
write.table(g1_cells, "/home/shared/SpatialCrosstalk/codes/scenic/g1_cells.txt",quote = F, row.names = F, col.names = F)
write.table(g2_cells, "/home/shared/SpatialCrosstalk/codes/scenic/g2_cells.txt",quote = F, row.names = F, col.names = F)
write.table(g3_cells, "/home/shared/SpatialCrosstalk/codes/scenic/g3_cells.txt",quote = F, row.names = F, col.names = F)
tfs <- XT_results$SSC %>% unique() %>% sort()
write.table(tfs, "/home/shared/SpatialCrosstalk/codes/scenic/tf_list.txt",quote = F, row.names = F, col.names = F)

# Prepare the data
allgenes <- rownames(exp_data)
keep <- !grepl("Rik$|^Gm", allgenes) & !grepl("^AC-", allgenes) & !grepl("^BC-", allgenes) & !grepl("^Mir-", allgenes) & !grepl("^mt-", allgenes) & !grepl("^Rps|^Rpl", allgenes)
exp_data <- exp_data[keep, ]
SeuratObj <- CreateSeuratObject(exp_data, meta.data = meta)
SeuratObj <- NormalizeData(SeuratObj)
Idents(SeuratObj) <- SeuratObj$annotation
markers <- FindMarkers(SeuratObj, ident.1 = target_type, min.pct = 0.25, only.pos = T)
markers <- markers %>% filter(p_val_adj < 0.01)
markers <- rownames(markers) %>% sort()

# Read the results
# subgroup1
SCENIC_results1 <- read.table("/home/shared/SpatialCrosstalk/results/Chen2022/output_1.csv", header = T, sep = ",")
g1_cells <- read.table("/home/shared/SpatialCrosstalk/codes/scenic/g1_cells.txt") %>% unlist() %>% as.vector()
XT_filtered <- filter(XT_results, Cell %in% g1_cells)
SCENIC_results1 <- filter(SCENIC_results1, TF %in% XT_filtered$SSC) %>% filter(target %in% markers)
# df1 <- filter(SCENIC_results1, TF %in% genes) %>% filter(target %in% markers)
head(SCENIC_results1)
topk <- 5

features <- c("Ddx5","Hif1a","Twist1","Mapk1","Mapk8","Myc")
SCENIC_filtered <- SCENIC_results1 %>% filter(TF %in% features)
df1 <- SCENIC_filtered %>% 
    group_by(TF) %>% 
    slice_max(order_by = importance, n = topk) %>%
    ungroup()
cells <- cluster_results$Group1$Cells
exp_used <- exp_data[,cells]
PlotXT_MultiCircularBar(df1[,c(2,1,3)])
filen <- paste0("./results/Chen2022/supp/SCENIC_EdgeBundle1.PDF")
ggsave(filename = filen, width = 16, height = 16)

# subgroup2
SCENIC_results2 <- read.table("/home/shared/SpatialCrosstalk/results/Chen2022/output_2.csv", header = T, sep = ",")
g2_cells <- read.table("/home/shared/SpatialCrosstalk/codes/scenic/g2_cells.txt") %>% unlist() %>% as.vector()
XT_filtered <- filter(XT_results, Cell %in% g2_cells)
SCENIC_results2 <- filter(SCENIC_results2, TF %in% XT_filtered$SSC) %>% filter(target %in% markers)
# df2 <- filter(SCENIC_results2, TF %in% genes) %>% filter(target %in% markers)
head(SCENIC_results2)
topk <- 5

features <- c("Ddx5","Hif1a","Twist1","Mapk1","Mapk8","Myc")
SCENIC_filtered <- SCENIC_results2 %>% filter(TF %in% features)
df2 <- SCENIC_filtered %>% 
    group_by(TF) %>% 
    slice_max(order_by = importance, n = topk) %>%
    ungroup()
cells <- cluster_results$Group2$Cells
exp_used <- exp_data[,cells]
PlotXT_MultiCircularBar(df2[,c(2,1,3)])
filen <- paste0("./results/Chen2022/supp/SCENIC_EdgeBundle2.PDF")
ggsave(filename = filen, width = 16, height = 16, device = cairo_pdf)

# subgroup3
SCENIC_results3 <- read.table("/home/shared/SpatialCrosstalk/results/Chen2022/output_3.csv", header = T, sep = ",")
g3_cells <- read.table("/home/shared/SpatialCrosstalk/codes/scenic/g3_cells.txt")%>% unlist() %>% as.vector()
XT_filtered <- filter(XT_results, Cell %in% g3_cells)
SCENIC_results3 <- filter(SCENIC_results3, TF %in% XT_filtered$SSC) %>% filter(target %in% markers)
# df3 <- filter(SCENIC_results3, TF %in% genes) %>% filter(target %in% markers)
head(SCENIC_results3)
topk <- 5
features <- c("Ddx5","Hif1a","Twist1","Mapk1","Mapk8","Myc")

SCENIC_filtered <- SCENIC_results3 %>% filter(TF %in% features)
df3 <- SCENIC_filtered %>% 
    group_by(TF) %>% 
    slice_max(order_by = importance, n = topk) %>%
    ungroup()
cells <- cluster_results$Group3$Cells
exp_used <- exp_data[,cells]
PlotXT_MultiCircularBar(df3[,c(2,1,3)])
filen <- paste0("./results/Chen2022/supp/SCENIC_EdgeBundle3.PDF")
ggsave(filename = filen, width = 16, height = 16, device = cairo_pdf)

### Plot
df1$group <- "Group1"
df2$group <- "Group2"
df3$group <- "Group3"
df_all <- bind_rows(df1, df2, df3)
group_levels <- c("Group1", "Group2", "Group3")
colnames(df_all) <- c("source","target","strength","group")
# complete all combinations, fill missing strengths with 0
all_sources <- sort(unique(df_all$source))
all_targets <- sort(unique(df_all$target))

df_plot <- df_all %>%
  mutate(
    source = as.character(source),
    target = as.character(target),
    group  = factor(group, levels = group_levels)
  ) %>%
  tidyr::complete(
    source = all_sources,
    target = all_targets,
    group  = group_levels,
    fill = list(strength = 0)
  )

# build x-axis levels: for each source, show 3 group columns
x_levels <- as.vector(sapply(all_sources, function(s) paste(s, group_levels, sep = "__")))

df_plot <- df_plot %>%
  mutate(
    x_col = factor(paste(source, group, sep = "__"), levels = x_levels),
    target = factor(target, levels = rev(all_targets))
  )

# only plot positive-strength links
df_nonzero <- df_plot %>%
  filter(strength > 0)
sep_pos <- seq(3.5, length(x_levels) - 0.5, by = 3)
df_nonzero[df_nonzero$target == "Basp1", ]
p <- ggplot(df_nonzero, aes(x = x_col, y = target)) +
  geom_point(aes(fill = strength), shape = 21, size = 6) +
  viridis::scale_fill_viridis(option = "H") +
  scale_x_discrete(
    drop = FALSE,   # <- critical
    labels = rep(group_levels, times = length(all_sources))
  ) +
  scale_y_discrete(drop = FALSE) +
  # coord_flip() +
  theme_bw(base_size = 14) +
  theme(text = element_text(size = 18)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_line(color = "grey90")
  ) +
  labs(x = NULL, y = NULL, fill = "Strength") +
  geom_vline(xintercept = sep_pos, color = "grey70")
p
p +
  annotate(
    "text",
    x = source_centers,
    y = length(all_targets) + 1,
    label = all_sources,
    size = 4,
    fontface = 2
  )
filen <- paste0("./results/Chen2022/figure0130/SCENIC_EdgeBundle.PDF")
ggsave(filename = filen, width = 10, height = 16)


