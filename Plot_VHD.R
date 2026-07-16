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

exp_data <- readRDS("./datasets/VisiumHD/expression_used.RDS")
meta <- read.table("./datasets/VisiumHD/meta_used.csv",sep = ",")
table(meta$annotation)
tumor_types <- c("Tumor I","Tumor II","Tumor III","Tumor IV","Tumor V")
meta$tumor <- ifelse(meta$annotation %in% tumor_types, "Tumor","Non-tumor")

# Data overview
ggplot() +
  geom_point(data = meta, aes(x = x, y = y, color = tumor),size = 6) +
  scale_color_manual(values = c("#dcdcdc", "#ff345d")) +
  theme_void() +
  theme(text = element_text(size = 24),legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 12)))

ggsave(filename = "./results/VHD/Figure0614/Data.PDF", width = 12, height = 12)

ggplot() +
  geom_point(data = meta, aes(x = x, y = y, color = annotation),size = 6) +
  theme_void() +
  theme(text = element_text(size = 24),legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 12)))

ggsave(filename = "./results/VHD/supp/Data_detail.PDF", width = 20, height = 12)

# Read results
filen <- "./results/VHD/VHD.RDS"
allresults <- readRDS(filen)
CCC_origin <- allresults$CCC
CCC_mat <- Aggr_CCC_results(CCC_origin, flag = "Sender") # sender*LR
CCC_df <- tidyfst::mat_df(CCC_mat)
GRN_filtered <- Process_GRN(allresults$GRN, beta_thres = 1e-5)
XT_results <- Aggr_XT_results(CCC_df, GRN_filtered, threshold = 1e-8)
XT_results <- XT_results %>% 
    group_by(Cell) %>% 
    mutate(Weight_new = Weight/sum(Weight)) %>% 
    ungroup()

receivers <- rownames(CCC_mat)
tumorcell <- rownames(meta[meta$tumor == "Tumor", ])
nontumorcell <- rownames(meta[meta$tumor == "Non-tumor", ])
meta <- meta[receivers, ]

# ITGAV expression 
data <- exp_data["NR4A1",receivers]
df <- meta
df$value <- 0
df[receivers,"value"] <- data
ggplot() +
    geom_point(
        data = df[df$value > 0,],
        aes(x = x, y = y, color = value),
        alpha = 0.9,
        size = 4
    ) +
    scale_color_gradientn(
        colors = c("#ffe65c","#0011ff"),
        name = "Expression"
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

# CCC
CCC_mat <- Aggr_CCC_results(CCC_origin[nontumorcell], flag = "Sender") # sender*LR
Ls <- sub("=.*", "", colnames(CCC_mat))
Rs <- sub(".*=", "", colnames(CCC_mat))
df <- data.frame(Ligand = Ls, Receptor = Rs, Weight = colSums(CCC_mat))
pdf("./results/VHD/Figure0614/DominantLR_nontumor.PDF", width = 8,height = 8)
SigXTalkR::PlotCCI_CirclePlot(df, topk = 15)
dev.off()

CCC_mat <- Aggr_CCC_results(CCC_origin[tumorcell], flag = "Sender") # sender*LR
Ls <- sub("=.*", "", colnames(CCC_mat))
Rs <- sub(".*=", "", colnames(CCC_mat))
df <- data.frame(Ligand = Ls, Receptor = Rs, Weight = colSums(CCC_mat))
pdf("./results/VHD/Figure0614/DominantLR_tumor.PDF", width = 8,height = 8)
SigXTalkR::PlotCCI_CirclePlot(df, topk = 15)
dev.off()

colors <- c("#d6d600","blue")
features <- c("COL1A1=ITGAV","COL1A1=SDC4","CEACAM5=CEACAM6")
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[1])
PlotSpa_SRVector(meta = meta, data = data,threshold = .5, colors = colors)
filen <- paste0("./results/VHD/Figure0614/Vector_",features[1],".PDF")
ggsave(filename = filen, device = cairo_pdf, width = 6, height = 6)

data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[2])
PlotSpa_SRVector(meta = meta, data = data,threshold = .5)
filen <- paste0("./results/VHD/Figure0614/Vector_",features[2],".PDF")
ggsave(filename = filen, device = cairo_pdf, width = 6, height = 6)

data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[3])
PlotSpa_SRVector(meta = meta, data = data,threshold = .75)
filen <- paste0("./results/VHD/Figure0614/Vector_",features[3],".PDF")
ggsave(filename = filen, device = cairo_pdf, width = 6, height = 6)

# Downstream
XT1 <- XT_results[XT_results$Cell %in% nontumorcell, ]
XT2 <- XT_results[XT_results$Cell %in% tumorcell, ]
XT2_aggr <- aggregate(Weight ~ SSC, data = XT2, sum)

# Receptor comparison 
Rec_aggr <- aggregate(Weight ~ Receptor + Cell, data = XT_results, mean)
Rec_aggr$tumor <- meta[Rec_aggr$Cell,"tumor"]
Rec_aggr$tumor <- factor(Rec_aggr$tumor, levels = c("Tumor","Non-tumor"))
unique(Rec_aggr$Receptor)
Receptors_use <- c("CD44", "CLDN3", "CXCR4", "IL2RG", "ITGA6", "ITGAV", 
"ITGB4", "ITGB5", "MMP14", "PDGFRB", "TYROBP")
Rec_aggr <- Rec_aggr[Rec_aggr$Receptor %in% Receptors_use, ]
ggplot(Rec_aggr, aes(x = Receptor, y = Weight, fill = tumor))+
  geom_boxplot(
    width = 0.3,
    outlier.size = 0.3,
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  stat_compare_means(
    aes(group = tumor),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE
  ) +
  ylim(0,3)+
  scale_fill_manual(values = c("Non-tumor" = "#9bcac8", "Tumor" = "#ca94c3")) +
  xlab("") +
  ylab("Average activation") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    text = element_text(size = 18)
  )
filen <- paste0("./results/VHD/Figure0614/Compare_receptor.PDF")
ggsave(filename = filen, width = 10, height = 6)

Tar_aggr <- aggregate(Weight ~ SSC + Cell, data = XT_results, sum) %>%
    group_by(Cell) %>% 
    mutate(Spe = Weight/sum(Weight)) %>%
    ungroup()
Tar_aggr$tumor <- meta[Tar_aggr$Cell,"tumor"]
Tar_aggr$tumor <- factor(Tar_aggr$tumor, levels = c("Tumor","Non-tumor"))
Tar_aggr$Cell <- as.character(Tar_aggr$Cell)

ggplot(Tar_aggr, aes(x = SSC, y = Weight, fill = tumor))+
  geom_boxplot(
    width = 0.5,
    outlier.size = 0.3,
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  stat_compare_means(
    aes(group = tumor),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE
  ) +
  ylim(0,1)+
  scale_fill_manual(values = c("Non-tumor" = "#9bcac8", "Tumor" = "#ca94c3")) +
  xlab("") +
  ylab("Average activation") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    text = element_text(size = 18)
  )
filen <- paste0("./results/VHD/Figure0614/Compare_tar.PDF")
ggsave(filename = filen, width = 10, height = 6)

# Receptor - Target
XT1_aggr <- aggregate(Weight ~ Receptor + SSC, data = XT1, sum)
XT1_aggr <- XT1_aggr %>% 
    arrange(Weight) %>% 
    slice_tail(n = 20)
Plot_EdgeBundle(XT1_aggr,text_size = 8, edge_colors = c("yellow","blue"),node_size = 7, edge_width = 2)
filen <- paste0("./results/VHD/Figure0614/edge_tar_nontumor.PDF")
ggsave(filename = filen, width = 10, height = 10)

XT2_aggr <- aggregate(Weight ~ Receptor + SSC, data = XT2, sum)
XT2_aggr <- XT2_aggr %>% 
    arrange(Weight) %>% 
    slice_tail(n = 20)
Plot_EdgeBundle(XT2_aggr,text_size = 8, edge_colors = c("yellow","blue"),node_size = 7, edge_width = 2)
filen <- paste0("./results/VHD/Figure0614/edge_tar_tumor.PDF")
ggsave(filename = filen, width = 10, height = 10)

# JUN
# Plot specificity
data <- Tar_aggr[Tar_aggr$SSC == "JUN", ]
wilcox.test(Spe ~ tumor, data = data)
df <- meta
df$value <- 0
df[data$Cell,"value"] <- data$Spe
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
    theme_void() +
    theme(text = element_text(size = 24),legend.position = "none") +
    guides(color = guide_colorbar(order = 1))
filen <- paste0("./results/VHD/Figure0614/Spe_JUN.PDF")
ggsave(filename = filen, width = 7.5, height = 7.5)

# Supp
# Upstream of JUN
XT_filtered <- XT_results[XT_results$SSC == "JUN", ] %>%
  group_by(Cell) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup()
df <- aggregate(Weight ~ Receptor + Ligand, data = XT_filtered, sum) %>% 
  arrange(desc(Weight)) %>%
  slice_head(n=50)
ggplot(df) +
  geom_point(aes(x = Ligand, y = Receptor, color = Weight),size = 8) +
  scale_color_gradientn(colors = c("white", "blue")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 24,angle = 90,vjust = 0.5)) +
  labs(x = "Ligand", y = "") +
  theme(text = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24))
filen <- paste0("./results/VHD/supp/upstream_JUN.PDF")
ggsave(filename = filen, width = 15, height = 7.5)

XT_rec <- aggregate(Fid ~ Cell + Receptor, data = XT_filtered, sum)

data <- XT_rec[XT_rec$Receptor == "ITGAV", ]
rownames(data) <- data$Cell
df <- meta
df$value <- 0
df[rownames(data),"value"] <- data$Fid
ggplot() +
    geom_point(
        data = df[df$value > 0,],
        aes(x = x, y = y, color = value),
        alpha = 0.9,
        size = 4
    ) +
    scale_color_gradientn(
        colors = c("yellow","blue"),
        name = "Fidelity"
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
    theme_void() +
    theme(text = element_text(size = 24),legend.position = "right") +
    guides(color = guide_colorbar(order = 1))
filen <- paste0("./results/VHD/supp/Fid_ITGAV.PDF")
ggsave(filename = filen, width = 7.5, height = 7.5)

data <- XT_rec[XT_rec$Receptor == "ITGB4", ]
rownames(data) <- data$Cell
df <- meta
df$value <- 0
df[rownames(data),"value"] <- data$Fid
ggplot() +
    geom_point(
        data = df[df$value > 0,],
        aes(x = x, y = y, color = value),
        alpha = 0.9,
        size = 4
    ) +
    scale_color_gradientn(
        colors = c("yellow","blue"),
        name = "Fidelity"
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
    theme_void() +
    theme(text = element_text(size = 24),legend.position = "right") +
    guides(color = guide_colorbar(order = 1))
filen <- paste0("./results/VHD/supp/Fid_ITGB4.PDF")
ggsave(filename = filen, width = 7.5, height = 7.5)

# upstream of NR4A1
XT_filtered <- XT_results[XT_results$SSC == "NR4A1", ] %>%
  group_by(Cell) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup()
XT1 <- XT_filtered[XT_filtered$Cell %in% nontumorcell, ]
XT2 <- XT_filtered[XT_filtered$Cell %in% tumorcell, ]
XT1_aggr <- aggregate(Weight ~ Ligand + Receptor, data = XT1, sum)
XT1_aggr <- XT1_aggr %>% 
    arrange(Weight) %>% 
    slice_tail(n = 20)
Plot_EdgeBundle(XT1_aggr,text_size = 8, edge_colors = c("yellow","blue"),node_size = 7, edge_width = 2)
filen <- paste0("./results/VHD/supp/edge_NR4A1_nontumor.PDF")
ggsave(filename = filen, width = 10, height = 10)

XT2_aggr <- aggregate(Weight ~ Ligand + Receptor, data = XT2, sum)
XT2_aggr <- XT2_aggr %>% 
    arrange(Weight) %>% 
    slice_tail(n = 20)
Plot_EdgeBundle(XT2_aggr,text_size = 8, edge_colors = c("yellow","blue"),node_size = 7, edge_width = 2)
filen <- paste0("./results/VHD/supp/edge_NR4A1_tumor.PDF")
ggsave(filename = filen, width = 10, height = 10)

# 
data <- exp_data["NR4A1", ]
df <- meta
df$value <- 0
df[names(data),"value"] <- data
ggplot() +
    geom_point(
        data = df[df$value > 0,],
        aes(x = x, y = y, color = value),
        alpha = 0.9,
        size = 4
    ) +
    scale_color_gradientn(
        colors = c("yellow","blue"),
        name = "Fidelity"
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
    theme_void() +
    theme(text = element_text(size = 24),legend.position = "right") +
    guides(color = guide_colorbar(order = 1))
filen <- paste0("./results/VHD/supp/exp_NR4A1.PDF")
ggsave(filename = filen, width = 7.5, height = 7.5)
