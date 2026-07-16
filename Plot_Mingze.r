closeAllConnections()
gc()
library(Seurat)
library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggnewscale)
library(ggforce)
library(ggpubr)
library(ggrepel)
library(scatterpie)
library(viridis)
library(Matrix)
library(ggalluvial)
library(RColorBrewer)
setwd("/home/shared/SpatialCrosstalk")
source("./codes/Process_results.R")
source("./codes/Visualization.R")

# Load the meta info
adata <- anndata::read_h5ad("./datasets/Mingze2025/Mingze_CosMx.h5ad")
exp_data_all <- adata$X %>% t()
rownames(exp_data_all) <- rownames(adata$var)
colnames(exp_data_all) <- rownames(adata$obs)

meta_all <- adata$obs
spatial_info <- adata$obsm$spatial
meta_all$x <- spatial_info[,1]; meta_all$y <- spatial_info[,2]
meta_all$annotation <- meta_all$celltype

# PD sample
CPID <- 31808
meta_PD <- meta_all[meta_all$CPID == CPID, ]
meta_PD <- meta_PD[,c("x","y","annotation")]
result_PD <- readRDS("./results/Mingze2025/Mingze_PD_new.RDS")
CCC_origin_PD <- result_PD$CCC
CCC_mat_PD <- Aggr_CCC_results(CCC_origin_PD, flag = "Sender") # sender*LR
CCC_df <- tidyfst::mat_df(CCC_mat_PD)
df_PD <- aggregate(value ~ col, data = CCC_df, sum)
df_PD$col <- as.character(df_PD$col)
df_PD$ligand <- sapply(strsplit(df_PD$col, split = "="), `[`, 1)
df_PD$receptor <- sapply(strsplit(df_PD$col, split = "="), `[`, 2)
df_PD <- df_PD[,c("ligand","receptor","value")]
colnames(df_PD) <- c("Ligand","Receptor","Weight")
pdf("./results/Mingze2025/supp/CCC_circle_PD.pdf", width = 8, height = 8)
SigXTalkR::PlotCCI_CirclePlot(df_PD, topk = 15)
dev.off()
GRN_filtered_PD <- Process_GRN(result_PD$GRN, beta_thres = 1e-5)
XT_results_PD <- Aggr_XT_results(CCC_df, GRN_filtered_PD, threshold = 1e-8) %>%
  group_by(Cell, Ligand) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()
Aggr_XT_PD <- aggregate(Weight ~ Ligand + Receptor + SSC, data = XT_results_PD, sum) # Cell group level
exp_PD <- exp_data_all[,rownames(meta_PD)]

# CR sample
CPID <- 31810
meta_CR <- meta_all[meta_all$CPID == CPID, ]
meta_CR <- meta_CR[,c("x","y","annotation")]
result_CR <- readRDS("./results/Mingze2025/Mingze_CR_new.RDS")
CCC_origin_CR <- result_CR$CCC
CCC_mat_CR <- Aggr_CCC_results(CCC_origin_CR, flag = "Sender") # sender*LR
CCC_df <- tidyfst::mat_df(CCC_mat_CR)
df_CR <- aggregate(value ~ col, data = CCC_df, sum)
df_CR$col <- as.character(df_CR$col)
df_CR$ligand <- sapply(strsplit(df_CR$col, split = "="), `[`, 1)
df_CR$receptor <- sapply(strsplit(df_CR$col, split = "="), `[`, 2)
df_CR <- df_CR[,c("ligand","receptor","value")]
colnames(df_CR) <- c("Ligand","Receptor","Weight")
pdf("./results/Mingze2025/supp/CCC_circle_CR.pdf", width = 8, height = 8)
SigXTalkR::PlotCCI_CirclePlot(df_CR, topk = 15)
dev.off()
GRN_filtered_CR <- Process_GRN(result_CR$GRN, beta_thres = 1e-5)
XT_results_CR <- Aggr_XT_results(CCC_df, GRN_filtered_CR, threshold = 1e-8) %>%
  group_by(Cell, Ligand) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()
Aggr_XT_CR <- aggregate(Weight ~ Ligand + Receptor + SSC, data = XT_results_CR, sum) # Cell group level
exp_CR <- exp_data_all[,rownames(meta_CR)]

############## FIGURES ##############

# 1. Overview of the data
ggplot() + 
  geom_point(data = meta_PD, aes(x = x, y = y, color = annotation),size = 4) +
  # scale_color_brewer()+
  theme_classic() +
  theme(text = element_text(size = 24))+
  guides(color = guide_legend(override.aes = list(size = 12))) +
  labs(color = "Cell type")
ggsave(filename = "./results/Mingze2025/supp/Data_PD.PDF", width = 12, height = 12)

ggplot() + 
  geom_point(data = meta_CR, aes(x = x, y = y, color = annotation),size = 4) +
  # scale_color_brewer()+
  theme_classic() +
  theme(text = element_text(size = 24))+
  guides(color = guide_legend(override.aes = list(size = 12))) +
  labs(color = "Cell type")
ggsave(filename = "./results/Mingze2025/supp/Data_CR.PDF", width = 12, height = 12)

# 2. Top-3 signals with the downstream targets
Ls <- sub("=.*", "", colnames(CCC_mat_CR))
Rs <- sub(".*=", "", colnames(CCC_mat_CR))
df <- data.frame(Ligand = Ls, Receptor = Rs, Weight = colSums(CCC_mat_CR))
pdf("./results/Mingze2025/supp/CCC_CR.PDF", width = 8,height = 8)
SigXTalkR::PlotCCI_CirclePlot(df, topk = 15)
dev.off()

Ls <- sub("=.*", "", colnames(CCC_mat_PD))
Rs <- sub(".*=", "", colnames(CCC_mat_PD))
df <- data.frame(Ligand = Ls, Receptor = Rs, Weight = colSums(CCC_mat_PD))
pdf("./results/Mingze2025/supp/CCC_PD.PDF", width = 8,height = 8)
SigXTalkR::PlotCCI_CirclePlot(df, topk = 15)
dev.off()

Plot_DominantLR(results = CCC_mat_PD)
Plot_DominantLR(results = CCC_mat_CR)


PD_tar <- aggregate(XT_results_PD$Weight, list(XT_results_PD$Ligand,XT_results_PD$SSC), sum)
colnames(PD_tar) <- c("Signal","Target","Value")
top3_signal_PD <- PD_tar %>%
  group_by(Signal) %>%
  summarise(SumValue = sum(Value), .groups = "drop") %>%
  arrange(desc(SumValue)) %>%
  slice_head(n = 3) %>% select(Signal) %>% unlist()

PD_tar <- PD_tar %>%
    group_by(Signal) %>%
    mutate(Spe = Value / sum(Value)) %>%
    ungroup()

PD_tar_filtered <- filter(PD_tar, Signal %in% top3_signal_PD)

CR_tar <- aggregate(XT_results_CR$Weight, list(XT_results_CR$Ligand,XT_results_CR$SSC), sum)
colnames(CR_tar) <- c("Signal","Target","Value")
top3_signal_CR <- CR_tar %>%
  group_by(Signal) %>%
  summarise(SumValue = sum(Value), .groups = "drop") %>%
  arrange(desc(SumValue)) %>%
  slice_head(n = 3) %>% select(Signal) %>% unlist()

intersect(top3_signal_PD, top3_signal_CR)
setdiff(top3_signal_PD, top3_signal_CR)
setdiff(top3_signal_CR, top3_signal_PD)

CR_tar <- CR_tar %>%
    group_by(Signal) %>%
    mutate(Spe = Value / sum(Value)) %>%
    ungroup()

CR_tar_filtered <- filter(CR_tar, Signal %in% top3_signal_CR)
CR_tar_filtered[CR_tar_filtered$Signal == "MIF","Signal"] <- "MIF_CR"
CR_tar_filtered[CR_tar_filtered$Signal == "FN1","Signal"] <- "FN1_CR"

seqs <- c("COL9A2","FN1","MIF","MIF_CR","FN1_CR","SPP1")
top3_XT <- rbind(PD_tar_filtered,CR_tar_filtered)
top3_XT$Signal <- factor(top3_XT$Signal,levels = seqs)
PlotXT_MultiCircularBar(df = top3_XT[,c(2,1,4)])
ggsave("./results/Mingze2025/figure0109/Dominant_signals.PDF", width = 12, height = 12)

mean(XT_results_CR$Weight)

# number/strength of CD44-pathways per cell: comparison
CD44_PD <- filter(XT_results_PD, Receptor == "CD44")
CD44_CR <- filter(XT_results_CR, Receptor == "CD44")

num_CD44_PD <- CD44_PD %>% group_by(Cell) %>% summarise(Num = n())
num_CD44_CR <- CD44_CR %>% group_by(Cell) %>% summarise(Num = n())
ggplot() +
  geom_density(data = num_CD44_PD, aes(x = Num), fill = "blue", alpha = 0.5) +
  geom_density(data = num_CD44_CR, aes(x = Num), fill = "red", alpha = 0.5) +
  labs(x = "Number of pathways", y = "Density") +
  theme_classic() +
  theme(text = element_text(size = 24)) +
  scale_fill_manual(values = c("blue", "red"), labels = c("PD", "CR")) +
  guides(fill = guide_legend(title = "Group"))
ggsave("./results/Mingze2025/supp/CD44_pathways_num.PDF", width = 12, height = 12)


# number/strength of ITGB1-pathways per cell: comparison
ITGB1_PD <- filter(XT_results_PD, Receptor == "ITGB1")
ITGB1_CR <- filter(XT_results_CR, Receptor == "ITGB1")
num_ITGB1_PD <- ITGB1_PD %>% group_by(Cell) %>% summarise(Num = n())
num_ITGB1_CR <- ITGB1_CR %>% group_by(Cell) %>% summarise(Num = n())
ggplot() +
  geom_density(data = num_ITGB1_PD, aes(x = Num), fill = "blue", alpha = 0.5) +
  geom_density(data = num_ITGB1_CR, aes(x = Num), fill = "red", alpha = 0.5) +
  labs(x = "Number of pathways", y = "Density") +
  theme_classic() +
  theme(text = element_text(size = 24)) +
  scale_fill_manual(values = c("blue", "red"), labels = c("PD", "CR")) +
  guides(fill = guide_legend(title = "Group"))
ggsave("./results/Mingze2025/supp/ITGB1_pathways_num.PDF", width = 12, height = 12)

# specificity distribution of CD44-pathways
targets <- c("CTNNB1","CCND1","MYC")
CD44_CR <- CD44_CR %>% 
  group_by(Cell) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()
df_CR <- aggregate(Spe ~ SSC + Cell, data = CD44_CR, sum) %>% filter(SSC %in% targets)
CD44_PD <- CD44_PD %>% 
  group_by(Cell) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()
df_PD <- aggregate(Spe ~ SSC + Cell, data = CD44_PD, sum) %>% filter(SSC %in% targets)
df_PD$Group <- "PD"
df_CR$Group <- "CR"
df_all <- rbind(df_PD, df_CR)
df_all$Group <- factor(df_all$Group, levels = c("CR","PD"))
ggplot(df_all, aes(x = SSC, y = Spe, fill = Group)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24)) +
  labs(x = "Target", y = "Specificity") +
  guides(fill = guide_legend(title = "Group"))

t.test(Spe ~ Group, data = df_all[df_all$SSC == "CTNNB1", ])
ggsave("./results/Mingze2025/supp/CD44_pathways_spe.PDF", width = 8, height = 12)

# specificity distribution of ITGB1-pathways
ITGB1_CR <- ITGB1_CR %>% 
  group_by(Cell) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()
df_CR <- aggregate(Spe ~ SSC + Cell, data = ITGB1_CR, sum) %>% filter(SSC %in% targets)
df_CR$Group <- "CR"

ITGB1_PD <- ITGB1_PD %>% 
  group_by(Cell) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()
df_PD <- aggregate(Spe ~ SSC + Cell, data = ITGB1_PD, sum) %>% filter(SSC %in% targets)
df_PD$Group <- "PD"

df_all <- rbind(df_PD, df_CR)
ggplot(df_all, aes(x = SSC, y = Spe, fill = Group)) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24)) +
  labs(x = "Target", y = "Specificity") +
  scale_fill_manual(values = c("blue", "red"), labels = c("PD", "CR")) +
  guides(fill = guide_legend(title = "Group"))
wilcox.test(Spe ~ Group, data = df_all %>% filter(SSC == "MYC"))
ggsave("./results/Mingze2025/supp/ITGB1_pathways_spe.PDF", width = 8, height = 12)

# FN1 distribution of most specific genes

XT_cell_PD <- filter(XT_results_PD, Ligand == "FN1")
XT_cell_PD$beta <- abs(XT_cell_PD$beta)
XT_cell_PD <- aggregate(XT_cell_PD$beta, list(XT_cell_PD$Cell, XT_cell_PD$SSC), sum)
colnames(XT_cell_PD) <- c("Cell","Target","Value")
XT_cell_PD <- tidyfst::df_mat(XT_cell_PD, row = Cell, col = Target, value = Value)
XT_cell_PD[is.na(XT_cell_PD)] <- 0
XT_cell_PD <- XT_cell_PD[rowSums(XT_cell_PD) > 1e-5,]
Dom_targets <- colnames(XT_cell_PD)[max.col(XT_cell_PD)]
topk <- 5
freq <- table(Dom_targets)
top5 <- names(sort(freq, decreasing = TRUE))[1:topk]
Dom_targets_collapsed <- ifelse(Dom_targets %in% top5, Dom_targets, "Others")
table(Dom_targets_collapsed)
names(Dom_targets_collapsed) <- names(Dom_targets) <- rownames(XT_cell_PD)
meta_PD$Spe_target <- "Others"
meta_PD[names(Dom_targets_collapsed),"Spe_target"] <- Dom_targets_collapsed
ggplot() + 
  geom_point(data = filter(meta_PD,Spe_target != "Others"), aes(x = x, y = y, color = Spe_target),size = 3) +
  geom_point(data = filter(meta_PD,Spe_target == "Others"), aes(x = x, y = y),size = 3,color = "grey") +
  # scale_color_brewer()+
  theme_void() +
  theme(text = element_text(family = "Arial",size = 24))+
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 12)))
ggsave("./results/Mingze2025/supp/Dominant_targets_PD.PDF", device = cairo_pdf,width = 12, height = 12)

XT_cell_CR <- filter(XT_results_CR, Ligand == "FN1")
XT_cell_CR$beta <- abs(XT_cell_CR$beta)
XT_cell_CR <- aggregate(XT_cell_CR$beta, list(XT_cell_CR$Cell, XT_cell_CR$SSC), sum)
colnames(XT_cell_CR) <- c("Cell","Target","Value")
XT_cell_CR <- tidyfst::df_mat(XT_cell_CR, row = Cell, col = Target, value = Value)
XT_cell_CR[is.na(XT_cell_CR)] <- 0
XT_cell_CR <- XT_cell_CR[rowSums(XT_cell_CR) > 1e-5,]
Dom_targets <- colnames(XT_cell_CR)[max.col(XT_cell_CR)]
topk <- 5
freq <- table(Dom_targets)
top5 <- names(sort(freq, decreasing = TRUE))[1:topk]
Dom_targets_collapsed <- ifelse(Dom_targets %in% top5, Dom_targets, "Others")
table(Dom_targets_collapsed)
names(Dom_targets_collapsed) <- names(Dom_targets) <- rownames(XT_cell_CR)
meta_CR$Spe_target <- "Others"
meta_CR[names(Dom_targets_collapsed),"Spe_target"] <- Dom_targets_collapsed
ggplot() + 
  geom_point(data = filter(meta_CR,Spe_target != "Others"), aes(x = x, y = y, color = Spe_target),size = 3) +
  geom_point(data = filter(meta_CR,Spe_target == "Others"), aes(x = x, y = y),size = 3,color = "grey") +
  # scale_color_brewer()+
  theme_void() +
  theme(text = element_text(family = "Arial",size = 24))+
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 12)))
ggsave("./results/Mingze2025/supp/Dominant_targets_CR.PDF", device = cairo_pdf,width = 12, height = 12)

# multi-layer
XT_filtered_PD <- aggregate(XT_results_PD$Weight, list(XT_results_PD$Ligand,XT_results_PD$Receptor,XT_results_PD$SSC), sum)
colnames(XT_filtered_PD) <- c("Ligand","Receptor","Target","Value")
XT_filtered_PD <- filter(XT_filtered_PD,Ligand %in% c("MIF","FN1")) %>% filter(Target %in% c("CTNNB1","CCND1","MYC"))
Plot_XTAlluvial(XT_filtered_PD,topk = 10, vertical = F)
ggsave("./results/Mingze2025/figure0109/Allu_PD.PDF", width = 12, height = 6)

XT_filtered_CR <- aggregate(XT_results_CR$Weight, list(XT_results_CR$Ligand,XT_results_CR$Receptor,XT_results_CR$SSC), sum)
colnames(XT_filtered_CR) <- c("Ligand","Receptor","Target","Value")
XT_filtered_CR <- filter(XT_filtered_CR,Ligand %in% c("MIF","FN1")) %>% filter(Target %in% c("CTNNB1","CCND1","MYC"))
Plot_XTAlluvial(XT_filtered_CR,topk = 10, vertical = F)
ggsave("./results/Mingze2025/figure0109/Allu_CR.PDF", width = 12, height = 6)

# Comparison: COL9A2-CD44, PD vs CR
XT_COL9A2_PD <- XT_results_PD %>%
  group_by(SSC) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  filter(Ligand == "COL9A2" & Receptor == "CD44")
XT_COL9A2_CR <- XT_results_CR %>%
  group_by(SSC) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  filter(Ligand == "COL9A2" & Receptor == "CD44")
wilcox.test(XT_COL9A2_PD$Fid, XT_COL9A2_CR$Fid)
df1 <- data.frame(Group = "PD", Fid = XT_COL9A2_PD$Fid)
df2 <- data.frame(Group = "CR", Fid = XT_COL9A2_CR$Fid)
df <- rbind(df1, df2)
ggplot(df, aes(x = Group, y = Fid, fill = Group)) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  theme(text = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24)) +
  labs(x = "", y = "Fidelity") +
  ylim(0, 0.00025) +
  scale_fill_manual(values = c("blue", "red"), labels = c("CR", "PD")) +
  theme(legend.position = "none")
ggsave("./results/Mingze2025/figure0109/COL9A2_CD44_Fid.PDF", width = 2.5, height = 3)

# Comparison: SPP1-CD44, PD vs CR
XT_SPP1_PD <- XT_results_PD %>%
  group_by(SSC) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  filter(Ligand == "SPP1" & Receptor == "CD44")
XT_SPP1_CR <- XT_results_CR %>%
  group_by(SSC) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  filter(Ligand == "SPP1" & Receptor == "CD44")
wilcox.test(XT_SPP1_PD$Fid, XT_SPP1_CR$Fid)
df1 <- data.frame(Group = "PD", Fid = XT_SPP1_PD$Fid)
df2 <- data.frame(Group = "CR", Fid = XT_SPP1_CR$Fid)
df <- rbind(df1, df2)
ggplot(df, aes(x = Group, y = Fid, fill = Group)) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  theme(text = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24)) +
  labs(x = "", y = "Fidelity") +
  ylim(0, 0.00025) +
  scale_fill_manual(values = c("blue", "red"), labels = c("CR", "PD")) +
  theme(legend.position = "none")
ggsave("./results/Mingze2025/figure0109/SPP1_CD44_Fid.PDF", width = 2.5, height = 3)

# comparison: TGFB1 and TNF, PD vs CR
# find the most dominant regulator signal for TGFB1 and TNF in PD and CR, respectively
XT_TT_PD <- XT_results_PD %>%
  group_by(SSC) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  filter(SSC %in% c("STAT3", "TWIST1"))
XT_TT_PD <- aggregate(Fid ~ SSC + Ligand, sum, data = XT_TT_PD) %>% 
  arrange(desc(Fid)) %>%
  slice_head(n = 30)
ligands <- unique(XT_TT_PD$Ligand)
ggplot(XT_TT_PD) +
  geom_point(aes(x = Ligand, y = SSC, color = Fid),size = 8) +
  scale_color_gradientn(colors = c("white", "blue"),limits = c(0,0.4)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 24,angle = 90,vjust = 0.5)) +
  labs(x = "Ligand", y = "") +
  theme(text = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24))
ggsave("./results/Mingze2025/figure0109/STAT3_TWIST1_PD.PDF", width = 12.5, height = 4)

XT_TT_CR <- XT_results_CR %>%
  group_by(SSC) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  filter(SSC %in% c("STAT3", "TWIST1")) %>%
  filter(Ligand %in% ligands)
XT_TT_CR <- aggregate(Fid ~ SSC + Ligand, sum, data = XT_TT_CR) 
ggplot(XT_TT_CR) +
  geom_point(aes(x = Ligand, y = SSC, color = Fid),size = 8) +
  scale_color_gradientn(colors = c("white", "blue"),limits = c(0,0.4)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 24,angle = 90,vjust = 0.5)) +
  labs(x = "Ligand", y = "") +
  theme(text = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24))
ggsave("./results/Mingze2025/figure0109/STAT3_TWIST1_CR.PDF", width = 12.5, height = 4)

# try alluvial package
XT_toplot <- XT_filtered_PD %>% arrange(desc(Value)) %>% slice_head(n = 10)
pdf("./results/Mingze2025/figure0109/Allu_PD.PDF", width = 8, height = 4)
alluvial::alluvial(XT_toplot[,1:3],
                  freq = XT_toplot$Value, border = NA, alpha = 0.5,
                  col = ifelse(XT_toplot$Ligand == "MIF", "pink", "skyblue"),
                  axis_labels = c("Ligand", "Receptor", "Target"))
dev.off()

XT_toplot <- XT_filtered_CR %>% arrange(desc(Value)) %>% slice_head(n = 10)
pdf("./results/Mingze2025/figure0109/Allu_CR.PDF", width = 8, height = 4)
alluvial::alluvial(XT_toplot[,1:3],
                  freq = XT_toplot$Value, border = NA, alpha = 0.6,
                  col = ifelse(XT_toplot$Ligand == "MIF", "blue", "yellow"),
                  axis_labels = c("Ligand", "Receptor", "Target"))
dev.off()

# 3. GRN-PATHWAY score
KEGG_paths <- readRDS("/home/shared/SpatialCrosstalk/databases/raw/KEGG_PPR_human.RDS")
KEGG_paths$Gene <- paste0(KEGG_paths$From,"=",KEGG_paths$To)
KEGG_paths <- KEGG_paths[,c("Gene","Pathway","Subcategory")]

# PD
GRN_KEGG_PD <- inner_join(GRN_filtered_PD,KEGG_paths,by = "Gene")
head(GRN_KEGG_PD)
GRN_KEGG_PD$beta <- abs(GRN_KEGG_PD$beta)
GRN_KEGG_PD <- GRN_KEGG_PD[,c("Cell","Gene","Pathway","Subcategory","beta")]
GRN_KEGG_PD <-  GRN_KEGG_PD %>%
  group_by(Cell) %>%
  mutate(Spe = beta/sum(beta)) %>%
  ungroup()
df <- aggregate(GRN_KEGG_PD$beta,list(GRN_KEGG_PD$Pathway), sum)
colnames(df) <- c("Pathway","Strength")
df$Strength <- df$Strength/sum(df$Strength)
df <- df[order(df$Strength, decreasing = T), ]
df <- df[1:20, ]
df <- df %>% mutate(Pathway = fct_reorder(Pathway, Strength))

ggplot(
    data = df,
    aes(
        x = Pathway, # First variable on the X-axis
        y = Strength, # Second variable on the X-axis
        fill = Pathway
    )) +
    coord_flip() +
    geom_bar(stat = "identity", color = "black") +
    theme_bw() +
    labs(y = "Relative activation", x = "Pathway") +
    theme(text = element_text(family = "Arial", size = 24)) +
    theme(legend.position = "none")
ggsave("./results/Mingze2025/figure0109/KEGG_PD.PDF", device = cairo_pdf,width = 18, height = 12)

# CR
GRN_KEGG_CR <- inner_join(GRN_filtered_CR,KEGG_paths,by = "Gene")
head(GRN_KEGG_CR)
GRN_KEGG_CR$beta <- abs(GRN_KEGG_CR$beta)
GRN_KEGG_CR <- GRN_KEGG_CR[,c("Cell","Gene","Pathway","Subcategory","beta")]
GRN_KEGG_CR <-  GRN_KEGG_CR %>%
  group_by(Cell) %>%
  mutate(Spe = beta/sum(beta)) %>%
  ungroup() 
df <- aggregate(GRN_KEGG_CR$beta,list(GRN_KEGG_CR$Pathway), sum)
colnames(df) <- c("Pathway","Strength")
df$Strength <- df$Strength/sum(df$Strength)
df <- df[order(df$Strength, decreasing = T), ]
df <- df[1:20, ]
df <- df %>% mutate(Pathway = fct_reorder(Pathway, Strength))

ggplot(
    data = df,
    aes(
        x = Pathway, # First variable on the X-axis
        y = Strength, # Second variable on the X-axis
        fill = Pathway
    )) +
    coord_flip() +
    geom_bar(stat = "identity", color = "black") +
    theme_bw() +
    labs(y = "Relative activation", x = "Pathway") +
    theme(text = element_text(family = "Arial", size = 24)) +
    theme(legend.position = "none")
ggsave("./results/Mingze2025/figure0109/KEGG_CR.PDF", device = cairo_pdf,width = 18, height = 12)

# Diff
df1 <- aggregate(GRN_KEGG_PD$Spe,list(GRN_KEGG_PD$Pathway), mean)
df2 <- aggregate(GRN_KEGG_CR$Spe,list(GRN_KEGG_CR$Pathway), mean)
colnames(df1) <- c("Pathway","Spe_PD")
colnames(df2) <- c("Pathway","Spe_CR")
df <- merge(df1,df2,by = "Pathway", all = T)
df[is.na(df)] <- 0
df$Diff <-  df$Spe_PD - df$Spe_CR
df <- df[order(df$Diff, decreasing = T), ]
df <- df[1:20, ]
df <- df %>% mutate(Pathway = fct_reorder(Pathway, Diff))

ggplot(
    data = df,
    aes(
        x = Pathway, # First variable on the X-axis
        y = Diff, # Second variable on the X-axis
        fill = Pathway
    )) +
    coord_flip() +
    geom_bar(stat = "identity", color = "black") +
    theme_bw() +
    labs(y = "Relative activation", x = "Pathway") +
    theme(text = element_text(size = 24)) +
    theme(legend.position = "none")
ggsave("./results/Mingze2025/figure0109/KEGG_diff_PD.PDF",width = 14, height = 12)

# Volcano plot of GRN
CCC_mat <- Aggr_CCC_results(CCC_origin_PD, flag = "Sender") # sender*LR
CCC_df_PD <- tidyfst::mat_df(CCC_mat)
colnames(CCC_df_PD) <- c("Cell","Gene","beta")
CCC_mat <- Aggr_CCC_results(CCC_origin_CR, flag = "Sender") # sender*LR
CCC_df_CR <- tidyfst::mat_df(CCC_mat)
colnames(CCC_df_CR) <- c("Cell","Gene","beta")
p <- Plot_GRNVolcano(
    CCC_df_PD,CCC_df_CR,
    sig_lfc_cut = 1,
    sig_p_cut = 0.05,
    label_lfc_cut = 5,
    label_fdr_cut = 0.005,
    cols = c(
        "Up-regulated" = "red",
        "Down-regulated" = "blue",
        "No change" = "grey80"
    )
)
p + ylim(0,300)
# CDH1-CDH1
ggsave("./results/Mingze2025/figure0109/CCC_volc.PDF", device = cairo_pdf,width = 12, height = 12)

features <- c("SPP1=ITGB1","COL9A2=CD44")
data <- Aggr_CCC_results(CCC_origin_PD, flag = "LR", LR_used = features[1])
p1 <- PlotSpa_SRVector(meta = meta_PD, data = data,threshold = .75)
# filen <- paste0("./results/Mingze2025/figure0109/PD_vec",features[1],".PDF")
# ggsave(filen, device = cairo_pdf,width = 8, height = 8)
data <- Aggr_CCC_results(CCC_origin_CR, flag = "LR", LR_used = features[1])
p2 <- PlotSpa_SRVector(meta = meta_CR, data = data,threshold = .75)
p1+p2
filen <- paste0("./results/Mingze2025/figure0109/vec",features[1],".PDF")
ggsave(filen, device = cairo_pdf,width = 16, height = 8)

data <- Aggr_CCC_results(CCC_origin_PD, flag = "LR", LR_used = features[2])
p1 <- PlotSpa_SRVector(meta = meta_PD, data = data,threshold = .75)
# filen <- paste0("./results/Mingze2025/figure0109/PD_vec",features[1],".PDF")
# ggsave(filen, device = cairo_pdf,width = 8, height = 8)
data <- Aggr_CCC_results(CCC_origin_CR, flag = "LR", LR_used = features[2])
p2 <- PlotSpa_SRVector(meta = meta_CR, data = data,threshold = .75)
p1+p2
filen <- paste0("./results/Mingze2025/figure0109/vec",features[2],".PDF")
ggsave(filen, device = cairo_pdf,width = 16, height = 8)

data <- Aggr_CCC_results(CCC_origin_PD, flag = "LR", LR_used = features[3])
p1 <- PlotSpa_SRVector(meta = meta_PD, data = data,threshold = .5)
# filen <- paste0("./results/Mingze2025/figure0109/PD_vec",features[1],".PDF")
# ggsave(filen, device = cairo_pdf,width = 8, height = 8)
data <- Aggr_CCC_results(CCC_origin_CR, flag = "LR", LR_used = features[3])
p2 <- PlotSpa_SRVector(meta = meta_CR, data = data,threshold = .5)
p1+p2
filen <- paste0("./results/Mingze2025/figure0109/vec",features[3],".PDF")
ggsave(filen, device = cairo_pdf,width = 16, height = 8)

data <- Aggr_CCC_results(CCC_origin_PD, flag = "LR", LR_used = features[4])
p1 <- PlotSpa_SRVector(meta = meta_PD, data = data,threshold = .75)
# filen <- paste0("./results/Mingze2025/figure0109/PD_vec",features[1],".PDF")
# ggsave(filen, device = cairo_pdf,width = 8, height = 8)
data <- Aggr_CCC_results(CCC_origin_CR, flag = "LR", LR_used = features[4])
p2 <- PlotSpa_SRVector(meta = meta_CR, data = data,threshold = .75)
p1+p2
filen <- paste0("./results/Mingze2025/figure0109/vec",features[4],".PDF")
ggsave(filen, device = cairo_pdf,width = 16, height = 8)

# differentially crosstalk
aggr_XT_PD <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_results_PD, sum) %>% 
  group_by(SSC) %>% 
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  group_by(Ligand) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup() %>%
  mutate(Pathway = paste0(Ligand,"_",Receptor,"_",SSC))
aggr_XT_CR <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_results_CR, sum) %>%
  group_by(SSC) %>% 
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup() %>%
  group_by(Ligand) %>%
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup() %>%
  mutate(Pathway = paste0(Ligand,"_",Receptor,"_",SSC))
diff_XT <- merge(aggr_XT_PD, aggr_XT_CR, by = "Pathway", all = T, suffixes = c("_PD","_CR")) %>%
  filter(!is.na(Weight_PD) | !is.na(Weight_CR)) %>%
  mutate(Weight_PD = ifelse(is.na(Weight_PD), 0, Weight_PD)) %>%
  mutate(Weight_CR = ifelse(is.na(Weight_CR), 0, Weight_CR)) %>%
  mutate(Diff_Weight = Fid_PD - Fid_CR)
diff_XT_used <- diff_XT[diff_XT$Spe_PD > 0.01 & diff_XT$Spe_CR > 0.01 & abs(diff_XT$Diff_Weight) > 0.1, ] %>%
  arrange(desc(abs(Diff_Weight))) 

lig_used <- c("COL9A2","SPP1","CDH1")
receptor_used <- c("CD44","ITGB1","CDH1")
SSC_used <- c("EGFR","TGFB1","TNF","ESR1")
XT_filtered_PD <- filter(XT_results_PD, Ligand %in% lig_used & Receptor %in% receptor_used & SSC %in% SSC_used)
XT_filtered_PD <- aggregate(list(XT_filtered_PD$Value, XT_filtered_PD$beta), 
  list(XT_filtered_PD$Ligand, XT_filtered_PD$Receptor, XT_filtered_PD$SSC), sum)
colnames(XT_filtered_PD) <- c("Ligand","Receptor","SSC","Value","beta")
XT_filtered_CR <- filter(XT_results_CR, Ligand %in% lig_used & Receptor %in% receptor_used & SSC %in% SSC_used)
XT_filtered_CR <- aggregate(list(XT_filtered_CR$Value, XT_filtered_CR$beta), 
  list(XT_filtered_CR$Ligand, XT_filtered_CR$Receptor, XT_filtered_CR$SSC), sum)
colnames(XT_filtered_CR) <- c("Ligand","Receptor","SSC","Value","beta")

color_use <- c("aquamarine","bisque","pink")
names(color_use) <- c("Ligand","Receptor","Target")

Plot_MultiLayer_Singlecell(XT_filtered_PD, start_trim = 0.3, end_trim = 0.3, arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,rect_height = 0.3, label_size = 12, layer_fill_colors = color_use)
ggsave(filename = "/home/shared/SpatialCrosstalk/results/Mingze2025/figure0109/multilayer_PD.PDF",width = 4200, height = 2400, units = "px")

Plot_MultiLayer_Singlecell(XT_filtered_CR, start_trim = 0.3, end_trim = 0.3, arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,rect_height = 0.3, label_size = 12, layer_fill_colors = color_use)
ggsave(filename = "/home/shared/SpatialCrosstalk/results/Mingze2025/figure0109/multilayer_CR.PDF",width = 3200, height = 2400, units = "px")

# COL9A2-CD44 vs SPP1-CD44

Gene <- "SPP1"
wilcox.test(exp_PD[Gene,],exp_CR[Gene,])
df <- meta_PD
df$value <- exp_PD[Gene,]
df <- df[order(df$value, decreasing = T), ]

ggplot()+
  geom_point(data = df, aes(x = x,y = y, color = value),size = 4)+
  scale_color_gradientn(colours = c("lightgray","blue"))+
  theme_void() +
  theme(text = element_text(family = "Arial",size = 24))+
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 12)))

df <- meta_CR
df$value <- exp_CR[Gene,]
df <- df[order(df$value, decreasing = T), ]

ggplot()+
  geom_point(data = df, aes(x = x,y = y, color = value),size = 4)+
  scale_color_gradientn(colours = c("white","blue"))+
  theme_void() +
  theme(text = element_text(family = "Arial",size = 24))+
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 12)))

path1 <- "SPP1=CD44"; path2 <- "COL9A2=CD44"
PlotSpa_Comparison2(CCC_mat_CR, path1, path2, meta_CR)
PlotSpa_Comparison2(CCC_mat_CR, path1, path2, meta_CR,mode = "other",pt_size = 3.5)

PlotSpa_Comparison2(CCC_mat_PD, path1, path2, meta_PD)
PlotSpa_Comparison2(CCC_mat_PD, path1, path2, meta_PD,mode = "other",pt_size = 3.5)

# Clustering of specificity
# Check the distribution of each target
Spe_PD <- aggregate(Weight ~ Cell + SSC, data = XT_results_PD, sum)
# Var_PD <- aggregate(Spe ~ SSC,data = Spe_PD,length)
Gene <- "CTNNB1"
data <- Spe_PD[Spe_PD$SSC == Gene, "Weight"] %>% as.vector() %>% unlist()
names(data) <- Spe_PD[Spe_PD$SSC == Gene, "Cell"] %>% as.vector() %>% unlist()
PlotSpa_Scatter(data, meta_PD)
data <- exp_PD[Gene, ]
PlotSpa_Scatter(data, meta_PD)

Plot_FidSpe(Aggr_XT_results = Aggr_XT_PD, feature = Gene, mode = "Fid")
Plot_FidSpe(Aggr_XT_results = Aggr_XT_CR, feature = Gene, mode = "Fid")

# Clustering & Regression
allgenes <- XT_results_PD$SSC %>% unique() %>% sort()
reg_results <- list()
Group_results_PD <- list()
for (Gene in allgenes){
  message(Gene)
  Group_result_PD <- Group_targets(XT_results_PD, feature = Gene, cluster_k = 4)
  reg_result <- Regression_target(XT_results = XT_results_PD, group_results = Group_result_PD, exp = exp_PD, feature = Gene)
  reg_results[[Gene]] <- reg_result
  Group_results_PD[[Gene]] <- Group_result_PD
}

# Use CR for regression
allgenes2 <- XT_results_CR$SSC %>% unique() %>% sort()
tg_spe <- aggregate(Spe ~ SSC+Cell, data = XT_results_CR, mean)
tg_spe <- tg_spe[!duplicated(tg_spe), ]
Z_new <- tg_spe$Spe
names(Z_new) <- tg_spe$SSC
intersect(allgenes, allgenes2)

pred_results <- list()
for (yf in allgenes2) {

  message(yf)
  spe_filtered <- tg_spe[tg_spe$SSC == yf, c("Cell","Spe")]
  cell_used <- spe_filtered$Cell %>% as.character()
  y_real <- exp_CR[yf, cell_used]
  z_val <- spe_filtered$Spe

  pred <- rep(NA_real_, length(y_real))
  chosen_interval <- rep(NA_integer_, length(y_real))


  # If all 4 models are NULL for this y, return NULL-like results
  if (all(sapply(reg_results[[yf]], is.null))) {
    pred_results[[yf]] <- data.frame(
      y_feature = yf,
      sample = cell_used,
      interval = NA_integer_,
      z = z_val,
      y_real = y_real,
      y_pred = NA_real_,
      error = NA_real_,
      abs_error = NA_real_,
      squared_error = NA_real_
    )
    next
  }

  for (s in seq_along(y_real)) {

    # intervals for this y
    ints <- Group_results_PD[[yf]]

    # first: find interval containing z
    k <- which(sapply(ints, function(int) {
      z_val[s] >= int$min_spe && z_val[s] <= int$max_spe
    }))

    if (length(k) > 0) {
      k <- k[1]
    } else {
      # if z is outside all intervals, choose closest interval by distance to interval range
      dists <- sapply(ints, function(int) {
        if (z_val[s] < int[["min_spe"]]) {
          int[["min_spe"]] - z_val[s]
        } else if (z_val[s] > int[["max_spe"]]) {
          z_val[s] - int[["max_spe"]]
        } else {
          0
        }
      })
      k <- which.min(dists)
    }

    # if selected model is NULL, choose closest interval among non-NULL models
    if (is.null(reg_results[[yf]][[k]])) {

      non_null_idx <- which(!sapply(reg_results[[yf]], is.null))

      dists_nonnull <- sapply(non_null_idx, function(kk) {
        int <- ints[[kk]]

        if (z_val[s] < int[["min_spe"]]) {
          int[["min_spe"]] - z_val[s]
        } else if (z_val[s] > int[["max_spe"]]) {
          z_val[s] - int[["max_spe"]]
        } else {
          0
        }
      })

      k <- non_null_idx[which.min(dists_nonnull)]
    }

    chosen_interval[s] <- k

    model <- reg_results[[yf]][[k]]
    receptors <- rownames(model$beta)
    X_new <- exp_CR[receptors, cell_used]
    x_s <- matrix(X_new[, s], nrow = 1)
    colnames(x_s) <- rownames(X_new)

    pred[s] <- as.numeric(
      predict(model, newx = x_s, s = "lambda.min")
    )
  }

  pred_results[[yf]] <- data.frame(
    y_feature = yf,
    sample = cell_used,
    interval = chosen_interval,
    z = z_val,
    y_real = y_real,
    y_pred = pred,
    error = y_real - pred,
    abs_error = abs(y_real - pred),
    squared_error = (y_real - pred)^2
  )
}

results_df <- do.call(rbind, pred_results)

summary_df <- results_df %>%
  group_by(y_feature) %>%
  summarise(
    n = sum(!is.na(y_pred)),
    RMSE = ifelse(n == 0, NA_real_, sqrt(mean(squared_error, na.rm = TRUE))),
    MAE  = ifelse(n == 0, NA_real_, mean(abs_error, na.rm = TRUE)),
    cor  = ifelse(n <= 1, NA_real_, cor(y_real, y_pred, use = "complete.obs")),
    .groups = "drop"
  )
