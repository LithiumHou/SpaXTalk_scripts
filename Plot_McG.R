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
library(scatterpie)
library(viridis)
library(Matrix)
library(ggalluvial)
library(RColorBrewer)
library(ggraph)
library(alluvial)
setwd("/home/shared/SpatialCrosstalk")
source("./codes/Process_results.R")
source("./codes/Visualization.R")

SeuratObj <- readRDS("./datasets/McGregor2025/McGregor_object.RDS")
meta_all <- SeuratObj@meta.data
meta_all <- meta_all[,c("x","y","annotation")]
clusters <- table(SeuratObj$CellAnnotation.Level0) %>% as.data.frame()
cluster_n <- clusters[clusters$Freq > 100, "Var1"] %>% as.character()
SeuratObj$annotation <- SeuratObj$CellAnnotation.Level0
table(SeuratObj$annotation)
target_type <- "Endothelium"

filen <- "./results/McGregor2025/McGregor_result_new.RDS"
allresults <- readRDS(filen)
CCC_origin <- allresults$CCC
GRN_filtered <- Process_GRN(allresults$GRN, beta_thres = 1e-5)
meta_all <- SeuratObj@meta.data[,c("x","y","annotation")]
CCC_mat <- Aggr_CCC_results(CCC_origin, flag = "Sender") # sender*LR
CCC_df <- tidyfst::mat_df(CCC_mat)
XT_results <- Aggr_XT_results(CCC_df, GRN_filtered, threshold = 0) 
XT_results <- XT_results %>% 
    group_by(Cell) %>% 
    mutate(Weight_new = Weight/sum(Weight)) %>% 
    ungroup() %>%
    group_by(Cell, SSC) %>%
    mutate(Fid = Weight/sum(Weight)) %>%
    ungroup()
table(XT_results$Receptor)

output_dir <- './results/McGregor2025/figure1124'
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
setwd(output_dir)

## Data overview
ggplot() +
  geom_point(data = meta_all, aes(x = x, y = y, color = annotation),size = 1) +
  # scale_color_brewer()+
  theme_classic() +
  theme(text = element_text(size = 24))+
  guides(color = guide_legend(override.aes = list(size = 12))) +
  labs(color = "Cell type")
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/data.PDF",width = 3600, height = 3000, units = "px")


#######################################
####### 1. Individual cell view #######
#######################################

# Select the cell with the most pathways
neighbor_annotation <- lapply(CCC_origin, function(x) {
    senders <- sub("=.*", "", colnames(x))
    temp <- meta_all[senders, "annotation"]
    temp
})
names(neighbor_annotation) <- names(CCC_origin)
str(neighbor_annotation)

homo_only <- which(sapply(neighbor_annotation, function(x) all(x == "Endothelium")))
mixed_low <- which(sapply(neighbor_annotation, function(x) length(unique(x)) > 1 & length(unique(x)) < 4 & any(x == "Endothelium")))
mixed_high <- which(sapply(neighbor_annotation, function(x) length(unique(x)) > 3 & any(x == "Endothelium")))
hetero_only <- which(sapply(neighbor_annotation, function(x) length(unique(x)) >= 1 & all(x != "Endothelium") ))

meta_all$category <- "Non_endothelial"
meta_all[names(homo_only),"category"] <- "Homo"
meta_all[names(mixed_low),"category"] <- "Mixed_low"
meta_all[names(mixed_high),"category"] <- "Mixed_high"
meta_all[names(hetero_only),"category"] <- "Hetero"

tempColors = c(brewer.pal(4, "Set3"),"#EEEEEE")

meta_all$category <- factor(meta_all$category, levels = c("Homo", "Hetero","Mixed_low", "Mixed_high", "Non_endothelial"))
ggplot(meta_all, aes(x, y)) +
  geom_point(
    data = subset(meta_all, category == "Non_endothelial"),
    color = "grey70",
    size = 0.6,
    alpha = 0.3
  ) +
  geom_point(
    data = subset(meta_all, category != "Non_endothelial"),
    aes(color = category),
    size = 1.5,
    alpha = 0.8
  ) +
  scale_color_brewer()+
  theme_bw() +
  scale_color_manual(values = tempColors) +
  theme(text = element_text(family = "Arial",size = 24))+
  guides(color = guide_legend(override.aes = list(size = 12)))
ggsave(filename = "./Category.PDF", device = cairo_pdf,width = 3200, height = 2400, units = "px")


data <- data.frame(category = c("Homo", "Hetero", "Mixed_low","Mixed_high"), count = c(length(homo_only), length(hetero_only), length(mixed_low), length(mixed_high)))
sum(data$count)

data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette="Set3") +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 18))+
  theme(legend.position = "none")
ggsave(filename = "./neighbor_info.PDF", device = cairo_pdf,width = 2000, height = 1800, units = "px")

dist_threshold <- 50
color_use <- c("aquamarine","bisque","pink")
names(color_use) <- c("Ligand","Receptor","Target")

# 1. homo_only
homo_used <- intersect(names(homo_only), XT_results$Cell)
cellnames <- (homo_used)[2]

CCC_mat <- CCC_origin[[cellnames]]
CCC_mat <- CCC_mat[rowSums(CCC_mat) > 0, ]

dim(CCC_mat)
XT_filtered <- XT_results[XT_results$Cell == cellnames, ]
XT_filtered <- XT_filtered[XT_filtered$Weight > 0.005, ]


Plot_MultiLayer_Singlecell(XT_filtered, start_trim = 0.3, end_trim = 0.3, 
  arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,
  rect_height = 0.3, label_size = 12,layer_fill_colors = color_use)
ggsave(filename = "./multilayer1.PDF", width = 4800, height = 2400, units = "px")

Plot_Signaling_Radial(
    CCC_mat,
    meta_all,
    receiver_name = cellnames,
    node_size    = 24,
    offset_step  = 0.03,
    edge_trim   = 0.075,
    arrow_len_mm = 3.5,
    width_range  = c(0.8, 3.5),
    node_palette = NULL
)
ggsave(filename = "./sender_arrow1.PDF", device = cairo_pdf,width = 3600, height = 3000, units = "px")

center <-  meta_all[cellnames, c("x", "y")]
dist <- sqrt((meta_all$x - center$x)^2 + (meta_all$y - center$y)^2)
neighbors <- meta_all[dist <= dist_threshold, ]
neighbors$iscenter <- rownames(neighbors) == cellnames
ggplot(neighbors, aes(x, y)) +
  geom_point(
    data = subset(neighbors, iscenter == F),
    color = "gray70",
    size = 5.6
  ) +
  geom_point(
    data = subset(neighbors, iscenter == T),
    color = "red",
    size = 9.5
  ) +
  scale_color_brewer()+
  theme_bw() +
  theme(text = element_text(size = 24),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 12)))
ggsave(filename = "./neighbor1.PDF", width = 6, height = 6)

# 2. lowly mixed
mixed_low_used <- intersect(names(mixed_low), XT_results$Cell)

cellnames <- (mixed_low_used)[2]
meta_all[cellnames,]
CCC_mat <- CCC_origin[[cellnames]]
CCC_mat <- CCC_mat[rowSums(CCC_mat) > 0, ]

dim(CCC_mat)
XT_filtered <- XT_results[XT_results$Cell == cellnames, ]
XT_filtered <- XT_filtered[XT_filtered$Weight > 0.0002, ]

Plot_MultiLayer_Singlecell(XT_filtered, start_trim = 0.3, end_trim = 0.3, arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,rect_height = 0.3, label_size = 12,layer_fill_colors = color_use)
ggsave(filename = "./multilayer2.PDF", width = 3600, height = 2400, units = "px")

Plot_Signaling_Radial(
    CCC_mat,
    meta_all,
    receiver_name = cellnames,
    node_size    = 24,
    offset_step  = 0.03,
    edge_trim   = 0.075,
    arrow_len_mm = 3.5,
    width_range  = c(0.8, 3.5),
    node_palette = NULL
)
ggsave(filename = "./sender_arrow2.PDF", device = cairo_pdf,width = 3600, height = 3000, units = "px")

center <-  meta_all[cellnames, c("x", "y")]
dist <- sqrt((meta_all$x - center$x)^2 + (meta_all$y - center$y)^2)
neighbors <- meta_all[dist <= dist_threshold, ]
neighbors$iscenter <- rownames(neighbors) == cellnames
ggplot(neighbors, aes(x, y)) +
  geom_point(
    data = subset(neighbors, iscenter == F),
    color = "gray70",
    size = 5.6
  ) +
  geom_point(
    data = subset(neighbors, iscenter == T),
    color = "red",
    size = 9.5
  ) +
  scale_color_brewer()+
  theme_bw() +
  theme(text = element_text(size = 24),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 12)))
ggsave(filename = "./neighbor2.PDF", width = 6, height = 6)

# 3. highly mixed
mixed_high_used <- intersect(names(mixed_high), XT_results$Cell)
cellnames <- (mixed_high_used)[3]
meta_all[cellnames,]
CCC_mat <- CCC_origin[[cellnames]]
CCC_mat <- CCC_mat[rowSums(CCC_mat) > 0, ]

dim(CCC_mat)
XT_filtered <- XT_results[XT_results$Cell == cellnames, ]
XT_filtered <- XT_filtered[XT_filtered$Weight > 0.07, ]
Plot_MultiLayer_Singlecell(XT_filtered, start_trim = 0.3, end_trim = 0.3, arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,rect_height = 0.3, label_size = 12,layer_fill_colors = color_use)
ggsave(filename = "./multilayer3.PDF", width = 3600, height = 2400, units = "px")

Plot_Signaling_Radial(
    CCC_mat,
    meta_all,
    receiver_name = cellnames,
    node_size    = 24,
    offset_step  = 0.03,
    edge_trim   = 0.1,
    arrow_len_mm = 3.5,
    width_range  = c(0.8, 3.5),
    node_palette = NULL
)
ggsave(filename = "./sender_arrow3.PDF", device = cairo_pdf,width = 3600, height = 3000, units = "px")

center <-  meta_all[cellnames, c("x", "y")]
dist <- sqrt((meta_all$x - center$x)^2 + (meta_all$y - center$y)^2)
neighbors <- meta_all[dist <= dist_threshold, ]
neighbors$iscenter <- rownames(neighbors) == cellnames
ggplot(neighbors, aes(x, y)) +
  geom_point(
    data = subset(neighbors, iscenter == F),
    color = "gray70",
    size = 5.6
  ) +
  geom_point(
    data = subset(neighbors, iscenter == T),
    color = "red",
    size = 9.5
  ) +
  scale_color_brewer()+
  theme_bw() +
  theme(text = element_text(size = 24),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 12)))
ggsave(filename = "./neighbor3.PDF", width = 6, height = 6)

##################################
####### 2. Cell group view #######
##################################

# Find the dominant signal and sender group

LR_mat <- Aggr_CCC_results(CCC_origin, flag = "Receiver")
senders_group <- meta_all[rownames(LR_mat),"annotation"]
LR_mat_group <- rowsum(as.matrix(LR_mat), group = senders_group)
LR_mat_group_top <- LR_mat_group[, order(colSums(LR_mat_group), decreasing = TRUE)[1:20]]
LR_mat_group_top <- LR_mat_group_top[rowSums(LR_mat_group_top) > 1, ]
LR_mat_group_top <- log10(LR_mat_group_top + 1)
df <- tidyfst::mat_df(LR_mat_group_top)
colnames(df) <- c("Sender","LR","Value")
df$LR <- gsub("="," - ",df$LR)
ggplot(df, aes(x = Sender, y = LR)) +
  geom_point(aes(color = Value),size = 5) +
  scale_color_gradientn(colours =  c("white", "#7bdbf8", "blue"),name = "CCC strength") +
  theme_bw() +
  theme(text = element_text(size = 24))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "bottom")+
  labs(x = "Sender cell type", y = "Ligand-receptor pair")
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/sender_LR.PDF",width = 16, height = 9)

CCC_mat <- Aggr_CCC_results(CCC_origin, flag = "Sender") # 
colnames(CCC_mat) <- gsub("="," - ",colnames(CCC_mat))
Plot_DominantLR(CCC_mat)
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/dominantLR2.PDF",width = 8, height = 11, device = cairo_pdf)

mat_fid <- CCC_mat/rowSums(CCC_mat)
colnames(mat_fid) <- gsub("="," \u2192 ",colnames(mat_fid))

mat_fid[is.na(mat_fid)] <- 0
max_fid_LR <- colnames(mat_fid)[max.col(mat_fid)]
top5 <- names(sort(table(max_fid_LR), decreasing = TRUE))[1:5]
max_fid_LR <- ifelse(max_fid_LR %in% top5,max_fid_LR,"Others")
names(max_fid_LR) <- rownames(mat_fid)
max_fid_LR <- factor(max_fid_LR, levels = c(top5,"Others"))
df <- meta_all
df$value <- "0"
df[names(max_fid_LR),"value"] <- max_fid_LR

tempColors = c(brewer.pal(5, "Set3"),"#EEEEEE")

ggplot(df, aes(x, y)) +
  geom_point(
    data = subset(df, value == "0"),
    color = "grey70",
    size = 0.6,
    alpha = 0.3
  ) +
  geom_point(
    data = subset(df, value != "0"),
    aes(color = value),
    size = 1.5,
    alpha = 0.8
  ) +
  scale_color_brewer()+
  theme_bw() +
  scale_color_manual(values = tempColors) +
  theme(text = element_text(family = "Arial",size = 24))+
  guides(color = guide_legend(override.aes = list(size = 12)))

# CCC
Ls <- sub("=.*", "", colnames(CCC_mat))
Rs <- sub(".*=", "", colnames(CCC_mat))
df <- data.frame(Ligand = Ls, Receptor = Rs, Weight = colSums(CCC_mat))

SigXTalkR::PlotCCI_CirclePlot(df, topk = 10)

# most specific target
GRN_used <- GRN_filtered[GRN_filtered$Receptor == "ITGA1", ]
XT_used <- XT_results[XT_results$Receptor == "PECAM1", ]
SSC_order <- XT_used %>%
  group_by(SSC) %>%
  summarise(mean_weight = median(Weight)) %>%
  arrange(mean_weight) %>%
  pull(SSC)
XT_used$SSC <- factor(XT_used$SSC, levels = SSC_order)
ggplot(XT_used, aes(x = SSC, y = Weight,fill = SSC)) +
  geom_boxplot() +
  theme_bw() +
  coord_flip() +
  theme(text = element_text(size = 24))+
  labs(x = "Target", y = "Pathway regulatory strength") +
  theme(legend.position = "none")
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/COL4A1_boxplot.PDF",width = 8, height = 14)



# 
XT_results <- XT_results %>%
  group_by(Cell, SSC) %>%
  mutate(Fid = Weight/sum(Weight)) %>%
  ungroup()
XT_aggr <- aggregate(Weight ~ Ligand + Receptor + SSC, data = XT_results, sum) 
XT_toplot <- XT_aggr %>% arrange(desc(Weight)) %>% slice_head(n=10)
uf <- unique(XT_toplot$Ligand) %>% sort()
cols <- paletteer::paletteer_c("grDevices::Set 3", length(uf))
col_map <- setNames(cols, uf)
feature_colors <- col_map[XT_toplot$Ligand]
pdf("./Allu_all.PDF", width = 6, height = 5.5)
alluvial::alluvial(XT_toplot[,1:3],gap.width = 0.02,xw = 0.2,cw = 0.2,
                  freq = XT_toplot$Weight, border = NA, alpha = 0.75,
                  col = feature_colors
                  )
dev.off()

# Plot_XTAlluvial(XT_aggr, topk = 15,vertical = F)
# ggsave(filename = "./Allu_all.PDF", device = cairo_pdf,width = 3600, height = 3000, units = "px")

# comparison
path1 <- "LAMC1=ITGB4"; path2 <- "THBS1=ITGB3"
PlotSpa_Comparison(XT_results, path1, path2, meta_all, path_flag = "LR",value_flag = "CCC")
filen <- paste0(path1,"_vs_",path2,".PDF")
ggsave(filename = filen, width = 3200, height = 3000, units = "px")

XT1 <- XT_results[XT_results$Ligand == "LAMC1" & XT_results$Receptor == "ITGB4" , ]
XT2 <- XT_results[XT_results$Ligand == "THBS1" & XT_results$Receptor == "ITGB3" , ]
mean(XT1$Fid);  mean(XT2$Fid)
wilcox.test(XT1$Fid, XT2$Fid)
df1 <- data.frame(Strength = XT1$Fid, group = "LAMC1")
df2 <- data.frame(Strength = XT2$Fid, group = "THBS1")
data <- rbind(df1,df2)
data %>% ggplot( aes(x = group, y=Strength, fill=group)) +
    geom_boxplot(alpha=0.6,outliers = F) +
    viridis::scale_fill_viridis(discrete=TRUE, option = "C") +
    ylim(c(0,1.1)) +
    theme_bw() +
    theme(text = element_text(size = 24))+
    theme(
      legend.position="none"
    )
ggsave(filename = "./Fid_comparison_LAMC1_vs_THBS1.PDF", width = 1000, height = 1000, units = "px")


# track the senders
CCC_sender <- Aggr_CCC_results(CCC_origin, flag = "Receiver")
CCC_group <- meta_all[rownames(CCC_sender),"annotation"]
CCC_sender <- rowsum(CCC_sender, group = CCC_group)
rowSums(CCC_sender)/sum(CCC_sender)*100
CCC_sender <- tidyfst::mat_df(CCC_sender)
colnames(CCC_sender) <- c("Type","LR","Value")
CCC_sender <- CCC_sender[CCC_sender$Value > 0, ]
CCC_sender$Ligand <- sub("=.*", "", CCC_sender$LR)
CCC_sender$Receptor <- sub(".*=", "", CCC_sender$LR)
Plot_XTAlluvial(CCC_sender[,c("Type","Ligand","Receptor","Value")], topk = 10,vertical = F)
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/CCC_sender.PDF",device = cairo_pdf,width = 3600, height = 3000, units = "px")

df <- aggregate(Value ~ Type, data = CCC_sender, sum)
df$Type <- as.character(df$Type)
data <- df %>%
  group_by(Type) %>%
  summarise(Value = sum(Value), .groups = "drop") %>%
  arrange(desc(Value)) %>%
  mutate(
    Type = ifelse(row_number() <= 4, Type, "Others")
  ) %>%
  group_by(Type) %>%
  summarise(Value = sum(Value), .groups = "drop") %>%
  arrange(desc(Value)) %>%
  mutate(
    fraction = Value / sum(Value),
    label = paste0(scales::percent(fraction))

  )
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$Type,"\n",data$label)

ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Type)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette="Pastel2") +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")
ggsave(filename = "./Sender_all.PDF", width = 2700, height = 2400, units = "px")

data$Type <- factor(data$Type, levels = c("Endothelium","Fibroblasts","Pericytes","Muscularis","Others") %>% rev())
ggplot(data) +
  geom_bar(stat="identity", aes(x = fraction, y = Type, fill = Type))+
  theme_classic() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 24))
ggsave(filename = "./Sender_all.PDF", width = 2700, height = 2400, units = "px")

# Supp figures
# 1. specificity of TFs
XT_aggr <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_results, sum) %>%
  group_by(Ligand) %>% 
  mutate(Spe = Weight/sum(Weight)) %>%
  ungroup()
df <- aggregate(Spe ~ SSC, data = XT_aggr, sum)
df %>% arrange(-Spe) %>%
  slice_head(n = 15) %>% 
  mutate(SSC = fct_reorder(SSC, (Spe))) %>%
  ggplot(aes(x = SSC, y = Spe, fill = SSC)) +
    geom_bar(stat="identity") +
    coord_flip() +
    theme_bw() +
    labs(y = "Total specificity") +
    theme(text = element_text(size = 24))+
    theme(
      legend.position="none"
    )
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/Spe_TF.PDF",width = 1800, height = 2400, units = "px")

# 2. Spe distribution of CTNNB1, ATF4
Tar_aggr <- aggregate(Weight ~ SSC + Cell, data = XT_results, sum) %>%
    group_by(Cell) %>% 
    mutate(Spe = Weight/sum(Weight)) %>%
    ungroup()
data <- Tar_aggr[Tar_aggr$SSC == "CTNNB1", ]
data$Cell <- data$Cell %>% as.character()
df <- meta_all
df$value <- 0
df[data$Cell,"value"] <- data$Spe
ggplot() +
    geom_point(
        data = df[df$value == 0,],
        aes(x = x, y = y),
        color = "lightgray",
        alpha = 0.4,
        size = 2
    ) +
    ggnewscale::new_scale("color") +   
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
    coord_fixed() + 
    # Theme
    theme_void() +
    theme(text = element_text(size = 24),legend.position = "right") +
    guides(color = guide_colorbar(order = 1))
filen <- "/home/shared/SpatialCrosstalk/results/McGregor2025/supp/Spe_CTNNB1.PDF"
ggsave(filename = filen, width = 15, height = 15)

# 3. PECAM1 downstream
XT_used <- XT_results[XT_results$Receptor == "PECAM1", ]
df <- aggregate(Weight ~ SSC, data = XT_used, sum)
df$Spe <- df$Weight/sum(df$Weight)
df %>% arrange(-Spe) %>%
  slice_head(n = 15) %>% 
  mutate(SSC = fct_reorder(SSC, (Spe))) %>%
  ggplot(aes(x = SSC, y = Spe, fill = SSC)) +
    geom_bar(stat="identity") +
    coord_flip() +
    theme_bw() +
    labs(y = "Total specificity") +
    theme(text = element_text(size = 24))+
    theme(
      legend.position="none"
    )
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/Spe_PECAM1.PDF",width = 1200, height = 1800, units = "px")

# 4. PECAM fidelity
XT1_aggr <- XT_aggr[XT_aggr$SSC == "AKT2", ]
XT1_aggr <- XT1_aggr %>% 
    arrange(Weight) %>% 
    slice_tail(n = 20)
XT1_aggr[XT1_aggr$Ligand == XT1_aggr$Receptor, "Receptor"] <- paste0(unlist(XT1_aggr[XT1_aggr$Ligand == XT1_aggr$Receptor, "Receptor"]), 1)
Plot_EdgeBundle(XT1_aggr[,c(1,2,4)],text_size = 8, edge_colors = c("yellow","blue"),node_size = 7, edge_width = 2)
ggsave("/home/shared/SpatialCrosstalk/results/McGregor2025/supp/Fid_AKT2.PDF",width = 10, height = 10)

