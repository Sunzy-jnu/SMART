setwd("D:/Rshuju/shiyan")
source("D:/Rshuju/shiyan/SMART/utils/utils.R")
source("D:/Rshuju/shiyan/SMART/utils/plot.R")
library(ggplot2)
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(ggpubr)
library(multcomp)
library(ggprism)
library(ggplotify)
library(scatterpie)
library(RColorBrewer)
library(showtext)
library(Seurat)



#################### Preprocess the original scRNA-seq data

### The original scRNA-seq data can be downloaded from https://www.spatialresearch.org/resources-published-datasets/doi-10-1016-j-cell-2019-11-025/.

sc_counts <- read.table("D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_scRNA-seq/all_cells_count_matrix_filtered.tsv", header = TRUE, row.names = 1, sep = "\t")
sc_metadata <- read.table("D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_scRNA-seq/all_cells_meta_data_filtered.tsv", header = TRUE, sep = "\n")

new_metadata <- c()
i <- 1
temp <- c()
while (i <= nrow(sc_metadata)){
  if (grepl(".1|.2", sc_metadata[i,])){
    new_metadata <- c(new_metadata, temp)
    temp <- sc_metadata[i,]
  } else{
    temp <- paste(temp, sc_metadata[i,], sep = "")
  }
  i <- i + 1
}
new_metadata <- c(new_metadata, temp)

new_metadata <- t(as.data.frame(strsplit(new_metadata, split = "\t")))
rownames(new_metadata) <- new_metadata[, 1]
new_metadata <- new_metadata[, -1]
colnames(new_metadata) <- c("nGene", "nUMI", "experiment", "Phase", "res.0.7", "celltype", "state")

identical(colnames(sc_counts), rownames(new_metadata))
identical(as.integer(colSums(sc_counts)), as.integer(new_metadata[, "nUMI"]))
identical(as.integer(colSums(sc_counts > 0)), as.integer(new_metadata[, "nGene"]))

table(new_metadata[, "celltype"])

sc_labels <- new_metadata[, "celltype"]

keep_cells <- names(sc_labels)[-which(sc_labels %in% c("Erythrocytes", "Immune cells"))]

sc_counts <- t(sc_counts[, keep_cells])
sc_labels <- sc_labels[keep_cells]

sc_labels[grep("Atrial cardiomyocytes", sc_labels)] <- "Atrial_cardiomyocytes"
sc_labels[grep("Capillary endothelium", sc_labels)] <- "Capillary_endothelium"
sc_labels[grep("Cardiac neural crest cells & Schwann progenitor cells", sc_labels)] <- "Cardiac_neural_crest_Schwann_progenitor"
sc_labels[grep("Endothelium / pericytes / adventitia", sc_labels)] <- "Endothelium_pericytes_adventitia"
sc_labels[grep("Epicardial cells", sc_labels)] <- "Epicardial_cells"
sc_labels[grep("Epicardium-derived cells", sc_labels)] <- "Epicardium_derived_cells"
sc_labels[grep("related to cardiac skeleton connective tissue", sc_labels)] <- "Fibroblast_like_1"
sc_labels[grep("related to larger vascular development", sc_labels)] <- "Fibroblast_like_2"
sc_labels[grep("related to smaller vascular development", sc_labels)] <- "Fibroblast_like_3"
sc_labels[grep("Myoz2-enriched cardiomyocytes", sc_labels)] <- "Myoz2_enriched_cardiomyocytes"
sc_labels[grep("Smooth muscle cells / fibroblast-like)", sc_labels)] <- "Smooth_muscle_cells_fibroblast_like"
sc_labels[grep("Ventricular cardiomyocytes", sc_labels)] <- "Ventricular_cardiomyocytes"


write.csv(sc_counts, file = "D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_scRNA-seq/human_heart/sc_data/sc_counts.csv")
sc_labels_m <- as.matrix(sc_labels)
colnames(sc_labels_m) <- "celltype"
write.csv(sc_labels_m, file = "D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_scRNA-seq/human_heart/sc_data/sc_labels.csv")

write.csv(new_metadata, file = "D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_scRNA-seq/human_heart/sc_data/sc_metadata.csv")



#################### Preprocess the original spatial transcriptomics data

### The original spatial transcriptomics data can be downloaded from https://www.spatialresearch.org/resources-published-datasets/doi-10-1016-j-cell-2019-11-025/.

st_counts <- read.table("D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_ST/ST/filtered_matrix.tsv", header = TRUE, row.names = 1, sep = "\t")

library(clusterProfiler)
library(org.Hs.eg.db)
gene_symbol <- bitr(substr(rownames(st_counts), 1, 15), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
gene_symbol <- gene_symbol[!duplicated(gene_symbol$ENSEMBL),]
gene_symbol <- gene_symbol[!duplicated(gene_symbol$SYMBOL),]

rownames(st_counts) <- substr(rownames(st_counts), 1, 15)
st_counts <- st_counts[gene_symbol$ENSEMBL,]
rownames(st_counts) <- gene_symbol$SYMBOL

st_location <- as.data.frame(Reduce(rbind, strsplit(substring(colnames(st_counts), 2), split = "x")))
rownames(st_location) <- colnames(st_counts)
colnames(st_location) <- c("id", "x", "y")

for (i in unique(st_location$id)){
  idx_i <- which(st_location$id == i)
  st_counts_i <- t(as.matrix(st_counts[, idx_i]))
  st_location_i <- st_location[idx_i,]
  save_dir <- paste("D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_ST/ST/st_data", paste("sample", i, sep = "_"), sep = "/")
  if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
  }
  write.csv(st_counts_i, file = paste(save_dir, "st_counts.csv", sep = "/"))
  write.csv(st_location_i, file = paste(save_dir, "st_location.csv", sep = "/"))
}



#################### Run SMART

### SMART package can be installed by install.packages("Code/SMART_1.1.0.tar.gz", repos = NULL, type = "source").
### SMART package is also available at https://github.com/Zhangxf-ccnu/SMART.

library(SMART)

sc_counts_path <- "D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_scRNA-seq/human_heart/sc_data/sc_counts.csv"
sc_labels_path <- "D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_scRNA-seq/human_heart/sc_data/sc_labels.csv"

sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels

sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10, n_cores = 6, save_dir = "D:/Rshuju/shiyan/SMART_results/human_heart/sc_results")


st_counts_all <- list()
st_location_all <- list()
st_results_all <- list()

for (i in list.files("D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_ST/ST/st_data")){
  st_counts_path <- paste(paste("D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_ST/ST/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
  st_location_path <- paste(paste("D:/Rshuju/shiyan/SMART_data/human_heart/heart_filtered_ST/ST/st_data", i, sep = "/"), "st_location.csv", sep = "/")
  st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
  st_location <- read.csv(st_location_path, header = TRUE, row.names = 1)
  spot_count <- nrow(st_counts)

  st_results <- st_train(st_counts, sc_results = sc_results, spatial_coords = as.matrix(st_location), lambda_spatial = 0.8, neighbor_radius = 5, save_dir = paste("D:/Rshuju/shiyan/SMART_results/human_heart/st_results", i, sep = "/"))
  st_counts_all[[i]] <- st_counts
  st_location_all[[i]] <- st_location
  st_results_all[[i]] <- st_results
}

save(sc_counts, sc_labels, sc_results, st_counts_all, st_location_all, st_results_all, file = "D:/Rshuju/shiyan/SMART_results/human_heart/st_results/HEART.Rdata")





#################### Analysis

output_dir <- "D:/Rshuju/shiyan/SMART_results/human_heart/Analysis"
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}


### Load results
load("D:/Rshuju/shiyan/SMART_results/human_heart/st_results/HEART.Rdata")

analysis_sections <- paste("sample", 1:19, sep = "_")
analysis_sections_5_PCW <- paste("sample", 1:4, sep = "_")
analysis_sections_6_PCW <- paste("sample", 5:13, sep = "_")
analysis_sections_9_PCW <- paste("sample", 14:19, sep = "_")
analysis_sections_partial <- paste("sample", c(1, 9, 15), sep = "_")

st_counts_all <- st_counts_all[analysis_sections]
st_location_all <- st_location_all[analysis_sections]
st_results_all <- st_results_all[analysis_sections]

for (i in 5:19){
  temp <- paste("sample", i, sep = "_")
  st_location_all[[temp]][, "y"] <- max(st_location_all[[temp]][, "y"]) + min(st_location_all[[temp]][, "y"]) - st_location_all[[temp]][, "y"]
}

st_location_lim <- c()
for (i in analysis_sections){
  st_location_i <- st_location_all[[i]]
  st_location_lim <- rbind(st_location_lim, c(min(st_location_i[,"x"]), max(st_location_i[,"x"]), min(st_location_i[,"y"]), max(st_location_i[,"y"])))
}
colnames(st_location_lim) <- c("x_min", "x_max", "y_min", "y_max")
rownames(st_location_lim) <- analysis_sections

x_range <- max(st_location_lim[, "x_max"] - st_location_lim[, "x_min"]) + 1
y_range <- max(st_location_lim[, "y_max"] - st_location_lim[, "y_min"]) + 1

x_center <- (1 + x_range) / 2
y_center <- (1 + y_range) / 2

for (i in analysis_sections){
  x_shift <- (st_location_lim[i, "x_min"] + st_location_lim[i, "x_max"]) / 2 - x_center
  y_shift <- (st_location_lim[i, "y_min"] + st_location_lim[i, "y_max"]) / 2 - y_center
  st_location_all[[i]][, "x"] <- st_location_all[[i]][, "x"] - x_shift
  st_location_all[[i]][, "y"] <- st_location_all[[i]][, "y"] - y_shift
}



########## Analysis of cell type proportions

proportion_all <- lapply(analysis_sections, function(x){st_results_all[[x]]$proportion})
names(proportion_all) <- analysis_sections
unique_labels <- sort(colnames(proportion_all[[1]]))
proportion_all <- lapply(proportion_all, function(x){x[, unique_labels]})
use_color <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_labels))
names(use_color) <- unique_labels

use_color_section <- scales::hue_pal()(length(analysis_sections))
names(use_color_section) <- analysis_sections

use_color_section_stage <- rep(scales::hue_pal()(3), times = c(4, 9, 6))
names(use_color_section_stage) <- analysis_sections

use_color_stage <- scales::hue_pal()(3)
names(use_color_stage) <- c("5_PCW", "6_PCW", "9_PCW")

sample_stage <- c(rep("5_PCW", 4), rep("6_PCW", 9), rep("9_PCW", 6))
names(sample_stage) <- analysis_sections


### Plot the scatterpie
plot_scatterpie_2(proportion_all, st_location_all, img_list = NULL, mname = analysis_sections, use_color = use_color, fname = "proportion_scatterpie_list", save_dir = output_dir, n_cols = 5, pie_r = 0.5, width = 15, height = 15)
plot_scatterpie_2(proportion_all[analysis_sections_partial], st_location_all[analysis_sections_partial], img_list = NULL, mname = analysis_sections_partial, use_color = use_color, fname = "proportion_scatterpie_list_partial", save_dir = output_dir, n_cols = 3, pie_r = 0.5, width = 10, height = 7)


### Plot the heatmap for each cell type
plot_scatter_heatmap_2(proportion_all[analysis_sections_5_PCW], st_location_all[analysis_sections_5_PCW], img_list = NULL, mname = analysis_sections_5_PCW, feature = unique_labels, scale = TRUE, fname = "scaled_proportion_heatmap_list_5_PCW", save_dir = output_dir, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 30, height = 20)
plot_scatter_heatmap_2(proportion_all[analysis_sections_6_PCW], st_location_all[analysis_sections_6_PCW], img_list = NULL, mname = analysis_sections_6_PCW, feature = unique_labels, scale = TRUE, fname = "scaled_proportion_heatmap_list_6_PCW", save_dir = output_dir, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 30, height = 40)
plot_scatter_heatmap_2(proportion_all[analysis_sections_9_PCW], st_location_all[analysis_sections_9_PCW], img_list = NULL, mname = analysis_sections_9_PCW, feature = unique_labels, scale = TRUE, fname = "scaled_proportion_heatmap_list_9_PCW", save_dir = output_dir, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 30, height = 30)

plot_scatter_heatmap_2(proportion_all[analysis_sections_partial], st_location_all[analysis_sections_partial], img_list = NULL, mname = analysis_sections_partial, feature = c("Ventricular_cardiomyocytes", "Atrial_cardiomyocytes", "Smooth_muscle_cells_fibroblast_like"), scale = TRUE, fname = "scaled_proportion_heatmap_list_partial", save_dir = output_dir, feat_by_sect = TRUE, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 11, height = 11)


### Find marker genes for each cell type based on scRNA-seq reference
marker_genes <- find_markers(t(sc_counts[, intersect(colnames(sc_counts), Reduce("intersect", lapply(st_counts_all, colnames)))]), sc_labels, n_markers = NULL)
marker_genes <- marker_genes[unique_labels]
saveRDS(marker_genes, file = paste(output_dir, "marker_genes.rds", sep = "/"))
write.csv(as.data.frame(lapply(marker_genes, function(x){x[1:5]})), file = paste(output_dir, "marker_genes.csv", sep = "/"))


### Plot the gene expression of markers for each cell type
st_counts_norm_all <- lapply(st_counts_all, function(x){x / rowSums(x)})
st_counts_norm_metagene_all <- lapply(st_counts_norm_all, function(x){
  temp <- sapply(marker_genes, function(y){log2(rowMeans(x[, y[1:5]]) * 1e4 + 1)})
  temp
})
plot_scatter_heatmap_2(st_counts_norm_metagene_all[analysis_sections_5_PCW], st_location_all[analysis_sections_5_PCW], img_list = NULL, mname = analysis_sections_5_PCW, feature = unique_labels, scale = TRUE, fname = "scaled_metagene_exp_heatmap_list_5_PCW", save_dir = output_dir, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 30, height = 20)
plot_scatter_heatmap_2(st_counts_norm_metagene_all[analysis_sections_6_PCW], st_location_all[analysis_sections_6_PCW], img_list = NULL, mname = analysis_sections_6_PCW, feature = unique_labels, scale = TRUE, fname = "scaled_metagene_exp_heatmap_list_6_PCW", save_dir = output_dir, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 30, height = 40)
plot_scatter_heatmap_2(st_counts_norm_metagene_all[analysis_sections_9_PCW], st_location_all[analysis_sections_9_PCW], img_list = NULL, mname = analysis_sections_9_PCW, feature = unique_labels, scale = TRUE, fname = "scaled_metagene_exp_heatmap_list_9_PCW", save_dir = output_dir, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 30, height = 30)
plot_scatter_heatmap_2(st_counts_norm_metagene_all[analysis_sections_partial], st_location_all[analysis_sections_partial], img_list = NULL, mname = analysis_sections_partial, feature = c("Ventricular_cardiomyocytes", "Atrial_cardiomyocytes", "Smooth_muscle_cells_fibroblast_like"), scale = TRUE, fname = "scaled_metagene_exp_heatmap_list_partial", save_dir = output_dir, feat_by_sect = TRUE, point_size = 3.8, x_lim = c(1, x_range), y_lim = c(1, y_range), width = 11, height = 11)



### Plot the colocalization of cell types
if (!dir.exists(paste(output_dir, "colocalization", sep = "/"))){
  dir.create(paste(output_dir, "colocalization", sep = "/"), recursive = TRUE)
}

colocalization_list <- list()
for (object in analysis_sections){
  colocalization <- cor(proportion_all[[object]], method = "pearson")
  fname <- paste(object, "colocalization_prop_cor", sep = "_")
  pdf(paste(paste(output_dir, "colocalization", sep = "/"), paste(fname, "pdf", sep = "."), sep = "/"), width = 8, height = 8)
  fig <- pheatmap::pheatmap(colocalization, cellwidth = 20, cellheight = 20, main = fname, cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 90, clustering_distance_rows = as.dist(1 - colocalization), clustering_distance_cols = as.dist(1 - colocalization))
  print(fig)
  dev.off()
  colocalization_list[[object]] <- colocalization
}


fig_list <- list()
for (i in seq(length(analysis_sections))){
  colocalization <- cor(proportion_all[[analysis_sections[i]]], method = "pearson")
  fig_list[[i]] <- as.grob(pheatmap::pheatmap(colocalization, cellwidth = 20, cellheight = 20, main = analysis_sections[i], cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 90, clustering_distance_rows = as.dist(1 - colocalization), clustering_distance_cols = as.dist(1 - colocalization), silent = TRUE))
}

fname <- "colocalization_prop_cor"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 28, height = 35)
fig <- ggarrange(plotlist = fig_list, ncol = 4, nrow = 5)
print(fig)
dev.off()

fname <- "colocalization_prop_cor_5_PCW"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 14, height = 14)
fig <- ggarrange(plotlist = fig_list[1:4], ncol = 2, nrow = 2)
print(fig)
dev.off()

fname <- "colocalization_prop_cor_6_PCW"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 14, height = 35)
fig <- ggarrange(plotlist = fig_list[5:13], ncol = 2, nrow = 5)
print(fig)
dev.off()

fname <- "colocalization_prop_cor_9_PCW"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 14, height = 21)
fig <- ggarrange(plotlist = fig_list[14:19], ncol = 2, nrow = 3)
print(fig)
dev.off()



colocalization_stage <- lapply(c("5_PCW", "6_PCW", "9_PCW"), function(x){Reduce("+", colocalization_list[get(paste("analysis_sections", x, sep = "_"))]) / length(get(paste("analysis_sections", x, sep = "_")))})
names(colocalization_stage) <- c("5_PCW", "6_PCW", "9_PCW")

fname <- "colocalization_prop_cor_stage"
fig_list <- list()
for (i in c("5_PCW", "6_PCW", "9_PCW")){
  colocalization <- colocalization_stage[[i]]
  fig_list[[i]] <- as.grob(pheatmap::pheatmap(colocalization, cellwidth = 20, cellheight = 20, main = i, cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 90, clustering_distance_rows = as.dist(1 - colocalization), clustering_distance_cols = as.dist(1 - colocalization), silent = TRUE))
}
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 14, height = 14)
fig <- ggarrange(plotlist = fig_list, ncol = 2, nrow = 2)
print(fig)
dev.off()





########## The difference of cell type proportions among three stages

### Plot the overall composition of cell types across different samples
overall_proportion_all <- t(sapply(proportion_all, colMeans))

plot_df <- as.data.frame(cbind(rownames(overall_proportion_all), overall_proportion_all))
colnames(plot_df) <- c("sample", colnames(overall_proportion_all))
plot_df <- reshape2::melt(plot_df, id.vars = c("sample"), variable.name = "label", value.name = "prop")
plot_df$prop <- as.numeric(plot_df$prop)
plot_df$sample <- factor(plot_df$sample, levels = analysis_sections)
plot_df$label <- factor(plot_df$label, levels = unique_labels)
plot_df$stage <- factor(sample_stage[plot_df$sample], levels = c("5_PCW", "6_PCW", "9_PCW"))

fname <- "overall_proportion_by_sample_barplot"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 12, height = 5)
fig <- ggbarplot(plot_df, x = "sample", y = "prop", fill = "label", color = "label", width = 0.7, ylim = c(0, 1)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "right") +
  labs(title = fname) +
  scale_fill_manual(values = use_color) +
  scale_color_manual(values = use_color)
print(fig)
dev.off()


fname <- "overall_proportion_by_stage_boxplot"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 12, height = 8)
fig <- ggplot(data = plot_df, aes(x = stage, y = as.numeric(prop), color = stage)) +
  geom_boxplot(size = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  labs(x = "", y = "", title = fname) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background  = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "bottom") +
  scale_color_manual(values = use_color_stage) +
  facet_wrap(~label, ncol = 4, scales = "free_y")
print(fig)
dev.off()


plot_df <- plot_df[which(plot_df$label %in% c("Ventricular_cardiomyocytes", "Epicardial_cells", "Epicardium_derived_cells")),]
plot_df$label <- factor(plot_df$label, levels = c("Ventricular_cardiomyocytes", "Epicardial_cells", "Epicardium_derived_cells"))

fname <- "overall_proportion_by_stage_boxplot_partial"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 8, height = 4)
fig <- ggplot(data = plot_df, aes(x = stage, y = as.numeric(prop), color = stage)) +
  geom_boxplot(size = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  labs(x = "", y = "", title = fname) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background  = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "bottom") +
  scale_color_manual(values = use_color_stage) +
  facet_wrap(~label, ncol = 4, scales = "free_y")
print(fig)
dev.off()



overall_proportion_all_stage_median <- t(sapply(c("5_PCW", "6_PCW", "9_PCW"), function(x){
  apply(overall_proportion_all[get(paste("analysis_sections", x, sep = "_")),], 2, median)
}))
plot_df <- as.data.frame(cbind(rownames(overall_proportion_all_stage_median), overall_proportion_all_stage_median))
colnames(plot_df) <- c("stage", colnames(overall_proportion_all_stage_median))

plot_df <- reshape2::melt(plot_df, id.vars = c("stage"), variable.name = "label", value.name = "prop")
plot_df$prop <- as.numeric(plot_df$prop)
plot_df$stage <- factor(plot_df$stage, levels = c("5_PCW", "6_PCW", "9_PCW"))
plot_df$label <- factor(plot_df$label, levels = unique_labels)
fname <- "overall_proportion_by_stage_lineplot"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 6, height = 4)
fig <- ggplot(data = plot_df, aes(x = stage, y = as.numeric(prop), color = label, group = label)) +
  geom_point(aes(shape = label), size = 2) +
  geom_line(linewidth = 0.5) +
  labs(x = "Stage", y = "Proportion", title = fname) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background  = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_shape_manual(values = c(16, 17, 15, 18, 1, 2, 0, 5, 4, 6, 7, 9)) +
  scale_color_manual(values = use_color)
print(fig)
dev.off()



########## Analysis of cell-type-specific gene expression

### Visualize cell type-specific gene expression and corresponding scRNA-seq reference
keep_genes <- colnames(sc_results$archetypes_list[[1]])
mu_list_all <- lapply(st_results_all, function(x){
  temp <- x$mu
  temp <- lapply(temp, function(y){y[, keep_genes]})
  temp
})


mu_reference_res <- plot_mu_reference_2(mu_list_all, proportion_all, sc_counts, sc_labels, archetypes_list = sc_results$archetypes_list, min_prop = 0.1, do_hvgs = FALSE, use_color = use_color_section_stage, f_name = "mu_reference_list", save_dir = output_dir, return_res = TRUE, n_cols = 3, width = 15, height = 15)
saveRDS(mu_reference_res, file = paste(output_dir, "mu_reference_res.rds", sep = "/"))


### Focus on Ventricular_cardiomyocytes

Ventricular_cardiomyocytes_mu_all <- lapply(analysis_sections, function(x){mu_list_all[[x]][["Ventricular_cardiomyocytes"]][which(proportion_all[[x]][, "Ventricular_cardiomyocytes"] > 0.1),]})
names(Ventricular_cardiomyocytes_mu_all) <- analysis_sections


Ventricular_cardiomyocytes_mu_stage <- lapply(c("5_PCW", "6_PCW", "9_PCW"), function(x){Reduce("rbind", Ventricular_cardiomyocytes_mu_all[get(paste("analysis_sections", x, sep = "_"))])})
names(Ventricular_cardiomyocytes_mu_stage) <- c("5_PCW", "6_PCW", "9_PCW")
Ventricular_cardiomyocytes_mu_stage <- lapply(Ventricular_cardiomyocytes_mu_stage, na.omit)


Ventricular_cardiomyocytes_mu_stage_mtx <- Reduce("rbind", Ventricular_cardiomyocytes_mu_stage)
Ventricular_cardiomyocytes_mu_stage_mtx_scale <- scale(Ventricular_cardiomyocytes_mu_stage_mtx[, colSums(Ventricular_cardiomyocytes_mu_stage_mtx) != 0], center = TRUE, scale = TRUE)
Ventricular_cardiomyocytes_mu_stage_scale <- lapply(names(Ventricular_cardiomyocytes_mu_stage), function(x){Ventricular_cardiomyocytes_mu_stage_mtx_scale[rownames(Ventricular_cardiomyocytes_mu_stage[[x]]),]})
names(Ventricular_cardiomyocytes_mu_stage_scale) <- names(Ventricular_cardiomyocytes_mu_stage)


genes_retain <- colnames(Ventricular_cardiomyocytes_mu_stage_mtx_scale)

mean_5 <- apply(Ventricular_cardiomyocytes_mu_stage[["5_PCW"]][, genes_retain], 2, mean)
mean_6 <- apply(Ventricular_cardiomyocytes_mu_stage[["6_PCW"]][, genes_retain], 2, mean)
mean_9 <- apply(Ventricular_cardiomyocytes_mu_stage[["9_PCW"]][, genes_retain], 2, mean)

scale_mean_5 <- apply(Ventricular_cardiomyocytes_mu_stage_scale[["5_PCW"]][, genes_retain], 2, mean)
scale_mean_6 <- apply(Ventricular_cardiomyocytes_mu_stage_scale[["6_PCW"]][, genes_retain], 2, mean)
scale_mean_9 <- apply(Ventricular_cardiomyocytes_mu_stage_scale[["9_PCW"]][, genes_retain], 2, mean)

pvalue_5_6_less <- sapply(genes_retain, function(x){
  t.test(Ventricular_cardiomyocytes_mu_stage_scale[["5_PCW"]][, x], Ventricular_cardiomyocytes_mu_stage_scale[["6_PCW"]][, x], alternative = "less")$p.value
})
pvalue_5_6_greater <- sapply(genes_retain, function(x){
  t.test(Ventricular_cardiomyocytes_mu_stage_scale[["5_PCW"]][, x], Ventricular_cardiomyocytes_mu_stage_scale[["6_PCW"]][, x], alternative = "greater")$p.value
})
pvalue_6_9_less <- sapply(genes_retain, function(x){
  t.test(Ventricular_cardiomyocytes_mu_stage_scale[["6_PCW"]][, x], Ventricular_cardiomyocytes_mu_stage_scale[["9_PCW"]][, x], alternative = "less")$p.value
})
pvalue_6_9_greater <- sapply(genes_retain, function(x){
  t.test(Ventricular_cardiomyocytes_mu_stage_scale[["6_PCW"]][, x], Ventricular_cardiomyocytes_mu_stage_scale[["9_PCW"]][, x], alternative = "greater")$p.value
})

temp <- (pvalue_5_6_less < 0.05) * 1 + (pvalue_5_6_greater < 0.05) * (-1) + (pvalue_6_9_less < 0.05) * 10 + (pvalue_6_9_greater < 0.05) * (-10)

genes_retain_res <- cbind(genes_retain, mean_5, mean_6, mean_9, scale_mean_5, scale_mean_6, scale_mean_9, pvalue_5_6_less, pvalue_5_6_greater, pvalue_6_9_less, pvalue_6_9_greater)
colnames(genes_retain_res) <- c("gene", "mean_5", "mean_6", "mean_9", "scale_mean_5", "scale_mean_6", "scale_mean_9", "pvalue_5_6_less", "pvalue_5_6_greater", "pvalue_6_9_less", "pvalue_6_9_greater")


all_patterns_Ventricular_cardiomyocytes <- list()
all_patterns_Ventricular_cardiomyocytes[["Increase - Increase"]] <- genes_retain_res[temp == 11, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Increase - Stable"]] <- genes_retain_res[temp == 1, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Increase - Decrease"]] <- genes_retain_res[temp == -9, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Stable - Increase"]] <- genes_retain_res[temp == 10, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Stable - Stable"]] <- genes_retain_res[temp == 0, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Stable - Decrease"]] <- genes_retain_res[temp == -10, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Decrease - Increase"]] <- genes_retain_res[temp == 9, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Decrease - Stable"]] <- genes_retain_res[temp == -1, , drop = FALSE]
all_patterns_Ventricular_cardiomyocytes[["Decrease - Decrease"]] <- genes_retain_res[temp == -11, , drop = FALSE]

if (!dir.exists(paste(output_dir, "all_patterns_Ventricular_cardiomyocytes", sep = "/"))){
  dir.create(paste(output_dir, "all_patterns_Ventricular_cardiomyocytes", sep = "/"), recursive = TRUE)
}

for (i in names(all_patterns_Ventricular_cardiomyocytes)){
  write.csv(all_patterns_Ventricular_cardiomyocytes[[i]], file = paste(paste(output_dir, "all_patterns_Ventricular_cardiomyocytes", sep = "/"), paste(i, "csv", sep = "."), sep = "/"))
}


if (!dir.exists(paste(output_dir, "enrichment_Ventricular_cardiomyocytes", sep = "/"))){
  dir.create(paste(output_dir, "enrichment_Ventricular_cardiomyocytes", sep = "/"), recursive = TRUE)
}

enrich_res <- list()
for (i in names(all_patterns_Ventricular_cardiomyocytes)){
  enrich_res[[i]] <- plot_enrichGO(all_patterns_Ventricular_cardiomyocytes[[i]][, "gene"], n_category = 12, fname = i, save_dir = paste(output_dir, "enrichment_Ventricular_cardiomyocytes", sep = "/"), width = 5, height = 5, return_results = TRUE)
}

fig_list <- lapply(enrich_res, function(x){x$fig})
fname <- "enrichment_Ventricular_cardiomyocytes"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 24)
fig <- ggarrange(plotlist = fig_list, ncol = 2, nrow = 4)
print(fig)
dev.off()


fig_list <- list()
fig_list[[1]] <- barplot(enrich_res[["Stable - Increase"]]$ego, showCategory = 6, drop = TRUE, label_format = 100) +
  labs(title = "Stable - Increase") +
  theme_prism(base_size = 10, base_fontface = "plain")
fig_list[[2]] <- barplot(enrich_res[["Decrease - Increase"]]$ego, showCategory = 6, drop = TRUE, label_format = 100) +
  labs(title = "Decrease - Increase") +
  theme_prism(base_size = 10, base_fontface = "plain")
fig_list[[3]] <- barplot(enrich_res[["Stable - Decrease"]]$ego, showCategory = 6, drop = TRUE, label_format = 100) +
  labs(title = "Stable - Decrease") +
  theme_prism(base_size = 10, base_fontface = "plain")
fig_list[[4]] <- barplot(enrich_res[["Increase - Decrease"]]$ego, showCategory = 6, drop = TRUE, label_format = 100) +
  labs(title = "Increase - Decrease") +
  theme_prism(base_size = 10, base_fontface = "plain")
fname <- "enrichment_Ventricular_cardiomyocytes_partial"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 24, height = 3)
fig <- ggarrange(plotlist = fig_list, ncol = 4, nrow = 1)
print(fig)
dev.off()



if (!dir.exists(paste(output_dir, "scale_lineplot_Ventricular_cardiomyocytes", sep = "/"))){
  dir.create(paste(output_dir, "scale_lineplot_Ventricular_cardiomyocytes", sep = "/"), recursive = TRUE)
}

fig_list <- list()
for (i in names(all_patterns_Ventricular_cardiomyocytes)){
  plot_genes <- all_patterns_Ventricular_cardiomyocytes[[i]]
  fname <- i
  plot_df <- as.data.frame(plot_genes[, c("gene", "scale_mean_5", "scale_mean_6", "scale_mean_9"), drop = FALSE])
  colnames(plot_df) <- c("gene", "4.5-5 PCW", "6.5 PCW", "9 PCW")
  plot_df <- reshape2::melt(plot_df, id.vars = c("gene"), variable.name = "stage", value.name = "value")
  pdf(paste(paste(output_dir, "scale_lineplot_Ventricular_cardiomyocytes", sep = "/"), paste(fname, "pdf", sep = "."), sep = "/"), width = 3, height = 2.5)
  fig_list[[i]] <- ggplot(data = plot_df, aes(x = stage, y = as.numeric(value), color = gene, group = gene)) +
    geom_point(alpha = 0.1) +
    geom_line(linewidth = 0.5, alpha = 0.1) +
    geom_smooth(aes(group = 1), method = "loess", se = FALSE) +
    labs(x = "", y = "", title = fname) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.background  = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  print(fig_list[[i]])
  dev.off()
}

fname <- "scale_lineplot_Ventricular_cardiomyocytes"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 7.5)
fig <- ggarrange(plotlist = fig_list, ncol = 3, nrow = 3)
print(fig)
dev.off()

fname <- "scale_lineplot_Ventricular_cardiomyocytes_partial"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 12, height = 2.5)
fig <- ggarrange(plotlist = fig_list[c("Stable - Increase", "Decrease - Increase", "Stable - Decrease", "Increase - Decrease")], ncol = 4, nrow = 1)
print(fig)
dev.off()

### Analysis of cell-type-specific gene expression with ANOVA
# Prepare Data for ANOVA Analysis
anv_data <- data.frame()
for (stage in names(Ventricular_cardiomyocytes_mu_stage_scale)) {
  temp_data <- as.data.frame(Ventricular_cardiomyocytes_mu_stage_scale[[stage]])
  temp_data$stage <- stage
  temp_data$spot <- rownames(temp_data)
  anv_data <- rbind(anv_data, temp_data)
}

anv_data_long <- anv_data %>%
  pivot_longer(cols = all_of(genes_retain), names_to = "gene", values_to = "expression")

anova_results <- data.frame(gene = genes_retain, p_value = NA, f_statistic = NA)

for (i in 1:length(genes_retain)) {
  gene <- genes_retain[i]
  gene_data <- anv_data_long %>% filter(gene == !!gene)
  if (nrow(gene_data) < 3 || length(unique(gene_data$stage)) < 3) {
    next
  }
  model <- aov(expression ~ stage, data = gene_data)
  aov_summary <- summary(model)
  anova_results[i, "p_value"] <- aov_summary[[1]][["Pr(>F)"]][1]
  anova_results[i, "f_statistic"] <- aov_summary[[1]][["F value"]][1]
}

# Multiple testing correction (FDR)
anova_results$padj <- p.adjust(anova_results$p_value, method = "fdr")

# Screen for significantly differentially expressed genes (FDR < 0.05)
sig_genes <- anova_results %>% 
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(padj)

# Add expression pattern classification for significant genes
sig_genes$pattern <- sapply(sig_genes$gene, function(g) {
  for (pattern in names(all_patterns_Ventricular_cardiomyocytes)) {
    if (g %in% all_patterns_Ventricular_cardiomyocytes[[pattern]][, "gene"]) {
      return(pattern)
    }
  }
  return("Unknown")
})

write.csv(anova_results, file.path(output_dir, "anova_results.csv"), row.names = FALSE)
write.csv(sig_genes, file.path(output_dir, "significant_genes_anova.csv"), row.names = FALSE)

### Visualization of results ANOVA
anova_dir <- file.path(output_dir, "anova_analysis")
if (!dir.exists(anova_dir)) {
  dir.create(anova_dir, recursive = TRUE)
}

pdf(file.path(anova_dir, "volcano_plot.pdf"), width = 8, height = 6)
ggplot(anova_results, aes(x = log2(f_statistic), y = -log10(padj))) + geom_point(aes(color = padj < 0.05), size = 2) + scale_color_manual(values = c("gray", "red")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") + labs(title = "Volcano Plot of ANOVA Results", x = "Log2(F-statistic)", y = "-Log10(Adjusted p-value)") + theme_minimal() + theme(legend.position = "bottom")
dev.off()

# The heat map shows the top 50 significant genes
if (nrow(sig_genes) > 0) {
  top_genes <- sig_genes$gene[1:min(50, nrow(sig_genes))]
  heatmap_data <- anv_data_long %>% 
    filter(gene %in% top_genes) %>%
    group_by(gene, stage) %>%
    summarise(mean_expr = mean(expression), .groups = "drop") %>%
    pivot_wider(names_from = stage, values_from = mean_expr)
  
  rownames(heatmap_data) <- heatmap_data$gene
  heatmap_data <- heatmap_data[, -1]

  pdf(file.path(anova_dir, "top_genes_heatmap.pdf"), width = 10, height = 12)
  pheatmap::pheatmap(as.matrix(heatmap_data), scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, main = "Top 50 Differentially Expressed Genes", fontsize = 8)
  dev.off()
}

# Create box plots for each gene expression pattern
if (nrow(sig_genes) > 0) {
  for (pattern in unique(sig_genes$pattern)) {
    pattern_genes <- sig_genes$gene[sig_genes$pattern == pattern]
    if (length(pattern_genes) > 0) {
      pattern_data <- anv_data_long %>% filter(gene %in% pattern_genes)
      pdf(file.path(anova_dir, paste0("boxplot_", gsub(" ", "_", pattern), ".pdf")), width = 10, height = 8)
      p <- ggplot(pattern_data, aes(x = stage, y = expression)) + geom_boxplot(aes(fill = stage)) + facet_wrap(~gene, scales = "free_y") + labs(title = paste("Gene Expression by Stage -", pattern), x = "Developmental Stage", y = "Normalized Expression") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
      dev.off()
    }
  }
}

# ANOVA post hoc examination (Tukey HSD)
if (nrow(sig_genes) > 0) {
  tukey_results <- list()
  for (gene in sig_genes$gene[1:min(20, nrow(sig_genes))]) { 
    gene_data <- anv_data_long %>% filter(gene == !!gene)
    model <- aov(expression ~ stage, data = gene_data)
    tukey <- TukeyHSD(model)
    tukey_results[[gene]] <- tukey$stage
  }
  pdf(file.path(anova_dir, "tukey_test_results.pdf"), width = 12, height = 10)
  for (gene in names(tukey_results)) {
    plot(tukey_results[[gene]], main = paste("Tukey HSD Test -", gene))
  }
  dev.off()
}

### Integrate the original analysis and ANOVA results
# Add ANOVA significance information to each mode in the original analysis
for (pattern in names(all_patterns_Ventricular_cardiomyocytes)) {
  genes <- all_patterns_Ventricular_cardiomyocytes[[pattern]][, "gene"]
  sig_genes_in_pattern <- genes[genes %in% sig_genes$gene]
  all_patterns_Ventricular_cardiomyocytes[[pattern]] <- 
    cbind(all_patterns_Ventricular_cardiomyocytes[[pattern]], anova_padj = ifelse(genes %in% sig_genes$gene, sig_genes$padj[match(genes, sig_genes$gene)], NA))
  write.csv(all_patterns_Ventricular_cardiomyocytes[[pattern]], file = paste(paste(output_dir, "all_patterns_Ventricular_cardiomyocytes", sep = "/"), paste(pattern, "csv", sep = "."), sep = "/"))
}

if (!dir.exists(paste(output_dir, "integrated_analysis", sep = "/"))) {
  dir.create(paste(output_dir, "integrated_analysis", sep = "/"), recursive = TRUE)
}

# Select the most significant genes for visualization
top_integrated_genes <- sig_genes$gene[1:min(10, nrow(sig_genes))]
integrated_data <- anv_data_long %>% 
  filter(gene %in% top_integrated_genes) %>%
  group_by(gene, stage) %>%
  summarise(mean_expr = mean(expression), .groups = "drop")
pdf(file.path(output_dir, "integrated_analysis", "integrated_expression_profiles.pdf"), width = 12, height = 8)
ggplot(integrated_data, aes(x = stage, y = mean_expr, group = gene)) + geom_line(aes(color = gene), size = 1.2) + geom_point(aes(color = gene), size = 3) + facet_wrap(~gene, scales = "free_y") + labs(title = "Integrated Expression Profiles of Top ANOVA Genes", x = "Developmental Stage", y = "Mean Normalized Expression") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Visualize the proportion of significant genes in each expression pattern
pattern_sig_counts <- sapply(names(all_patterns_Ventricular_cardiomyocytes), function(pattern) {
  genes <- all_patterns_Ventricular_cardiomyocytes[[pattern]][, "gene"]
  sum(genes %in% sig_genes$gene)
})

pattern_total_counts <- sapply(names(all_patterns_Ventricular_cardiomyocytes), function(pattern) {
  nrow(all_patterns_Ventricular_cardiomyocytes[[pattern]])
})

pattern_sig_ratio <- pattern_sig_counts / pattern_total_counts
pattern_sig_data <- data.frame(Pattern = names(all_patterns_Ventricular_cardiomyocytes), Significant = pattern_sig_counts, Total = pattern_total_counts, Ratio = pattern_sig_ratio)

pdf(file.path(output_dir, "integrated_analysis", "pattern_significance_ratios.pdf"), width = 10, height = 6)
ggplot(pattern_sig_data, aes(x = reorder(Pattern, -Ratio), y = Ratio)) + geom_bar(stat = "identity", fill = "skyblue") + geom_text(aes(label = paste0(Significant, "/", Total)), vjust = -0.5) + labs(title = "Proportion of Significant Genes in Each Expression Pattern", x = "Expression Pattern", y = "Proportion of Significant Genes") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
