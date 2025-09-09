setwd("D:/Rshuju/shiyan")
source("D:/Rshuju/shiyan/SMART/utils/utils.R")
source("D:/Rshuju/shiyan/SMART/utils/plot.R")
library(ggplot2)
library(ggpubr)
library(ggprism)
library(scatterpie)


#################### Preprocess the original scRNA-seq data

### The original scRNA-seq data can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056.
sc_counts <- read.table("D:/Rshuju/shiyan/SMART_data/melanoma/GSE72056_melanoma_single_cell_revised_v2.txt", header = TRUE, sep = "\t", quote = "")
sc_counts <- sc_counts[!duplicated(sc_counts[, 1]),]
rownames(sc_counts) <- sc_counts[, 1]
sc_counts <- sc_counts[, -1]
idx_malignant <- which(sc_counts[2,] == 2)        
idx_nonmalignant <- which(sc_counts[2,] == 1)    
idx_unresolved <- which(sc_counts[2,] == 0)       

sc_labels <- as.numeric(sc_counts[3,])
names(sc_labels) <- colnames(sc_counts)

keep_cells <- c(intersect(idx_malignant, as.numeric(which(sc_labels == 0))), as.numeric(which(sc_labels != 0)))
keep_cells <- keep_cells[-which(keep_cells %in% intersect(idx_malignant, as.numeric(which(sc_labels != 0))))]

sc_labels <- sc_labels[keep_cells]

sc_labels[which(sc_labels == 0)] <- "Malignant"
sc_labels[which(sc_labels == 1)] <- "T"
sc_labels[which(sc_labels == 2)] <- "B"
sc_labels[which(sc_labels == 3)] <- "Macro"
sc_labels[which(sc_labels == 4)] <- "Endo"
sc_labels[which(sc_labels == 5)] <- "CAF"
sc_labels[which(sc_labels == 6)] <- "NK"

sc_counts <- sc_counts[-c(1:3), keep_cells]
sc_counts <- t(as.matrix(sc_counts))
sc_counts <- sc_counts[, -grep("_", colnames(sc_counts))]

sc_counts <- round(2^sc_counts - 1)

write.csv(sc_counts, file = "D:/Rshuju/shiyan/SMART_data/melanoma/sc_data/sc_counts.csv")
sc_labels_m <- as.matrix(sc_labels)
colnames(sc_labels_m) <- "celltype"
write.csv(sc_labels_m, file = "D:/Rshuju/shiyan/SMART_data/melanoma/sc_data/sc_labels.csv")



#################### Preprocess the original spatial transcriptomics data

### The original spatial transcriptomics data can be downloaded from https://www.spatialresearch.org/resources-published-datasets/doi-10-1158-0008-5472-can-18-0747/.
st_counts <- read.table("D:/Rshuju/shiyan/SMART_data/melanoma/ST-Melanoma-Datasets_1/ST_mel1_rep2_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
st_location <- as.data.frame(Reduce(rbind, strsplit(substring(colnames(st_counts), 2), split = "x")))

rownames(st_location) <- colnames(st_counts)
colnames(st_location) <- c("x", "y")

st_counts <- t(as.matrix(st_counts))

colnames(st_counts) <- as.vector(sapply(colnames(st_counts), function(x){strsplit(x, split = " ")[[1]][1]}))

idx_keep <- which(rowSums(st_counts) > 100)
st_counts <- st_counts[idx_keep,]

st_location <- st_location[rownames(st_counts),]

save_dir <- "D:/Rshuju/shiyan/SMART_data/melanoma/st_data"
if (!dir.exists(save_dir)){
  dir.create(save_dir, recursive = TRUE)
}

write.csv(st_counts, file = paste(save_dir, "st_counts.csv", sep = "/"))
write.csv(st_location, file = paste(save_dir, "st_location.csv", sep = "/"))



#################### Run SMART
##install.packages('remotes') remotes::install_version(package = 'Seurat', version = package_version('4.3.0'))
##install.packages("devtools") devtools::install_github("Zhangxf-ccnu/SMART")
### SMART package can be installed by install.packages("SMART/SMART_1.1.0.tar.gz", repos = NULL, type = "source").
### SMART package is also available at https://github.com/Zhangxf-ccnu/SMART.

library(SMART)
sc_counts_path <- "D:/Rshuju/shiyan/SMART_data/melanoma/sc_data/sc_counts.csv"
sc_labels_path <- "D:/Rshuju/shiyan/SMART_data/melanoma/sc_data/sc_labels.csv"

sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)

names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels

sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10, n_cores = 6, save_dir = "D:/Rshuju/shiyan/SMART_results/melanoma/sc_results")

st_counts_path <- "D:/Rshuju/shiyan/SMART_data/melanoma/st_data/st_counts.csv"
st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
st_coords <- as.matrix(dist(st_location))

st_results <- st_train(st_counts = st_counts, sc_results = sc_results, spatial_coords = st_coords, lambda_spatial = 0.8, neighbor_radius = 5, k_neighbors = 10, adaptive_sigma = 1.2, save_mu_csv = TRUE, save_dir = paste("D:/Rshuju/shiyan/SMART_results/melanoma/st_results"))

################## Analysis

output_dir <- "D:/Rshuju/shiyan/SMART_results/melanoma/Analysis"
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}


### Load results

sc_results <- readRDS("D:/Rshuju/shiyan/SMART_results/melanoma/sc_results/sc_results.rds")
st_results <- readRDS("D:/Rshuju/shiyan/SMART_results/melanoma/st_results/st_results.rds")

sc_counts <- read.csv("D:/Rshuju/shiyan/SMART_data/melanoma/sc_data/sc_counts.csv", header = TRUE, row.names = 1)
sc_labels <- read.csv("D:/Rshuju/shiyan/SMART_data/melanoma/sc_data/sc_labels.csv", header = TRUE, row.names = 1)
names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels

st_counts <- read.csv("D:/Rshuju/shiyan/SMART_data/melanoma/st_data/st_counts.csv", header = TRUE, row.names = 1)
st_location <- read.csv("D:/Rshuju/shiyan/SMART_data/melanoma/st_data/st_location.csv", header = TRUE, row.names = 1)
colnames(st_location) <- c("y", "x")

region_labels <- read.csv("D:/Rshuju/shiyan/SMART_data/melanoma/st_data/region_labels.csv", header = TRUE, row.names = 1)
names_region_labels <- rownames(region_labels)
region_labels <- region_labels$region_labels
names(region_labels) <- names_region_labels

lymphoid_spots <- names(region_labels)[region_labels == "Lymphoid"]
melanoma_spots <- names(region_labels)[region_labels == "Melanoma"]
stromal_spots <- names(region_labels)[region_labels == "Stromal"]
plot_scatter_label(as.data.frame(region_labels), st_location, img = NULL, use_color = c("Others" = "grey90", "Lymphoid" = "#F8766D", "Melanoma" = "#00BFC4", "Stromal" = "#7081DA"), fname = "region_labels", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 1, point_size = 5, width = 5.8, height = 5)


########## Analysis of cell type proportions

proportion <- st_results$proportion
unique_labels <- sort(colnames(proportion))
use_color <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_labels))
names(use_color) <- unique_labels


### Plot the scatterpie
plot_scatterpie(proportion, st_location, img = NULL, use_color = use_color, fname = "proportion_scatterpie", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, pie_r = 0.5, width = 7, height = 7)


### Plot the heatmap for each cell type
plot_scatter_heatmap(proportion, st_location, img = NULL, feature = unique_labels, scale = TRUE, fname = "scaled_proportion_heatmap", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 4, point_size = 3, width = 12, height = 8, legend_position = "bottom")
plot_scatter_heatmap(proportion, st_location, img = NULL, feature = c("Malignant", "B", "T", "CAF", "Endo"), scale = TRUE, fname = "scaled_proportion_heatmap_partial", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 5, point_size = 2, width = 10, height = 3, legend_position = "bottom")


### Find marker genes for each cell type based on scRNA-seq reference
marker_genes <- find_markers(t(sc_counts[, intersect(colnames(sc_counts), colnames(st_counts))]), sc_labels, n_markers = NULL)
marker_genes <- marker_genes[unique_labels]
saveRDS(marker_genes, file = paste(output_dir, "marker_genes.rds", sep = "/"))
write.csv(as.data.frame(lapply(marker_genes, function(x){x[1:5]})), file = paste(output_dir, "marker_genes.csv", sep = "/"))


### Plot the gene expression of markers for each cell type
st_counts_norm <- st_counts / rowSums(st_counts)
st_counts_norm_metagene <- sapply(marker_genes, function(x){log2(rowMeans(st_counts_norm[, x[1:5]]) * 1e4 + 1)})
plot_scatter_heatmap(st_counts_norm_metagene, st_location, img = NULL, feature = unique_labels, scale = TRUE, fname = "scaled_metagene_exp_heatmap", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 4, point_size = 3, width = 12, height = 8, legend_position = "bottom")
plot_scatter_heatmap(st_counts_norm_metagene, st_location, img = NULL, feature = c("Malignant", "B", "T", "CAF", "Endo"), scale = TRUE, fname = "scaled_metagene_exp_heatmap_partial", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 5, point_size = 2, width = 10, height = 3, legend_position = "bottom")

# Calculate the Pearson correlation coefficient for metagene expression and the corresponding ratio for each cell type
proportion <- st_results$proportion
metagene_expr <- st_counts_norm_metagene
corr_results <- data.frame(CellType = character(), Correlation = numeric(), PValue = numeric(), AdjustedPValue = numeric(), stringsAsFactors = FALSE)
for (cell_type in colnames(metagene_expr)) {
  if (cell_type %in% colnames(proportion)) {
    expr <- metagene_expr[, cell_type]
    prop <- proportion[, cell_type]
    corr_test <- cor.test(expr, prop, method = "pearson")
    corr_results <- rbind(corr_results, data.frame(CellType = cell_type, Correlation = corr_test$estimate, PValue = corr_test$p.value, AdjustedPValue = p.adjust(corr_test$p.value, method = "fdr"), stringsAsFactors = FALSE))
  } else {warning(paste("Cell type", cell_type, "not found in proportion matrix."))}
}
write.csv(corr_results, file.path(output_dir, "metagene_proportion_correlation.csv"), row.names = FALSE)

### Region-based enrichment of cell types
enrich_res <- enrichment_celltype(proportion, region_labels, n_shuffle = 10000)
saveRDS(enrich_res, file = paste(output_dir, "enrich_res.rds", sep = "/"))
# enrich_res <- readRDS(paste(output_dir, "enrich_res.rds", sep = "/"))
write.csv(as.data.frame(enrich_res$diff_mean), file = paste(output_dir, "enrichment_celltype_region.csv", sep = "/"))

plot_df <- enrich_res$diff_mean[c("Lymphoid", "Melanoma", "Stromal"),]
plot_df <- as.data.frame(cbind(rownames(plot_df), plot_df))
colnames(plot_df) <- c("Region", colnames(enrich_res$diff_mean))
plot_df <- reshape2::melt(plot_df, id.vars = c("Region"), variable.name = "Label", value.name = "value")
plot_df$Region <- factor(plot_df$Region, levels = c("Lymphoid", "Melanoma", "Stromal"))
plot_df$Label <- factor(plot_df$Label, levels = unique_labels)

fname <- "enrichment_celltype_region"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 4.5, height = 4)
fig <- ggplot() +
  geom_point(data = plot_df, aes(x = Region, y = Label, size = abs(as.numeric(value)), color = as.character(sign(as.numeric(value))))) +
  theme_bw() +
  scale_size(range = c(1, 8)) +
  scale_color_manual(values = c("-1" = "#D13854", "1" = "#43BF55")) +
  scale_y_discrete(limits = rev(levels(plot_df$Label))) +
  labs(title = fname, color = "Direction", size = "Effect size")
print(fig)
dev.off()



########## Analysis of cell-type-specific gene expression

### Visualize cell type-specific gene expression and corresponding scRNA-seq reference
mu_list <- st_results$mu
keep_genes <- colnames(sc_results$archetypes_list[[1]])
mu_list <- lapply(mu_list, function(x){x[, keep_genes]})

mu_reference_res <- plot_mu_reference(mu_list, proportion, sc_counts, sc_labels, archetypes_list = sc_results$archetypes_list, min_prop = 0.1, f_name = "mu_reference", save_dir = output_dir, return_res = TRUE, n_cols = 4, width = 12, height = 6)
saveRDS(mu_reference_res, file = paste(output_dir, "mu_reference_res.rds", sep = "/"))


### Focus on CAF cells
keep_spots <- which(proportion[, "CAF"] > 0.1)
CAF_mu <- mu_list[["CAF"]][keep_spots,]

### Plot the heatmap of alpha for each archetypes
CAF_alpha <- st_results$alpha[keep_spots, grep("CAF", colnames(st_results$alpha))]

unique_alpha <- colnames(CAF_alpha)
use_color_alpha <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_alpha))
names(use_color_alpha) <- unique_alpha

plot_scatter_heatmap(CAF_alpha, st_location, img = NULL, feature = colnames(CAF_alpha), scale = FALSE, fname = "CAF_alpha_heatmap", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 5, point_size = 2.5, width = 12, height = 6, legend_position = "bottom")

### Plot the scatterpie
plot_scatterpie(CAF_alpha, st_location, img = NULL, use_color = use_color_alpha, fname = "CAF_alpha_scatterpie", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, pie_r = 0.45, point_size = 6.3, width = 7, height = 7)

### Segmentation of lymphoid region and melanoma region on the pie plot
lymphoid_df <- data.frame(x = c(20.5, 20.5, 19.5, 19.5, 20.5, 20.5, 28.5, 28.5, 22.5, 20.5),
                          y = c(27.5, 23.5, 23.5, 22.5, 22.5, 20.5, 20.5, 21.5, 27.5, 27.5))
lymphoid_df$xend <- c(lymphoid_df$x[2:nrow(lymphoid_df)], lymphoid_df$x[1])
lymphoid_df$yend <- c(lymphoid_df$y[2:nrow(lymphoid_df)], lymphoid_df$y[1])

melanoma_df <- data.frame(x = c(8.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5, 14.5, 14.5, 15.5,
                                15.5, 16.5, 16.5, 17.5, 17.5, 18.5, 18.5, 16.5, 16.5, 14.5,
                                14.5, 12.5, 12.5, 11.5, 11.5, 10.5, 10.5, 6.5, 6.5, 4.5,
                                4.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5, 5.5, 5.5, 8.5),
                          y = c(25.5, 25.5, 24.5, 24.5, 23.5, 23.5, 22.5, 22.5, 21.5, 21.5,
                                19.5, 19.5, 9.5, 9.5, 8.5, 8.5, 7.5, 7.5, 6.5, 6.5,
                                7.5, 7.5, 10.5, 10.5, 12.5, 12.5, 13.5, 13.5, 14.5, 14.5,
                                15.5, 15.5, 16.5, 16.5, 19.5, 19.5, 20.5, 20.5, 22.5, 25.5))
melanoma_df$xend <- c(melanoma_df$x[2:nrow(melanoma_df)], melanoma_df$x[1])
melanoma_df$yend <- c(melanoma_df$y[2:nrow(melanoma_df)], melanoma_df$y[1])


plot_df <- as.data.frame(apply(cbind(st_location, proportion), 2, as.numeric))
rownames(plot_df) <- rownames(proportion)

fname <- "segmentation_lymphoid_melanoma"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 7, height = 7)
fig <- ggplot() +
  geom_scatterpie(data = plot_df, aes(x = x, y = y, r = 0.5), cols = unique_labels, color = NA) +
  coord_fixed() +
  scale_fill_manual(values = use_color) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background  = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size = 15, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "right") +
  labs(title = fname) +
  geom_segment(data = lymphoid_df, aes(x = x, y = y, xend = xend, yend = yend), lty = 1, color = "purple", linewidth = 1) +
  geom_segment(data = melanoma_df, aes(x = x, y = y, xend = xend, yend = yend), lty = 1, color = "black", linewidth = 1)
print(fig)
dev.off()


### Find differentially expressed genes in CAF cells between stromal region and melanoma region
CAF_stromal_spots <- intersect(rownames(CAF_mu), stromal_spots)
CAF_melanoma_spots <- intersect(rownames(CAF_mu), melanoma_spots)

marker_genes_CAF <- find_markers(t(CAF_mu[c(CAF_stromal_spots, CAF_melanoma_spots), ] * 1e4), c(rep("stromal", length(CAF_stromal_spots)), rep("Melanoma", length(CAF_melanoma_spots))), n_markers = NULL)
saveRDS(marker_genes_CAF, file = paste(output_dir, "marker_genes_CAF.rds", sep = "/"))

marker_genes_CAF_use <- c(marker_genes_CAF[["stromal"]][1:30], marker_genes_CAF[["Melanoma"]][1:30])
plot_df <- scale(log2(CAF_mu[c(CAF_stromal_spots, CAF_melanoma_spots), marker_genes_CAF_use] * 1e4 + 1), center = TRUE, scale = TRUE)

annotation_row <- data.frame(Region = c(rep("stromal", length(CAF_stromal_spots)), rep("Melanoma", length(CAF_melanoma_spots))))
rownames(annotation_row) <- rownames(plot_df)
annotation_col <- data.frame(Region = c(rep("stromal", length(marker_genes_CAF[["stromal"]][1:30])), rep("Melanoma", length(marker_genes_CAF[["Melanoma"]][1:30]))))
rownames(annotation_col) <- colnames(plot_df)
annotation_colors <- list(Region = c("stromal" = "#F8766D", "Melanoma" = "#00CAFFC4"))

CAF_stromal_spots_order <- order(rowMeans(length(CAF_stromal_spots) + 1 - apply(plot_df[1:length(CAF_stromal_spots), 1:length(marker_genes_CAF[["stromal"]][1:30])], 2, rank)))
CAF_melanoma_spots_order <- order(rowMeans(apply(plot_df[(length(CAF_stromal_spots) + 1):nrow(plot_df), (length(marker_genes_CAF[["stromal"]][1:30]) + 1):ncol(plot_df)], 2, rank)))
plot_df <- plot_df[c(CAF_stromal_spots[CAF_stromal_spots_order], CAF_melanoma_spots[CAF_melanoma_spots_order]),]

fname <- "marker_genes_CAF"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 15)
fig <- pheatmap::pheatmap(t(plot_df), cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = annotation_col, annotation_col = annotation_row, annotation_colors = annotation_colors, show_rownames = TRUE, show_colnames = FALSE)
print(fig)
dev.off()


CAF_mu_markers_exp <- sapply(marker_genes_CAF, function(x){log2(rowMeans(CAF_mu[c(CAF_stromal_spots, CAF_melanoma_spots), x, drop = FALSE]) * 1e4 + 1)})
plot_scatter_heatmap(CAF_mu_markers_exp, st_location, img = NULL, feature = NULL, scale = TRUE, fname = "scaled_CAF_mu_markers_heatmap", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 2, point_size = 2.1, width = 5, height = 4)


if (!dir.exists(paste(output_dir, "enrichGO_CAF_stromal_melanoma", sep = "/"))){
  dir.create(paste(output_dir, "enrichGO_CAF_stromal_melanoma", sep = "/"), recursive = TRUE)
}

enrichGO_CAF_stromal <- plot_enrichGO(marker_genes_CAF[["stromal"]], n_category = 12, fname = "enrichGO_CAF_stromal", save_dir = paste(output_dir, "enrichGO_CAF_stromal_melanoma", sep = "/"), width = 5, height = 6, return_results = TRUE)
enrichGO_CAF_melanoma <- plot_enrichGO(marker_genes_CAF[["Melanoma"]], n_category = 12, fname = "enrichGO_CAF_melanoma", save_dir = paste(output_dir, "enrichGO_CAF_stromal_melanoma", sep = "/"), width = 5, height = 6, return_results = TRUE)

fig_list <- list()
fig_list[[1]] <- barplot(enrichGO_CAF_stromal$ego, showCategory = 6, drop = TRUE, label_format = 100) +
  labs(title = "enrichGO_CAF_stromal") +
  theme_prism(base_size = 10, base_fontface = "plain")
fig_list[[2]] <- barplot(enrichGO_CAF_melanoma$ego, showCategory = 6, drop = TRUE, label_format = 100) +
  labs(title = "enrichGO_CAF_melanoma") +
  theme_prism(base_size = 10, base_fontface = "plain")
fname <- "enrichGO_CAF_stromal_melanoma"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 12, height = 3)
fig <- ggarrange(plotlist = fig_list, ncol = 2, nrow = 1)
print(fig)
dev.off()

### Find genes in malignant cells varying with the proportion of immune cells in the melanoma area

keep_spots <- intersect(rownames(proportion)[which(proportion[, "Malignant"] > 0.1)], melanoma_spots)
Malignant_mu <- mu_list[["Malignant"]][keep_spots,]

test_genes <- colnames(Malignant_mu)[colSums(Malignant_mu) != 0]
test_proportion <- rowSums(proportion[keep_spots, c("B", "T", "Macro", "NK")])

temp_prop <- cbind(proportion[keep_spots, "Malignant"], test_proportion, rowSums(proportion[keep_spots, c("CAF", "Endo")]))
colnames(temp_prop) <- c("Malignant", "Immune", "Stromal")
rownames(temp_prop) <- keep_spots
use_color_temp <- c("#F8766D", "#00BFC4", "grey90")
names(use_color_temp) <- c("Immune", "Malignant", "Stromal")

### Plot the scatterpie
plot_scatterpie(temp_prop, st_location, img = NULL, use_color = use_color_temp, fname = "Malignant_Immune_proportion_scatterpie", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, pie_r = 0.5, point_size = 6.75, width = 7, height = 7)

plot_scatter_heatmap(temp_prop[, "Immune", drop = FALSE], st_location, img = NULL, feature = "Immune", scale = TRUE, fname = "Malignant_scaled_Immune_proportion_heatmap", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 1, point_size = 2.2, width = 3, height = 3, legend_position = "right")


test_cor_pvalue <- c()
for (i in seq(length(test_genes))){
  temp <- cor.test(Malignant_mu[, test_genes[i]], test_proportion, method = "pearson")
  test_cor_pvalue <- rbind(test_cor_pvalue, c(temp$estimate, temp$p.value))
}
test_cor_pvalue <- cbind(test_cor_pvalue, p.adjust(test_cor_pvalue[,2], method = "BH"))
rownames(test_cor_pvalue) <- test_genes
colnames(test_cor_pvalue) <- c("correlation", "p_value", "p_adjust")

test_cor_pvalue <- test_cor_pvalue[order(test_cor_pvalue[,1], decreasing = TRUE),]
write.csv(as.matrix(test_cor_pvalue), file = paste(output_dir, "Malignant_Immune_test_cor_pvalue.csv", sep = "/"))

sig_positive_genes <- rownames(test_cor_pvalue)[(test_cor_pvalue[, 3] < 0.05) & (test_cor_pvalue[, 1] > 0)]
sig_negative_genes <- rev(rownames(test_cor_pvalue)[(test_cor_pvalue[, 3] < 0.05) & (test_cor_pvalue[, 1] < 0)])
write.csv(as.matrix(sig_positive_genes), file = paste(output_dir, "Malignant_Immune_sig_positive_genes.csv", sep = "/"))
write.csv(as.matrix(sig_negative_genes), file = paste(output_dir, "Malignant_Immune_sig_negative_genes.csv", sep = "/"))


plot_df <- data.frame(x = seq(nrow(test_cor_pvalue)), y = rev(test_cor_pvalue[, 1]), text = rev(rownames(test_cor_pvalue)), logp = rev(- log10(test_cor_pvalue[, 3])))

f_name <- "Malignant_Immune_test_cor_pvalue"
pdf(paste(output_dir, paste(f_name, "pdf", sep = "."), sep = "/"), width = 4, height = 4.5)
fig <- ggplot() +
  geom_point(data = plot_df, aes(x = x, y = y, color = logp), size = 1, shape = 21, stroke = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Malignant_Immune_test_cor_pvalue", x = "Index", y = "Correlation") +
  scale_color_gradientn(colors = c("#5E4FA2", "#3E96B6", "#99D5A4", "#FDFEBD", "#FDDB88", "#F67948", "#9E0142")) +
  theme_prism(base_size = 10, base_fontface = "plain")
print(fig)
dev.off()



### Plot the average gene expression for each gene cluster
Malignant_mu_Immune_gene <- sapply(list(sig_positive_genes, sig_negative_genes), function(x){log2(rowMeans(Malignant_mu[, x, drop = FALSE]) * 1e4 + 1)})
colnames(Malignant_mu_Immune_gene) <- c("positive", "negative")
plot_scatter_heatmap(Malignant_mu_Immune_gene, st_location, img = NULL, feature = c("positive", "negative"), scale = TRUE, fname = "scaled_Malignant_mu_Immune_gene_exp_heatmap", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 2, point_size = 2.1, width = 5, height = 4, legend_position = "right")



### NicheNet
library(nichenetr)
library(tidyverse)
library(circlize)

# The data "lr_network_human_21122021.rds" and "ligand_target_matrix_nsga2r_final.rds" can be downloaded from https://zenodo.org/records/7074291.

lr_network <- readRDS("D:/Rshuju/shiyan/SMART_data/melanoma/melanoma_NicheNet/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("D:/Rshuju/shiyan/SMART_data/melanoma/melanoma_NicheNet/ligand_target_matrix_nsga2r_final.rds")

keep_spots_immune <- keep_spots[apply(proportion[keep_spots, c("B", "T", "Macro")], 1, max) > 0.1]     # 87 / 150

communication_spot_labels <- rep("Others", nrow(st_location))
names(communication_spot_labels) <- rownames(st_location)
communication_spot_labels[keep_spots] <- "Non-communication spots"
communication_spot_labels[keep_spots_immune] <- "Communication spots"

plot_scatter_label(as.data.frame(communication_spot_labels), st_location, img = NULL, use_color = c("Others" = "grey90", "Communication spots" = "#F8766D", "Non-communication spots" = "#00BFC4"), fname = "communication_spot_labels", save_dir = output_dir, reverse_x = FALSE, reverse_y = FALSE, n_cols = 1, point_size = 5, width = 6.7, height = 5)

keep_spots_B <- intersect(keep_spots_immune, rownames(proportion)[which(proportion[, "B"] > 0.1)])
keep_spots_T <- intersect(keep_spots_immune, rownames(proportion)[which(proportion[, "T"] > 0.1)])
keep_spots_Macro <- intersect(keep_spots_immune, rownames(proportion)[which(proportion[, "Macro"] > 0.1)])

keep_genes <- convert_alias_to_symbols(keep_genes, "human", verbose = FALSE)
mu_list <- lapply(mu_list, function(x){
  temp <- x
  colnames(temp) <- keep_genes
  temp
})

expressed_genes_B <- keep_genes[log2(colMeans(mu_list[["B"]][keep_spots_B, , drop = FALSE] * 1e6) + 1) >= 4]                 # 366 / 2015
expressed_genes_T <- keep_genes[log2(colMeans(mu_list[["T"]][keep_spots_T, , drop = FALSE] * 1e6) + 1) >= 4]                 # 433 / 2015
expressed_genes_Macro <- keep_genes[log2(colMeans(mu_list[["Macro"]][keep_spots_Macro, , drop = FALSE] * 1e6) + 1) >= 4]     # 441 / 2015

expressed_genes_Malignant <- keep_genes[log2(colMeans(mu_list[["Malignant"]][keep_spots_immune, , drop = FALSE] * 1e6) + 1) >= 4]     # 347 / 2015

geneset_oi_pos <- intersect(convert_alias_to_symbols(sig_positive_genes, "human", verbose = FALSE), rownames(ligand_target_matrix))     # 846 / 896
geneset_oi_neg <- intersect(convert_alias_to_symbols(sig_negative_genes, "human", verbose = FALSE), rownames(ligand_target_matrix))     # 421 / 447
geneset_oi <- c(geneset_oi_pos, geneset_oi_neg)     # 1267
target_type_df <- tibble(
  target = geneset_oi,
  target_type = c(rep("Positive", length(geneset_oi_pos)), rep("Negative", length(geneset_oi_neg)))
)
background_expressed_genes <- intersect(convert_alias_to_symbols(test_genes, "human", verbose = FALSE), rownames(ligand_target_matrix))     # 1813

ligands <- lr_network %>% pull(from) %>% unique()     # 1226
expressed_ligands_B <- intersect(ligands, expressed_genes_B)             # 41 / 366
expressed_ligands_T <- intersect(ligands, expressed_genes_T)             # 69 / 433
expressed_ligands_Macro <- intersect(ligands, expressed_genes_Macro)     # 106 / 441
expressed_ligands <- Reduce("union", list(expressed_ligands_B, expressed_ligands_T, expressed_ligands_Macro))     # 140

receptors <- lr_network %>% pull(to) %>% unique()     # 1067
expressed_receptors <- intersect(receptors, expressed_genes_Malignant)     # 29 / 347

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()     # 18

ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)     # 18 x 5

ligand_activities %>% pull(auroc) %>% mean()              #  0.5072625
ligand_activities %>% pull(aupr) %>% mean()               # 0.6462429
ligand_activities %>% pull(aupr_corrected) %>% mean()     # 0.007459588
ligand_activities %>% pull(pearson) %>% mean()            # 0.02440453

best_upstream_ligands <- ligand_activities %>% arrange(-aupr_corrected) %>% pull(test_ligand)

intersect(best_upstream_ligands, expressed_ligands_B)         # 2
intersect(best_upstream_ligands, expressed_ligands_T)         # 15
intersect(best_upstream_ligands, expressed_ligands_Macro)     # 12

ligand_expression_tbl <- tibble(
  ligand = best_upstream_ligands,
  B = log2(colMeans(mu_list[["B"]][keep_spots_B, best_upstream_ligands, drop = FALSE] * 1e6) + 1),
  T = log2(colMeans(mu_list[["T"]][keep_spots_T, best_upstream_ligands, drop = FALSE] * 1e6) + 1),
  Macro = log2(colMeans(mu_list[["Macro"]][keep_spots_Macro, best_upstream_ligands, drop = FALSE] * 1e6) + 1)
)

B_specific_ligands <- ligand_expression_tbl %>% filter(B >= 4 & T < 4 & Macro < 4) %>% pull(ligand)         # 2
T_specific_ligands <- ligand_expression_tbl %>% filter(B < 4 & T >= 4 & Macro < 4) %>% pull(ligand)         # 2
Macro_specific_ligands <- ligand_expression_tbl %>% filter(B < 4 & T < 4 & Macro >= 4) %>% pull(ligand)     # 6
B_T_ligands <- ligand_expression_tbl %>% filter(B >= 4 & T >= 4 & Macro < 4) %>% pull(ligand)               # 1
B_Macro_ligands <- ligand_expression_tbl %>% filter(B >= 4 & T < 4 & Macro >= 4) %>% pull(ligand)           # 2
T_Macro_ligands <- ligand_expression_tbl %>% filter(B < 4 & T >= 4 & Macro >= 4) %>% pull(ligand)           # 2
B_T_Macro_ligands <- ligand_expression_tbl %>% filter(B >= 4 & T >= 4 & Macro >= 4) %>% pull(ligand)        # 3

ligand_type_indication_df <- tibble(
  ligand_type = c(rep("B-specific", times = length(B_specific_ligands)),
                  rep("T-specific", times = length(T_specific_ligands)),
                  rep("Macro-specific", times = length(Macro_specific_ligands)),
                  rep("B-T", times = length(B_T_ligands)),
                  rep("B-Macro", times = length(B_Macro_ligands)),
                  rep("T-Macro", times = length(T_Macro_ligands)),
                  rep("B-T-Macro", times = length(B_T_Macro_ligands))),
  ligand = c(B_specific_ligands, T_specific_ligands, Macro_specific_ligands, B_T_ligands, B_Macro_ligands, T_Macro_ligands, B_T_Macro_ligands)
)

active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>% bind_rows()
active_ligand_target_links_df <- active_ligand_target_links_df %>% inner_join(target_type_df) %>% inner_join(ligand_type_indication_df)

write.csv(as.matrix(active_ligand_target_links_df), file = paste(output_dir, "active_ligand_target_links_df.csv", sep = "/"))

cutoff_include_all_ligands <- active_ligand_target_links_df$weight %>% quantile(0.80)     # 0.06475703
active_ligand_target_links_df_circos <- active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)     # 100 x 5
ligands_to_remove <- setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())     # 4
targets_to_remove <- setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())     # 45

circos_links <- active_ligand_target_links_df %>% filter(!target %in% targets_to_remove & !ligand %in% ligands_to_remove)     # 250 x 5

grid_col_ligand <- c("B-specific" = "#E41A1C", "B-T" = "#BF5A5B", "T-specific" = "#999999", "T-Macro" = "#CC8C4D", "Macro-specific" = "#FF7F00", "B-Macro" = "#F24D0E", "B-T-Macro" = "#D4663C")
grid_col_target <- c("Positive" = "#E1C62F", "Negative" = "#00CBFF")

grid_col_tbl_ligand <- tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target <- tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links <- circos_links %>% mutate(ligand = paste(ligand, " "))
circos_links <- circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle <- circos_links %>% dplyr::select(ligand, target, weight)

ligand_color <- circos_links %>% distinct(ligand, color_ligand_type)
grid_ligand_color <- ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color <- circos_links %>% distinct(target, color_target_type)
grid_target_color <- target_color$color_target_type %>% set_names(target_color$target)

grid_col <- c(grid_ligand_color, grid_target_color)

transparency <- circos_links %>% mutate(weight = (weight - min(weight)) / (max(weight) - min(weight))) %>% mutate(transparency = 1 - weight) %>% .$transparency

active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

ligand_order <- c(B_specific_ligands, B_T_ligands, T_specific_ligands, T_Macro_ligands, Macro_specific_ligands, B_Macro_ligands, B_T_Macro_ligands) %>% c(paste(., " ")) %>% intersect(circos_links$ligand)


dist_targets_pos = dist(active_ligand_target_links[intersect(geneset_oi_pos, circos_links$target),], method = "binary")
hclust_targets_pos = hclust(dist_targets_pos, method = "ward.D2")
order_targets_pos = hclust_targets_pos$labels[hclust_targets_pos$order]     # 56
dist_targets_neg = dist(active_ligand_target_links[intersect(geneset_oi_neg, circos_links$target),], method = "binary")
hclust_targets_neg = hclust(dist_targets_neg, method = "ward.D2")
order_targets_neg = hclust_targets_neg$labels[hclust_targets_neg$order]     # 14
order_targets <- c(order_targets_pos, order_targets_neg)                    # 70


order <- c(ligand_order, order_targets)

width_same_cell_same_ligand_type <- 0.5
width_different_cell <- 6
width_ligand_target <- 15
width_same_cell_same_target_type <- 0.5

gaps <- c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-specific") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-T") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "T-specific") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "T-Macro") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macro-specific") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-Macro") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-T-Macro") %>% distinct(ligand) %>% nrow() - 1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Positive") %>% distinct(target) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Negative") %>% distinct(target) %>% nrow() - 1)),
  width_ligand_target
)

circos.clear()
pdf(paste(output_dir, "ligand_target_circos.pdf", sep = "/"), width = 15, height = 15)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1, order = order, link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col, transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)
circos.clear()
dev.off()


library(cowplot)
library(ggpubr)
library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)

vis_ligand_immune_expression <- ligand_expression_tbl %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(paste(ligand_expression_tbl$ligand, " ")) %>% .[ligand_order %>% rev(),]
p_ligand_immune_expression <- vis_ligand_immune_expression %>% make_heatmap_ggplot("Ligands", "Immune cells", color = color[100], legend_position = "top", x_axis_position = "top", legend_title = "Ligand expression\nin immune cells") + theme(axis.text.y = element_text(face = "italic"))

vis_ligand_target <- active_ligand_target_links %>% magrittr::set_colnames(paste(colnames(.), " ")) %>% .[order_targets, ligand_order %>% rev()] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Ligands","Target genes", color = "purple", legend_position = "top", x_axis_position = "top", legend_title = "Ligand-target\nregulatory potential") + scale_fill_gradient2(low = "whitesmoke", high = "purple", breaks = c(0, 0.1, 0.2)) + theme(axis.text.x = element_text(face = "italic"))

vis_target_expression <- rbind(log2(colMeans(mu_list[["Malignant"]][setdiff(keep_spots, keep_spots_immune), order_targets, drop = FALSE] * 1e6) + 1), log2(colMeans(mu_list[["Malignant"]][keep_spots_immune, order_targets, drop = FALSE] * 1e6) + 1)) %>% as.matrix() %>% magrittr::set_rownames(c("Non-communication\nspots", "Communication\nspots")) %>% scale(scale = FALSE)
p_target_expression = vis_target_expression  %>% make_threecolor_heatmap_ggplot("Expression", "Target genes", low_color = color[1], mid_color = color[50], mid = 0, high_color = color[100], legend_position = "top", x_axis_position = "top", legend_title = "Target gene expression\nin malignant cells") + theme(axis.text.x = element_text(face = "italic"))

rownames(test_cor_pvalue) <- convert_alias_to_symbols(rownames(test_cor_pvalue), "human", verbose = FALSE)
vis_target_immune_cor <- test_cor_pvalue[order_targets, 1] %>% as.matrix() %>% t()
p_target_immune_cor = vis_target_immune_cor  %>% make_threecolor_heatmap_ggplot("Correlation", "Target genes", low_color = color[1], mid_color = color[50], mid = 0, high_color = color[100], legend_position = "top", x_axis_position = "top", legend_title = "Correlation") + theme(axis.text.x = element_text(face = "italic"), axis.text.y = element_blank())

pdf(paste(output_dir, "ligand_target_matrix.pdf", sep = "/"), width = 12, height = 13)
figures_without_legend <- plot_grid(
  p_ligand_immune_expression + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  NULL,
  p_target_expression + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  NULL,
  p_target_immune_cor + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  align = "hv",
  nrow = 3,
  rel_widths = c(ncol(vis_ligand_immune_expression) + 10, ncol(vis_ligand_target)),
  rel_heights = c(nrow(vis_ligand_immune_expression), nrow(vis_target_expression) + 2, nrow(vis_target_immune_cor) + 2)
)
legends <- plot_grid(
  as_ggplot(get_legend(p_ligand_immune_expression)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_target_expression)),
  as_ggplot(get_legend(p_target_immune_cor)),
  nrow = 2,
  align = "h"
)
fig <- plot_grid(figures_without_legend, legends, rel_heights = c(10, 2), nrow = 2, align = "hv")
print(fig)
dev.off()



# KEGG and GO enrichment analysis

active_target_pos <- active_ligand_target_links_df %>% filter(target_type == "Positive") %>% pull(target) %>% unique()     # 86 / 846
active_target_neg <- active_ligand_target_links_df %>% filter(target_type == "Negative") %>% pull(target) %>% unique()     # 29 / 421

if (!dir.exists(paste(output_dir, "enrichKEGG_target", sep = "/"))){
  dir.create(paste(output_dir, "enrichKEGG_target", sep = "/"), recursive = TRUE)
}

enrichKEGG_active_target_pos <- plot_enrichKEGG(active_target_pos, n_category = 12, object = "human", organism = "hsa", fname = "enrichKEGG_active_target_pos", save_dir = paste(output_dir, "enrichKEGG_target", sep = "/"), width = 5, height = 6, return_results = TRUE)
enrichKEGG_active_target_neg <- plot_enrichKEGG(active_target_neg, n_category = 12, object = "human", organism = "hsa", fname = "enrichKEGG_active_target_neg", save_dir = paste(output_dir, "enrichKEGG_target", sep = "/"), width = 5, height = 6, return_results = TRUE)

fig_list <- list()
fig_list[["positive"]] <- enrichKEGG_active_target_pos$fig
fig_list[["negative"]] <- enrichKEGG_active_target_neg$fig

fname <- "enrichKEGG_target"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 6)
fig <- ggarrange(plotlist = fig_list, ncol = 2, nrow = 1)
print(fig)
dev.off()

if (!dir.exists(paste(output_dir, "enrichGO_target", sep = "/"))){
  dir.create(paste(output_dir, "enrichGO_target", sep = "/"), recursive = TRUE)
}

enrichGO_active_target_pos <- plot_enrichGO(active_target_pos, n_category = 12, object = "human", fname = "enrichGO_active_target_pos", save_dir = paste(output_dir, "enrichGO_target", sep = "/"), width = 5, height = 6, return_results = TRUE)
enrichGO_active_target_neg <- plot_enrichGO(active_target_neg, n_category = 12, object = "human", fname = "enrichGO_active_target_neg", save_dir = paste(output_dir, "enrichGO_target", sep = "/"), width = 5, height = 6, return_results = TRUE)

fig_list <- list()
fig_list[["positive"]] <- enrichGO_active_target_pos$fig
fig_list[["negative"]] <- enrichGO_active_target_neg$fig

fname <- "enrichGO_target"
pdf(paste(output_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 6)
fig <- ggarrange(plotlist = fig_list, ncol = 2, nrow = 1)
print(fig)
dev.off()


######## ligand_target_circos_partial

cutoff_include_all_ligands <- active_ligand_target_links_df$weight %>% quantile(0.90)     # 0.07338175
active_ligand_target_links_df_circos <- active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove <- setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove <- setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links <- active_ligand_target_links_df %>% filter(!target %in% targets_to_remove & !ligand %in% ligands_to_remove)

circos_links <- circos_links %>% mutate(ligand = paste(ligand, " "))
circos_links <- circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle <- circos_links %>% dplyr::select(ligand, target, weight)

ligand_color <- circos_links %>% distinct(ligand, color_ligand_type)
grid_ligand_color <- ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color <- circos_links %>% distinct(target, color_target_type)
grid_target_color <- target_color$color_target_type %>% set_names(target_color$target)

grid_col <- c(grid_ligand_color, grid_target_color)

transparency <- circos_links %>% mutate(weight = (weight - min(weight)) / (max(weight) - min(weight))) %>% mutate(transparency = 1 - weight) %>% .$transparency

ligand_order <- intersect(ligand_order, circos_links$ligand)
order_targets <- intersect(order_targets, circos_links$target)

order <- c(ligand_order, order_targets)

width_same_cell_same_ligand_type <- 0.6
width_different_cell <- 3
width_ligand_target <- 10
width_same_cell_same_target_type <- 0.6

gaps <- c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-T") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "T-specific") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "T-Macro") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macro-specific") %>% distinct(ligand) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "B-Macro") %>% distinct(ligand) %>% nrow() - 1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Positive") %>% distinct(target) %>% nrow() - 1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Negative") %>% distinct(target) %>% nrow() - 1)),
  width_ligand_target
)

circos.clear()
pdf(paste(output_dir, "ligand_target_circos_partial.pdf", sep = "/"), width = 6.4, height = 6.8)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1, order = order, link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col, transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA)
circos.clear()
dev.off()

write.csv(as.matrix(active_ligand_target_links_df_circos), file = paste(output_dir, "active_ligand_target_links_df_circos_partial.csv", sep = "/"))
write.csv(as.matrix(order_targets), file = paste(output_dir, "order_targets_partial.csv", sep = "/"))
write.csv(as.matrix(ligand_order), file = paste(output_dir, "ligand_order_partial.csv", sep = "/"))


