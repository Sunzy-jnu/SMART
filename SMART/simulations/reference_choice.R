setwd("/home/wangmg/Documents/deconvolution/Data_and_code_for_reproduction")
source("Code/utils/generate.R")
source("Code/utils/metrics.R")

library(ggplot2)
library(reshape2)


#################### Construct datasets

### The original data can be downloaded from https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data.
### The code for dividing the raw data into separate datasets according to experiments and sequencing platforms is in file "Code/utils/pbmc_preprocess.R".
### We start this experiment based on the separate datasets.


dataset_path <- "Data/raw_data/simulations/unpaired_scenario/pbmc/pbmc_datasets"
dataset_names <- c("pbmc1_10x_v2", "pbmc1_10x_v3", "pbmc1_Drop", "pbmc1_inDrops", "pbmc2_10x_v2", "pbmc2_Drop", "pbmc2_inDrops")

counts_list <- list()
labels_list <- list()

for (i in dataset_names){
    counts_list[[i]] <- read.csv(paste(dataset_path, paste(i, "counts.csv", sep = "_"), sep = "/"), header = TRUE, row.names = 1)
    labels_list[[i]] <- read.csv(paste(dataset_path, paste(i, "labels.csv", sep = "_"), sep = "/"), header = TRUE, row.names = 1)
}

keep_labels <- c("B", "CD14_monocyte", "CD16_monocyte", "CD4_T", "Cytotoxic_T")

for (i in dataset_names){
    keep_idx_i <- labels_list[[i]][, "celltype"] %in% keep_labels
    labels_list[[i]] <- labels_list[[i]][keep_idx_i, , drop = FALSE]
    counts_list[[i]] <- counts_list[[i]][keep_idx_i,]
}


save_dir <- "Data/processed_data/simulations/reference_choice"
if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
}

dataset_generation <- "pbmc1_10x_v2"
dataset_reference <- c("pbmc1_10x_v3", "pbmc1_Drop", "pbmc1_inDrops", "pbmc2_10x_v2", "pbmc2_Drop", "pbmc2_inDrops")


################### Generate pseudo spots using pbmc1_10x_v2

counts_generation <- counts_list[[dataset_generation]]
labels_generation <- labels_list[[dataset_generation]]
names_labels_generation <- rownames(labels_generation)
labels_generation <- labels_generation$celltype
names(labels_generation) <- names_labels_generation

generate_pseudo_st(counts_generation, labels_generation, n_spots = 1000, min_cells = 10, max_cells = 30, min_types = 1, max_types = NULL, temperature = 1000, do_plot = TRUE, save_dir = paste(save_dir, "st_data", sep = "/"))


################### Save other datasets as references

sc_save_dir <- paste(save_dir, "sc_data", sep = "/")
if (!dir.exists(sc_save_dir)){
    dir.create(sc_save_dir, recursive = TRUE)
}

for (i in dataset_reference){
    write.csv(counts_list[[i]], file = paste(sc_save_dir, paste(i, "counts.csv", sep = "_"), sep = "/"))
    write.csv(labels_list[[i]], file = paste(sc_save_dir, paste(i, "labels.csv", sep = "_"), sep = "/"))
}



#################### Run PANDA

### PANDA package can be installed by install.packages("Code/PANDA_1.1.0.tar.gz", repos = NULL, type = "source").
### PANDA package is also available at https://github.com/Zhangxf-ccnu/PANDA.

library(PANDA)

dataset_reference <- c("pbmc1_10x_v3", "pbmc1_Drop", "pbmc1_inDrops", "pbmc2_10x_v2", "pbmc2_Drop", "pbmc2_inDrops")

for (i in dataset_reference){
    sc_counts_path <- paste("Data/processed_data/simulations/reference_choice/sc_data", paste(i, "counts.csv", sep = "_"), sep = "/")
    sc_labels_path <- paste("Data/processed_data/simulations/reference_choice/sc_data", paste(i, "labels.csv", sep = "_"), sep = "/")
    sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
    sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
    names_sc_labels <- rownames(sc_labels)
    sc_labels <- sc_labels$celltype
    names(sc_labels) <- names_sc_labels
    sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10, n_cores = 6, save_dir = paste(c("Results/simulations/reference_choice/PANDA_results", i, "sc_results"), collapse = "/"))
    st_counts_path <- "Data/processed_data/simulations/reference_choice/st_data/st_counts.csv"
    st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
    st_results <- st_train(st_counts, sc_results = sc_results, save_dir = paste(c("Results/simulations/reference_choice/PANDA_results", i, "st_results"), collapse = "/"))
}



#################### Evaluation

dataset_reference <- c("pbmc1_10x_v3", "pbmc1_Drop", "pbmc1_inDrops", "pbmc2_10x_v2", "pbmc2_Drop", "pbmc2_inDrops")
results_dir <- "Results/simulations/reference_choice/PANDA_results"

save_dir <- "Analysis/simulations/reference_choice"
if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
}


cal_log_likelihood <- function(Y, N, omega, phi, gamma, genes){
    eps <- 1e-12
    n_spots <- dim(Y)[1]
    n_genes <- length(genes)
    Y <- Y[, genes]
    phi <- phi[, genes]
    gamma <- gamma[genes,]
    temp <- omega %*% phi
    results <- sum(diag(t(Y) %*% log((N %*% t(rep(1, n_genes))) * temp * (rep(1, n_spots) %*% t(exp(gamma))) + eps))) -
                t(N) %*% temp %*% exp(gamma)
    results <- results / (n_spots * n_genes)
    return(as.numeric(results))
}


st_counts_path <- "Data/processed_data/simulations/reference_choice/st_data/st_counts.csv"
st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
st_counts <- as.matrix(st_counts)
N <- rowSums(st_counts)
Y <- st_counts


truth_dir <- "Data/processed_data/simulations/reference_choice/st_data"
type_prop_truth <- read.csv(list.files(truth_dir, pattern = "st_type_prop.csv", full.names = TRUE), header = TRUE, row.names = 1)
unique_labels <- colnames(type_prop_truth)
type_exp_truth_list <- lapply(unique_labels, function(x){read.csv(list.files(truth_dir, pattern = paste(x, "exp.csv", sep = "_"), full.names = TRUE), header = TRUE, row.names = 1)})
names(type_exp_truth_list) <- unique_labels
type_signature_truth_list <- lapply(type_exp_truth_list, function(x){x / (rowSums(x) + 1e-12)})


keep_genes_signature_list <- list()
keep_genes_gamma_list <- list()
for (i in dataset_reference){
    sc_results <- readRDS(paste(paste(c(results_dir, i, "sc_results"), collapse = "/"), "sc_results.rds", sep = "/"))
    st_results <- readRDS(paste(paste(c(results_dir, i, "st_results"), collapse = "/"), "st_results.rds", sep = "/"))
    keep_genes_signature_list[[i]] <- colnames(sc_results$archetypes_list[[1]])
    keep_genes_gamma_list[[i]] <- rownames(st_results$gamma)
}
keep_genes_signature <- Reduce("intersect", keep_genes_signature_list)
keep_genes_gamma <- Reduce("intersect", keep_genes_gamma_list)
length(keep_genes_signature)   # 129
length(keep_genes_gamma)       # 62


results_all <- c()
results_median <- c()
log_likelihood_all <- c()

for (i in dataset_reference){
    st_results <- readRDS(paste(paste(c(results_dir, i, "st_results"), collapse = "/"), "st_results.rds", sep = "/"))
    log_likelihood <- cal_log_likelihood(Y, N, omega = st_results$omega, phi = st_results$phi, gamma = st_results$gamma, genes = keep_genes_gamma)
    log_likelihood_all <- c(log_likelihood_all, log_likelihood)
    type_prop_pred <- st_results$proportion
    type_signature_pred <- st_results$mu
    prop_pearson_spot <- cal_prop_pearson(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
    prop_RMSE_spot <- cal_prop_RMSE(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
    prop_JSD_spot <- cal_prop_JSD(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
    signature_pearson_spot <- cal_signature_pearson(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes_signature, level = "spot")
    signature_RMSE_spot <- cal_signature_RMSE(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes_signature, level = "spot")
    signature_JSD_spot <- cal_signature_JSD(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes_signature, level = "spot")
    results_all <- rbind(results_all, cbind(i, prop_pearson_spot, prop_RMSE_spot, prop_JSD_spot, signature_pearson_spot, signature_RMSE_spot, signature_JSD_spot))
    results_median <- rbind(results_median, c(i, median(prop_pearson_spot), median(prop_RMSE_spot), median(prop_JSD_spot), median(signature_pearson_spot), median(signature_RMSE_spot), median(signature_JSD_spot)))
}

colnames(results_all) <- c("reference", "prop_pearson", "prop_RMSE", "prop_JSD", "signature_pearson", "signature_RMSE", "signature_JSD")
colnames(results_median) <- c("reference", "prop_pearson", "prop_RMSE", "prop_JSD", "signature_pearson", "signature_RMSE", "signature_JSD")
names(log_likelihood_all) <- dataset_reference

saveRDS(results_all, file = paste(save_dir, "results_all.rds", sep = "/"))
saveRDS(results_median, file = paste(save_dir, "results_median.rds", sep = "/"))
saveRDS(log_likelihood_all, file = paste(save_dir, "log_likelihood_all.rds", sep = "/"))


#################### Plot

# evaluation_prop
plot_df <- as.data.frame(results_all[, c("reference", "prop_pearson", "prop_RMSE", "prop_JSD")])
plot_df <- melt(plot_df, id.vars = "reference", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$reference <- factor(plot_df$reference, levels = dataset_reference)

fname <- "evaluation_prop"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 5)
fig <- ggplot(data = plot_df, aes(x = reference, y = as.numeric(value), fill = reference)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Reference", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# evaluation_signature
plot_df <- as.data.frame(results_all[, c("reference", "signature_pearson", "signature_RMSE", "signature_JSD")])
plot_df <- melt(plot_df, id.vars = "reference", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$reference <- factor(plot_df$reference, levels = dataset_reference)

fname <- "evaluation_signature"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 5)
fig <- ggplot(data = plot_df, aes(x = reference, y = as.numeric(value), fill = reference)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Reference", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()



# likelihood_prop
plot_df <- as.data.frame(cbind(dataset_reference, log_likelihood_all, results_median[, c("prop_pearson", "prop_RMSE", "prop_JSD")]))
colnames(plot_df) <- c("reference", "log_likelihood", "prop_pearson", "prop_RMSE", "prop_JSD")
plot_df <- melt(plot_df, id.vars = c("reference", "log_likelihood"), variable.name = "metric", value.name = "value")
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$reference <- factor(plot_df$reference, levels = dataset_reference)

fname <- "likelihood_prop"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 4.8)
fig <- ggplot(data = plot_df, aes(x = as.numeric(value), y = as.numeric(log_likelihood))) + 
        geom_point(aes(shape = reference, color = reference), size = 4) +
        labs(x = "", y = "Average log likelihood", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
              panel.grid.major.x = element_line(color = "grey", linewidth = 0.3),
              panel.grid.major.y = element_line(color = "grey", linewidth = 0.3),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12, color = "black"),
              axis.text.y = element_text(size = 12, color = "black"),
              axis.title.x = element_text(size = 13),
              axis.title.y = element_text(size = 13),
              strip.text.x = element_text(size = 13),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") +
        scale_shape_manual(values = c(16, 17, 15, 18, 1, 2, 0, 5, 4)) +
        facet_wrap(~metric, ncol = 3, scales="free_x")
print(fig)
dev.off()



# likelihood_signature
plot_df <- as.data.frame(cbind(dataset_reference, log_likelihood_all, results_median[, c("signature_pearson", "signature_RMSE", "signature_JSD")]))
colnames(plot_df) <- c("reference", "log_likelihood", "signature_pearson", "signature_RMSE", "signature_JSD")
plot_df <- melt(plot_df, id.vars = c("reference", "log_likelihood"), variable.name = "metric", value.name = "value")
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$reference <- factor(plot_df$reference, levels = dataset_reference)

fname <- "likelihood_signature"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 10, height = 4.8)
fig <- ggplot(data = plot_df, aes(x = as.numeric(value), y = as.numeric(log_likelihood))) + 
        geom_point(aes(shape = reference, color = reference), size = 4) +
        labs(x = "", y = "Average log likelihood", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
              panel.grid.major.x = element_line(color = "grey", linewidth = 0.3),
              panel.grid.major.y = element_line(color = "grey", linewidth = 0.3),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12, color = "black"),
              axis.text.y = element_text(size = 12, color = "black"),
              axis.title.x = element_text(size = 13),
              axis.title.y = element_text(size = 13),
              strip.text.x = element_text(size = 13),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") +
        scale_shape_manual(values = c(16, 17, 15, 18, 1, 2, 0, 5, 4)) +
        facet_wrap(~metric, ncol = 3, scales="free_x")
print(fig)
dev.off()


