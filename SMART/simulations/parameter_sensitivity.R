setwd("/home/wangmg/Documents/deconvolution/Data_and_code_for_reproduction")
source("Code/utils/metrics.R")

library(ggplot2)
library(reshape2)



n_archetypes_list <- c(5, 10, 15, 20)
sigma_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)


#################### Paired

### PANDA package can be installed by install.packages("Code/PANDA_1.1.0.tar.gz", repos = NULL, type = "source").
### PANDA package is also available at https://github.com/Zhangxf-ccnu/PANDA.

library(PANDA)

sc_counts_path <- "Data/processed_data/simulations/paired_scenario/sc_data/counts_validation.csv"
sc_labels_path <- "Data/processed_data/simulations/paired_scenario/sc_data/labels_validation.csv"

sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels


for (n_archetypes in n_archetypes_list){
    sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = n_archetypes, n_cores = 6, save_dir = paste("Results/simulations/parameter_sensitivity/paired/PANDA_results/sc_results", paste("n_archetypes", n_archetypes, sep = "_"), sep = "/"))
    st_counts_path <- paste(paste("Data/processed_data/simulations/paired_scenario/st_data", "uniform_ST", sep = "/"), "st_counts.csv", sep = "/")
    st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
    for (sigma in sigma_list){
        st_results <- st_train(st_counts, sc_results = sc_results, sigma = sigma, save_dir = paste("Results/simulations/parameter_sensitivity/paired/PANDA_results/st_results", paste(c("n_archetypes", n_archetypes, "sigma", substring(sigma, 3)), collapse = "_"), sep = "/"))
    }
}



#################### Unpaired

### PANDA package can be installed by install.packages("Code/PANDA_1.1.0.tar.gz", repos = NULL, type = "source").
### PANDA package is also available at https://github.com/Zhangxf-ccnu/PANDA.

library(PANDA)

sc_counts_path <- "Data/processed_data/simulations/unpaired_scenario/sc_data/pbmc2_inDrops_counts.csv"
sc_labels_path <- "Data/processed_data/simulations/unpaired_scenario/sc_data/pbmc2_inDrops_labels.csv"

sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels


for (n_archetypes in n_archetypes_list){
    sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = n_archetypes, n_cores = 6, save_dir = paste("Results/simulations/parameter_sensitivity/unpaired/PANDA_results/sc_results", paste("n_archetypes", n_archetypes, sep = "_"), sep = "/"))
    st_counts_path <- paste(paste("Data/processed_data/simulations/unpaired_scenario/st_data", "uniform_ST", sep = "/"), "st_counts.csv", sep = "/")
    st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
    for (sigma in sigma_list){
        st_results <- st_train(st_counts, sc_results = sc_results, sigma = sigma, save_dir = paste("Results/simulations/parameter_sensitivity/unpaired/PANDA_results/st_results", paste(c("n_archetypes", n_archetypes, "sigma", substring(sigma, 3)), collapse = "_"), sep = "/"))
    }
}



###################################################################################################
# Evaluation
# Paired

save_dir <- "Analysis/simulations/parameter_sensitivity/paired"
if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
}

paired_results_all <- c()
paired_results_median <- c()

truth_dir <- paste("Data/processed_data/simulations/paired_scenario/st_data", "uniform_ST", sep = "/")
type_prop_truth <- read.csv(list.files(truth_dir, pattern = "st_type_prop.csv", full.names = TRUE), header = TRUE, row.names = 1)
unique_labels <- colnames(type_prop_truth)
type_exp_truth_list <- lapply(unique_labels, function(x){read.csv(list.files(truth_dir, pattern = paste(x, "exp.csv", sep = "_"), full.names = TRUE), header = TRUE, row.names = 1)})
names(type_exp_truth_list) <- unique_labels
type_signature_truth_list <- lapply(type_exp_truth_list, function(x){x / (rowSums(x) + 1e-12)})

keep_genes <- readRDS("Analysis/simulations/paired_scenario/comparison_results/uniform_ST/keep_genes.rds")

for (n_archetypes in n_archetypes_list){
    sc_results <- readRDS(paste(paste("Results/simulations/parameter_sensitivity/paired/PANDA_results/sc_results", paste("n_archetypes", n_archetypes, sep = "_"), sep = "/"), "sc_results.rds", sep = "/"))
    for (sigma in sigma_list){
        st_results <- readRDS(paste(paste("Results/simulations/parameter_sensitivity/paired/PANDA_results/st_results", paste(c("n_archetypes", n_archetypes, "sigma", substring(sigma, 3)), collapse = "_"), sep = "/"), "st_results.rds", sep = "/"))
        type_prop_pred <- st_results$proportion
        type_signature_pred <- st_results$mu
        prop_pearson_spot <- cal_prop_pearson(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
        prop_RMSE_spot <- cal_prop_RMSE(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
        prop_JSD_spot <- cal_prop_JSD(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
        signature_pearson_spot <- cal_signature_pearson(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")
        signature_RMSE_spot <- cal_signature_RMSE(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")
        signature_JSD_spot <- cal_signature_JSD(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")
        paired_results_all <- rbind(paired_results_all, cbind(n_archetypes, sigma, prop_pearson_spot, prop_RMSE_spot, prop_JSD_spot, signature_pearson_spot, signature_RMSE_spot, signature_JSD_spot))
        paired_results_median <- rbind(paired_results_median, c(n_archetypes, sigma, median(prop_pearson_spot), median(prop_RMSE_spot), median(prop_JSD_spot), median(signature_pearson_spot), median(signature_RMSE_spot), median(signature_JSD_spot), sc_results$run_time, st_results$run_time))
    }
}

colnames(paired_results_all) <- c("n_archetypes", "sigma", "prop_pearson", "prop_RMSE", "prop_JSD", "signature_pearson", "signature_RMSE", "signature_JSD")
colnames(paired_results_median) <- c("n_archetypes", "sigma", "prop_pearson", "prop_RMSE", "prop_JSD", "signature_pearson", "signature_RMSE", "signature_JSD", "sc_time", "st_time")

saveRDS(paired_results_all, file = paste(save_dir, "paired_results_all.rds", sep = "/"))
saveRDS(paired_results_median, file = paste(save_dir, "paired_results_median.rds", sep = "/"))


# Unpaired

save_dir <- "Analysis/simulations/parameter_sensitivity/unpaired"
if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
}

unpaired_results_all <- c()
unpaired_results_median <- c()

truth_dir <- paste("Data/processed_data/simulations/unpaired_scenario/st_data", "uniform_ST", sep = "/")
type_prop_truth <- read.csv(list.files(truth_dir, pattern = "st_type_prop.csv", full.names = TRUE), header = TRUE, row.names = 1)
unique_labels <- colnames(type_prop_truth)
type_exp_truth_list <- lapply(unique_labels, function(x){read.csv(list.files(truth_dir, pattern = paste(x, "exp.csv", sep = "_"), full.names = TRUE), header = TRUE, row.names = 1)})
names(type_exp_truth_list) <- unique_labels
type_signature_truth_list <- lapply(type_exp_truth_list, function(x){x / (rowSums(x) + 1e-12)})

keep_genes <- readRDS("Analysis/simulations/unpaired_scenario/comparison_results/uniform_ST/keep_genes.rds")

for (n_archetypes in n_archetypes_list){
    sc_results <- readRDS(paste(paste("Results/simulations/parameter_sensitivity/unpaired/PANDA_results/sc_results", paste("n_archetypes", n_archetypes, sep = "_"), sep = "/"), "sc_results.rds", sep = "/"))
    for (sigma in sigma_list){
        st_results <- readRDS(paste(paste("Results/simulations/parameter_sensitivity/unpaired/PANDA_results/st_results", paste(c("n_archetypes", n_archetypes, "sigma", substring(sigma, 3)), collapse = "_"), sep = "/"), "st_results.rds", sep = "/"))
        type_prop_pred <- st_results$proportion
        type_signature_pred <- st_results$mu
        prop_pearson_spot <- cal_prop_pearson(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
        prop_RMSE_spot <- cal_prop_RMSE(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
        prop_JSD_spot <- cal_prop_JSD(prediction = type_prop_pred, truth = type_prop_truth, level = "spot")
        signature_pearson_spot <- cal_signature_pearson(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")
        signature_RMSE_spot <- cal_signature_RMSE(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")
        signature_JSD_spot <- cal_signature_JSD(prediction_list = type_signature_pred, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")
        unpaired_results_all <- rbind(unpaired_results_all, cbind(n_archetypes, sigma, prop_pearson_spot, prop_RMSE_spot, prop_JSD_spot, signature_pearson_spot, signature_RMSE_spot, signature_JSD_spot))
        unpaired_results_median <- rbind(unpaired_results_median, c(n_archetypes, sigma, median(prop_pearson_spot), median(prop_RMSE_spot), median(prop_JSD_spot), median(signature_pearson_spot), median(signature_RMSE_spot), median(signature_JSD_spot), sc_results$run_time, st_results$run_time))
    }
}

colnames(unpaired_results_all) <- c("n_archetypes", "sigma", "prop_pearson", "prop_RMSE", "prop_JSD", "signature_pearson", "signature_RMSE", "signature_JSD")
colnames(unpaired_results_median) <- c("n_archetypes", "sigma", "prop_pearson", "prop_RMSE", "prop_JSD", "signature_pearson", "signature_RMSE", "signature_JSD", "sc_time", "st_time")

saveRDS(unpaired_results_all, file = paste(save_dir, "unpaired_results_all.rds", sep = "/"))
saveRDS(unpaired_results_median, file = paste(save_dir, "unpaired_results_median.rds", sep = "/"))



###################################################################################################
# Plot
# Paired
save_dir <- "Analysis/simulations/parameter_sensitivity/paired"

# paired_prop_sigma_3
plot_df <- as.data.frame(paired_results_all[paired_results_all[, "sigma"] == 0.3, c("n_archetypes", "prop_pearson", "prop_RMSE", "prop_JSD")])
plot_df <- melt(plot_df, id.vars = "n_archetypes", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)

fname <- "paired_prop_sigma_3"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# paired_signature_sigma_3
plot_df <- as.data.frame(paired_results_all[paired_results_all[, "sigma"] == 0.3, c("n_archetypes", "signature_pearson", "signature_RMSE", "signature_JSD")])
plot_df <- melt(plot_df, id.vars = "n_archetypes", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)

fname <- "paired_signature_sigma_3"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# paired_prop_n_archetypes_10
plot_df <- as.data.frame(paired_results_all[paired_results_all[, "n_archetypes"] == 10, c("sigma", "prop_pearson", "prop_RMSE", "prop_JSD")])
plot_df <- melt(plot_df, id.vars = "sigma", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "paired_prop_n_archetypes_10"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = sigma, y = as.numeric(value), fill = sigma)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Sigma", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# paired_signature_n_archetypes_10
plot_df <- as.data.frame(paired_results_all[paired_results_all[, "n_archetypes"] == 10, c("sigma", "signature_pearson", "signature_RMSE", "signature_JSD")])
plot_df <- melt(plot_df, id.vars = "sigma", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "paired_signature_n_archetypes_10"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = sigma, y = as.numeric(value), fill = sigma)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Sigma", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# paired_integrated_prop
plot_df <- as.data.frame(paired_results_all[, c("n_archetypes", "sigma", "prop_pearson", "prop_RMSE", "prop_JSD")])
plot_df <- melt(plot_df, id.vars = c("n_archetypes", "sigma"), variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "paired_integrated_prop"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 12)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.text.y = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_grid(sigma~metric)
print(fig)
dev.off()


# paired_integrated_signature
plot_df <- as.data.frame(paired_results_all[, c("n_archetypes", "sigma", "signature_pearson", "signature_RMSE", "signature_JSD")])
plot_df <- melt(plot_df, id.vars = c("n_archetypes", "sigma"), variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "paired_integrated_signature"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 12)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.text.y = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_grid(sigma~metric)
print(fig)
dev.off()



# Unpaired
save_dir <- "Analysis/simulations/parameter_sensitivity/unpaired"

# unpaired_prop_sigma_3
plot_df <- as.data.frame(unpaired_results_all[unpaired_results_all[, "sigma"] == 0.3, c("n_archetypes", "prop_pearson", "prop_RMSE", "prop_JSD")])
plot_df <- melt(plot_df, id.vars = "n_archetypes", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)

fname <- "unpaired_prop_sigma_3"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# unpaired_signature_sigma_3
plot_df <- as.data.frame(unpaired_results_all[unpaired_results_all[, "sigma"] == 0.3, c("n_archetypes", "signature_pearson", "signature_RMSE", "signature_JSD")])
plot_df <- melt(plot_df, id.vars = "n_archetypes", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)

fname <- "unpaired_signature_sigma_3"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# unpaired_prop_n_archetypes_10
plot_df <- as.data.frame(unpaired_results_all[unpaired_results_all[, "n_archetypes"] == 10, c("sigma", "prop_pearson", "prop_RMSE", "prop_JSD")])
plot_df <- melt(plot_df, id.vars = "sigma", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "unpaired_prop_n_archetypes_10"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = sigma, y = as.numeric(value), fill = sigma)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Sigma", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# unpaired_signature_n_archetypes_10
plot_df <- as.data.frame(unpaired_results_all[unpaired_results_all[, "n_archetypes"] == 10, c("sigma", "signature_pearson", "signature_RMSE", "signature_JSD")])
plot_df <- melt(plot_df, id.vars = "sigma", variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "unpaired_signature_n_archetypes_10"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 4)
fig <- ggplot(data = plot_df, aes(x = sigma, y = as.numeric(value), fill = sigma)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Sigma", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_wrap(~metric, ncol = 3, scales = "free")
print(fig)
dev.off()


# unpaired_integrated_prop
plot_df <- as.data.frame(unpaired_results_all[, c("n_archetypes", "sigma", "prop_pearson", "prop_RMSE", "prop_JSD")])
plot_df <- melt(plot_df, id.vars = c("n_archetypes", "sigma"), variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("prop_pearson", "prop_RMSE", "prop_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "unpaired_integrated_prop"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 12)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.text.y = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_grid(sigma~metric)
print(fig)
dev.off()


# unpaired_integrated_signature
plot_df <- as.data.frame(unpaired_results_all[, c("n_archetypes", "sigma", "signature_pearson", "signature_RMSE", "signature_JSD")])
plot_df <- melt(plot_df, id.vars = c("n_archetypes", "sigma"), variable.name = "metric", value.name = "value")
rownames(plot_df) <- NULL
plot_df$metric <- factor(plot_df$metric, levels = c("signature_pearson", "signature_RMSE", "signature_JSD"), labels = c("PCC", "RMSE", "JSD"))
plot_df$n_archetypes <- factor(plot_df$n_archetypes, levels = n_archetypes_list)
plot_df$sigma <- factor(plot_df$sigma, levels = sigma_list)

fname <- "unpaired_integrated_signature"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 9, height = 12)
fig <- ggplot(data = plot_df, aes(x = n_archetypes, y = as.numeric(value), fill = n_archetypes)) + 
        geom_boxplot(size = 0.2, outlier.size = 0.01) + 
        labs(x = "Number of archetypes", y = "", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
              plot.background  = element_blank(),
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12), 
              strip.text.x = element_text(size = 12),
              strip.text.y = element_text(size = 12),
              strip.background = element_rect(fill = "grey90"),
              legend.position = "bottom") + 
        facet_grid(sigma~metric)
print(fig)
dev.off()


