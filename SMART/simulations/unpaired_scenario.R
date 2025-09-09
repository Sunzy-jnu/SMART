setwd("D:/Rshuju/shiyan")
source("SMART/utils/generate.R")
source("SMART/utils/metrics.R")
source("SMART/utils/utils.R")
source("SMART/utils/plot.R")


#### Unpaired scenario


### The original data can be downloaded from https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data.
### The code for dividing the raw data into separate datasets according to experiments and sequencing platforms is in file "SMART/utils/pbmc_preprocess.R".


#################### Visualize the generation dataset and the validation dataset

sc_counts_generation_path = "SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_counts.csv"
sc_labels_generation_path = "SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_labels.csv"
sc_counts_validation_path = "SMART_data/unpaired_scenario/sc_data/pbmc2_inDrops_counts.csv"
sc_labels_validation_path = "SMART_data/unpaired_scenario/sc_data/pbmc2_inDrops_labels.csv"

sc_counts_generation <- read.csv(sc_counts_generation_path, header = TRUE, row.names = 1)
sc_labels_generation <- read.csv(sc_labels_generation_path, header = TRUE, row.names = 1)
names_sc_labels_generation <- rownames(sc_labels_generation)
sc_labels_generation <- sc_labels_generation$celltype
names(sc_labels_generation) <- names_sc_labels_generation

sc_counts_validation <- read.csv(sc_counts_validation_path, header = TRUE, row.names = 1)
sc_labels_validation <- read.csv(sc_labels_validation_path, header = TRUE, row.names = 1)
names_sc_labels_validation <- rownames(sc_labels_validation)
sc_labels_validation <- sc_labels_validation$celltype
names(sc_labels_validation) <- names_sc_labels_validation


sc_counts <- rbind(sc_counts_generation, sc_counts_validation)
sc_labels <- c(sc_labels_generation, sc_labels_validation)
sc_domain <- c(rep("Generation", length(sc_labels_generation)), rep("Validation", length(sc_labels_validation)))

sc_hvgs <- find_hvgs(t(sc_counts), n_hvgs = 2000)
sc_counts_norm <- sc_counts / rowSums(sc_counts)

plot_df <- as.data.frame(prcomp(sc_counts_norm[, sc_hvgs], center = TRUE, scale. = TRUE)$x[, 1:min(30, dim(sc_counts_norm)[1])])
umap_settings <- umap::umap.defaults
umap_settings$n_neighbors <- min(15, dim(sc_counts_norm)[1])
plot_df <- as.data.frame(umap::umap(as.matrix(plot_df), config = umap_settings)$layout)
plot_df <- as.data.frame(cbind(plot_df, sc_labels, sc_domain))
colnames(plot_df) <- c("dim1", "dim2", "label", "domain")
saveRDS(plot_df, file = "SMART_results/simulation/unpaired_scenario/visualization/visualization_plot_df.rds")

library(ggplot2)

pdf(paste("SMART_results/simulation/unpaired_scenario/visualization", paste("UMAP_labels", "pdf", sep = "."), sep = "/"), width = 6, height = 4)
fig <- ggplot(data = plot_df, aes(x = dim1, y = dim2, color = as.factor(label))) +
        geom_point(size = 0.2) +
        # coord_equal() +
        coord_fixed(ratio = 1.091815) +
        labs(title = "UMAP_labels", x = "UMAP 1", y = "UMAP 2") +
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
              plot.background  = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_text(size = 10),
              legend.title = element_text(size = 12),
              legend.key = element_rect(color = "transparent", fill = "white"),
              legend.key.size = unit(0.45, "cm"),
              legend.position = "right") +
        guides(color = guide_legend(title = "Cell type", override.aes = list(size=2)))
print(fig)
dev.off()

pdf(paste("SMART_results/simulation/unpaired_scenario/visualization", paste("UMAP_domain", "pdf", sep = "."), sep = "/"), width = 6, height = 4)
fig <- ggplot(data = plot_df, aes(x = dim1, y = dim2, color = as.factor(domain))) +
        geom_point(size = 0.2) +
        # coord_equal() +
        coord_fixed(ratio = 1.091815) +
        labs(title = "UMAP_domain", x = "UMAP 1", y = "UMAP 2") +
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
              plot.background  = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_text(size = 10),
              legend.title = element_text(size = 12),
              legend.key = element_rect(color = "transparent", fill = "white"),
              legend.key.size = unit(0.45, "cm"),
              legend.position = "right") +
        guides(color = guide_legend(title = "Domain", override.aes = list(size=2)))
print(fig)
dev.off()







################### Generate pseudo spots

sc_counts_path <- "SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_counts.csv"
sc_labels_path <- "SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_labels.csv"

sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels

simulation_settings <- list(ST_uniform = list(tech_type = "ST", coord_type = "uniform", min_cells = 10, max_cells = 30, temperature = 1000 ),
  ST_boundary = list(tech_type = "ST", coord_type = "boundary", min_cells = 10, max_cells = 30, temperature = 0.1),
  Visium_uniform = list(tech_type = "Visium", coord_type = "uniform", min_cells = 1, max_cells = 10, temperature = 1000),
  Visium_boundary = list(tech_type = "Visium", coord_type = "boundary", min_cells = 1, max_cells = 10, temperature = 0.1)
)

for (setting in names(simulation_settings)) {
  coords <- generate_coords(1000, interval = 100, type = simulation_settings[[setting]]$coord_type)
  d <- as.matrix(dist(coords, method = "euclidean"))
  save_dir <- file.path("SMART_data/unpaired_scenario/st_data", setting)
  generate_pseudo_st(sc_counts = sc_counts, sc_labels = sc_labels, coords = coords, n_spots = 1000, min_cells = simulation_settings[[setting]]$min_cells, max_cells = simulation_settings[[setting]]$max_cells, temperature = simulation_settings[[setting]]$temperature, tech_type = simulation_settings[[setting]]$tech_type, save_dir = save_dir)
}

#################### Run SMART

### SMART package can be installed by install.packages("SMART/SMART/R/SMART_1.1.0.tar.gz", repos = NULL, type = "source").
### SMART package is also available at https://github.com/Zhangxf-ccnu/PANDA.

library(SMART)

sc_counts_path <- "SMART_data/unpaired_scenario/sc_data/pbmc2_inDrops_counts.csv"
sc_labels_path <- "SMART_data/unpaired_scenario/sc_data/pbmc2_inDrops_labels.csv"

sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels

sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10, n_cores = 6, save_dir = "SMART_results/simulation/unpaired_scenario/SMART_results/sc_results")

for (i in list.files("SMART_data/unpaired_scenario/st_data")){
  st_counts_path <- file.path("SMART_data", "unpaired_scenario", "st_data", i, "st_counts.csv")
  st_coords_path <- file.path("SMART_data", "unpaired_scenario", "st_data", i, "st_coords.csv")

  st_counts <- as.matrix(read.csv(st_counts_path, header = TRUE, row.names = 1))
  spatial_coords <- as.matrix(read.csv(st_coords_path, header = TRUE, row.names = 1))

  sc_genes <- Reduce(intersect, lapply(sc_results$archetypes_list, colnames))
  common_genes <- intersect(colnames(st_counts), sc_genes)
  st_counts <- st_counts[, common_genes]

  st_results <- st_train(
    st_counts = st_counts,
    sc_results = sc_results,
    spatial_coords = spatial_coords,
    lambda_spatial = 0.8,
    neighbor_radius = 5,
    save_mu_csv = TRUE,
    save_dir = paste("SMART_results/simulation/unpaired_scenario/SMART_results/st_results", i, sep = "/"))
}

#####Run PANDA
library(PANDA)

sc_counts_path <- "SMART_data/unpaired_scenario/sc_data/pbmc2_inDrops_counts.csv"
sc_labels_path <- "SMART_data/unpaired_scenario/sc_data/pbmc2_inDrops_labels.csv"

sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
names_sc_labels <- rownames(sc_labels)
sc_labels <- sc_labels$celltype
names(sc_labels) <- names_sc_labels

sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10, n_cores = 6, save_dir = "SMART_results/simulation/unpaired_scenario/comparison_methods/PANDA_results/sc_results")

for (i in list.files("SMART_data/unpaired_scenario/st_data")){
  st_counts_path <- paste(paste("SMART_data/unpaired_scenario/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
  st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
  st_results <- st_train(st_counts, sc_results = sc_results, save_mu_csv = TRUE, save_dir = paste("SMART_results/simulation/unpaired_scenario/comparison_methods/PANDA_results/st_results", i, sep = "/"))
}


#################### Compare with other methods
all_methods <- c("SMART", "PANDA", "cell2location", "DestVI", "RCTD", "SPOTlight", "stereoscope")

use_color <- scales::hue_pal()(length(all_methods))
names(use_color) = all_methods

suffix_all <- c("ST_boundary", "ST_uniform", "Visium_boundary", "Visium_uniform")

integrated_results_list <- list()

for (suffix in suffix_all){
  cat(suffix, "\n")

  truth_dir <- paste("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data", suffix, sep = "/")
  results_dir <- list(SMART = paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/SMART_results/st_results", suffix, sep = "/"),
                      PANDA = paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/PANDA_results/st_results", suffix, sep = "/"),
                      cell2location = paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/cell2location_results", suffix, sep = "/"),
                      DestVI = paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/DestVI_results", suffix, sep = "/"),
                      RCTD = paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/RCTD_results", suffix, sep = "/"),
                      SPOTlight = paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/SPOTlight_results", suffix, sep = "/"),
                      stereoscope = paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/stereoscope_scvi_results", suffix, sep = "/"))
  output_dir <- paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_results", suffix, sep = "/")

  if (!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }

    ### Cell type proportion

    type_prop_truth <- read.csv(list.files(truth_dir, pattern = "st_type_prop.csv", full.names = TRUE), header = TRUE, row.names = 1)
    type_prop_pred_list <- lapply(all_methods, function(x){read.csv(list.files(results_dir[[x]], pattern = "prop.csv", full.names = TRUE), header = TRUE, row.names = 1)})
    names(type_prop_pred_list) <- all_methods

    saveRDS(type_prop_pred_list, file = paste(output_dir, "type_prop_pred_list.rds", sep = "/"))

    prop_pearson_spot_list <- lapply(type_prop_pred_list, function(x){cal_prop_pearson(prediction = x, truth = type_prop_truth, level = "spot")})
    prop_RMSE_spot_list <- lapply(type_prop_pred_list, function(x){cal_prop_RMSE(prediction = x, truth = type_prop_truth, level = "spot")})
    prop_JSD_spot_list <- lapply(type_prop_pred_list, function(x){cal_prop_JSD(prediction = x, truth = type_prop_truth, level = "spot")})

    plot_evaluation(prop_pearson_spot_list, prop_RMSE_spot_list, prop_JSD_spot_list, cname = c("PCC", "RMSE", "JSD"), mname = all_methods, use_color = use_color, average = FALSE, fname = "prop_evaluation_spot_box", save_dir = output_dir, n_cols = 3, width = 8, height = 3.6)


    ### Cell type signature

    unique_labels <- colnames(type_prop_truth)
    type_exp_truth_list <- lapply(unique_labels, function(x){read.csv(list.files(truth_dir, pattern = paste(x, "exp.csv", sep = "_"), full.names = TRUE), header = TRUE, row.names = 1)})
    names(type_exp_truth_list) <- unique_labels

    type_signature_truth_list <- lapply(type_exp_truth_list, function(x){x / (rowSums(x) + 1e-12)})

    type_signature_pred_list <- list()
    type_signature_pred_list[["SMART"]] <- lapply(unique_labels, function(x){read.csv(list.files(results_dir[["SMART"]], pattern = paste(x, "mu.csv", sep = "_"), full.names = TRUE), header = TRUE, row.names = 1)})
    names(type_signature_pred_list[["SMART"]]) <- unique_labels
    type_signature_pred_list[["PANDA"]] <- lapply(unique_labels, function(x){read.csv(list.files(results_dir[["PANDA"]], pattern = paste(x, "mu.csv", sep = "_"), full.names = TRUE), header = TRUE, row.names = 1)})
    names(type_signature_pred_list[["PANDA"]]) <- unique_labels
    type_signature_pred_list[["DestVI"]] <- lapply(unique_labels, function(x){read.csv(list.files(results_dir[["DestVI"]], pattern = paste(x, "signature.csv", sep = "_"), full.names = TRUE), header = TRUE, row.names = 1)})
    names(type_signature_pred_list[["DestVI"]]) <- unique_labels
    for (i in c("cell2location", "RCTD", "SPOTlight", "stereoscope")){
      temp <- as.matrix(read.csv(list.files(results_dir[[i]], pattern = "signature.csv", full.names = TRUE), header = TRUE, row.names = 1))
      temp <- lapply(unique_labels, function(x){
        temp2 <- rep(1, nrow(type_signature_truth_list[[x]])) %*% t(temp[x,])
        rownames(temp2) <- rownames(type_signature_truth_list[[x]])
        temp2
      })
      names(temp) <- unique_labels
      type_signature_pred_list[[i]] <- temp
    }
    type_signature_pred_list <- type_signature_pred_list[all_methods]

    saveRDS(type_signature_pred_list, file = paste(output_dir, "type_signature_pred_list.rds", sep = "/"))
    # type_signature_pred_list <- readRDS(paste(output_dir, "type_signature_pred_list.rds", sep = "/"))

    keep_genes <- Reduce("intersect", lapply(type_signature_pred_list, function(x){colnames(x[[1]])}))
    saveRDS(keep_genes, file = paste(output_dir, "keep_genes.rds", sep = "/"))

    signature_pearson_spot_list <- lapply(type_signature_pred_list, function(x){cal_signature_pearson(prediction_list = x, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")})
    signature_RMSE_spot_list <- lapply(type_signature_pred_list, function(x){cal_signature_RMSE(prediction_list = x, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")})
    signature_JSD_spot_list <- lapply(type_signature_pred_list, function(x){cal_signature_JSD(prediction_list = x, truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0.1, keep_genes = keep_genes, level = "spot")})

    plot_evaluation(signature_pearson_spot_list, signature_RMSE_spot_list, signature_JSD_spot_list, cname = c("PCC", "RMSE", "JSD"), mname = all_methods, use_color = use_color, average = FALSE, fname = "signature_evaluation_spot_box", save_dir = output_dir, n_cols = 3, width = 8, height = 3.6)


    prop_signature_pearson <- cal_signature_pearson_individual(prediction_list = type_signature_pred_list[["SMART"]], truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0, keep_genes = keep_genes)
    prop_signature_RMSE <- cal_signature_RMSE_individual(prediction_list = type_signature_pred_list[["SMART"]], truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0, keep_genes = keep_genes)
    prop_signature_JSD <- cal_signature_JSD_individual(prediction_list = type_signature_pred_list[["SMART"]], truth_list = type_signature_truth_list, prop_truth = type_prop_truth, min_prop = 0, keep_genes = keep_genes)

    plot_prop_signature_acc(prop_signature_pearson, fname = "prop_signature_pearson", save_dir = output_dir, size = 0.5, color = "deepskyblue3", n_cols = 3, x_breaks = seq(0, 1, by = 0.1), width = 10, height = 11)
    plot_prop_signature_acc(prop_signature_RMSE, fname = "prop_signature_RMSE", save_dir = output_dir, size = 0.5, color = "deepskyblue3", n_cols = 3, x_breaks = seq(0, 1, by = 0.1), width = 10, height = 11)
    plot_prop_signature_acc(prop_signature_JSD, fname = "prop_signature_JSD", save_dir = output_dir, size = 0.5, color = "deepskyblue3", n_cols = 3, x_breaks = seq(0, 1, by = 0.1), width = 10, height = 11)


    rm(type_signature_truth_list)
    rm(type_signature_pred_list)
    gc()



    ### Save results

    saveRDS(list(prop_pearson_spot_list = prop_pearson_spot_list,
                 prop_RMSE_spot_list = prop_RMSE_spot_list,
                 prop_JSD_spot_list = prop_JSD_spot_list,
                 signature_pearson_spot_list = signature_pearson_spot_list,
                 signature_RMSE_spot_list = signature_RMSE_spot_list,
                 signature_JSD_spot_list = signature_JSD_spot_list), file = paste(output_dir, "integrated_results.rds", sep = "/"))

    all_metrics <- c("prop_pearson_spot",
                     "prop_RMSE_spot",
                     "prop_JSD_spot",
                     "signature_pearson_spot",
                     "signature_RMSE_spot",
                     "signature_JSD_spot")

    integrated_results <- matrix(0, nrow = length(all_metrics), ncol = length(all_methods), dimnames = list(all_metrics, all_methods))
    for (metric in all_metrics){
        temp <- sapply(get(paste(metric, "list", sep = "_")), function(x){median(x)})
        for (method in names(temp)){
            integrated_results[metric, method] <- temp[method]
        }
    }
    write.csv(integrated_results, file = paste(output_dir, "integrated_results.csv", sep = "/"))
    integrated_results_list[[suffix]] <- integrated_results
}

saveRDS(integrated_results_list, file = "D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_results/integrated_results_list.rds")

#integrated_results_list <- readRDS("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_results/integrated_results_list.rds")


prop_metrics <- c("prop_pearson_spot", "prop_RMSE_spot", "prop_JSD_spot")
signature_metrics <- c("signature_pearson_spot", "signature_RMSE_spot", "signature_JSD_spot")


rank_matrix_prop <- c()
for (i in suffix_all){
    for (j in prop_metrics){
        if (j == "prop_pearson_spot"){
            temp <- length(integrated_results_list[[i]][j,]) + 1 - rank(integrated_results_list[[i]][j,])
        } else{
            temp <- rank(integrated_results_list[[i]][j,])
        }
        rank_matrix_prop <- rbind(rank_matrix_prop, temp)
    }
}
rownames(rank_matrix_prop) <- paste(rep(suffix_all, each = length(prop_metrics)), rep(prop_metrics, times = length(suffix_all)), sep = "_")

rank_matrix_signature <- c()
for (i in suffix_all){
    for (j in signature_metrics){
        if (j == "signature_pearson_spot"){
            temp <- length(integrated_results_list[[i]][j,]) + 1 - rank(integrated_results_list[[i]][j,])
        } else{
            temp <- rank(integrated_results_list[[i]][j,])
        }
        rank_matrix_signature <- rbind(rank_matrix_signature, temp)
    }
}
rownames(rank_matrix_signature) <- paste(rep(suffix_all, each = length(signature_metrics)), rep(signature_metrics, times = length(suffix_all)), sep = "_")


write.csv(rank_matrix_prop, file = "SMART_results/simulation/unpaired_scenario/comparison_results/rank_matrix_prop.csv")
write.csv(rank_matrix_signature, file = "SMART_results/simulation/unpaired_scenario/comparison_results/rank_matrix_signature.csv")

mean_rank_prop <- colMeans(rank_matrix_prop)
mean_rank_signature <- colMeans(rank_matrix_signature)


plot_rank(mean_rank_prop, mname = all_methods, use_color = use_color, fname = "mean_rank_prop", save_dir = "SMART_results/simulation/unpaired_scenario/comparison_results", width = 7.5, height = 4)
plot_rank(mean_rank_signature, mname = all_methods, use_color = use_color, fname = "mean_rank_signature", save_dir = "SMART_results/simulation/unpaired_scenario/comparison_results", width = 7.5, height = 4)



metrics = c("PCC", "RMSE", "JSD")
library(ggplot2)
for (i in suffix_all){
    plot_df <- Reduce("rbind", lapply(seq(3), function(j){
        temp <- as.data.frame(cbind(t(integrated_results_list[[i]][c(prop_metrics[j], signature_metrics[j]), all_methods]), all_methods, metrics[j]))
        colnames(temp) <- c("x", "y", "method", "metric")
        rownames(temp) <- NULL
        temp
    }))
    plot_df$method <- factor(plot_df$method, levels = all_methods)
    plot_df$metric <- factor(plot_df$metric, levels = metrics)
    fname <- i
    pdf(paste("SMART_results/simulation/unpaired_scenario/comparison_results", paste(fname, "pdf", sep = "."), sep = "/"), width = 7.5, height = 3.8)
    fig <- ggplot(data = plot_df, aes(x = as.numeric(x), y = as.numeric(y), color = method, shape = method)) +
            geom_point(size = 2.5) +
            labs(x = "Cell type proportion", y = "Cell type signature", title = fname) +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
                  plot.background  = element_blank(),
                  axis.text.x = element_text(size = 12, color = "black"),
                  axis.text.y = element_text(size = 12, color = "black"),
                  axis.title.x = element_text(size = 13),
                  axis.title.y = element_text(size = 13),
                  strip.text.x = element_text(size = 13),
                  strip.background = element_rect(fill = "grey90"),
                  legend.position = "bottom") +
            scale_shape_manual(values = c(16, 17, 15, 18, 1, 2, 0, 5, 4)) +
            scale_color_manual(values = use_color) +
            facet_wrap(~metric, ncol = 3, scales = "free")
    print(fig)
    dev.off()
}


