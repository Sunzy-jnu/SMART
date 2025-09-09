setwd("/home/wangmg/Documents/deconvolution/Data_and_code_for_reproduction")

library(ggplot2)
library(reshape2)



#################### Run PANDA

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

sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10, do_plot = FALSE, do_parallel = FALSE, n_cores = 1, save_dir = "Results/simulations/time_memory/PANDA_results/sc_results")

st_counts_path <- paste(paste("Data/processed_data/simulations/paired_scenario/st_data", "uniform_ST", sep = "/"), "st_counts.csv", sep = "/")
st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
st_results <- st_train(st_counts, sc_results = sc_results, save_dir = "Results/simulations/time_memory/PANDA_results/st_results")


write.csv(t(as.matrix(c("sc_elapsed_time" = sc_results$run_time, "st_elapsed_time" = st_results$run_time, "sc_peak" = sc_results$peak_ram, "st_peak" = st_results$peak_ram))), file = paste("Results/simulations/time_memory", "PANDA_time_memory.csv", sep = "/"))



#################### Compare with other methods
save_dir <- "Analysis/simulations/time_memory"
if (!dir.exists(save_dir)){
    dir.create(save_dir)
}

all_methods <- c("PANDA", "cell2location", "DestVI", "RCTD", "SpatialDWLS", "SPOTlight", "stereoscope", "STRIDE")

use_color <- scales::hue_pal()(length(all_methods))
names(use_color) = all_methods


results_path <- list(PANDA = "Results/simulations/time_memory/PANDA_time_memory.csv",
                     cell2location = "Results/simulations/time_memory/cell2location_time_memory.csv",
                     DestVI = "Results/simulations/time_memory/DestVI_time_memory.csv",
                     RCTD = "Results/simulations/time_memory/RCTD_time_memory.csv",
                     SpatialDWLS = "Results/simulations/time_memory/SpatialDWLS_time_memory.csv",
                     SPOTlight = "Results/simulations/time_memory/SPOTlight_time_memory.csv",
                     stereoscope = "Results/simulations/time_memory/stereoscope_scvi_time_memory.csv",
                     STRIDE = "Results/simulations/time_memory/STRIDE_time_memory.csv")


sc_elapsed_time <- c()
st_elapsed_time <- c()
sc_peak <- c()
st_peak <- c()

for (i in all_methods){
    time_memory_i <- read.csv(results_path[[i]], header = TRUE, row.names = 1)
    if ("sc_elapsed_time" %in% colnames(time_memory_i)){
        sc_elapsed_time <- c(sc_elapsed_time, time_memory_i[, "sc_elapsed_time"])
        st_elapsed_time <- c(st_elapsed_time, time_memory_i[, "st_elapsed_time"])
        sc_peak <- c(sc_peak, time_memory_i[, "sc_peak"])
        st_peak <- c(st_peak, time_memory_i[, "st_peak"])
    } else{
        sc_elapsed_time <- c(sc_elapsed_time, NA)
        st_elapsed_time <- c(st_elapsed_time, time_memory_i[, "elapsed_time"])
        sc_peak <- c(sc_peak, NA)
        st_peak <- c(st_peak, time_memory_i[, "peak"])
    }
}


sc_elapsed_time <- sc_elapsed_time / 60
st_elapsed_time <- st_elapsed_time / 60
sc_peak <- sc_peak / 1024
st_peak <- st_peak / 1024


# comparison_time
plot_df <- as.data.frame(cbind(all_methods, sc_elapsed_time, st_elapsed_time))
colnames(plot_df) <- c("method", "On scRNA-seq reference", "On spatial transcriptomics data")
plot_df <- melt(plot_df, id.vars = "method", variable.name = "object", value.name = "value")
plot_df$method <- factor(plot_df$method, levels = all_methods)

fname <- "comparison_time"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 6, height = 5)
fig <- ggplot(data = plot_df, aes(x = method, y = as.numeric(value), fill = object)) + 
        geom_bar(position = "dodge", stat = "identity") +
        labs(x = "Method", y = "Time (min)", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
              plot.background  = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
              axis.text.y = element_text(size = 12, color = "black"),
              axis.title.x = element_text(size = 13),
              axis.title.y = element_text(size = 13),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              legend.position = "bottom")
print(fig)
dev.off()



# comparison_memory
plot_df <- as.data.frame(cbind(all_methods, sc_peak, st_peak))
colnames(plot_df) <- c("method", "On scRNA-seq reference", "On spatial transcriptomics data")
plot_df <- melt(plot_df, id.vars = "method", variable.name = "object", value.name = "value")
plot_df$method <- factor(plot_df$method, levels = all_methods)

fname <- "comparison_memory"
pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = 6, height = 5)
fig <- ggplot(data = plot_df, aes(x = method, y = as.numeric(value), fill = object)) + 
        geom_bar(position = "dodge", stat = "identity") +
        labs(x = "Method", y = "Memory (GiB)", title = fname) + 
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
              plot.background  = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
              axis.text.y = element_text(size = 12, color = "black"),
              axis.title.x = element_text(size = 13),
              axis.title.y = element_text(size = 13),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              legend.position = "bottom")
print(fig)
dev.off()


