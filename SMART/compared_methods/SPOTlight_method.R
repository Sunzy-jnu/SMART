run_SPOTlight <- function(st_counts_path, sc_counts_path, sc_labels_path, output_dir = "."){

    library(SPOTlight)
    library(SingleCellExperiment)
    library(scran)
    library(peakRAM)

    if (!dir.exists(output_dir)){
        dir.create(output_dir, recursive = TRUE)
    }

    sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
    sc_counts <- t(sc_counts)

    sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
    use_labels <- sc_labels

    st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
    st_counts <- t(st_counts)


    # sc
    sc_start_time <- Sys.time()

    sc_peak <- peakRAM(
        sce <- SingleCellExperiment(assays = list(counts = sc_counts), colData = use_labels),
        sce <- logNormCounts(sce),
        genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce)),
        dec <- modelGeneVar(sce, subset.row = genes),
        hvg <- getTopHVGs(dec, n = 3000),
        colLabels(sce) <- colData(sce)$celltype,
        mgs <- scoreMarkers(sce, subset.row = genes),
        mgs_fil <- lapply(names(mgs), function(i) {
            x <- mgs[[i]]
            x <- x[x$mean.AUC > 0.8, ]
            x <- x[order(x$mean.AUC, decreasing = TRUE), ]
            x$gene <- rownames(x)
            x$cluster <- i
            data.frame(x)
        }),
        mgs_df <- do.call(rbind, mgs_fil),
        idx <- split(seq(ncol(sce)), sce$celltype),
        n_cells <- 100,
        cs_keep <- lapply(idx, function(i) {
            n <- length(i)
            if (n < n_cells)
                n_cells <- n
            sample(i, n_cells)
        }),
        sce <- sce[, unlist(cs_keep)],

        # res <- SPOTlight(x = sce, y = st_counts, groups = as.character(sce$celltype), mgs = mgs_df, hvg = hvg, weight_id = "mean.AUC", group_id = "cluster", gene_id = "gene")
        mod_ls <- trainNMF(x = sce, y = st_counts, groups = as.character(sce$celltype), mgs = mgs_df, hvg = hvg, weight_id = "mean.AUC", group_id = "cluster", gene_id = "gene")
    )

    sc_end_time <- Sys.time()
    sc_elapsed_time <- as.numeric(difftime(sc_end_time, sc_start_time, units = "secs"))
    sc_peak <- max(sc_peak$Peak_RAM_Used_MiB)


    # st
    st_start_time <- Sys.time()

    st_peak <- peakRAM(
        res <- runDeconvolution(x = st_counts, mod = mod_ls[["mod"]], ref = mod_ls[["topic"]]),
        weights <- res$mat
    )

    st_end_time <- Sys.time()
    st_elapsed_time <- as.numeric(difftime(st_end_time, st_start_time, units = "secs"))
    st_peak <- max(st_peak$Peak_RAM_Used_MiB)

    signature <- t(NMF::basis(mod_ls[["mod"]]) %*% t(mod_ls[["topic"]]))


    write.csv(as.matrix(weights), file = paste(output_dir, "SPOTlight_prop.csv", sep = "/"))
    write.csv(as.matrix(signature), file = paste(output_dir, "SPOTlight_signature.csv", sep = "/"))
    write.csv(t(as.matrix(c("sc_elapsed_time" = sc_elapsed_time, "st_elapsed_time" = st_elapsed_time, "sc_peak" = sc_peak, "st_peak" = st_peak))), file = paste(output_dir, "SPOTlight_time_memory.csv", sep = "/"))

}





#################### Run SPOTlight for paired scenario
for (i in list.files("D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data")){
    st_counts_path <- paste(paste("D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
    output_dir <- paste("D:/Rshuju/shiyan/SMART_results/simulation/paired_scenario/comparison_methods/SPOTlight_results", i, sep = "/")
    run_SPOTlight(st_counts_path = st_counts_path,
                  sc_counts_path = "D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/counts_validation.csv",
                  sc_labels_path = "D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/labels_validation.csv",
                  output_dir = output_dir)
}


#################### Run SPOTlight for unpaired scenario
for (i in list.files("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data")){
    st_counts_path <- paste(paste("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
    output_dir <- paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/SPOTlight_results", i, sep = "/")
    run_SPOTlight(st_counts_path = st_counts_path,
                  sc_counts_path = "D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_counts.csv",
                  sc_labels_path = "D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_labels.csv",
                  output_dir = output_dir)
}


