run_RCTD <- function(st_counts_path, sc_counts_path, sc_labels_path, output_dir = ".", doublet_mode = "full", max_cores = 10){

    library(spacexr)
    library(Matrix)
    library(peakRAM)

    if (!dir.exists(output_dir)){
        dir.create(output_dir, recursive = TRUE)
    }

    sc_counts <- read.csv(sc_counts_path, header = TRUE, row.names = 1)
    sc_counts <- t(sc_counts)

    sc_labels <- read.csv(sc_labels_path, header = TRUE, row.names = 1)
    names_sc_labels <- rownames(sc_labels)
    sc_labels <- sc_labels$celltype
    names(sc_labels) <- names_sc_labels

    use_labels <- sc_labels

    st_counts <- read.csv(st_counts_path, header = TRUE, row.names = 1)
    st_counts <- t(st_counts)

    start_time <- Sys.time()
    peak <- peakRAM(
        reference <- Reference(sc_counts, as.factor(use_labels)),
        puck <- SpatialRNA(coords = NULL, counts = st_counts, use_fake_coords = TRUE),

        myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, UMI_min = 50),
        myRCTD <- run.RCTD(myRCTD, doublet_mode = doublet_mode),

        signature <- t(myRCTD@cell_type_info$info[[1]]),
        proportion <- normalize_weights(myRCTD@results$weights)
    )

    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    peak <- max(peak$Peak_RAM_Used_MiB)

    write.csv(as.matrix(proportion), file = paste(output_dir, "RCTD_prop.csv", sep = "/"))
    write.csv(as.matrix(signature), file = paste(output_dir, "RCTD_signature.csv", sep = "/"))
    write.csv(t(as.matrix(c("elapsed_time" = elapsed_time, "peak" = peak))), file = paste(output_dir, "RCTD_time_memory.csv", sep = "/"))

}





#################### Run RCTD for paired scenario
for (i in list.files("D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data")){
    st_counts_path <- paste(paste("D:/Rshuju/shiyan/SMART_data/paired_scenario/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
    output_dir <- paste("D:/Rshuju/shiyan/SMART_results/simulation/paired_scenario/comparison_methods/RCTD_results", i, sep = "/")
    run_RCTD(st_counts_path = st_counts_path,
             sc_counts_path = "D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/counts_validation.csv",
             sc_labels_path = "D:/Rshuju/shiyan/SMART_data/paired_scenario/sc_data/labels_validation.csv",
             output_dir = output_dir,
             doublet_mode = "full",
             max_cores = 3)
}


#################### Run RCTD for unpaired scenario
for (i in list.files("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data")){
    st_counts_path <- paste(paste("D:/Rshuju/shiyan/SMART_data/unpaired_scenario/st_data", i, sep = "/"), "st_counts.csv", sep = "/")
    output_dir <- paste("D:/Rshuju/shiyan/SMART_results/simulation/unpaired_scenario/comparison_methods/RCTD_results", i, sep = "/")
    run_RCTD(st_counts_path = st_counts_path,
             sc_counts_path = "D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_counts.csv",
             sc_labels_path = "D:/Rshuju/shiyan/SMART_data/unpaired_scenario/sc_data/pbmc1_10x_v2_labels.csv",
             output_dir = output_dir,
             doublet_mode = "full",
             max_cores = 3)
}

