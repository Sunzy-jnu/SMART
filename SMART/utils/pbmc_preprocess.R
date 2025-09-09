counts <- Matrix::readMM("Data/raw_data/simulations/unpaired_scenario/pbmc/counts.umi.txt")
cells <- as.vector(read.table("Data/raw_data/simulations/unpaired_scenario/pbmc/cells.umi.new.txt")$V1)
genes <- as.vector(read.table("Data/raw_data/simulations/unpaired_scenario/pbmc/genes.umi.txt")$V1)
genes <- as.vector(sapply(genes, function(x){strsplit(x, split = "_")[[1]][1]}))
meta <- read.table("Data/raw_data/simulations/unpaired_scenario/pbmc/meta.txt", header = TRUE, row.names = 1, quote = "", sep = "\t")
meta <- meta[-1,]

dim(counts)     # 33694x44433
length(cells)   # 44433
length(genes)   # 33694

rownames(counts) <- genes
colnames(counts) <- cells

keep_cells <- intersect(rownames(meta), cells)

counts <- counts[,keep_cells]
meta <- meta[keep_cells,]

unique_experiment <- unique(meta$Experiment)
unique_method <- unique(meta$Method)
unique_labels <- unique(meta$CellType)

unique_method_abbre <- list("CEL-Seq2" = "Celseq2", "10x Chromium (v2) A" = "10x_v2", "10x Chromium (v2) B" = "10x_v2", "10x Chromium (v3)" = "10x_v3", "Drop-seq" = "Drop", "Seq-Well" = "Seqwell", "inDrops" = "inDrops", "10x Chromium (v2)" = "10x_v2")

unique_labels_abbre <- list("CD4+ T cell" = "CD4_T", "Cytotoxic T cell" = "Cytotoxic_T", "Natural killer cell" = "NK", "CD16+ monocyte" = "CD16_monocyte", "CD14+ monocyte" = "CD14_monocyte", "Megakaryocyte" = "Megakaryocyte", "B cell" = "B", "Dendritic cell" = "Dendritic", "Plasmacytoid dendritic cell" = "Plasmacytoid_dendritic", "Unassigned" = "Unassigned")

all_method <- meta$Method
all_method <- sapply(all_method, function(x){unique_method_abbre[[x]]})
names(all_method) <- NULL

all_labels <- meta$CellType
all_labels <- sapply(all_labels, function(x){unique_labels_abbre[[x]]})
names(all_labels) <- NULL

unique_method <- unique(all_method)

info_list <- list()
unique_labels_list <- list()
tables_labels_list <- list()
save_dir <- "Data/raw_data/simulations/unpaired_scenario/pbmc/pbmc_datasets"
if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
}


for (i in unique_experiment){
    for (j in unique_method){
        idx_ij <- which((meta$Experiment == i) & (all_method == j) & (all_labels != "Unassigned"))
        counts_ij <- counts[,idx_ij]
        labels_ij <- all_labels[idx_ij]
        labels_ij <- as.matrix(labels_ij)
        colnames(labels_ij) <- "celltype"
        rownames(labels_ij) <- colnames(counts_ij)
        info_list[[paste(i, j, sep = "_")]] <- dim(counts_ij)
        unique_labels_list[[paste(i, j, sep = "_")]] <- sort(unique(labels_ij))
        tables_labels_list[[paste(i, j, sep = "_")]] <- table(labels_ij)
        write.csv(t(as.matrix(counts_ij)), file = paste(save_dir, paste(c(i, j, "counts.csv"), collapse = "_"), sep = "/"))
        write.csv(labels_ij, file = paste(save_dir, paste(c(i, j, "labels.csv"), collapse = "_"), sep = "/"))
    }
}


identical(unique_labels_list[["pbmc1_10x_v2"]], unique_labels_list[["pbmc2_inDrops"]])
