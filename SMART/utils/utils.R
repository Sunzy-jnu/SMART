TISCH_preprocess <- function(files_path, n_cells_per_type = NULL, output_dir = NULL, do_return = FALSE, prefix = ""){
    library(Seurat)
    library(readr)
    expression <- Read10X_h5(list.files(files_path, pattern = "expression.h5", full.names = TRUE))
    # expression <- round(as.matrix(exp(expression) - 1))
    meta <- read_tsv(list.files(files_path, pattern = "CellMetainfo_table.tsv", full.names = TRUE))
    meta <- as.data.frame(meta)
    rownames(meta) <- meta$Cell
    meta <- meta[colnames(expression),]
    labels_df <- data.frame(celltype = gsub("/", "", meta[,"Celltype (major-lineage)"]), row.names = rownames(meta))

    if (!is.null(n_cells_per_type)){
        idx <- split(seq(dim(labels_df)[1]), labels_df$celltype)
        keep_idx <- Reduce(c, lapply(idx, function(x){
            n_x <- length(x)
            if (n_x < n_cells_per_type){
                n_cells_per_type <- n_x
            }
            sample(x, n_cells_per_type)
        }))
        expression <- expression[,keep_idx]
        labels_df <- data.frame(celltype = labels_df$celltype[keep_idx], row.names = rownames(labels_df)[keep_idx])
    }

    expression <- round(as.matrix(exp(expression) - 1))

    if (!is.null(output_dir)){
        if (!dir.exists(output_dir)){
            dir.create(output_dir, recursive = TRUE)
        }
        write.csv(t(expression), file = paste(output_dir, paste(prefix, "counts.csv", sep = ""), sep = "/"))
        write.csv(labels_df, file = paste(output_dir, paste(prefix, "labels.csv", sep = ""), sep = "/"))
    }

    if (do_return){
        return(list(counts = expression, labels = labels_df))
    }
}



find_hvgs <- function(counts, n_hvgs = 2000){
    library(Seurat)
    library(dplyr)
    dat <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 0)
    dat <- NormalizeData(dat)
    dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = n_hvgs)
    hvgs <- VariableFeatures(dat)
    return(hvgs)
}



find_markers <- function(counts, labels, n_markers = NULL){
    library(Seurat)
    library(dplyr)
    dat <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 0)
    dat <- NormalizeData(dat)
    dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(dat)
    dat <- ScaleData(dat, features = all.genes)
    Idents(dat) <- labels
    markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    unique_labels <- unique(labels)
    markers_list <- lapply(unique_labels, function(x){
        temp <- markers[markers$cluster == x,]
        temp <- temp[order(temp$avg_log2FC, decreasing = TRUE),]
        temp})
    names(markers_list) <- unique_labels
    if (!is.null(n_markers)){
        marker_genes <- lapply(markers_list, function(x){
            temp <- x$gene
            n_markers <- min(n_markers, length(temp))
            temp[1:n_markers]})
    } else{
        marker_genes <- lapply(markers_list, function(x){x$gene})
    }
    return(marker_genes)
}



find_clusters <- function(counts, n_hvgs = 2000, resolution = 0.5){
    library(Seurat)
    library(dplyr)
    dat <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 0)
    dat <- NormalizeData(dat)
    dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = n_hvgs)
    all.genes <- rownames(dat)
    dat <- ScaleData(dat, features = all.genes)
    dat <- RunPCA(dat, features = VariableFeatures(object = dat))
    dat <- FindNeighbors(dat, dims = 1:10)
    dat <- FindClusters(dat, resolution = resolution)
    return(Idents(dat))
}



cosine_similarity <- function(data_1, data_2){
    data_1 <- as.matrix(data_1)
    data_2 <- as.matrix(data_2)
    data <- rbind(data_1, data_2)
    data_scale <- scale(data[, colSums(data) != 0], center = TRUE, scale = TRUE)
    data_1 <- data_scale[1:nrow(data_1),]
    data_2 <- data_scale[(nrow(data_1) + 1):nrow(data_scale),]
    data_1 <- data_1 / sqrt(rowSums(data_1^2))
    data_2 <- data_2 / sqrt(rowSums(data_2^2))
    res <- data_1 %*% t(data_2)
    return(res)
}



cosine_distance <- function(data){
    data <- as.matrix(data)
    temp <- data / sqrt(rowSums(data^2))
    res <- 1 - temp %*% t(temp)
    diag(res) <- 0
    return(res)
}

# 计算余弦相似度
cosine_similarity <- function(x, y) {
  dot_product <- x %*% t(y)
  norm_x <- sqrt(rowSums(x^2))
  norm_y <- sqrt(rowSums(y^2))
  return(dot_product / outer(norm_x, norm_y))
}

max_min_sampling <- function(data, n_points = 5){
    n_data <- nrow(data)
    data_scale <- scale(data[, colSums(data) != 0], center = TRUE, scale = TRUE)
    data_dist <- cosine_distance(data_scale)
    all_points <- which.max(rowSums(data_scale^2))
    if (n_points > 1){
        for (i in 2:n_points){
            all_points <- c(all_points, which.max(apply(as.matrix(data_dist[, all_points]), 1, min)))
        }
    }
    return(all_points)
}



enrichment_celltype <- function(proportion, regions, n_shuffle = 10000){
    unique_regions <- sort(unique(regions))
    true_mean <- t(sapply(unique_regions, function(x){colMeans(proportion[regions == x,])}))
    null_mean <- sapply(seq(n_shuffle), function(i){
        permute_i <- sample(seq(length(regions)), length(regions))
        temp <- t(sapply(unique_regions, function(x){colMeans(proportion[permute_i[regions == x],])}))
        temp
    }, simplify = "array")
    diff_mean <- sapply(seq(n_shuffle), function(i){true_mean - null_mean[,,i]}, simplify = "array")
    std <- apply(diff_mean, 1:2, sd)
    diff_mean <- apply(diff_mean, 1:2, mean)
    diff_mean <- diff_mean / std
    return(list(diff_mean = diff_mean, true_mean = true_mean, null_mean = null_mean))
}
