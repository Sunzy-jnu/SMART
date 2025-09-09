source("D:/Rshuju/shiyan/NPANDA_1.1.0/utils/generate.R")

split_dataset <- function(sc_counts, sc_labels, min_cells_per_type = 30, prefix = "", save_dir = "."){

    if (!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }

    unique_labels <- unique(sc_labels)
    keep_types <- unique_labels[table(sc_labels)[unique_labels] > min_cells_per_type * 2]
    keep_idx <- sc_labels %in% keep_types

    sc_counts <- sc_counts[keep_idx,]
    sc_labels <- sc_labels[keep_idx]

    idx_generation <- c()
    idx_validation <- c()

    for (i in keep_types){
        idx_i <- which(sc_labels == i)
        n_generation <- as.numeric(round(length(idx_i) / 2))
        idx_i_generation <- sample(seq(length(idx_i)), n_generation)
        idx_generation <- c(idx_generation, idx_i[idx_i_generation])
        idx_validation <- c(idx_validation, idx_i[-idx_i_generation])
    }

    idx_generation <- sort(idx_generation)
    idx_validation <- sort(idx_validation)

    counts_generation <- sc_counts[idx_generation,]
    counts_validation <- sc_counts[idx_validation,]

    labels_generation <- sc_labels[idx_generation]
    labels_validation <- sc_labels[idx_validation]

    labels_generation <- as.matrix(labels_generation)
    colnames(labels_generation) <- "celltype"
    labels_validation <- as.matrix(labels_validation)
    colnames(labels_validation) <- "celltype"

    write.csv(counts_generation, file = paste(save_dir, paste(prefix, "counts_generation.csv", sep = ""), sep = "/"))
    write.csv(counts_validation, file = paste(save_dir, paste(prefix, "counts_validation.csv", sep = ""), sep = "/"))
    write.csv(labels_generation, file = paste(save_dir, paste(prefix, "labels_generation.csv", sep = ""), sep = "/"))
    write.csv(labels_validation, file = paste(save_dir, paste(prefix, "labels_validation.csv", sep = ""), sep = "/"))
}

generate_coords <- function(total_spots, interval = 100, type = "uniform") {
  if (type == "uniform") {
    coords <- data.frame(
      x = runif(total_spots, 0, interval),
      y = runif(total_spots, 0, interval)
    )
  } else if (type == "boundary") {
    coords <- data.frame(
      x = c(runif(total_spots/2, 0, interval/5), runif(total_spots/2, 4*interval/5, interval)),
      y = c(runif(total_spots/2, 0, interval/5), runif(total_spots/2, 4*interval/5, interval))
    )
  } else {
    stop("Invalid coordinate type")
  }
  return(coords)
}

# 生成伪空间数据函数
generate_pseudo_st <- function(sc_counts, sc_labels, coords, n_spots = 1000,
                               min_cells = 10, max_cells = 30, temperature = 1e-1,
                               tech_type = "ST", save_dir = ".") {
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  sc_counts <- as.matrix(sc_counts)
  Z <- sc_counts / rowSums(sc_counts)
  n_genes <- ncol(sc_counts)

  unique_labels <- unique(sc_labels)
  n_unique_labels <- length(unique_labels)

  # 初始化数据结构
  st_counts <- matrix(0, nrow = n_spots, ncol = n_genes, dimnames = list(1:n_spots, colnames(sc_counts)))
  st_type_memb <- matrix(0, nrow = n_spots, ncol = n_unique_labels, dimnames = list(1:n_spots, unique_labels))
  st_type_prop <- matrix(0, nrow = n_spots, ncol = n_unique_labels, dimnames = list(1:n_spots, unique_labels))
  st_type_exp_list <- lapply(unique_labels, function(x) {
    matrix(0, nrow = n_spots, ncol = n_genes, dimnames = list(1:n_spots, colnames(sc_counts)))
  })
  names(st_type_exp_list) <- unique_labels

  # 预计算每个细胞类型的信息
  idx_type_list <- list()
  prob_type_list <- list()
  knn_type_list <- list()
  for (cell_type in unique_labels) {
    idx <- which(sc_labels == cell_type)
    mean_expr <- colSums(Z[idx, ])
    dist <- sapply(idx, function(i) 1 - sum(Z[i, ] * mean_expr) / sqrt(sum(Z[i, ]^2) * sum(mean_expr^2)))
    prob <- exp(dist / temperature)
    prob <- prob / sum(prob)
    knn_result <- knn.covertree::find_knn(Z[idx, ], k = max_cells, distance = "cosine")
    knn_indices <- apply(knn_result$index, 2, function(x) idx[x])
    rownames(knn_indices) <- rownames(sc_counts)[idx]

    idx_type_list[[cell_type]] <- idx
    prob_type_list[[cell_type]] <- prob
    knn_type_list[[cell_type]] <- knn_indices
  }

  # 生成每个空间位点
  for (s in 1:n_spots) {
    s_counts <- rep(0, n_genes)
    s_type_memb <- rep(0, n_unique_labels)  # 确保长度与 st_type_prop 的列数一致

    n_cells <- round(runif(1, min = min_cells, max = max_cells))
    n_types <- round(runif(1, min = 1, max = n_unique_labels))
    selected_types <- sample(unique_labels, n_types)

    memb <- rep(0, n_types)
    while (sum(memb) < min_cells) {
      prop <- gtools::rdirichlet(1, rep(1, n_types))
      memb <- round(n_cells * prop)
    }
    s_type_memb[match(selected_types, unique_labels)] <- memb  # 使用 match 确保正确索引
    s_type_prop_row <- s_type_memb / sum(s_type_memb)

    st_type_prop[s, ] <- s_type_prop_row

    for (cell_type in unique_labels) {
      if (s_type_memb[match(cell_type, unique_labels)] == 0) next

      idx <- idx_type_list[[cell_type]]
      prob <- prob_type_list[[cell_type]]
      knn <- knn_type_list[[cell_type]]

      if (s_type_memb[match(cell_type, unique_labels)] == 1) {
        selected_cell <- sample(idx, size = 1, prob = prob)
        st_type_exp_list[[cell_type]][s, ] <- sc_counts[selected_cell, ]
      } else {
        selected_cell <- sample(idx, size = 1, prob = prob)
        neighbors <- knn[rownames(sc_counts)[selected_cell], 1:(s_type_memb[match(cell_type, unique_labels)] - 1)]
        st_type_exp_list[[cell_type]][s, ] <- colSums(sc_counts[c(selected_cell, neighbors), ])
      }

      s_counts <- s_counts + st_type_exp_list[[cell_type]][s, ]
    }

    st_counts[s, ] <- s_counts
    st_type_memb[s, ] <- s_type_memb
  }

  # 保存数据
  write.csv(coords[1:n_spots, ], file = file.path(save_dir, "st_coords.csv"))
  write.csv(st_counts, file = file.path(save_dir, "st_counts.csv"))
  write.csv(st_type_memb, file = file.path(save_dir, "st_type_memb.csv"))
  write.csv(st_type_prop, file = file.path(save_dir, "st_type_prop.csv"))
  saveRDS(st_type_exp_list, file = file.path(save_dir, "st_type_exp_list.rds"))

  # 生成可视化
  if (tech_type %in% c("ST", "Visium")) {
    save_subdir <- file.path(save_dir, "figs")
    if (!dir.exists(save_subdir)) {
      dir.create(save_subdir)
    }

    for (cell_type in unique_labels) {
      write.csv(as.matrix(st_type_exp_list[[cell_type]]), file = paste(save_dir, paste(c("st", cell_type, "exp.csv"), collapse = "_"), sep = "/"))
      expr <- st_type_exp_list[[cell_type]][rowSums(st_type_exp_list[[cell_type]]) > 0, ]

      plot_simulation(sc_counts[sc_labels == cell_type, ], expr,
                      fname = paste0("simulation_", cell_type),
                      save_dir = save_subdir)
    }
  }
}
