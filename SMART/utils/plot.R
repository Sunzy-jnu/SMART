plot_simulation <- function(counts, simulations, counts_labels = NULL, simulations_labels = NULL, fname = "simulation", save_dir = ".", width = 7, height = 7){

    library(ggplot2)

    if (!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }

    if (is.null(counts_labels)){
        counts_labels <- rep("truth", dim(counts)[1])
    }

    if (is.null(simulations_labels)){
        simulations_labels <- rep("simulation", dim(simulations)[1])
    }

    gene_intersect <- intersect(colnames(counts), colnames(simulations))
    counts <- counts / rowSums(counts)
    counts <- counts[, gene_intersect]
    simulations <- simulations / rowSums(simulations)
    simulations <- simulations[, gene_intersect]

    all_data <- rbind(counts, simulations)
    all_labels <- c(counts_labels, simulations_labels)
    all_domains <- c(rep("truth", length(counts_labels)), rep("simulation", length(simulations_labels)))

    # 检查合并后的数据矩阵是否存在非数值元素
    if (any(is.na(all_data)) || any(is.nan(all_data)) || any(is.infinite(all_data))) {
      # 记录非数值元素位置
      na_rows <- which(rowSums(is.na(all_data)) > 0)
      nan_rows <- which(rowSums(is.nan(all_data)) > 0)
      inf_rows <- which(rowSums(is.infinite(all_data)) > 0)

      if (length(na_rows) > 0) {
        warning(paste("数据中存在 NA 值，位于行:", paste(na_rows, collapse = ", ")))
      }
      if (length(nan_rows) > 0) {
        warning(paste("数据中存在 NaN 值，位于行:", paste(nan_rows, collapse = ", ")))
      }
      if (length(inf_rows) > 0) {
        warning(paste("数据中存在 Inf 值，位于行:", paste(inf_rows, collapse = ", ")))
      }

      # 删除包含非数值元素的行
      all_data <- na.omit(all_data)
      all_labels <- all_labels[complete.cases(all_data)]
      all_domains <- all_domains[complete.cases(all_data)]
    }

    plot_df <- as.data.frame(prcomp(all_data[, colSums(all_data) != 0], center = TRUE, scale. = TRUE)$x[,1:30])
    plot_df <- as.data.frame(umap::umap(as.matrix(plot_df))$layout)

    plot_df <- cbind(plot_df, all_labels, all_domains)
    colnames(plot_df) <- c("dim1", "dim2", "labels", "domains")

    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot(data = plot_df, aes(x = dim1, y = dim2, color = as.factor(labels), shape = as.factor(domains))) + geom_point() + labs(x = "UMAP_1", y = "UMAP_2")
    print(fig)
    dev.off()
}



plot_archetypes <- function(Z, phi, f_name = "archetypes", save_dir = ".", width = 7, height = 7){

    library(ggplot2)
    library(ggrepel)

    if (!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }

    gene_intersect <- intersect(colnames(Z), colnames(phi))
    Z <- Z[, gene_intersect]
    phi <- phi[, gene_intersect]

    idx_gene <- which(colSums(Z) != 0)
    pca_res <- prcomp(Z[, idx_gene], center = TRUE, scale. = TRUE)
    ndims <- min(30, dim(Z)[1])
    Z_pca <- as.data.frame(pca_res$x[, 1:ndims])
    umap_settings <- umap::umap.defaults
    umap_settings$n_neighbors <- min(15, dim(Z)[1])
    umap_res <- umap::umap(as.matrix(Z_pca), config = umap_settings)
    Z_umap <- as.data.frame(umap_res$layout)
    colnames(Z_umap) <- c("dim1", "dim2")
    phi_pca <- as.data.frame(predict(pca_res, phi[, idx_gene])[, 1:ndims])
    phi_umap <- as.data.frame(predict(umap_res, phi_pca))
    phi_umap <- cbind(phi_umap, 1:nrow(phi))
    colnames(phi_umap) <- c("dim1", "dim2", "name")

    pdf(paste(save_dir, paste(f_name, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot() +
            geom_point(aes(x = dim1, y = dim2), color = "#70B2DE", shape = 16, size = 1, data = Z_umap) +
            geom_point(aes(x = dim1, y = dim2), color = "#F0525F", shape = 17, size = 2, data = phi_umap) +
            geom_text_repel(aes(x = dim1, y = dim2, label = name), data = phi_umap, fontface = "bold", size = 2) +
            labs(title = f_name, x = "UMAP_1", y = "UMAP_2") +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  legend.position = "none")
    print(fig)
    dev.off()
}



plot_evaluation <- function(..., cname, mname, use_color = NULL, average = TRUE, fname = "evaluation", save_dir = ".", n_cols = 4, width = 7, height = 7){
    library(ggplot2)
    library(ggpubr)
    plot_df <- as.data.frame(cbind(...))
    colnames(plot_df) <- cname
    if(average){
        plot_df <- as.data.frame(cbind(as.vector(rep(colnames(plot_df), each = dim(plot_df)[1])), as.vector(sapply(Reduce(c, plot_df), mean)), as.vector(rep(rownames(plot_df), times = dim(plot_df)[2]))))
        colnames(plot_df) <- c("metric", "value", "method")
        rownames(plot_df) <- NULL
        plot_df$metric <- factor(plot_df$metric, levels = cname)
        plot_df$method <- factor(plot_df$method, levels = mname)
        pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
        fig <- ggplot(data = plot_df, aes(x = method, y = as.numeric(value), fill = method)) +
                geom_bar(stat = "identity") +
                geom_text(aes(label = round(as.numeric(value), 3)), position = position_dodge(0.9), vjust = -0.3, size = 2) +
                labs(x = "", y = "", title = fname) +
                theme(panel.background = element_blank(),
                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
                      plot.background  = element_blank(),
                      axis.text.x = element_blank(),
                    #   axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                      axis.text.y = element_text(size = 12),
                      strip.text.x = element_text(size = 15),
                      legend.position = "bottom") +
                facet_wrap(~metric, ncol = n_cols, scales = "free")
        if (!is.null(use_color)){
            fig <- fig + scale_fill_manual(values = use_color)
        }
        print(fig)
        dev.off()
    }else{
        plot_df <- as.data.frame(cbind(as.vector(rep(colnames(plot_df), each = dim(plot_df)[1] * length(plot_df[1,1][[1]]))), as.vector(unlist(Reduce(c, Reduce(c, plot_df)))), as.vector(rep(rep(rownames(plot_df), each = length(plot_df[1,1][[1]])), times = dim(plot_df)[2]))))
        colnames(plot_df) <- c("metric", "value", "method")
        rownames(plot_df) <- NULL
        plot_df$metric <- factor(plot_df$metric, levels = cname)
        plot_df$method <- factor(plot_df$method, levels = mname)
        pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
        fig <- ggplot(data = plot_df, aes(x = method, y = as.numeric(value), fill = method)) +
                geom_boxplot(size = 0.2, outlier.size = 0.01) +
                stat_compare_means(aes(label = after_stat(p.signif)), method = "wilcox.test", ref.group = "SMART", vjust = -0.3, size = 2) +
                labs(x = "", y = "", title = fname) +
                theme(panel.background = element_blank(),
                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
                      plot.background  = element_blank(),
                      axis.text.x = element_blank(),
                    #   axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                      axis.ticks.x = element_blank(),
                      axis.text.y = element_text(size = 12),
                      strip.text.x = element_text(size = 12),
                      strip.background = element_rect(fill = "grey90"),
                      legend.position = "bottom") +
                facet_wrap(~metric, ncol = n_cols, scales = "free")
        if (!is.null(use_color)){
            fig <- fig + scale_fill_manual(values = use_color)
        }
        print(fig)
        dev.off()
    }
}



plot_prop_signature_acc <- function(data, fname = "prop_signature_acc", save_dir = ".", size = 0.5, color = "deepskyblue3", y_lab = "Accuracy of cell type signature estimation", n_cols = 3, x_breaks = seq(0, 1, by = 0.1), width = 7, height = 7){
    library(ggplot2)
    plot_df <- Reduce("rbind", lapply(names(data), function(x){
        temp <- as.data.frame(cbind(data[[x]], x))
        colnames(temp) <- c("prop", "sig_acc", "label")
        temp
    }))
    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot(data = plot_df, aes(x = as.numeric(prop), y = as.numeric(sig_acc))) +
            geom_point(size = size, color = color) +
            # geom_smooth(aes(group = 1), method = "loess", se = FALSE) +
            labs(x = "Ground truth cell type proportion", y = y_lab, title = fname) +
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
            scale_x_continuous(breaks = x_breaks) +
            facet_wrap(~label, ncol = n_cols)
    print(fig)
    dev.off()
}



plot_rank <- function(data, mname, use_color = NULL, fname = "rank", save_dir = ".", width = 7, height = 7){
    library(ggplot2)
    plot_df <- as.data.frame(cbind(names(data), data))
    colnames(plot_df) <- c("method", "value")
    plot_df$method <- factor(plot_df$method, levels = mname)
    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot() +
            geom_bar(data = plot_df, aes(x = as.numeric(value), y = reorder(method, -as.numeric(value)), fill = method), stat = "identity", width = 0.7) +
            labs(x = "Rank", y = "Method", title = fname) +
            theme(panel.background = element_blank(),
                  panel.border = element_blank(),
                  plot.background  = element_blank(),
                  axis.line.x = element_line(color = "black"),
                  axis.line.y = element_line(color = "black"),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_text(size = 10)) +
            scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))
    if (!is.null(use_color)){
        fig <- fig + scale_fill_manual(values = use_color)
    }
    print(fig)
    dev.off()
}



plot_image <- function(img, location = NULL, fname = "image", save_dir = ".", reverse_x = FALSE, reverse_y = TRUE, point_size = 1, width = 7, height = 7){
    library(ggplot2)
    red <- img[,,1]
    green <- img[,,2]
    blue <- img[,,3]
    plot_df <- data.frame(x = rep(1:ncol(img), each = nrow(img)), y = rep(1:nrow(img), ncol(img)), red = as.vector(red), green = as.vector(green), blue = as.vector(blue))
    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot() +
            geom_raster(data = plot_df, aes(x = x, y = y, fill = rgb(red, green, blue))) +
            coord_fixed() +
            scale_fill_identity() +
            theme_void() +
            labs(title = fname)
    if (reverse_y){
        fig <- fig + scale_y_continuous(trans = scales::reverse_trans())
    }
    if (reverse_x){
        fig <- fig + scale_x_continuous(trans = scales::reverse_trans())
    }
    if (!is.null(location)){
        fig <- fig +
                geom_point(data = location, aes(x = x, y = y), color = "#5E4FA2", size = point_size)
    }
    print(fig)
    dev.off()
}



plot_scatter_label <- function(data, location, img = NULL, use_color = NULL, fname = "scatter_label", save_dir = ".", reverse_x = FALSE, reverse_y = TRUE, n_cols = 1, point_size = 1, width = 7, height = 7){
    library(ggplot2)
    library(reshape2)
    library(ggnewscale)
    plot_df <- as.data.frame(cbind(location[rownames(data), c("x", "y")], data))
    plot_df <- melt(plot_df, id.vars = c("x", "y"), variable.name = "feature", value.name = "value")

    if (is.null(use_color)){
        if (length(unique(c(data))) <= 9){
            use_color <- RColorBrewer::brewer.pal(length(unique(c(data))), "Set1")
        } else{
            use_color <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(c(data))))
        }
    }

    if (!is.null(img)){
        red <- img[,,1]
        green <- img[,,2]
        blue <- img[,,3]
        img_df <- data.frame(x = rep(1:ncol(img), each = nrow(img)), y = rep(1:nrow(img), ncol(img)), red = as.vector(red), green = as.vector(green), blue = as.vector(blue))
    }

    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot()
    if (!is.null(img)){
        fig <- fig +
                geom_raster(data = img_df, aes(x = x, y = y, fill = rgb(red, green, blue))) +
                scale_fill_identity() +
                new_scale_fill()
    }
    fig <- fig +
            geom_point(data = location, aes(x = x, y = y), fill = "grey80", size = point_size, colour = "black", shape = 21, stroke = 0.1) +
            geom_point(data = plot_df, aes(x = x, y = y, fill = as.factor(value)), size = point_size, colour = "black", shape = 21, stroke = 0.1) +
            coord_fixed() +
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
            scale_fill_manual(values = use_color) +
            facet_wrap(~feature, ncol = n_cols) +
            labs(title = fname)
    if (reverse_y){
        fig <- fig + scale_y_continuous(trans = scales::reverse_trans())
    }
    if (reverse_x){
        fig <- fig + scale_x_continuous(trans = scales::reverse_trans())
    }
    print(fig)
    dev.off()
}



plot_scatterpie <- function(data, location, img = NULL, use_color = NULL, fname = "scatterpie", save_dir = ".", reverse_x = FALSE, reverse_y = TRUE, pie_r = 2, point_size = 2, width = 7, height = 7){
    library(ggplot2)
    library(scatterpie)
    library(RColorBrewer)
    library(ggpubr)
    library(ggnewscale)
    plot_df <- as.data.frame(apply(cbind(location[rownames(data),], data), 2, as.numeric))
    rownames(plot_df) <- rownames(data)
    unique_labels <- sort(colnames(data))

    if (is.null(use_color)){
        if (length(unique_labels) <= 9){
            use_color <- RColorBrewer::brewer.pal(length(unique_labels), "Set1")
        } else{
            use_color <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_labels))
        }
    }

    if (!is.null(img)){
        red <- img[,,1]
        green <- img[,,2]
        blue <- img[,,3]
        img_df <- data.frame(x = rep(1:ncol(img), each = nrow(img)), y = rep(1:nrow(img), ncol(img)), red = as.vector(red), green = as.vector(green), blue = as.vector(blue))
    }

    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot()
    if (!is.null(img)){
        fig <- fig +
                geom_raster(data = img_df, aes(x = x, y = y, fill = rgb(red, green, blue)), alpha = 0.5) +
                scale_fill_identity() +
                new_scale_fill()
    }
    fig <- fig +
            geom_point(data = location, aes(x = x, y = y), fill = "grey90", size = point_size, colour = "grey90", shape = 21, stroke = 0.01) +
            geom_scatterpie(data = plot_df, aes(x = x, y = y, r = pie_r), cols = unique_labels, color = NA) +
            coord_fixed() +
            scale_fill_manual(values = use_color) +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                #   panel.border = element_blank(),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text(size = 15, face = "bold"),
                  strip.background = element_blank(),
                  legend.title = element_blank(),
                  legend.position = "right") +
            labs(title = fname)
    if (reverse_y){
        fig <- fig + scale_y_continuous(trans = scales::reverse_trans())
    }
    if (reverse_x){
        fig <- fig + scale_x_continuous(trans = scales::reverse_trans())
    }
    print(fig)
    dev.off()
}



plot_scatterpie_2 <- function(data_list, location_list, img_list = NULL, mname = NULL, use_color = NULL, fname = "scatterpie_list", save_dir = ".", reverse_x = FALSE, reverse_y = TRUE, n_cols = 3, pie_r = 2, width = 7, height = 7){
    library(ggplot2)
    library(scatterpie)
    library(RColorBrewer)
    library(ggpubr)
    library(ggnewscale)
    if (is.null(mname)){
        mname <- names(data_list)
    }
    if (is.null(mname)){
        mname <- seq(length(data_list))
    }
    unique_labels <- sort(Reduce("intersect", lapply(data_list, colnames)))
    plot_df <- Reduce("rbind", lapply(mname, function(x){
        temp <- as.data.frame(apply(cbind(location_list[[x]][rownames(data_list[[x]]),], data_list[[x]][, unique_labels]), 2, as.numeric))
        temp <- cbind(temp, x)
        colnames(temp) <- c(colnames(location_list[[x]]), unique_labels, "section")
        temp
    }))
    plot_df$section <- factor(plot_df$section, levels = mname)

    if (is.null(use_color)){
        if (length(unique_labels) <= 9){
            use_color <- RColorBrewer::brewer.pal(length(unique_labels), "Set1")
        } else{
            use_color <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique_labels))
        }
    }

    if (!is.null(img_list)){
        img_df <- Reduce("rbind", lapply(mname, function(x){
            img <- img_list[[x]]
            red <- img[,,1]
            green <- img[,,2]
            blue <- img[,,3]
            data.frame(x = rep(1:ncol(img), each = nrow(img)), y = rep(1:nrow(img), ncol(img)), red = as.vector(red), green = as.vector(green), blue = as.vector(blue), section = x)
        }))
    }

    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot()
    if (!is.null(img_list)){
        fig <- fig +
                geom_raster(data = img_df, aes(x = x, y = y, fill = rgb(red, green, blue))) +
                scale_fill_identity() +
                new_scale_fill()
    }
    fig <- fig +
            geom_scatterpie(data = plot_df, aes(x = x, y = y, r = pie_r), cols = unique_labels, color = NA) +
            coord_fixed() +
            scale_fill_manual(values = use_color) +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text(size = 15, face = "bold"),
                  strip.background = element_blank(),
                  legend.title = element_blank(),
                  legend.position = "bottom") +
            labs(title = fname) +
            facet_wrap(~section, ncol = n_cols)
    if (reverse_y){
        fig <- fig + scale_y_continuous(trans = scales::reverse_trans())
    }
    if (reverse_x){
        fig <- fig + scale_x_continuous(trans = scales::reverse_trans())
    }
    print(fig)
    dev.off()
}



plot_scatter_heatmap <- function(data, location, img = NULL, feature = NULL, scale = FALSE, fname = "scatter_heatmap", save_dir = ".", reverse_x = FALSE, reverse_y = TRUE, n_cols = 4, point_size = 1.5, width = 7, height = 7, legend_position = "right", panel_border = TRUE){
    library(ggplot2)
    library(reshape2)
    library(ggnewscale)
    if (is.null(feature)){
        feature <- colnames(data)
    }
    if (scale){
        data <- apply(data[, feature, drop = FALSE], 2, function(y){(y - min(y)) / (max(y) - min(y))})
    }
    plot_df <- as.data.frame(apply(cbind(location[rownames(data), c("x", "y")], data[, feature, drop = FALSE]), 2, as.numeric))
    plot_df <- melt(plot_df, id.vars = c("x", "y"), variable.name = "feature", value.name = "value")
    plot_df$feature <- factor(plot_df$feature, levels = feature)

    if (!is.null(img)){
        red <- img[,,1]
        green <- img[,,2]
        blue <- img[,,3]
        img_df <- data.frame(x = rep(1:ncol(img), each = nrow(img)), y = rep(1:nrow(img), ncol(img)), red = as.vector(red), green = as.vector(green), blue = as.vector(blue))
    }

    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot()
    if (!is.null(img)){
        fig <- fig +
                geom_raster(data = img_df, aes(x = x, y = y, fill = rgb(red, green, blue)), alpha = 1) +
                scale_fill_identity() +
                new_scale_fill()
    }
    fig <- fig +
            geom_point(data = location, aes(x = x, y = y), fill = "grey90", size = point_size, colour = "black", shape = 21, stroke = 0.05) +
            geom_point(data = plot_df, aes(x = x, y = y, fill = value), size = point_size, colour = "black", shape = 21, stroke = 0.05) +
            coord_fixed() +
            scale_fill_gradientn(colors = c("#5E4FA2", "#3E96B6", "#99D5A4", "#FDFEBD", "#FDDB88", "#F67948", "#9E0142")) +
            # scale_fill_gradient2(low = "#989797", high = "red") +
            facet_wrap(~feature, ncol = n_cols) +
            labs(title = fname)
    if (panel_border){
        fig <- fig +
                theme(panel.background = element_blank(),
                    #   panel.border = element_blank(),
                    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                    plot.background  = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    strip.text = element_text(size = 15, face = "bold"),
                    strip.background = element_blank(),
                    legend.title = element_blank(),
                    legend.position = legend_position)
    } else{
        fig <- fig +
                theme(panel.background = element_blank(),
                      panel.border = element_blank(),
                    # panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                    plot.background  = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    strip.text = element_text(size = 15, face = "bold"),
                    strip.background = element_blank(),
                    legend.title = element_blank(),
                    legend.position = legend_position)
    }
    if (reverse_y){
        fig <- fig + scale_y_continuous(trans = scales::reverse_trans())
    }
    if (reverse_x){
        fig <- fig + scale_x_continuous(trans = scales::reverse_trans())
    }
    print(fig)
    dev.off()
}



plot_scatter_heatmap_2 <- function(data_list, location_list, img_list = NULL, mname = NULL, feature = NULL, scale = FALSE, fname = "scatter_heatmap_list", save_dir = ".", reverse_x = FALSE, reverse_y = TRUE, feat_by_sect = FALSE, point_size = 1.5, x_lim = NULL, y_lim = NULL, width = 7, height = 7, legend_position = "bottom"){
    library(ggplot2)
    library(reshape2)
    library(ggnewscale)
    if (is.null(feature)){
        feature <- sort(Reduce("intersect", lapply(data_list, colnames)))
    }
    if (scale){
        data_list <- lapply(data_list, function(x){apply(x[, feature, drop = FALSE], 2, function(y){(y - min(y)) / (max(y) - min(y))})})
    }
    if (is.null(mname)){
        mname <- names(data_list)
    }
    if (is.null(mname)){
        mname <- seq(length(data_list))
    }
    plot_df <- Reduce("rbind", lapply(mname, function(x){
        temp <- cbind(location_list[[x]][rownames(data_list[[x]]), c("x", "y"), drop = FALSE], data_list[[x]][, feature, drop = FALSE])
        temp <- cbind(temp, x)
        colnames(temp) <- c(c("x", "y"), feature, "section")
        temp
    }))
    for (i in c(c("x", "y"), feature)){
        plot_df[, i] <- as.numeric(plot_df[, i])
    }
    plot_df <- melt(plot_df, id.vars = c("x", "y", "section"), variable.name = "feature", value.name = "value")
    plot_df$feature <- factor(plot_df$feature, levels = feature)
    plot_df$section <- factor(plot_df$section, levels = mname)

    location <- Reduce("rbind", lapply(mname, function(x){
        temp <- as.data.frame(apply(location_list[[x]][, c("x", "y")], 2, as.numeric))
        temp <- cbind(temp, x)
        colnames(temp) <- c(c("x", "y"), "section")
        temp
    }))
    location$section <- factor(location$section, levels = mname)

    if (!is.null(img_list)){
        img_df <- Reduce("rbind", lapply(mname, function(x){
            img <- img_list[[x]]
            red <- img[,,1]
            green <- img[,,2]
            blue <- img[,,3]
            data.frame(x = rep(1:ncol(img), each = nrow(img)), y = rep(1:nrow(img), ncol(img)), red = as.vector(red), green = as.vector(green), blue = as.vector(blue), section = x)
        }))
    }

    pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot()
    if (!is.null(img_list)){
        fig <- fig +
                geom_raster(data = img_df, aes(x = x, y = y, fill = rgb(red, green, blue))) +
                scale_fill_identity() +
                new_scale_fill()
    }
    fig <- fig +
            geom_point(data = location, aes(x = x, y = y), fill = "grey90", size = point_size, colour = "black", shape = 21, stroke = 0.05) +
            geom_point(data = plot_df, aes(x = x, y = y, fill = value), size = point_size, colour = "black", shape = 21, stroke = 0.05) +
            coord_fixed() +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text(size = 15, face = "bold"),
                  strip.background = element_blank(),
                  legend.title = element_blank(),
                  legend.position = legend_position) +
            scale_fill_gradientn(colors = c("#5E4FA2", "#3E96B6", "#99D5A4", "#FDFEBD", "#FDDB88", "#F67948", "#9E0142")) +
            labs(title = fname)
    if (reverse_y & (is.null(y_lim))){
        fig <- fig + scale_y_continuous(trans = scales::reverse_trans())
    } else if (reverse_y & (!is.null(y_lim))) {
        fig <- fig + scale_y_continuous(trans = scales::reverse_trans(), limits = c(y_lim[2], y_lim[1]))
    } else if ((!reverse_y) & (!is.null(y_lim))) {
        fig <- fig + scale_y_continuous(limits = y_lim)
    }
    if (reverse_x & (is.null(x_lim))){
        fig <- fig + scale_x_continuous(trans = scales::reverse_trans())
    } else if (reverse_x & (!is.null(x_lim))){
        fig <- fig + scale_x_continuous(trans = scales::reverse_trans(), limits = c(x_lim[2], x_lim[1]))
    } else if ((!reverse_x) & (!is.null(x_lim))){
        fig <- fig + scale_x_continuous(limits = x_lim)
    }
    if (feat_by_sect){
        fig <- fig + facet_grid(feature~section)
    } else{
        fig <- fig + facet_grid(section~feature)
    }
    print(fig)
    dev.off()
}



plot_mu_reference <- function(mu_list, proportion, sc_counts, sc_labels, archetypes_list = NULL, min_prop = 0.2, f_name = "mu_reference", save_dir = ".", return_res = FALSE, n_cols = 4, width = 7, height = 7){

    library(ggplot2)
    library(ggrepel)

    if (!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }

    Z <- sc_counts / rowSums(sc_counts)

    gene_intersect <- intersect(colnames(Z), colnames(mu_list[[1]]))
    Z <- Z[, gene_intersect]
    mu_list <- lapply(mu_list, function(x){x[, gene_intersect]})

    unique_labels <- names(mu_list)
    pca_res_list <- list()
    umap_res_list <- list()
    plot_df <- Reduce(rbind, lapply(unique_labels, function(x){
        idx_sc <- which(sc_labels == x)
        idx_gene <- which(colSums(Z[idx_sc,]) != 0)
        pca_res <- prcomp(Z[idx_sc, idx_gene], center = TRUE, scale. = TRUE)
        pca_res_list[[x]] <- pca_res
        ndims <- min(30, length(idx_sc))
        Z_pca <- as.data.frame(pca_res$x[, 1:ndims])
        umap_settings <- umap::umap.defaults
        umap_settings$n_neighbors <- min(15, length(idx_sc))
        umap_res <- umap::umap(as.matrix(Z_pca), config = umap_settings)
        umap_res_list[[x]] <- umap_res
        Z_umap <- as.data.frame(umap_res$layout)
        idx_mu <- which(proportion[, x] > min_prop)
        if (length(idx_mu) > 0){
            mu <- mu_list[[x]][idx_mu, idx_gene, drop = FALSE]
            mu_pca <- as.data.frame(predict(pca_res, mu)[, 1:ndims, drop = FALSE])
            mu_umap <- as.data.frame(predict(umap_res, mu_pca))
            temp <- rbind(Z_umap, mu_umap)
            domain <- c(rep("sc", length(idx_sc)), rep("st", length(idx_mu)))
            prop <- c(rep(0, length(idx_sc)), proportion[idx_mu, x])
            name <- c(rownames(Z[idx_sc,]), rownames(mu))
            temp <- as.data.frame(cbind(temp, domain, prop, x, name))
            colnames(temp) <- c("dim1", "dim2", "domain", "prop", "label", "name")
        } else{
            temp <- Z_umap
            domain <- rep("sc", length(idx_sc))
            prop <- rep(0, length(idx_sc))
            name <- rownames(Z[idx_sc,])
            temp <- as.data.frame(cbind(temp, domain, prop, x, name))
            colnames(temp) <- c("dim1", "dim2", "domain", "prop", "label", "name")
        }
        if (!is.null(archetypes_list)){
            arche <- archetypes_list[[x]][, gene_intersect[idx_gene]]
            arche_pca <- as.data.frame(predict(pca_res, arche)[, 1:ndims])
            arche_umap <- as.data.frame(predict(umap_res, arche_pca))
            arche_temp <- as.data.frame(cbind(arche_umap, "archetype", 0, x, substring(rownames(arche), nchar(x) + 2)))
            colnames(arche_temp) <- c("dim1", "dim2", "domain", "prop", "label", "name")
            temp <- rbind(temp, arche_temp)
        }
        temp
    }))

    pdf(paste(save_dir, paste(f_name, "reference_only.pdf", sep = "_"), sep = "/"), width = width, height = height)
    fig <- ggplot() +
            geom_point(aes(x = dim1, y = dim2, color = as.factor(domain), shape = as.factor(domain), size = as.factor(domain)), data = plot_df[(plot_df$domain == "sc") | (plot_df$domain == "archetype"),]) +
            geom_text_repel(aes(x = dim1, y = dim2, label = name), data = plot_df[plot_df$domain == "archetype",], fontface = "bold") +
            labs(title = paste(f_name, "reference_only", sep = "_"), x = "UMAP_1", y = "UMAP_2") +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text(size = 10, face = "bold"),
                  legend.position = "right") +
            scale_color_manual(values = c(sc = "#989797", archetype = "#00A3D5")) +
            scale_shape_manual(values = c(sc = 16, archetype = 17)) +
            scale_size_manual(values = c(sc = 2, archetype = 2.5)) +
            facet_wrap(~label, scales = "free", ncol = n_cols)
    print(fig)
    dev.off()

    pdf(paste(save_dir, paste(f_name, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot() +
            geom_point(aes(x = dim1, y = dim2), color = "#989797", shape = 16, size = 2, data = plot_df[plot_df$domain == "sc",]) +
            geom_point(aes(x = dim1, y = dim2, color = prop), shape = 16, size = 2, data = plot_df[plot_df$domain == "st",]) +
            labs(title = f_name, x = "UMAP_1", y = "UMAP_2") +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text(size = 10, face = "bold"),
                  legend.position = "right") +
            scale_color_gradient(low = "#fae1df", high =  "#f7061a") +
            facet_wrap(~label, scales = "free", ncol = n_cols)
    print(fig)
    dev.off()

    if (return_res){
        return(list(plot_df = plot_df, pca_res_list = pca_res_list, umap_res_list = umap_res_list))
    }
}



plot_mu_reference_2 <- function(mu_list_all, proportion_list, sc_counts, sc_labels, archetypes_list = NULL, min_prop = 0.2, do_hvgs = TRUE, use_color = NULL, f_name = "mu_reference_list", save_dir = ".", return_res = FALSE, n_cols = 4, width = 7, height = 7){

    library(ggplot2)
    library(ggrepel)
    library(ggnewscale)

    if (!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }

    if (is.null(use_color)){
        use_color <- scales::hue_pal()(length(mu_list_all))
    }

    Z <- sc_counts / rowSums(sc_counts)

    if (do_hvgs){
        hvgs <- find_hvgs(counts = t(sc_counts), n_hvgs = 2000)
        Z <- Z[, hvgs]
    }

    gene_intersect <- intersect(colnames(Z), Reduce("intersect", lapply(mu_list_all, function(x){colnames(x[[1]])})))
    Z <- Z[, gene_intersect]
    mu_list_all <- lapply(mu_list_all, function(x){lapply(x, function(y){y[, gene_intersect]})})

    unique_labels <- sort(unique(sc_labels))
    gene_list <- list()
    pca_res_list <- list()
    umap_res_list <- list()
    plot_df <- Reduce(rbind, lapply(unique_labels, function(x){
        idx_sc <- which(sc_labels == x)
        idx_gene <- which(colSums(Z[idx_sc,]) != 0)
        gene_list[[x]] <- gene_intersect[idx_gene]
        pca_res <- prcomp(Z[idx_sc, idx_gene], center = TRUE, scale. = TRUE)
        pca_res_list[[x]] <- pca_res
        ndims <- min(30, length(idx_sc))
        Z_pca <- as.data.frame(pca_res$x[, 1:ndims])
        umap_settings <- umap::umap.defaults
        umap_settings$n_neighbors <- min(15, length(idx_sc))
        umap_res <- umap::umap(as.matrix(Z_pca), config = umap_settings)
        umap_res_list[[x]] <- umap_res
        Z_umap <- as.data.frame(umap_res$layout)
        temp <- Z_umap
        domain <- rep("sc", length(idx_sc))
        prop <- rep(0, length(idx_sc))
        name <- rownames(Z[idx_sc,])
        for (i in names(mu_list_all)){
            idx_mu <- which(proportion_list[[i]][, x] > min_prop)
            if (length(idx_mu) > 0){
                mu <- mu_list_all[[i]][[x]][idx_mu, idx_gene, drop = FALSE]
                mu_pca <- as.data.frame(predict(pca_res, mu)[, 1:ndims, drop = FALSE])
                mu_umap <- as.data.frame(predict(umap_res, mu_pca))
                rownames(mu_umap) <- paste(i, rownames(mu_umap), sep = "_")
                temp <- rbind(temp, mu_umap)
                domain <- c(domain, rep(i, length(idx_mu)))
                prop <- c(prop, proportion_list[[i]][idx_mu, x])
                name <- c(name, rownames(mu))
            }
        }
        temp <- cbind(temp, domain, prop, x, name)
        colnames(temp) <- c("dim1", "dim2", "domain", "prop", "label", "name")
        if (!is.null(archetypes_list)){
            arche <- archetypes_list[[x]][, gene_intersect[idx_gene]]
            arche_pca <- as.data.frame(predict(pca_res, arche)[, 1:ndims])
            arche_umap <- as.data.frame(predict(umap_res, arche_pca))
            arche_temp <- as.data.frame(cbind(arche_umap, "archetype", rep(0, dim(arche)[1]), x, substring(rownames(arche), nchar(x) + 2)))
            colnames(arche_temp) <- c("dim1", "dim2", "domain", "prop", "label", "name")
            temp <- rbind(temp, arche_temp)
        }
        temp
    }))

    pdf(paste(save_dir, paste(f_name, "reference_only.pdf", sep = "_"), sep = "/"), width = width, height = height)
    fig <- ggplot() +
            # geom_point(aes(x = dim1, y = dim2), color = "#989797", data = plot_df[plot_df$domain == "sc",]) +
            geom_point(aes(x = dim1, y = dim2, color = as.factor(domain), shape = as.factor(domain), size = as.factor(domain)), data = plot_df[(plot_df$domain == "sc") | (plot_df$domain == "archetype"),]) +
            geom_text_repel(aes(x = dim1, y = dim2, label = name), data = plot_df[plot_df$domain == "archetype",], fontface = "bold") +
            labs(title = paste(f_name, "reference_only", sep = "_"), x = "UMAP_1", y = "UMAP_2") +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text(size = 10, face = "bold"),
                #   legend.title = element_blank(),
                #   legend.key = element_rect(color = "transparent", fill = "white"),
                #   legend.key.size = unit(0.45, "cm"),
                  legend.position = "right") +
            scale_color_manual(values = c(sc = "#989797", archetype = "#00A3D5")) +
            scale_shape_manual(values = c(sc = 16, archetype = 17)) +
            scale_size_manual(values = c(sc = 2, archetype = 2.5)) +
            facet_wrap(~label, scales = "free", ncol = n_cols)
    print(fig)
    dev.off()

    pdf(paste(save_dir, paste(f_name, "pdf", sep = "."), sep = "/"), width = width, height = height)
    fig <- ggplot() +
            geom_point(aes(x = dim1, y = dim2), color = "#989797", shape = 16, size = 2, data = plot_df[plot_df$domain == "sc",]) +
            geom_point(aes(x = dim1, y = dim2, color = as.factor(domain)), shape = 16, size = 2, data = plot_df[(plot_df$domain != "sc") & (plot_df$domain != "archetype"),]) +
            labs(title = f_name, x = "UMAP_1", y = "UMAP_2") +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
                  plot.background  = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text(size = 10, face = "bold"),
                #   legend.title = element_blank(),
                #   legend.key = element_rect(color = "transparent", fill = "white"),
                #   legend.key.size = unit(0.45, "cm"),
                  legend.position = "right") +
            scale_color_manual(values = use_color) +
            facet_wrap(~label, scales = "free", ncol = n_cols)
    print(fig)
    dev.off()

    if (return_res){
        return(list(plot_df = plot_df, pca_res_list = pca_res_list, umap_res_list = umap_res_list, gene_list = gene_list))
    }
}



plot_enrichGO <- function(genes, n_category = 15, object = "human", fname = "enrichGO", save_dir = ".", width = 7, height = 7, return_results = FALSE){
    library(clusterProfiler)
    library(ggplot2)
    if (object == "mouse"){
        library(org.Mm.eg.db)
        OrgDb <- "org.Mm.eg.db"
    } else{
        library(org.Hs.eg.db)
        OrgDb <- "org.Hs.eg.db"
    }
    ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
    ego <- enrichGO(ids$ENTREZID, OrgDb = OrgDb, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
    if (dim(ego)[1] > 0){
        write.csv(ego, file = paste(save_dir, paste(fname, "csv", sep = "."), sep = "/"))
        pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
        fig <- barplot(ego, showCategory = n_category, drop = TRUE) +
                labs(title = fname) +
                theme(axis.text.y = element_text(size = 10))
        print(fig)
        dev.off()
    } else{
        fig <- NULL
    }
    if (return_results){
        return(list(ego = ego, fig = fig))
    }
}



plot_enrichKEGG <- function(genes, n_category = 15, object = "human", organism = "hsa", fname = "enrichKEGG", save_dir = ".", width = 7, height = 7, return_results = FALSE){
    library(clusterProfiler)
    library(ggplot2)
    if (object == "mouse"){
        library(org.Mm.eg.db)
        OrgDb <- "org.Mm.eg.db"
    } else{
        library(org.Hs.eg.db)
        OrgDb <- "org.Hs.eg.db"
    }
    ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
    kk <- enrichKEGG(ids$ENTREZID, organism = organism, keyType = "kegg", pAdjustMethod = "BH")
    if (dim(kk)[1] > 0){
        write.csv(kk, file = paste(save_dir, paste(fname, "csv", sep = "."), sep = "/"))
        pdf(paste(save_dir, paste(fname, "pdf", sep = "."), sep = "/"), width = width, height = height)
        fig <- barplot(kk, showCategory = n_category, drop = TRUE) +
                labs(title = fname) +
                theme(axis.text.y = element_text(size = 10))
        print(fig)
        dev.off()
    } else{
        fig <- NULL
    }
    if (return_results){
        return(list(kk = kk, fig = fig))
    }
}

