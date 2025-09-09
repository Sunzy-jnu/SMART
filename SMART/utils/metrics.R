cal_prop_pearson <- function(prediction, truth, level = "all"){
    prediction <- as.matrix(prediction)
    truth <- as.matrix(truth)
    prediction <- prediction[rownames(truth), colnames(truth)]
    if (level == "spot"){
        prediction <- prediction / rowSums(prediction)
        truth <- truth / rowSums(truth)
        results <- diag(cor(t(prediction), t(truth), method = "pearson"))
    } else if (level == "type") {
        prediction <- t(t(prediction) / colSums(prediction))
        truth <- t(t(truth) / colSums(truth))
        results <- diag(cor(prediction, truth, method = "pearson"))
    } else{
        prediction <- c(prediction)
        truth <- c(truth)
        prediction <- prediction / sum(prediction)
        truth <- truth / sum(truth)
        results <- cor(prediction, truth, method = "pearson")
    }
    return(results)
}



cal_prop_RMSE <- function(prediction, truth, level = "all"){
    prediction <- as.matrix(prediction)
    truth <- as.matrix(truth)
    prediction <- prediction[rownames(truth), colnames(truth)]
    if (level == "spot"){
        prediction <- prediction / rowSums(prediction)
        truth <- truth / rowSums(truth)
        results <- sapply(seq(dim(prediction)[1]), function(x){sqrt(mean((prediction[x,] - truth[x,])^2))})
        names(results) <- rownames(prediction)
    } else if (level == "type") {
        prediction <- t(t(prediction) / colSums(prediction))
        truth <- t(t(truth) / colSums(truth))
        results <- sapply(seq(dim(prediction)[2]), function(x){sqrt(mean((prediction[,x] - truth[,x])^2))})
        names(results) <- colnames(prediction)
    } else{
        prediction <- c(prediction)
        truth <- c(truth)
        prediction <- prediction / sum(prediction)
        truth <- truth / sum(truth)
        results <- sqrt(mean((prediction - truth)^2))
    }
    return(results)
}



cal_prop_JSD <- function(prediction, truth, level = "all"){
    prediction <- as.matrix(prediction)
    truth <- as.matrix(truth)
    prediction <- prediction[rownames(truth), colnames(truth)]
    if (level == "spot"){
        prediction <- prediction / rowSums(prediction)
        truth <- truth / rowSums(truth)
        results <- sapply(seq(dim(prediction)[1]), function(x){cal_JSD(prediction[x,], truth[x,])})
        names(results) <- rownames(prediction)
    } else if (level == "type") {
        prediction <- t(t(prediction) / colSums(prediction))
        truth <- t(t(truth) / colSums(truth))
        results <- sapply(seq(dim(prediction)[2]), function(x){cal_JSD(prediction[,x], truth[,x])})
        names(results) <- colnames(prediction)
    } else{
        prediction <- c(prediction)
        truth <- c(truth)
        prediction <- prediction / sum(prediction)
        truth <- truth / sum(truth)
        results <- cal_JSD(prediction, truth)
    }
    return(results)
}



cal_JSD <- function(a, b){
    a <- pmax(a, 0)
    b <- pmax(b, 0)
    a <- a / sum(a)
    b <- b / sum(b)
    ab <- (a + b) / 2
    KL_1 <- sum(a * log2(a / (ab + 1e-12) + 1e-12))
    KL_2 <- sum(b * log2(b / (ab + 1e-12) + 1e-12))
    return((KL_1 + KL_2) / 2)
}



cal_signature_pearson <- function(prediction_list, truth_list, prop_truth, min_prop = 0, keep_genes = NULL, level = "all"){

    if (is.null(keep_genes)){
        keep_genes <- intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)})))
    }
    keep_genes <- intersect(keep_genes, intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)}))))

    prediction_list <- lapply(prediction_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    truth_list <- lapply(truth_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    prediction_list <- prediction_list[names(truth_list)]

    if (level == "spot"){
        n_spots <- dim(prop_truth)[1]
        results <- sapply(seq(n_spots), function(x){
            keep_labels <- colnames(prop_truth)[which(prop_truth[x,] > min_prop)]
            prediction_x <- as.numeric(Reduce(c, lapply(prediction_list[keep_labels], function(y){y[x,]})))
            truth_x <- as.numeric(Reduce(c, lapply(truth_list[keep_labels], function(y){y[x,]})))
            cor(prediction_x, truth_x, method = "pearson")
        })
        names(results) <- rownames(prop_truth)
    } else if (level == "type"){
        unique_labels <- colnames(prop_truth)
        results <- sapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            prediction_x <- as.numeric(c(as.matrix(prediction_list[[x]][keep_spots,])))
            truth_x <- as.numeric(c(as.matrix(truth_list[[x]][keep_spots,])))
            cor(prediction_x, truth_x, method = "pearson")
        })
        names(results) <- unique_labels
    } else{
        unique_labels <- colnames(prop_truth)
        prediction <- as.numeric(Reduce(c, lapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            c(as.matrix(prediction_list[[x]][keep_spots,]))
        })))
        truth <- as.numeric(Reduce(c, lapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            c(as.matrix(truth_list[[x]][keep_spots,]))
        })))
        results <- cor(prediction, truth, method = "pearson")
    }

    return(results)
}



cal_signature_RMSE <- function(prediction_list, truth_list, prop_truth, min_prop = 0, keep_genes = NULL, level = "all"){

    if (is.null(keep_genes)){
        keep_genes <- intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)})))
    }
    keep_genes <- intersect(keep_genes, intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)}))))

    prediction_list <- lapply(prediction_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    truth_list <- lapply(truth_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    prediction_list <- prediction_list[names(truth_list)]

    if (level == "spot"){
        n_spots <- dim(prop_truth)[1]
        results <- sapply(seq(n_spots), function(x){
            keep_labels <- colnames(prop_truth)[which(prop_truth[x,] > min_prop)]
            prediction_x <- as.numeric(Reduce(c, lapply(prediction_list[keep_labels], function(y){y[x,]})))
            truth_x <- as.numeric(Reduce(c, lapply(truth_list[keep_labels], function(y){y[x,]})))
            sqrt(mean((prediction_x - truth_x)^2))
        })
        names(results) <- rownames(prop_truth)
    } else if (level == "type"){
        unique_labels <- colnames(prop_truth)
        results <- sapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            prediction_x <- as.numeric(c(as.matrix(prediction_list[[x]][keep_spots,])))
            truth_x <- as.numeric(c(as.matrix(truth_list[[x]][keep_spots,])))
            sqrt(mean((prediction_x - truth_x)^2))
        })
        names(results) <- unique_labels
    } else{
        unique_labels <- colnames(prop_truth)
        prediction <- as.numeric(Reduce(c, lapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            c(as.matrix(prediction_list[[x]][keep_spots,]))
        })))
        truth <- as.numeric(Reduce(c, lapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            c(as.matrix(truth_list[[x]][keep_spots,]))
        })))
        results <- sqrt(mean((prediction - truth)^2))
    }

    return(results)
}



cal_signature_JSD <- function(prediction_list, truth_list, prop_truth, min_prop = 0, keep_genes = NULL, level = "all"){

    if (is.null(keep_genes)){
        keep_genes <- intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)})))
    }
    keep_genes <- intersect(keep_genes, intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)}))))

    prediction_list <- lapply(prediction_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    truth_list <- lapply(truth_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    prediction_list <- prediction_list[names(truth_list)]

    if (level == "spot"){
        n_spots <- dim(prop_truth)[1]
        results <- sapply(seq(n_spots), function(x){
            keep_labels <- colnames(prop_truth)[which(prop_truth[x,] > min_prop)]
            prediction_x <- as.numeric(Reduce(c, lapply(prediction_list[keep_labels], function(y){y[x,]})))
            truth_x <- as.numeric(Reduce(c, lapply(truth_list[keep_labels], function(y){y[x,]})))
            cal_JSD(prediction_x, truth_x)
        })
        names(results) <- rownames(prop_truth)
    } else if (level == "type"){
        unique_labels <- colnames(prop_truth)
        results <- sapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            prediction_x <- as.numeric(c(as.matrix(prediction_list[[x]][keep_spots,])))
            truth_x <- as.numeric(c(as.matrix(truth_list[[x]][keep_spots,])))
            cal_JSD(prediction_x, truth_x)
        })
        names(results) <- unique_labels
    } else{
        unique_labels <- colnames(prop_truth)
        prediction <- as.numeric(Reduce(c, lapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            c(as.matrix(prediction_list[[x]][keep_spots,]))
        })))
        truth <- as.numeric(Reduce(c, lapply(unique_labels, function(x){
            keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
            c(as.matrix(truth_list[[x]][keep_spots,]))
        })))
        results <- cal_JSD(prediction, truth)
    }

    return(results)
}



cal_signature_pearson_individual <- function(prediction_list, truth_list, prop_truth, min_prop = 0, keep_genes = NULL){

    if (is.null(keep_genes)){
        keep_genes <- intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)})))
    }
    keep_genes <- intersect(keep_genes, intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)}))))

    prediction_list <- lapply(prediction_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    truth_list <- lapply(truth_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    prediction_list <- prediction_list[names(truth_list)]

    unique_labels <- colnames(prop_truth)
    results <- lapply(unique_labels, function(x){
        keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
        temp <- sapply(keep_spots, function(y){cor(as.numeric(prediction_list[[x]][y,]), as.numeric(truth_list[[x]][y,]), method = "pearson")})
        temp <- cbind(prop_truth[keep_spots, x], temp)
        colnames(temp) <- c("proportion", "signature_PCC")
        temp
    })
    names(results) <- unique_labels

    return(results)
}



cal_signature_RMSE_individual <- function(prediction_list, truth_list, prop_truth, min_prop = 0, keep_genes = NULL){

    if (is.null(keep_genes)){
        keep_genes <- intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)})))
    }
    keep_genes <- intersect(keep_genes, intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)}))))

    prediction_list <- lapply(prediction_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    truth_list <- lapply(truth_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    prediction_list <- prediction_list[names(truth_list)]

    unique_labels <- colnames(prop_truth)
    results <- lapply(unique_labels, function(x){
        keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
        temp <- sapply(keep_spots, function(y){sqrt(mean((as.numeric(prediction_list[[x]][y,]) - as.numeric(truth_list[[x]][y,]))^2))})
        temp <- cbind(prop_truth[keep_spots, x], temp)
        colnames(temp) <- c("proportion", "signature_RMSE")
        temp
    })
    names(results) <- unique_labels

    return(results)
}



cal_signature_JSD_individual <- function(prediction_list, truth_list, prop_truth, min_prop = 0, keep_genes = NULL){

    if (is.null(keep_genes)){
        keep_genes <- intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)})))
    }
    keep_genes <- intersect(keep_genes, intersect(Reduce(intersect, lapply(prediction_list, function(x){colnames(x)})), Reduce(intersect, lapply(truth_list, function(x){colnames(x)}))))

    prediction_list <- lapply(prediction_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    truth_list <- lapply(truth_list, function(x){
        temp <- x[, keep_genes]
        temp <- temp / rowSums(temp)
        temp
    })
    prediction_list <- prediction_list[names(truth_list)]

    unique_labels <- colnames(prop_truth)
    results <- lapply(unique_labels, function(x){
        keep_spots <- rownames(prop_truth)[which(prop_truth[, x] > min_prop)]
        temp <- sapply(keep_spots, function(y){cal_JSD(as.numeric(prediction_list[[x]][y,]), as.numeric(truth_list[[x]][y,]))})
        temp <- cbind(prop_truth[keep_spots, x], temp)
        colnames(temp) <- c("proportion", "signature_JSD")
        temp
    })
    names(results) <- unique_labels

    return(results)
}


