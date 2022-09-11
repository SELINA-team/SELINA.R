#' Fine-tuning and predict for the query data
#'
#' @param queryObj Seurat object for query data.
#' @param model The pre-trained model.
#' @param path_out File path of the output files.
#' @param outprefix Prefix of the output files. DEFAULT: query.
#' @param disease Depend on your data is in some disease condition or not, choose from 'TRUE' or 'FALSE'. DEFAULT: 'FALSE'.
#' @param prob_cutoff Cutoff for prediction probability. DEFAULT: 0.9.
#'
#' @return A list which includes prediction results and corresponding probability for query data.
#'
#' @importFrom stats rnorm
#' @importFrom utils read.csv write.table
#' @importFrom torch nn_relu cuda_is_available torch_device torch_tensor torch_unsqueeze dataset nn_module dataloader nn_cross_entropy_loss nn_mse_loss optim_adam torch_float torch_load torch_save autograd_function with_no_grad
#' @importFrom Seurat FindMarkers GetAssayData DimPlot
#' @importFrom presto wilcoxauc sumGroups
#' @importFrom dplyr filter arrange desc
#' @importFrom ggplot2 ggsave
#' @export
#'
#'
query_predict <- function(queryObj, model, path_out, outprefix = "query", disease = FALSE, prob_cutoff = 0.9) {
  device <- if (torch::cuda_is_available()) torch::torch_device("cuda:0") else "cpu"
  params_tune1 <- c(0.0005, 50, 128)
  params_tune2 <- c(0.0001, 10, 128)
  message("Loading data")
  query_expr <- queryObj@assays$RNA@counts
  meta <- model$meta
  if (disease) {
    genes <- meta$genes
    ct_dic <- meta$celltypes
    disease_dic <- meta$cellsources
  } else {
    genes <- meta$genes
    ct_dic <- meta$celltypes
  }
  query_expr <- merge_query(genes, query_expr)
  nfeatures <- length(genes)
  nct <- length(ct_dic)
  network <- model$network
  if (disease) {
    network <- Disease_Classifier(network)$to(device = device)
    test_res <- test(query_expr, network, ct_dic, disease_dic, disease, device)
    pred_labels <- test_res$pred_labels
    pred_prob <- test_res$pred_prob
    disease_labels <- test_res$disease_labels
    disease_labels <- data.frame(Cell = colnames(query_expr), Prediction = disease_labels, stringsAsFactors = F)
    write.table(disease_labels, file.path(path_out, paste0(outprefix, "_cellsources.txt")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  } else {
    network <- Autoencoder(network, nfeatures, nct)
    message("Fine-tuning1")
    network <- tune1(query_expr, network, params_tune1, device)
    message("Fine-tuning2")
    network <- tune2(query_expr, network, params_tune2, device)
    network <- Normal_Classifier(network)$to(device = device)
    test_res <- test(query_expr, network, ct_dic, NULL, disease, device)
  }
  pred_labels <- test_res$pred_labels
  pred_prob <- test_res$pred_prob
  if (disease) {
    disease_labels <- test_res$disease_labels
  }

  write.table(data.frame(pred_prob, stringsAsFactors = F), file.path(path_out, paste0(outprefix, "_probability.txt")), col.names = TRUE, quote = FALSE, sep = "\t")
  message("Finish prediction")

  message("Begin downstream analysis")
  predict_downstream(
    queryObj, prob_cutoff, path_out,
    outprefix, pred_labels, pred_prob
  )
  message("Finish downstream analysis")
  if (disease) {
    return(list(pred_labels = pred_labels, pred_prob = pred_prob, disease_labels = disease_labels))
  } else {
    return(list(pred_labels = pred_labels, pred_prob = pred_prob))
  }
}

merge_query <- function(genes, query_expr) {
  model_expr <- data.frame(genes = rnorm(length(genes)), stringsAsFactors = F)
  rownames(model_expr) <- genes
  query_expr <- merge(model_expr, query_expr, by = "row.names", all.x = T)
  rownames(query_expr) <- query_expr$Row.names
  query_expr <- query_expr[, -match("genes", colnames(query_expr))]
  query_expr <- query_expr[, -1]
  query_expr[is.na(query_expr)] <- 0
  query_expr <- query_expr / colSums(query_expr) * 10000 # or matrix
  query_expr <- log2(query_expr + 1)
  return(query_expr)
}

Datasets <- torch::dataset(
  initialize = function(data) {
    self$expr <- data
  },
  .getitem = function(index) {
    torch::torch_tensor(self$expr[, index])
  },
  .length = function() {
    ncol(self$expr)
  }
)
# Autoencoder(network, nfeatures, nct)
Autoencoder <- torch::nn_module(
  "class_Autoencoder",
  initialize = function(network, nfeature, nct) {
    encoder <- c(network$feature$children, network$class_classifier$children)
    encoder_index <- c(1, 2, 4, 5, 7)
    self$encoder <- do.call(torch::nn_sequential, c(unlist(sapply(encoder_index, function(x) encoder[[x]]), use.names = F), torch::nn_relu()))
    self$decoder <- torch::nn_sequential(
      torch::nn_linear(in_features = nct, out_features = 50), torch::nn_relu(),
      torch::nn_linear(in_features = 50, out_features = 100), torch::nn_relu(),
      torch::nn_linear(in_features = 100, out_features = nfeature)
    )
  },
  forward = function(input_data) {
    output <- self$decoder(self$encoder(input_data))
    return(output)
  }
)


Normal_Classifier <- torch::nn_module(
  "Normal_Classifier",
  initialize = function(network) {
    self$classifier <- do.call(torch::nn_sequential, c(unlist(network$encoder$children[-length(network$encoder$children)], use.names = F), torch::nn_softmax(dim = 2)))
  },
  forward = function(input_data) {
    output <- self$classifier(input_data)
    return(output)
  }
)

Disease_Classifier <- torch::nn_module(
  "Disease_Classifier",
  initialize = function(network) {
    self$feature <- do.call(torch::nn_sequential, c(unlist(sapply(1:2, function(x) network$feature$children[[x]]), use.names = F)))
    self$celltype <- do.call(torch::nn_sequential, c(unlist(sapply(c(1, 2, 4), function(x) network$class_classifier$children[[x]]), use.names = F), torch::nn_softmax(dim = 2)))
    self$disease <- do.call(torch::nn_sequential, c(unlist(network$disease_classifier$children, use.names = F), torch::nn_softmax(dim = 2)))
  },
  forward = function(input_data) {
    celltype <- self$celltype(self$feature(input_data))
    disease <- self$disease(self$feature(input_data))
    return(list(celltype = celltype, disease = disease))
  }
)

tune1 <- function(test_df, network, params, device) {
  test_dat <- Datasets(test_df)
  lr <- params[1]
  n_epoch <- params[2]
  batch_size <- params[3]
  optimizer <- torch::optim_adam(network$parameters, lr = lr)
  loss <- torch::nn_mse_loss()
  loss <- loss$to(device = device)
  test_loader <- torch::dataloader(dataset = test_dat, batch_size = batch_size, shuffle = TRUE)
  for (i in 1:length(network$encoder$named_parameters())) {
    res <- network$encoder$named_parameters()[[i]]
    res$requires_grad <- FALSE
  }
  network <- network$to(device = device)
  for (epoch in 1:n_epoch) {
    coro::loop(for (batch in test_loader) {
      expr <- batch
      expr <- expr$to(dtype = torch_float())
      expr <- expr$to(device = device)
      output <- network(expr)
      err <- loss(output, expr)
      optimizer$zero_grad()
      err$backward()
      optimizer$step()
    })
  }
  message("Finish tuning1")
  return(network)
}


tune2 <- function(test_df, network, params, device) {
  test_dat <- Datasets(test_df)
  lr <- params[1]
  n_epoch <- params[2]
  batch_size <- params[3]
  optimizer <- torch::optim_adam(network$parameters, lr = lr)
  loss <- torch::nn_mse_loss()
  loss <- loss$to(device = device)
  test_loader <- torch::dataloader(dataset = test_dat, batch_size = batch_size, shuffle = TRUE)
  for (i in 1:length(network$encoder$named_parameters())) {
    res <- network$encoder$named_parameters()[[i]]
    res$requires_grad <- TRUE
  }
  for (i in 1:length(network$decoder$named_parameters())) {
    res <- network$encoder$named_parameters()[[i]]
    res$requires_grad <- FALSE
  }
  network <- network$to(device = device)
  for (epoch in 1:n_epoch) {
    coro::loop(for (batch in test_loader) {
      expr <- batch
      expr <- expr$to(dtype = torch_float())
      expr <- expr$to(device = device)
      output <- network(expr)
      err <- loss(output, expr)
      optimizer$zero_grad()
      err$backward()
      optimizer$step()
    })
  }
  message("Finish tuning2")
  return(network)
}

test <- function(test_df, network, ct_dic, disease_dic, disease, device) {
  test_dat <- Datasets(test_df)
  ct_dic_rev <- split(rep(names(ct_dic), sapply(ct_dic, length)), unlist(ct_dic))
  if (disease) {
    disease_dic_rev <- split(rep(names(disease_dic), sapply(disease_dic, length)), unlist(disease_dic))
  }
  test_loader <- torch::dataloader(dataset = test_dat, batch_size = ncol(test_df), shuffle = FALSE)
  pred_labels <- c()
  disease_labels <- c()
  pred_prob <- list()
  with_no_grad(
    coro::loop(for (batch in test_loader) {
      expr <- batch
      expr <- expr$to(dtype = torch_float())
      expr <- expr$to(device = device)
      if (disease) {
        class_output <- network(expr)[[1]]
        disease_output <- network(expr)[[2]]
        disease_labels <- c(disease_labels, as.numeric(disease_output$argmax(dim = 2)$cpu()))
      } else {
        class_output <- network(expr)
      }
      pred_labels <- c(pred_labels, as.numeric(class_output$argmax(dim = 2)$cpu()))
      pred_prob <- append(pred_prob, list(as.matrix(class_output$cpu())))
    })
  )
  pred_labels <- sapply(pred_labels, function(x) ct_dic_rev[[x]])
  pred_prob <- do.call("cbind", pred_prob)
  rownames(pred_prob) <- colnames(test_df)
  colnames(pred_prob) <- names(ct_dic)
  if (disease) {
    disease_labels <- sapply(disease_labels, function(x) disease_dic_rev[[x]])
    return(list(pred_labels = pred_labels, pred_prob = pred_prob, disease_labels = disease_labels))
  } else {
    return(list(pred_labels = pred_labels, pred_prob = pred_prob))
  }
}

predict_downstream <- function(queryObj, prob_cutoff, path_out, outprefix, pred_labels, pred_prob) {
  # main step
  # Filter cells with low prediction score
  filtered_label <- pred_filter(pred_labels, pred_prob, prob_cutoff)

  # Find differentially expressed genes for each cell type
  if (length(unique(filtered_label[filtered_label != "Unknown"])) > 1) {
    cluster.genes <- FindMarkers(object = queryObj[, filtered_label != "Unknown"], celltypes = filtered_label[filtered_label != "Unknown"])
    write.table(cluster.genes, file.path(path_out, paste0(outprefix, "_DiffGenes.tsv")), quote = F, sep = "\t", row.names = FALSE)
  }
  # Output umap plot with predicted cell type labels
  queryObj$pred <- filtered_label
  p <- DimPlot(object = queryObj[, filtered_label != "Unknown"], label = TRUE, pt.size = 0.2, repel = TRUE, group.by = "pred")
  ggsave(file.path(path_out, paste0(outprefix, "_pred.png")), p, width = 7, height = 5)
  write.table(data.frame(Cell = colnames(queryObj), Prediction = filtered_label, Cluster = queryObj$seurat_clusters), file.path(path_out, paste0(outprefix, "_predictions.txt")), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}

# Find differentially expressed genes
FindMarkers <- function(object, celltypes, features = NULL, min.pct = 0.1, logfc.threshold = 0.25,
                        only.pos = FALSE, return.thresh = 1e-2,
                        slot = "data") {
  matrix <- GetAssayData(object, slot = slot)
  features <- rownames(matrix)
  y <- celltypes
  y <- factor(y)
  test.res <- wilcoxauc(matrix, y)

  # Calculate logFC
  if (slot != "scale.data") {
    if (slot == "data") {
      X <- expm1(matrix)
    }
    group_sums <- sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      log(group_means[, g] + 1) - log(((cs - group_sums[g, ]) / (length(y) - gs[g])) + 1)
    }))
  } else {
    group_sums <- sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      group_means[, g] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
    }))
  }

  test.res$avg_logFC <- as.vector(lfc)
  res <- test.res[, c("pval", "avg_logFC", "pct_in", "pct_out", "padj", "group", "feature")]
  res[, c("pct_in", "pct_out")] <- round(res[, c("pct_in", "pct_out")] / 100, digits = 3)
  colnames(res) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "celltype", "gene")
  res <- res %>% dplyr::filter(.data$p_val < return.thresh &
    abs(.data$avg_logFC) > logfc.threshold &
    (.data$pct.1 > min.pct |
      .data$pct.2 > min.pct) &
    .data$gene %in% features)
  if (only.pos) {
    res <- res %>% dplyr::filter(.data$avg_logFC > 0)
  }
  res <- res %>% dplyr::arrange(.data$celltype, .data$p_val, dplyr::desc(.data$avg_logFC))
  return(res)
}

# Filter cells with low prediction probability
pred_filter <- function(pred_label, pred_prob, prob_cutoff) {
  pred_prob <- apply(pred_prob, MARGIN = 1, max)
  filtered_label <- pred_label
  for (cell_index in 1:length(pred_label)) {
    prob <- pred_prob[cell_index]
    if (prob < prob_cutoff) {
      filtered_label[cell_index] <- "Unknown"
    }
  }
  return(filtered_label)
}
