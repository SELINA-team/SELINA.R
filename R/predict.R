#' Fine-tuning and predict for the query data
#'
#' @param queryObj Seurat object for query data.
#' @param model The pre-trained model.
#' @param path_out File path of the output files.
#' @param outprefix Prefix of the output files. DEFAULT: query.
#' @param disease Depend on your data is in some disease condition or not, choose from 'TRUE' or 'FALSE'. DEFAULT: 'FALSE'.
#' @param cell_cutoff Cutoff for number of cells with the same cell type in 10 nearest neighbor cells(only used when the input is single-cell level expression matrix). DEFAULT: 5.
#' @param prob_cutoff Cutoff for prediction probability. DEFAULT: 0.9.
#'
#' @return A list which includes prediction results and corresponding probability for query data. This step will generate eight files. 
#'
#' @importFrom stats rnorm
#' @importFrom utils read.csv write.table
#' @importFrom torch nn_relu cuda_is_available torch_device torch_tensor torch_unsqueeze dataset nn_module dataloader nn_cross_entropy_loss nn_mse_loss optim_adam torch_float torch_load torch_save autograd_function with_no_grad
#' @import Seurat
#' @import ggplot2
#' @import dbscan
#' @import gtools
#' @import presto
#' @export
#'
query_predict <- function(queryObj, model, path_out, outprefix='query', disease = FALSE, cell_cutoff = 5, prob_cutoff = 0.9) {
  device <- if (torch::cuda_is_available()) torch::torch_device("cuda:0") else "cpu"
  params_tune1 <- c(0.0005, 50, 128)
  params_tune2 <- c(0.0001, 10, 128)

  message("Loading data")
  query_expr <- queryObj@assays$RNA@counts
  meta <- model$meta
  genes <- meta$genes
  ct_dic <- meta$celltypes
  query_expr <- merge_query(genes, query_expr)
  nfeatures <- length(genes)
  nct <- length(ct_dic)
  network <- model$network
  network <- Autoencoder(network, nfeatures, nct)
  if (!disease) {
    message("Fine-tuning1")
    network <- tune1(query_expr, network, params_tune1, device)
    message("Fine-tuning2")
    network <- tune2(query_expr, network, params_tune2, device)
  }
  network <- Classifier(network)$to(device = device)
  test_res <- test(query_expr, network, ct_dic, device)
  pred_labels <- test_res$pred_labels
  pred_prob <- test_res$pred_prob
  
  write.table(data.frame(pred_prob, stringsAsFactors = F), file.path(path_out, paste0(outprefix, "_probability.txt")), col.names = TRUE, quote = FALSE, sep = "\t")

  message("Finish Prediction")

  message("Begin downstream analysis")
  predict_downstream(
    queryObj, cell_cutoff, prob_cutoff, path_out,
    outprefix, pred_labels, pred_prob
  )
  message("Finish downstream analysis")
  return(list(pred_labels=pred_labels, pred_prob=pred_prob))
}

merge_query <- function(genes, query_expr) {
  model_expr <- data.frame(genes = rnorm(length(genes)), stringsAsFactors = F)
  rownames(model_expr) <- genes
  query_expr <- merge(model_expr, query_expr, by = "row.names", all.x=T)
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


Classifier <- torch::nn_module(
  "class_Classifier",
  initialize = function(network) {
    self$classifier <- do.call(torch::nn_sequential, c(unlist(network$encoder$children[-length(network$encoder$children)], use.names = F), torch::nn_softmax(dim = 2)))
  },
  forward = function(input_data) {
    output <- self$classifier(input_data)
    return(output)
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
  message("Finish Tuning1")
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
  message("Finish Tuning2")
  return(network)
}


test <- function(test_df, network, ct_dic, device) {
  test_dat <- Datasets(test_df)
  ct_dic_rev <- split(rep(names(ct_dic), sapply(ct_dic, length)), unlist(ct_dic))
  test_loader <- torch::dataloader(dataset = test_dat, batch_size = ncol(test_df), shuffle = FALSE)

  with_no_grad(
    coro::loop(for (batch in test_loader) {
      expr <- batch
      expr <- expr$to(dtype = torch_float())
      expr <- expr$to(device = device)
      class_output <- network(expr)
      pred_labels <- as.numeric(class_output$argmax(dim = 2)$cpu())
      pred_prob <- class_output
    })
  )
  pred_labels <- sapply(pred_labels, function(x) ct_dic_rev[[x]])
  pred_prob <- as.matrix(pred_prob$cpu())
  rownames(pred_prob) <- colnames(test_df)
  colnames(pred_prob) <- names(ct_dic)
  return(list(pred_labels = pred_labels, pred_prob = pred_prob))
}


predict_downstream <- function(queryObj, cell_cutoff, prob_cutoff, path_out, outprefix, pred_labels, pred_prob) {
  # main step
  # Filter cells with low prediction score
  filtered_label <- pred_filter(queryObj, pred_labels, pred_prob, cell_cutoff, prob_cutoff)

  # Generate plots and files indicating the prediction quality for each cluster
  cluster_quality(queryObj, filtered_label, pred_prob, path_out, outprefix)

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
  res <- res %>% dplyr::arrange(.data$celltype, .data$p_val, desc(.data$avg_logFC))
  return(res)
}

# Filter cells with low prediction probability
pred_filter <- function(queryObj, pred_labels, pred_prob, cell_cutoff, prob_cutoff) {
  neighbor_id <- kNN(queryObj@reductions[["umap"]]@cell.embeddings, k = 10)$id
  pred_prob <- apply(pred_prob, MARGIN = 1, max)
  filtered_label <- pred_labels
  for (cell_index in 1:length(pred_labels)) {
    neighbor_cts <- pred_labels[neighbor_id[cell_index, ]]
    ct <- pred_labels[cell_index]
    prob <- pred_prob[cell_index]
    nsame <- length(which(neighbor_cts == ct))
    if (nsame < cell_cutoff | prob < prob_cutoff) {
      filtered_label[cell_index] <- "Unknown"
    }
  }
  return(filtered_label)
}

theme_box <- function(...) {
  require(grid)
  theme(
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "black"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(family = "Helvetica", size = 15, vjust = 0.5, hjust = 0.5, margin = margin(t = 10, b = 10)),
    axis.text.x = element_text(family = "Helvetica", size = 15, hjust = 0.5, vjust = 0.5, margin = margin(t = 5, b = 5), colour = "black"),
    axis.text.y = element_text(family = "Helvetica", size = 15, margin = margin(l = 5, r = 5), colour = "black"),
    axis.title.x = element_text(family = "Helvetica", size = 14, colour = "black"),
    axis.title.y = element_text(family = "Helvetica", size = 14, colour = "black"),
    axis.ticks.y = element_line(size = 0.1),
    legend.title = element_text(family = "Helvetica", size = 13, colour = "black"),
    legend.text = element_text(family = "Helvetica", size = 13, colour = "black"),
    legend.key.size = unit(20, "pt"),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.key = element_blank(),
  )
}

# Generate plots and files indicating the prediction quality for each cluster
cluster_quality <- function(queryObj, filtered_label, pred_prob, path_out, outprefix) {
  pred_prob <- apply(pred_prob, MARGIN = 1, max)
  unknown_percent <- c()
  for (cluster in 1:length(unique(queryObj$seurat_clusters))) {
    cluster_index <- which(queryObj$seurat_clusters == as.character(cluster - 1))
    unknown_percent <- c(unknown_percent, length(which(filtered_label[cluster_index] == "Unknown")) / length(cluster_index))
    names(unknown_percent)[cluster] <- as.character(cluster - 1)
  }
  pred_prob <- data.frame(prob = pred_prob, cluster = as.character(queryObj$seurat_clusters))
  unknown_percent <- data.frame(unknown_percent = unknown_percent, cluster = names(unknown_percent))
  pred_prob$cluster <- factor(pred_prob$cluster, levels <- mixedsort(unique(pred_prob$cluster)))
  unknown_percent$cluster <- factor(unknown_percent$cluster, levels <- mixedsort(unique(pred_prob$cluster)))
  png(file.path(path_out, paste0(outprefix, "_cluster_prob.png")), width = 300 + 120 * length(unique(queryObj$seurat_clusters)), height = 1500, res = 300)
  print(
    ggplot(pred_prob) +
      geom_boxplot(aes(x = cluster, y = prob), colour = "#4DBBD5CC", width = 0.6, outlier.shape = NA, lwd = 0.2) +
      labs(title = "", y = "Prob.", x = "Cluster") +
      theme_box()
  )
  dev.off()
  png(file.path(path_out, paste0(outprefix, "_unknown_percent.png")), width = 300 + 120 * length(unique(queryObj$seurat_clusters)), height = 1500, res = 300)
  print(
    ggplot(unknown_percent) +
      geom_col(aes(x = cluster, y = unknown_percent), colour = "#91D1C2CC", fill = "#91D1C2CC", position = position_dodge(0.7), width = 0.7) +
      labs(title = "", y = "Unknown percent.", x = "Cluster") +
      theme_classic() +
      theme(
        plot.title = element_text(family = "Helvetica", size = 13, vjust = 0.5, hjust = 0.5, margin = margin(t = 10, b = 10)),
        axis.text.x = element_text(family = "Helvetica", size = 10, hjust = 0, vjust = 0.5, angle = -90, margin = margin(t = 5, b = 5)),
        axis.text.y = element_text(family = "Helvetica", size = 10, margin = margin(l = 5, r = 5)),
        axis.title.x = element_text(family = "Helvetica", size = 12),
        axis.title.y = element_text(family = "Helvetica", size = 12),
        legend.title = element_text(family = "Helvetica", size = 11),
        legend.text = element_text(family = "Helvetica", size = 10),
        legend.key.size = unit(12, "pt")
      )
  )
  dev.off()
  write.table(pred_prob, file.path(path_out, paste0(outprefix, "_prob.txt")), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(unknown_percent, file.path(path_out, paste0(outprefix, "_unknown_percent.txt")), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}

