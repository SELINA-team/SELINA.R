#' Train model based on your own data
#'
#' @param seuratlist A list of seurat object for training.
#'
#' @return An MADA model.
#' @importFrom stats runif
#' @importFrom utils read.csv write.table
#' @importFrom grDevices dev.off png
#' @importFrom utils read.table
#' @importFrom torch nn_sequential cuda_is_available torch_device torch_tensor torch_unsqueeze dataset nn_module dataloader nn_cross_entropy_loss nn_mse_loss optim_adam torch_float torch_load torch_save autograd_function with_no_grad

#' @export
#'
#' @examples model <- train_model(list(rds1, rds2))
train_model <- function(seuratlist) {
  if(class(seuratlist)!='list') { stop("Input should be a list")}
  params_train <- c(0.0001, 50, 128)
  device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"

  res_pre <- preprocessing(seuratlist)
  train_data <- res_pre$train_sets
  celltypes <- res_pre$celltypes
  platforms <- res_pre$platforms
  genes <- res_pre$genes

  ct_dic <- label2dic(celltypes)
  plat_dic <- label2dic(platforms)
  nfeatures <- nrow(train_data)
  nct <- length(ct_dic)
  nplat <- length(plat_dic)

  network <- train(
    train_data, params_train, celltypes, platforms, nfeatures,
    nct, nplat, ct_dic, plat_dic, device
  )
  model_meta <- list(genes = genes, celltypes = ct_dic)
  
  message("All done")
  return(list(network=network, meta=model_meta))
}

se_SMOTE <- function(X, target, target_name, K = 5, add_size = 1000) {
  ncD <- ncol(X)
  n_target <- table(target)
  multiP <- add_size / n_target[target_name]
  if (multiP <= 1) {
    sample_num <- add_size
    dup_size <- 1
  } else {
    dup_size <- round(multiP) + 1
    sample_num <- round(add_size / dup_size)
  }

  P_set <- subset(X, target == target_name)[sample(n_target[target_name], sample_num), ]
  N_set <- X[setdiff(rownames(X), rownames(P_set)), ]

  P_class <- rep(target_name, nrow(P_set))
  N_class <- target[which(!(rownames(X) %in% rownames(P_set)))]

  sizeP <- nrow(P_set)
  sizeN <- nrow(N_set)
  knear <- smotefamily::knearest(P_set, P_set, K)
  sum_dup <- smotefamily::n_dup_max(sizeP + sizeN, sizeP, sizeN, dup_size)
  syn_dat <- NULL
  for (i in 1:sizeP) {
    if (is.matrix(knear)) {
      pair_idx <- knear[i, ceiling(runif(sum_dup) * K)]
    } else {
      pair_idx <- rep(knear[i], sum_dup)
    }
    g <- runif(sum_dup)
    P_i <- matrix(unlist(P_set[i, ]), sum_dup, ncD, byrow = TRUE)
    Q_i <- as.matrix(P_set[pair_idx, ])
    syn_i <- P_i + g * (Q_i - P_i)
    syn_dat <- rbind(syn_dat, syn_i)
  }
  P_set[, ncD + 1] <- P_class
  colnames(P_set) <- c(colnames(X), "class")
  N_set[, ncD + 1] <- N_class
  colnames(N_set) <- c(colnames(X), "class")
  rownames(syn_dat) <- NULL
  syn_dat <- data.frame(syn_dat, check.names = F)
  syn_dat[, ncD + 1] <- rep(target_name, nrow(syn_dat))
  colnames(syn_dat) <- c(colnames(X), "class")
  NewD <- rbind(P_set, syn_dat, N_set)
  rownames(NewD) <- NULL
  D_result <- list(
    data = NewD, syn_data = syn_dat, orig_N = N_set,
    orig_P = P_set, K = K, K_all = NULL, dup_size = sum_dup,
    outcast = NULL, eps = NULL, method = "SMOTE"
  )
  class(D_result) <- "smote_data"
  return(D_result)
}

label2dic <- function(label) {
  label_set <- unique(label)
  label_list <- list()
  for (i in 1:length(label_set)) {
    label_list[label_set[i]] <- i
  }
  return(label_list)
}

preprocessing <- function(seuratlist) {
  message("Loading data")
  train_sets = sapply(seuratlist, function(rds) rds@assays$RNA@data, USE.NAMES = F)
  celltypes <- list()
  platforms <- list()
  genes <- list()
  for (i in 1:length(seuratlist)) {
    celltypes[[i]] <- seuratlist[[i]]$Celltype
    platforms[[i]] <- seuratlist[[i]]$Platform
    genes[[i]] <- rownames(seuratlist[[i]])
  }
  genes <- Reduce(intersect, genes)
  ct_freqs <- table(unlist(celltypes, use.names = F))
  max_n <- max(ct_freqs)
  rct_freqs <- {}
  if (max_n < 500) {
    sample_n <- 100
  } else if (max_n < 1000) {
    sample_n <- 500
  } else {
    sample_n <- 1000
  }
  for (i in 1:length(ct_freqs)) {
    ct <- names(ct_freqs[i])
    if (ct_freqs[i] < sample_n) {
      rct_freqs[ct] <- ct_freqs[i]
    }
  }
  if(length(rct_freqs)>0){
    message('Start SMOTE')
    for (i in 1:length(seuratlist)) {
      sample_ct_freq <- list()
      ct_freq <- table(celltypes[i])
      if (length(ct_freq) > 1) {
        for (j in 1:length(rct_freqs)) {
          ct <- names(rct_freqs[j])
          freq <- rct_freqs[j]
          if ((ct %in% names(ct_freq)) & (ct_freq[ct] > 6)) {
            sample_ct_freq[ct] <- round(sample_n * ct_freq[ct] / freq)
          }
        }
        tmp <- data.frame(t(as.matrix(train_sets[[i]])), stringsAsFactors = F, check.names = F)
        tmp$class <- celltypes[[i]]
        tmp_N <- subset(tmp, !(tmp$class %in% names(sample_ct_freq)))
        celltype_N <- subset(celltypes[[i]], !(celltypes[[i]] %in% names(sample_ct_freq)))
        tmp_smote <- tmp_N
        celltype_smote <- celltype_N
        ### smote for rare cell-types
        for (t_name in names(sample_ct_freq)) {
          tmp_target <- subset(tmp, tmp$class == t_name)
          res <- se_SMOTE(tmp_target[, -ncol(tmp_target)], target = tmp_target$class, target_name = t_name, K = 5, add_size = (sample_ct_freq[[t_name]]-ct_freq[t_name]))
          tmp_smote <- rbind(tmp_smote, res$data)
          celltype_smote <- c(celltype_smote, rep(t_name, nrow(res$data)))
        }
        train_sets[[i]] <- t(tmp_smote[, -ncol(tmp_smote)])
        celltypes[[i]] <- celltype_smote
        platforms[[i]] <- rep(unique(platforms[[i]]), ncol(train_sets[[i]]))
      }
    }
    message('SMOTE Done')
  }else{
    for (i in 1:length(seuratlist)) {
      train_sets[[i]] <- data.frame(as.matrix(train_sets[[i]]), stringsAsFactors = F, check.names = F)
    }
  }
  platforms <- unlist(platforms, use.names = F)
  celltypes <- unlist(celltypes, use.names = F)
  for (i in 1:length(seuratlist)) {
    train_sets[[i]] <- train_sets[[i]][match(genes, rownames(train_sets[[i]])), ]
    train_sets[[i]][is.na(train_sets[[i]])] <- 0
    train_sets[[i]] <- train_sets[[i]] / colSums(train_sets[[i]]) * 10000 # or matrix
    train_sets[[i]] <- log2(train_sets[[i]] + 1)
  }
  train_sets <- do.call("cbind", train_sets)
  return(list(train_sets = train_sets, celltypes = celltypes, platforms = platforms, genes = genes))
}

trainDatasets <- torch::dataset(
  name = "trainDatasets",
  initialize = function(data, celltypes, platforms, ct_dic, plat_dic) {
    class_labels <- unlist(sapply(celltypes, function(x) ct_dic[x]), use.names = F)
    domain_labels <- unlist(sapply(platforms, function(x) plat_dic[x]), use.names = F)

    self$class_labels <- torch::torch_tensor(class_labels)
    self$domain_labels <- torch::torch_tensor(domain_labels)
    self$expr <- torch::torch_tensor(data)
  },
  .getitem = function(index) {
    return(list(
      torch::torch_tensor(self$expr[, index]),
      self$class_labels[index],
      self$domain_labels[index]
    ))
  },
  .length = function() {
    return(length(self$class_labels))
  }
)


GRL <- torch::autograd_function(
  forward = function(ctx, x, alpha) {
    ctx$save_for_backward(alpha = alpha)
    x
  },
  backward = function(ctx, grad_output) {
    list(
      x = grad_output * ctx$saved_variables$alpha
    )
  }
)

MADA <- torch::nn_module(
  "class_MADA",
  initialize = function(nfeatures, nct, nplat) {
    self$nct <- nct
    self$feature <- torch::nn_sequential(torch::nn_linear(nfeatures, 100), torch::nn_relu(), torch::nn_dropout(p = 0.5, inplace = FALSE))
    self$class_classifier <- torch::nn_sequential(torch::nn_linear(100, 50), torch::nn_relu(), torch::nn_dropout(p = 0.5, inplace = FALSE), torch::nn_linear(50, nct))
    self$domain_classifier <- torch::nn_module_list(lapply(1:nct, function(x) {
      torch::nn_sequential(torch::nn_linear(100, 25), torch::nn_relu(), torch::nn_linear(25, nplat))
    }))
  },
  forward = function(input_data, alpha, nct) {
    features <- self$feature(input_data)
    class_logits <- self$class_classifier(features)
    class_predictions <- torch::nn_softmax(dim = 2)
    reverse_features <- GRL(features, alpha)
    domain_logits <- list()
    for (class_idx in 1:nct) {
      wrf <- torch::torch_unsqueeze(class_predictions(class_logits)[, class_idx], 2) * reverse_features
      domain_logits[[length(domain_logits) + 1]] <- self$domain_classifier[[class_idx]](wrf)
    }
    return(list(class_logits = class_logits, domain_logits = domain_logits))
  }
)

train <- function(train_data, params, celltypes, platforms, nfeatures, nct, nplat,
                  ct_dic, plat_dic, device) {
  network <- MADA(nfeatures, nct, nplat)$train()
  lr <- params[1]
  n_epoch <- params[2]
  batch_size <- params[3]
  train_data <- trainDatasets(data = train_data, celltypes, platforms, ct_dic, plat_dic)
  optimizer <- torch::optim_adam(network$parameters, lr = lr)
  loss_class <- torch::nn_cross_entropy_loss()
  loss_domain <- torch::nn_cross_entropy_loss()
  network <- network$to(device = device)
  loss_class <- loss_class$to(device = device)
  loss_domain <- loss_domain$to(device = device)
  train_loader <- torch::dataloader(dataset = train_data, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
  len_train_loader <- length(train_loader)
  message("Begin training")
  for (epoch in 1:n_epoch) {
    i <- 1
    coro::loop(for (l_n in train_loader) {
      p <- as.double(i + epoch * len_train_loader) / n_epoch / len_train_loader
      i <- i + 1
      alpha <- 2. / (1. + exp(-10 * p)) - 1
      expr <- l_n[[1]]
      class_label <- l_n[[2]]
      domain_label <- l_n[[3]]
      expr <- expr$to(device = device)
      expr <- expr$to(dtype = torch_float())

      class_label <- class_label$to(device = device)
      domain_label <- domain_label$to(device = device)
      res_net <- network(
        input_data = expr,
        alpha = alpha,
        nct = nct
      )
      class_output <- res_net$class_logits
      domain_output <- res_net$domain_logits
      err_class <- loss_class(class_output, class_label)
      err_domain <- sapply(1:nct, function(class_idx) {
        loss_domain(domain_output[[class_idx]], domain_label)
      })
      loss_total <- (1 - alpha) * Reduce("+", err_domain) / nct + alpha * err_class

      optimizer$zero_grad()
      loss_total$backward()
      optimizer$step()
    })
    cat(sprintf("Loss at epoch %d: %3f\n", epoch, loss_total$item()))
  }
  message("Finish Training")
  return(network)
}