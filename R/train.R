#' Train model based on your own data
#'
#' @param path_in Folder includes expression and meta files.
#' @param disease Depend on your data is in some disease condition or not, choose from 'TRUE' or 'FALSE'. DEFAULT: 'FALSE'.
#'
#' @return An MADA model.
#' @importFrom stats runif
#' @importFrom reticulate source_python
#' @importFrom utils read.csv write.table
#' @importFrom grDevices dev.off png
#' @importFrom utils read.table
#' @importFrom torch nn_sequential cuda_is_available torch_device torch_tensor torch_unsqueeze dataset nn_module dataloader nn_cross_entropy_loss nn_mse_loss optim_adam torch_float torch_load torch_save autograd_function with_no_grad
#' @export
#'
train_model <- function(path_in, disease) {
  params_train <- c(0.0001, 50, 128)
  device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
  reticulate::source_python(system.file("python", "selina.py", package = "SELINA"))
  if (disease) {
    res_pre <- preprocessing(path_in, disease)
    train_data <- as.matrix(res_pre[[1]])
    celltypes <- res_pre[[2]]
    diseases <- res_pre[[3]]
    platforms <- res_pre[[4]]
    genes <- res_pre[[5]]
    ct_dic <- label2dic(celltypes)
    disease_dic <- label2dic(diseases)
    plat_dic <- label2dic(platforms)
    nfeatures <- nrow(train_data)
    nct <- length(ct_dic)
    ndis <- length(disease_dic)
    nplat <- length(plat_dic)

    network <- train(
      train_data, params_train, celltypes, diseases, platforms, nfeatures,
      nct, ndis, nplat, ct_dic, disease_dic, plat_dic, device, disease
    )
    model_meta <- list(genes = genes, cellsources = disease_dic, celltypes = ct_dic)

    message("All done")
    return(list(network = network, meta = model_meta))
  } else {
    res_pre <- preprocessing(path_in, disease)
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
      train_data, params_train, celltypes, NULL, platforms, nfeatures,
      nct, NULL, nplat, ct_dic, NULL, plat_dic, device, disease
    )
    model_meta <- list(genes = genes, celltypes = ct_dic)
    message("All done")
    return(list(network = network, meta = model_meta))
  }
}
label2dic <- function(label) {
  label_set <- unique(label)
  label_list <- list()
  for (i in 1:length(label_set)) {
    label_list[label_set[i]] <- i
  }
  return(label_list)
}

delete_multiple_element <- function(list_object, indices) {
  indices <- sort(indices, decreasing = TRUE)
  for (idx in indices) {
    if (idx < length(list_object)) {
      list_object[[idx]] <- NULL
    }
  }
  return(list_object)
}

trainDatasets <- torch::dataset(
  name = "trainDatasets",
  initialize = function(data, celltypes, diseases, platforms, ct_dic, disease_dic, plat_dic, disease) {
    class_labels <- unlist(sapply(celltypes, function(x) ct_dic[x]), use.names = F)
    domain_labels <- unlist(sapply(platforms, function(x) plat_dic[x]), use.names = F)
    self$class_labels <- torch::torch_tensor(class_labels)
    self$domain_labels <- torch::torch_tensor(domain_labels)
    if (disease) {
      self$disease_flag <- TRUE
      disease_labels <- unlist(sapply(diseases, function(x) disease_dic[x]), use.names = F)
      self$disease_labels <- torch::torch_tensor(disease_labels)
    } else {
      self$disease_flag <- FALSE
    }
    self$expr <- torch::torch_tensor(data)
  },
  .getitem = function(index) {
    if (self$disease_flag) {
      return(list(
        torch::torch_tensor(self$expr[, index]),
        self$class_labels[index],
        self$disease_labels[index],
        self$domain_labels[index]
      ))
    } else {
      return(list(
        torch::torch_tensor(self$expr[, index]),
        self$class_labels[index],
        self$domain_labels[index]
      ))
    }
  },
  .length = function() {
    return(length(self$class_labels))
  }
)

GRL <- torch::autograd_function(
  forward = function(ctx, x, alpha) {
    ctx$save_for_backward(alpha = alpha)
    return(x$view_as(x))
  },
  backward = function(ctx, grad_output) {
    output <- grad_output$neg() * ctx$saved_variables$alpha
    return(list(output, NULL))
  }
)

MADA <- torch::nn_module(
  "MADA",
  initialize = function(nfeatures, nct, ndis, nplat, disease) {
    if (disease) {
      self$disease_flag <- TRUE
      self$nct <- nct
      self$feature <- torch::nn_sequential(torch::nn_linear(nfeatures, 500), torch::nn_relu(), torch::nn_dropout(p = 0.5, inplace = FALSE))
      self$disease_classifier <- torch::nn_sequential(torch::nn_linear(500, 50), torch::nn_relu(), torch::nn_linear(50, ndis))
      self$class_classifier <- torch::nn_sequential(torch::nn_linear(500, 50), torch::nn_relu(), torch::nn_dropout(p = 0.5, inplace = FALSE), torch::nn_linear(50, nct))
      self$domain_classifier <- torch::nn_module_list(lapply(1:(nct + ndis), function(x) {
        torch::nn_sequential(torch::nn_linear(500, 25), torch::nn_relu(), torch::nn_linear(25, nplat))
      }))
    } else {
      self$disease_flag <- FALSE
      self$nct <- nct
      self$feature <- torch::nn_sequential(torch::nn_linear(nfeatures, 100), torch::nn_relu(), torch::nn_dropout(p = 0.5, inplace = FALSE))
      self$class_classifier <- torch::nn_sequential(torch::nn_linear(100, 50), torch::nn_relu(), torch::nn_dropout(p = 0.5, inplace = FALSE), torch::nn_linear(50, nct))
      self$domain_classifier <- torch::nn_module_list(lapply(1:nct, function(x) {
        torch::nn_sequential(torch::nn_linear(100, 25), torch::nn_relu(), torch::nn_linear(25, nplat))
      }))
    }
  },
  forward = function(input_data, alpha, nct, ndis) {
    features <- self$feature(input_data)
    class_logits <- self$class_classifier(features)
    class_predictions <- torch::nn_softmax(dim = 2)
    reverse_features <- GRL(features, alpha)
    domain_logits <- list()
    for (class_idx in 1:nct) {
      wrf <- torch::torch_unsqueeze(class_predictions(class_logits)[, class_idx], 2) * reverse_features
      domain_logits[[length(domain_logits) + 1]] <- self$domain_classifier[[class_idx]](wrf)
    }
    if (self$disease_flag) {
      disease_logits <- self$disease_classifier(features)
      disease_predictions <- torch::nn_softmax(dim = 2)
      for (dis_idx in 1:ndis) {
        wrf <- torch::torch_unsqueeze(disease_predictions(disease_logits)[, dis_idx], 2) * reverse_features
        domain_logits[[length(domain_logits) + 1]] <- self$domain_classifier[[dis_idx]](wrf)
      }
      return(list(class_logits = class_logits, disease_logits = disease_logits, domain_logits = domain_logits))
    } else {
      return(list(class_logits = class_logits, domain_logits = domain_logits))
    }
  }
)
# train(
#   train_data, params_train, celltypes, diseases, platforms, nfeatures,
#   nct, ndis, nplat, ct_dic, disease_dic, plat_dic, device, disease
# )
train <- function(train_data, params, celltypes, diseases, platforms, nfeatures, nct, ndis, nplat,
                  ct_dic, disease_dic, plat_dic, device, disease) {
  if (disease) {
    network <- MADA(nfeatures, nct, ndis, nplat, disease)$train()
    train_data <- trainDatasets(train_data, celltypes, diseases, platforms, ct_dic, disease_dic, plat_dic, disease)
  } else {
    network <- MADA(nfeatures, nct, ndis, nplat, disease)$train()
    train_data <- trainDatasets(train_data, celltypes, NULL, platforms, ct_dic, NULL, plat_dic, disease)
  }
  lr <- params[1]
  n_epoch <- params[2]
  batch_size <- params[3]
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
      if (disease) {
        expr <- l_n[[1]]
        class_label <- l_n[[2]]
        disease_label <- l_n[[3]]
        domain_label <- l_n[[4]]
        expr <- expr$to(device = device)
        expr <- expr$to(dtype = torch_float())
        class_label <- class_label$to(device = device)
        disease_label <- disease_label$to(device = device)
        domain_label <- domain_label$to(device = device)
        res_net <- network(
          input_data = expr,
          alpha = alpha,
          nct = nct,
          ndis = ndis
        )
        class_output <- res_net$class_logits
        disease_output <- res_net$disease_logits
        domain_output <- res_net$domain_logits
        err_class <- loss_class(class_output, class_label)
        err_disease <- loss_class(disease_output, disease_label)
        err_domain_ct <- sapply(1:nct, function(class_idx) {
          loss_domain(domain_output[[class_idx]], domain_label)
        })
        err_domain_dis <- sapply((nct + 1):(nct + ndis), function(dis_idx) {
          loss_domain(domain_output[[dis_idx]], domain_label)
        })
        err_domain <- c(err_domain_ct, err_domain_dis)
        loss_total <- (1 - alpha) / 2 * (1 - alpha) * (Reduce("+", err_domain_dis) / ndis) + (alpha * 2) * (1 - alpha) * (Reduce("+", err_domain_ct) / nct) + alpha * (err_class + err_disease)
      } else {
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
          nct = nct,
          ndis = NULL
        )
        class_output <- res_net$class_logits
        domain_output <- res_net$domain_logits
        err_class <- loss_class(class_output, class_label)
        err_domain <- sapply(1:nct, function(class_idx) {
          loss_domain(domain_output[[class_idx]], domain_label)
        })
        loss_total <- (1 - alpha) * Reduce("+", err_domain) / nct + alpha * err_class
      }
      optimizer$zero_grad()
      loss_total$backward()
      optimizer$step()
    })
    cat(sprintf("Loss at epoch %d: %3f\n", epoch, loss_total$item()))
  }
  message("Finish training")
  return(network)
}
