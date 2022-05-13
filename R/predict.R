#' Fine-tuning and predict for the query data
#'
#' @param query_expr File path of the query data matrix
#' @param model File path of the pre-trained model
#' @param path_out File path of the output files.
#' @param outprefix Prefix of the output files. DEFAULT: query.
#' @param disease Depend on your data is in some disease condition or not, choose from 'TRUE' or 'FALSE'. DEFAULT: 'FALSE'.
#' @param mode Single-cell level input or cluster level input, choose from 'single' or 'cluster'. DEFAULT: 'single'.
#'
#' @return Prediction results and corresponding probability for each query cell.
#' 
#' @importFrom stats rnorm 
#' @importFrom utils read.csv write.table
#' @importFrom torch cuda_is_available torch_device torch_tensor torch_unsqueeze dataset nn_module dataloader nn_cross_entropy_loss nn_mse_loss optim_adam torch_float torch_load torch_save autograd_function with_no_grad
#' @export
#'
#' @examples
query_predict <- function(query_expr, model, path_out, outprefix, disease, mode) {
  print("Loading data")
  query_expr <- read_expr(query_expr)
  meta <- readRDS(gsub("pt", "rds", (gsub("params", "meta", model))))
  genes <- meta$genes
  ct_dic <- meta$celltypes
  query_expr <- merge_query(genes, query_expr)
  nfeatures <- length(genes)
  nct <- length(ct_dic)
  network <- torch_load(model, device = device)
  network <- Autoencoder(network, nfeatures, nct)
  if ((!disease) & (mode != "cluster")) {
    print("Fine-tuning1")
    network <- tune1(query_expr, network, params_tune1)
    print("Fine-tuning2")
    network <- tune2(query_expr, network, params_tune2)
  }
  network <- Classifier(network)$to(device = device)
  test_res <- test(query_expr, network, ct_dic)
  pred_labels <- test_res$pred_labels
  pred_prob <- test_res$pred_prob
  
  ## make directory path_out
  write.table(data.frame(pred_labels, stringsAsFactors = F),
              paste0(path_out, "/", outprefix, "_predictions.txt"),
              row.names = F, sep = "\t", quote = F
  )
  write.table(data.frame(pred_prob, stringsAsFactors = F),
              paste0(path_out, "/", outprefix, "_probability.txt"),
              sep = "\t", quote = F
  )
  print("Finish Prediction")
}

device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
params_tune1 <- c(0.0005, 50, 128)
params_tune2 <- c(0.0001, 10, 128)

read_expr <- function(path) {
  expr <- data.table::fread(path, sep = "\t")
  expr <- data.frame(expr, stringsAsFactors = F, check.names = F)
  rownames(expr) <- expr$Gene
  expr <- expr[, -match("Gene", colnames(expr))]
  return(expr)
}

merge_query <- function(genes, query_expr) {
  model_expr <- data.frame(genes = rnorm(length(genes)), stringsAsFactors = F)
  rownames(model_expr) <- genes
  query_expr <- merge(model_expr, query_expr, by = "row.names")
  rownames(query_expr) <- query_expr$Row.names
  query_expr <- query_expr[, -match("genes", colnames(query_expr))]
  query_expr <- query_expr[, -1]
  query_expr[is.na(query_expr)] <- 0
  query_expr <- query_expr / colSums(query_expr) * 10000 # or matrix
  query_expr <- log2(query_expr + 1)
  return(query_expr)
}


Datasets <- dataset(
  initialize = function(data) {
    self$expr <- data
  },
  .getitem = function(index) {
    torch_tensor(self$expr[, index])
  },
  .length = function() {
    ncol(self$expr)
  }
)

Autoencoder <- nn_module(
  "class_Autoencoder",
  initialize = function(network, nfeature, nct) {
    encoder <- c(network$feature$children, network$class_classifier$children)
    encoder_index <- c(1, 2, 4, 5, 7)
    self$encoder <- do.call(nn_sequential, c(unlist(sapply(encoder_index, function(x) encoder[[x]]), use.names = F), nn_relu()))
    # self$encoder = nn_sequential(encoder[[1]],encoder[[2]],encoder[[4]],encoder[[5]],encoder[[7]],nn_relu() ) #
    self$decoder <- nn_sequential(
      nn_linear(in_features = nct, out_features = 50), nn_relu(),
      nn_linear(in_features = 50, out_features = 100), nn_relu(),
      nn_linear(in_features = 100, out_features = nfeature)
    )
  },
  forward = function(input_data) {
    output <- self$decoder(self$encoder(input_data))
    return(output)
  }
)


Classifier <- nn_module(
  "class_Classifier",
  initialize = function(network) {
    self$classifier <- do.call(nn_sequential, c(unlist(network$encoder$children[-length(network$encoder$children)], use.names = F), nn_softmax(dim = 1)))
  },
  forward = function(input_data) {
    output <- self$classifier(input_data)
    return(output)
  }
)
# test_df=query_expr
# params=params_tune1
tune1 <- function(test_df, network, params) {
  test_dat <- Datasets(test_df)
  lr <- params[1]
  n_epoch <- params[2]
  batch_size <- params[3]
  optimizer <- optim_adam(network$parameters, lr = lr)
  loss <- nn_mse_loss()
  loss <- loss$to(device = device)
  test_loader <- dataloader(dataset = test_dat, batch_size = batch_size, shuffle = TRUE)
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
  print("Finish Tuning1")
  return(network)
}


tune2 <- function(test_df, network, params) {
  test_dat <- Datasets(test_df)
  lr <- params[1]
  n_epoch <- params[2]
  batch_size <- params[3]
  optimizer <- optim_adam(network$parameters, lr = lr)
  loss <- nn_mse_loss()
  loss <- loss$to(device = device)
  test_loader <- dataloader(dataset = test_dat, batch_size = batch_size, shuffle = TRUE)
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
  print("Finish Tuning2")
  return(network)
}


test <- function(test_df, network, ct_dic) {
  test_dat <- Datasets(test_df)
  ct_dic_rev <- split(rep(names(ct_dic), sapply(ct_dic, length)), unlist(ct_dic))
  test_loader <- dataloader(dataset = test_dat, batch_size = ncol(test_df), shuffle = FALSE)
  i <- 1 # why i not increase
  with_no_grad(
    coro::loop(for (batch in test_loader) {
      i <- i + 1
      expr <- batch
      expr <- expr$to(dtype = torch_float())
      expr <- expr$to(device = device)
      class_output <- network(expr)
      pred_labels <- as.numeric(class_output$argmax(dim = 2))
      pred_prob <- class_output
    })
  )
  pred_labels <- sapply(pred_labels, function(x) ct_dic_rev[[x]])
  pred_prob <- as.matrix(pred_prob)
  rownames(pred_prob) <- colnames(test_df)
  colnames(pred_prob) <- names(ct_dic)
  return(list(pred_labels = pred_labels, pred_prob = pred_prob))
}