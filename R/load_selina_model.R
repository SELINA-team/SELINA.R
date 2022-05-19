#' load_selina_model
#'
#' @param path_model Path for loading SELINA model. Please put '.pth' and 'rds' in the same folder.
#'
#' @return A trained MADA model.
#' @importFrom torch torch_load cuda_is_available load_state_dict
#' @export
#'
load_selina_model <- function(path_model) {
  meta <- readRDS(gsub("_params.pth", "_meta.rds", path_model))
  model <- MADA(meta$node[1], meta$node[2], meta$node[3])
  model$load_state_dict(torch::load_state_dict(path_model))
  model <- list(network = model, meta = meta[c(1, 2)])
}


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
