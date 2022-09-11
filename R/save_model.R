#' save_model
#'
#' @param model Output of train_model
#' @param path_out Directory for saving model.
#' @param prefix Prefix of the output files. DEFAULT: train
#'
#' @return None
#' @importFrom torch torch_save
#'
#' @export
save_model <- function(model, path_out, prefix = "train") {
  torch_save(model$network, file.path(path_out, paste0(prefix, "_params.pt")))
  saveRDS(model$meta, file.path(path_out, paste0(prefix, "_meta.rds")))
}
