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
#'
#' @examples save_model(model, path_out, prefix)
save_model <- function(model, path_out, prefix='train') {
  torch_save(model$network, paste0(path_out, '/', prefix, '_params.pt'))
  saveRDS(model$meta, paste0(path_out, '/', prefix, '_meta.rds'))
  }