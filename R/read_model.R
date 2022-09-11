#' Read in model trained by yourself
#'
#' @param path_model Path for model file suffixed with '.pt'.
#'
#' @return A trained MADA model.
#' @importFrom torch torch_load cuda_is_available
#' @export
#'
read_model <- function(path_model) {
  device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
  network <- torch_load(path_model, device = device)
  meta <- readRDS(gsub("_params.pt", "_meta.rds", path_model))
  return(list(network = network, meta = meta))
}
