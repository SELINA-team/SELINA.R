#' read_model
#'
#' @param path_model Path for loading model.
#'
#' @return A MADA object.
#' @importFrom torch torch_load cuda_is_available
#' @export
#'
#' @examples model <- read_model(path_model)
read_model <- function(path_model) {
  device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
  model <- torch_load(path_model, device = device)
  return(model)
}