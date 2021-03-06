#' Load model trained by SELINA
#'
#' @param tissue Choose which tissue model you want to load, eg: Pancreas. Please see documentation for tissue name.
#'
#' @return A trained MADA model.
#' @importFrom torch torch_load cuda_is_available load_state_dict
#' @export
#'
load_selina_model <- function(tissue) {
  model <- read_model(system.file("extdata", paste0(tissue,"_params.pt"), package = "SELINA"))
  return(model)
}