#' @useDynLib gearcalib
.onUnload <- function (lib) {
  library.dynam.unload("gearcalib", lib)
}
