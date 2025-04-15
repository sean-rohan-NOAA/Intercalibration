#' @useDynLib gearcalib
.onUnload <- function (lib) {
  library.dynam.unload("gearcalib", lib)
  library.dynam.unload("gearcalib_nb", lib)
}
