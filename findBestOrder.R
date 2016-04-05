source("gearcalib.R") ## Get fitted model
obj <- fits[[1]]$obj
opt <- fits[[1]]$opt

## RW order:
e <- expand.grid(spectrum=1, gear=1:3)
fun <- function(i){
    obj$env$data$rw_order <- as.numeric(e[i,])
    obj$retape()
    obj$env$value.best <- Inf
    start <- opt$par
    start["logsdGearRW"] <- -1
    nlminb(start,obj$fn,obj$gr)
}
li <- lapply(1:nrow(e),fun)
e$nll <- sapply(li, function(x)x$objective)
xtabs(nll ~ spectrum + gear, data=e)
