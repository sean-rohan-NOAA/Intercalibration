source("gearcalib.R") ## Get fitted model

obj <- fits[[1]]$obj
getSample <- function(obj, as.list=TRUE){
    tmp <- obj$env$MC(n=1, keep=TRUE, antithetic=FALSE)
    samp <- attr(tmp,"samples")
    if(as.list){
        par <- obj$env$last.par.best
        par[obj$env$random] <- samp
        samp <- obj$env$parList(par=par)
    }
    samp
}

pl <- getSample(obj)
