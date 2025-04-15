##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Fit relative catch effiency by length group between to gears
##' @param d A list with four elements: N (matrix of integers), SweptArea(numeric vector), group(factor vector), and Gear(factor vector). 
##' @param fit0 If TRUE a Chisq-test of no size structure in gear effect is performed.
##' @param linearEffort If TRUE, catch is assumed proportional to swept area, otherwise proportional to a power function of swept area.
##' @return A list
gearcalibFitNB <- function(d,fit0=FALSE,logsd=NA,phi=NA,logsdnug=NA,logsdres=NA,logsdGearRW=NA,logalpha=0,logtheta=NA)
{
  check.gearcalib.data(d)
  nsize <- ncol(d$N)
  ngroup <- nlevels(d$group)
  nhaul <- nrow(d$N)
  ngear <- nlevels(d$Gear)
  
  ## Set default RW order
  rw_order <- c(1,1)
  
  data <- list(
    N=d$N,
    SweptArea=d$SweptArea,
    group=d$group,
    Gear=d$Gear,
    huge=10,
    tiny=0.01,
    rw_order=rw_order
  )
  
  parameters=list(
    logspectrum=matrix(0,ngroup,nsize),
    nugget=matrix(0,nhaul,nsize),
    residual=matrix(0,nhaul,nsize),
    loggear=numeric(nsize),
    logsd=-1,
    phi=0.9,
    logsdnug=-1,
    logsdres=-1,
    logsdGearRW=-1,
    logalpha = 0,
    logtheta = -1
  )
  
  random <- c("logspectrum","residual","loggear","nugget")
  
  map <- list()
  
  ## If any parameters have been specified in the call, do not estimate those
  ## parameters but fix them to the specified value
  setparameter <- function(name) {
    var <- get(name)
    if(!is.na(var))
    {
      map[[name]] <<- factor(NA)
      parameters[[name]] <<- var
    }
  }
  parameternames <- c("logsd","phi","logsdnug","logsdres","logsdGearRW","logalpha","logtheta")
  sapply(parameternames,setparameter)
  
  
  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    map = map, 
    random = random,
    DLL="gearcalib_nb"
  )
  
  
  lower <- 0*obj$par-Inf
  upper <- 0*obj$par+Inf
  if(any("phi" == names(obj$par)))
  {
    lower["phi"] <- 0
    upper["phi"] <- 0.99
  }
  
  
  obj$env$tracepar <- TRUE
  
  system.time( opt <- nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper) )
  Pvalue <- NA
  
  if(fit0)
  {
    parameters0 <- parameters
    parameters0$logsdGearRW <- -10
    
    obj0 <- MakeADFun(
      data = data,
      parameters = parameters0,
      map = c(map,list(logsdGearRW=factor(NA))),
      random = random,
      DLL="gearcalib_nb"
    )
    
    obj0$env$tracepar <- TRUE
    
    system.time( opt0 <- nlminb(obj0$par,obj0$fn,obj0$gr,lower=lower,upper=upper) )
    
    Pvalue <- pchisq(2*(opt0$objective-opt$objective),df=1)
    
    print(paste("Chisq-test of no size structure in gear effect:",Pvalue))
    
    rm(obj0)
  }
  
  rep <- sdreport(obj)
  repsum <- summary(rep,"random")
  s <- repsum[grep("gear",rownames(repsum)),]
  est <- 2*s[,1]
  sd <- 2*s[,2]
  
  ret <- list(d=d,
              Pvalue=Pvalue,rep=rep,opt=opt,obj=obj,est=est,sd=sd)
  
  class(ret)<-"gearcalibFit" 
  
  return(ret)
  
}
