##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Checks the a gear calib data set has the correct format
##' @param d A 'gear calib list'
##' @return nothing
check.gearcalib.data<-function(d){
    stopifnot(class(d)=="list")
    stopifnot( names(d)==c("N","SweptArea","group","Gear") )
    stopifnot(class(d$N)=="matrix")
    stopifnot(class(d$SweptArea)=="numeric")
    stopifnot(class(d$group)=="factor")
    stopifnot(class(d$Gear)=="factor")
    no.hauls <- nrow(d$N)
    stopifnot(dim(d$N)==dim(d$SweptArea))
    stopifnot(nlevels(d$group)==no.hauls/2)
    stopifnot(nlevels(d$Gear)==2)
    stopifnot(all(is.finite(d$SweptArea)))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Fit relative catch effiency by length group between to gears
##' @param d A list with four elements: N (matrix of integers), SweptArea(numeric vector), group(factor vector), and Gear(factor vector). 
##' @param fit0 If TRUE a Chisq-test of no size structure in gear effect is performed.
##' @return A list
gearcalibFit <- function(d,fit0=FALSE)
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
            phi=0.99,
            logsdnug=-1,
            logsdres=-1,
            logsdGearRW=-1,
            logalpha = 0
            )

        random <- c("logspectrum","residual","loggear","nugget")

        map <- list()
        
        obj <- MakeADFun(
            data = data,
            parameters = parameters,
            map = map, # phi=factor(NA),logsdnug=factor(NA),logsdGearRW=factor(NA)),
            random = random,
            DLL="gearcalib"
        )

        lower <- 0*obj$par-Inf
        upper <- 0*obj$par+Inf
        if(any("phi" == names(obj$par)))
            {
                lower["phi"] <- 0
                upper["phi"] <- 0.9999
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
                    DLL="gearcalib"
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
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Empirical estimate with bootstrap confidence regions
##' @param d gearcalib data
##' @param quantiles report these quantiles
##' @param Nboot number of bootstrap replicates
##' @return a list
boot <- function(d,quantiles = c(0.025,0.16,0.5,0.84,0.975),Nboot=1000)
    {
        GearNames <- levels(d$Gear)

        ## Raw data: Total catches
        I1 <- d$Gear==GearNames[1]
        I2 <- d$Gear==GearNames[2]

        totcatch1 <- apply(d$N[I1,],2,sum)
        totcatch2 <- apply(d$N[I2,],2,sum)

        Density1 <- totcatch1 / sum(d$SweptArea[I1])
        Density2 <- totcatch2 / sum(d$SweptArea[I2])

        tiny <- 1e-6
        
        RawEstimate <- (Density2 + tiny) / (Density1 + tiny)

        BootEstimate <- array(NA,c(Nboot,length(RawEstimate)))

        for(i in 1:Nboot)
            {
                ## Select random groups with replacement
                GroupBoot <- sample(1:length(levels(d$group)),length(levels(d$group)),replace=TRUE)
                
                Haul1 <- as.numeric(sapply(GroupBoot,function(g) which( (as.numeric(d$group) == g) & I1)))
                Haul2 <- as.numeric(sapply(GroupBoot,function(g) which( (as.numeric(d$group) == g) & I2)))

                Nb1 <- d$N[Haul1,]
                Nb2 <- d$N[Haul2,]

                Db1 <- apply(Nb1,2,sum) / sum(d$SweptArea[Haul1])
                Db2 <- apply(Nb2,2,sum) / sum(d$SweptArea[Haul2])

                BootEstimate[i,] <- (Db2 + tiny) / (Db1 + tiny) 
            }

        BootQuantiles <- apply(BootEstimate,2,function(x)quantile(x,quantiles))

        return(list(RawEstimate=RawEstimate,
                    BootEstimate=BootEstimate,
                    BootQuantiles=BootQuantiles,
                    Density1=Density1,
                    Density2=Density2))

    }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plots a gear calibration fit
##' @param fit an onject of class 'gearcalibFit' as fitted by gearcalib()
##' @param select a vector specifying which plots are wanted. Can be "relsel" for "Relative selectivity" or "density".
##' @param boot (optional) list with bootstrap estimates as produced by boot()
##' @return nothing
plot.gearcalibFit <- function(fit,select=c("relsel","density"),boot=NULL)
{
    
        Lvec <- 1:(ncol(fit$d$N))
        Lmin <- min(Lvec)
        Lmax <- max(Lvec)
        
        with(fit,{

            if("relsel" %in% select){
                 plot(range(Lvec),c(0,max(exp(est + sd %o% c(0,-2,2)))),type="n",
                      ylim=c(0,2),
                      xlim=range(Lvec),
                      xlab="Length group",
                      ylab=paste(levels(fit$d$Gear)[2]," vs. ",levels(fit$d$Gear)[1]),
                      main="Relative selectivity")
                 
                 lines(range(Lvec),rep(1,2),col="darkgrey",lwd=3,lty="dashed")
                 
                 polygon(c(Lvec,rev(Lvec)),c(exp(est-2*sd),rev(exp(est+2*sd))),
                         col="grey",border=NA)
                 lines(Lvec,exp(est),lwd=3)
                 if(!is.null(boot)){
                     points(Lvec,b$RawEstimate)
                     apply(b$BootQuantiles,1,function(x)lines(Lvec,x))
                 }
                 grid()
            }
                 

            if(!is.null(boot) && "density" %in% select){
                with(boot,{
                plot(Lvec,log10(1+Density1),
                      ylim=log10(1+range(c(Density1,Density2))),
                      xlim=c(Lmin,Lmax),
                      type="l",lty="dashed",
                      xlab="Length [cm]",ylab="Density (log10(N/A+1))")
                 points(Lvec,log10(1+Density1),pch="o")
                 lines(Lvec,log10(1+Density2))
                 points(Lvec,log10(1+Density2),pch="+")
                 legend("topright",legend=GearNames,lty=c("dashed","solid"),pch=c("o","+"))
                 
                 grid()
                })
            }
      })
    }
