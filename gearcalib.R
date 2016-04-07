library(TMB)
compile("gearcalib.cpp")
if(is.loaded("EvalDoubleFunObject")) dyn.unload("gearcalib.so")
dyn.load("gearcalib.so")

## Specify range of length classes to include in the model
Lmin <- 10
Lmax <- 80

Lvec <- seq(Lmin,Lmax)

Do.Intercal <- TRUE
Do.Verification <- FALSE

## Empirical estimate with bootstrap confidence regions
boot <- function(d,quantiles = c(0.025,0.16,0.5,0.84,0.975))
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

        Nboot <- 1000

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
                    BootQuantiles=BootQuantiles))

    }


ChooseGear <- function(d,GearNames)
    {
        ## Construct table of which gears are used in which groups
        tt <- table(d$group,d$Gear)

        ## Identify those groups where the two gears are used
        UseGroups <- tt[,GearNames[1]] & tt[,GearNames[2]] 

        ## Identify those hauls which are part of these groups
        I <- sapply(d$group,function(g)UseGroups[g])

        ## Use only if SweptArea is finite
        I <- I & is.finite(d$SweptArea)
        
        ## Select only these records
        d$SweptArea <- d$SweptArea[I]
        d$N <- d$N[I,Lmin:Lmax]
        d$group <- factor(d$group[I],levels=unique(d$group[I]))
        d$Gear <- factor(d$Gear[I],levels=GearNames)

        return(d)
    }


fitmodel <- function(Species,GearNames,d,fit0=FALSE)
    {
        d <- ChooseGear(d,GearNames)
        
        nsize <- ncol(d$N)
        ngroup <- nlevels(d$group)
        nhaul <- nrow(d$N)
        ngear <- nlevels(d$Gear)

        ## Set default RW order
        rw_order <- c(1,1)

        obj <- MakeADFun(
            data=list(
                N=d$N,
                SweptArea=d$SweptArea,
                group=d$group,
                Gear=d$Gear,
                huge=10,
                tiny=0.01,
                rw_order=rw_order
            ),
            parameters=list(
                logspectrum=matrix(0,ngroup,nsize),
                nugget=matrix(0,nhaul,nsize),
                residual=matrix(0,nhaul,nsize),
                loggear=numeric(nsize),
                logsd=-1,
                phi=.5,
                logsdnug=-1,
                logsdres=-1,
                logsdGearRW=-1,
                alpha = 1
            ),
            map = list(alpha=factor(NA)), # phi=factor(NA),logsdnug=factor(NA),logsdGearRW=factor(NA)),
            random=c("logspectrum","residual","loggear","nugget"),
            DLL="gearcalib"
        )


        lower <- 0*obj$par-Inf
        upper <- 0*obj$par+Inf
        lower["phi"] <- 0
        upper["phi"] <- 0.9999

        obj$env$tracepar <- TRUE

        system.time( opt <- nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper) )
        Pvalue <- NA

        if(fit0)
            {
                obj0 <- MakeADFun(
                    data=list(
                        N=d$N,
                        SweptArea=d$SweptArea,
                        group=d$group,
                        Gear=d$Gear,
                        huge=10,
                        tiny=0.01,
                        rw_order=rw_order
                        ),
                    parameters=list(
                        logspectrum=matrix(0,ngroup,nsize),
                        nugget=matrix(0,nhaul,nsize),
                        residual=matrix(0,nhaul,nsize),
                        loggear=numeric(nsize),
                        logsd=-1,
                        phi=.5,
                        logsdnug=-1,
                        logsdres=-1,
                        logsdGearRW=-10,
                        alpha = 1
                        ),
                    ## map = list(phi=factor(NA),logsdnug=factor(NA),logsdGearRW=factor(NA)),    
                    map = list(alpha=factor(NA),logsdGearRW=factor(NA)),    

                    random=c("logspectrum","residual","loggear","nugget"),
                    DLL="gearcalib"
                    )

                obj0$env$tracepar <- TRUE
                
                system.time( opt0 <- nlminb(obj0$par,obj0$fn,obj0$gr,lower=lower,upper=upper) )

                Pvalue <- pchisq(2*(opt0$objective-opt$objective),df=1)

                print(paste("Chisq-test of no size structure in gear effect:",Pvalue))

                rm(obj0)
            }

        save.image(file=paste("Image-",Species,"-",GearNames[2],"-vs-",GearNames[1],".Rdata",sep=""))

        rep <- sdreport(obj)

        repsum <- summary(rep)
        s <- repsum[grep("gear",rownames(repsum)),]
        est <- 2*s[,1]
        sd <- 2*s[,2]

        emp <- boot(d)

        b <- boot(d)
        
        ## Plots
        X11()

        plot(c(Lmin,Lmax),c(0,max(exp(est + sd %o% c(0,-2,2)))),type="n",
             ylim=c(0,2),
             xlim=c(Lmin,Lmax),
             xlab="Length [cm]",
             ylab=paste(GearNames[2]," vs. ",GearNames[1]),
             main=paste("Relative selectivity for", Species))

        lines(c(Lmin,Lmax),rep(1,2),col="darkgrey",lwd=3,lty="dashed")

        polygon(c(Lmin:Lmax,Lmax:Lmin),c(exp(est-2*sd),rev(exp(est+2*sd))),
                col="grey",border=NA)
        lines(Lmin:Lmax,exp(est),lwd=3)

        points(Lmin:Lmax,b$RawEstimate)
        apply(b$BootQuantiles,1,function(x)lines(Lmin:Lmax,x))

        grid()

        dev.copy2pdf(file=paste("Relsel-",Species,"-",GearNames[2],"-vs-",GearNames[1],".pdf",sep=""))

        ## plot(Lmin:Lmax,log10(1+Density1),
        ##      ylim=log10(1+range(c(Density1,Density2))),
        ##      xlim=c(Lmin,Lmax),
        ##      type="l",lty="dashed",
        ##      xlab="Length [cm]",ylab="Density (log10(N/A+1))")
        ## points(Lmin:Lmax,log10(1+Density1),pch="o")
        ## lines(Lmin:Lmax,log10(1+Density2))
        ## points(Lmin:Lmax,log10(1+Density2),pch="+")
        ## legend("topright",legend=GearNames,lty=c("dashed","solid"),pch=c("o","+"))

        ## grid()

        ## dev.copy2pdf(file=paste("Density-",Species,"-",GearNames[2],"-vs-",GearNames[1],".pdf",sep=""))

        ## #############################################################################################
        ## Plot example of a size spectrum
        pl <- obj$env$parList()
        spec <- pl$logspectrum

        ## Take the station witht the largest catch
        ii <- which.max(apply(spec,1,sum))

        iii <- which(d$group == levels(d$group)[ii])

        nug1 <- pl$residual[iii,]

        plot(Lmin:Lmax,exp(spec[ii,])*d$SweptArea[ii],log="y",type="l",ylim=c(1e-1,max(d$N[iii,])),
             xlim=c(Lmin,Lmax),
             xlab="Length [cm]",ylab="Size spectrum of population [#/cm]",lwd=3,
             main=paste(Species,"at station",d$group[ii]))
        grid()

        lines(Lmin:Lmax,exp(spec[ii,] + nug1[1,])*d$SweptArea[ii],lty="solid")
        lines(Lmin:Lmax,exp(spec[ii,] + nug1[2,])*d$SweptArea[ii],lty="dashed")

        points(Lmin:Lmax,d$N[iii[1],],pch="+")
        points(Lmin:Lmax,d$N[iii[2],],pch="o")

        legend("bottom",legend=c("Size structure at station",
                            paste(GearNames[1],", observed and modelled"),
                            paste(GearNames[2],", observed and modelled")),
               lwd=c(3,1,1),lty=c("solid","solid","dashed"),pch=c(NA,"+","o"))

        grid()

        dev.copy2pdf(file=paste("Sample-",Species,"-",GearNames[2],"-vs-",GearNames[1],".pdf",sep=""))

        res <- data.frame(L=Lmin:Lmax,est=est,sd=sd)

        write.table(res,file=paste("Results-",Species,"-",GearNames[2],"-vs-",GearNames[1],".csv",sep=""))
        write.table(res,file=paste("Estimates-",Species,"-",GearNames[2],"-vs-",GearNames[1],".csv",sep=""))

        return(list(Species=Species,GearNames=GearNames,d=d,res=res,Pvalue=Pvalue,rep=rep,opt=opt,obj=obj))
    }

fits <- list()

if(Do.Intercal)
    {
        Dir <- "Repaired-Trawl-Pairs"
        filename <- paste(Dir,"InterCal_Merluccius.RData",sep="/")
        load(filename)

        for(Species in c("Paradoxus"))
            for(gear in c("Afr_Old"))
                {
                    d <- eval(parse(text=paste("listIntercalData_",Species,sep="")))

                    d$group <- factor(d$group)
                    
                    ## Change order of gear such that Gisund is reference
                    d$Gear <- factor(d$Gear,levels=rev(unique(d$Gear)))

                    d$SweptArea <- d$SweptArea[,1]


                    ## Compare these two gear types
                    GearNames <- c("Gisund",gear)

                    fits[[length(fits)+1]] <- fitmodel(Species,GearNames,d)
                }
    }

if(Do.Verification)
    {
        
        Dir <- "Verification"
        filename <- paste(Dir,"InterCal_Merluccius_1998_1999.RData",sep="/")
        load(filename)


        for(Species in c("Capensis","Paradoxus"))
                {
                    d <- eval(parse(text=paste("listIntercalData_",Species,sep="")))

                    ## Pick only data from the 1990's
                    I <- grep("199",d$group)
                    d$Gear <- d$Gear[I]
                    d$N <- d$N[I,]
                    d$group <- d$group[I]
                    d$SweptArea <- d$SweptArea[I,]
                    
                    d$Gear[d$Gear == "Blue Sea"] <- "Blue_Sea"
                    d$group <- factor(d$group)

                    ## Change order of gear such that Gisund is reference
                    d$Gear <- factor(d$Gear,levels=rev(unique(d$Gear)))

                    d$SweptArea <- d$SweptArea[,1]


                    ## Compare these two gear types
                    GearNames <- c("Nansen","Blue_Sea")

                    fits[[length(fits)+1]] <- fitmodel(Species,GearNames,d)
                }
    }

save(fits,file="fits.Rdata")

casetab <- lapply(fits,function(f)return(
    list(Species=f$Species,
      Gear1=f$GearNames[1],
      Gear2=f$GearNames[2])))

casetab <- do.call(rbind,casetab)
partab <- t(sapply(fits,function(f)f$opt$par))
sdtab <- t(sapply(fits,function(f)sqrt(diag(f$rep$cov.fixed))))

colnames(sdtab) <- paste(colnames(sdtab),".sd",sep="")
partab <- cbind(casetab,partab,sdtab)

