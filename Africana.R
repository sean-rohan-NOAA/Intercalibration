for(Species in c("Capensis","Paradoxus"))
    {
        New <- read.table(paste("Results",Species,"Afr_New-vs-Gisund.csv",sep="-"))
        Old <- read.table(paste("Results",Species,"Afr_Old-vs-Gisund.csv",sep="-"))

        NewVsOld <- New$est - Old$est
        sd <- sqrt(New$sd^2 + Old$sd^2)

        xmin <- 10
        nsize <- length(NewVsOld)

        pdf(file=paste("Afr_New_Vs_Old_",Species,".pdf",sep=""))

        plot(NewVsOld,type="n",
             xlim=c(xmin,length(NewVsOld)),ylim=c(0,exp(max(NewVsOld+2*sd))),
             xlab="Length [cm]",ylab="Relative selectivity, Afr_New vs Afr_Old",
             main=Species)

        polygon(c(1:nsize,nsize:1),exp(c(NewVsOld+2*sd,rev(NewVsOld-2*sd))),
                col="grey",border=NA)

        lines(exp(NewVsOld),lwd=3)
        lines(rep(1,length(NewVsOld)))

        grid()

        dev.off()

    }
