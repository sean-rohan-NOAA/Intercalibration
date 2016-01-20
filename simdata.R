Dir <- "Repaired-Trawl-Pairs"
filename <- paste(Dir,"InterCal_Merluccius.RData",sep="/")
load(filename)

Species <- "Paradoxus"
gear <- "Afr_Old"

rm(d)

d <- eval(parse(text=paste("listIntercalData_",Species,sep="")))

d$group <- factor(d$group)

## Change order of gear such that Gisund is reference
d$Gear <- factor(d$Gear,levels=rev(unique(d$Gear)))

d$SweptArea <- d$SweptArea[,1]

## Construct table of which gears are used in which groups
tt <- table(d$group,d$Gear)

## Compare these two gear types
GearNames <- c("Gisund",gear)

## Identify those groups where the two gears are used
UseGroups <- tt[,GearNames[1]] & tt[,GearNames[2]] 

I <- sapply(d$group,function(g)UseGroups[g])

## Use only if SweptArea is finite
I <- I & is.finite(d$SweptArea)

## Select only these records
d$SweptArea <- d$SweptArea[I]
d$group <- factor(d$group[I],levels=unique(d$group[I]))
d$Gear <- factor(d$Gear[I],levels=GearNames)
d$N <- d$N[I,]

## Generate an "arbitrary" relative selectivity curve
relsel <- 0.5 + exp(-seq(0,1,length=ncol(d$N)))

Nsim <- array(NA,dim(d$N))

grouplevels <- levels(d$group)

Nstat <- 0

for(i in 1:length(grouplevels))
    {
        Stations <- which(d$group==grouplevels[i])
              
        Ntot <- sum(d$N[Stations,])

        PopSpec <- log(1+apply(d$N[Stations,]/d$SweptArea[Stations],2,mean))

        for(j in Stations)
            {
                ElogN <- PopSpec - 0.5*relsel
                if(d$Gear[j]==(levels(d$Gear)[1]))
                    {
                        ElogN <- ElogN + relsel
                    }

                nuggetphi <- 0.8
                nuggetsd <- 0.0
                nugget <- nuggetsd*(1-nuggetphi^2)*
                    arima.sim(model=list(ar=nuggetphi),n=ncol(d$N))

                phi <- phi + nugget

                Nsim[j,] <- rpois(n=ncol(Nsim),lambda=exp(phi)*d$SweptArea[j])
                                       
            }

        dsim <- d
        dsim$N <- Nsim

    }

graphics.off()
matplot(t(d$N[Stations,]))
X11()
matplot(t(dsim$N[Stations,]))

# simfit <- fitmodel(Species,GearNames,dsim)


# fits[[length(fits)+1]] <- fitmodel(Species,GearNames,d)


