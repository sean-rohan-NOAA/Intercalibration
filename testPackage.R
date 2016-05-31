library(gearcalib)

## Specify range of length classes to include in the model
Lmin <- 10
Lmax <- 80
dL <- 1
Lvec <- seq(Lmin,Lmax,dL)

### Choose and combine length groups
### N must be a matrix where columns represent different length classes.
### The length classes are in L0 and defaults to the column numbers
### We subsample from N the length classes in Lvec
ChooseLengthGroups <- function(N,Lvec)
{

    A <- array(0,c(length(Lvec),ncol(N)))
    for(i in 1:(length(Lvec)-1)) A[i,Lvec[i]:(Lvec[i+1]-1)] <- 1
    Lend <- tail(Lvec,2)
    A[length(Lvec),Lend[2]:(2*Lend[2]-Lend[1]-1)] <- 1
    
    No <- N %*% t(A)
    colnames(No) <- Lvec
    return(No)
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
        d$N <- d$N[I,]
        d$group <- factor(d$group[I],levels=unique(d$group[I]))
        d$Gear <- factor(d$Gear[I],levels=GearNames)

        return(d)
    }

Dir <- "Repaired-Trawl-Pairs"
filename <- paste(Dir,"InterCal_Merluccius.RData",sep="/")
load(filename)

d=listIntercalData_Capensis
d$SweptArea=d$SweptArea[,1]

check.gearcalib.data(d)
d$group=factor(d$group)

d$Gear=factor(d$Gear)

GearNames=c("Gisund","Afr_Old")
d <- ChooseGear(d,GearNames)
d$N <- ChooseLengthGroups(d$N,Lvec)

fit <- gearcalibFit(d)

b <- boot(d)

plot(fit,select="relsel")

plot(fit,select="density",boot=b)
