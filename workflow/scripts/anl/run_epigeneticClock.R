## Load or install necessary libraries: cgageR, wateRmelon, parallel

if (!require("cgageR", character.only=TRUE, quietly=TRUE))
{
  options(install.packages.compile.from.source="always")
  install.packages("remotes", dependencies=TRUE, repos="https://cloud.r-project.org/")
  remotes::install_github("metamaden/cgageR")
  library(cgageR)
}

library(wateRmelon)
library(parallel)

# Load arguements
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
threads <- as.numeric(args[2])

methylation <- read.table(gzfile(filename), sep="\t", header=FALSE, stringsAsFactors=FALSE)

removeDuplicates <- function(df)
{
  checkDuplicates <- table(df[,5])
  for(n in names(checkDuplicates[which(checkDuplicates > 1)]))
  {
    where <- which(df[,5] == n)
    for(i in 2:length(where)) df[where[i],5] <- NA
  }
  df <- df[!is.na(df[,5]),]
  return(df)
}

getClockOverlap <- function(sites)
{
  # is the tagged prove name in the clock list?
  truthTable <- unlist(mclapply(methylation[,5], function(x) { sum(x == sites) > 0 }, mc.cores=threads))
  sitesToReturn <- methylation[truthTable,1:5]
  sitesToReturn <- removeDuplicates(sitesToReturn)
  rownames(sitesToReturn) <- sitesToReturn[,5]
  sitesToReturn <- sitesToReturn[,1:4]
  sitesToReturn <- sitesToReturn[sites,]
  rownames(sitesToReturn) <- sites
  sitesToReturn[is.na(sitesToReturn[,4]),4] <- 0.5
  return(sitesToReturn)
}

# run EpiTOC clock
EpiToc_table <- getClockOverlap(unique(c(
                                         EpiTOCcpgs,
                                         as.character(hannumModel[,1]),
                                         as.character(HorvathLongCGlist[,1]
                                        ))))

predictedAges <- getAgeR(EpiToc_table[,c(4,4)], epitoc=TRUE, horvath=TRUE, hannum=TRUE)
predictedAges <- c( EpiToc=as.vector(unlist(predictedAges$EpiTOC.output$EpiTOC.Est)[1]),
                    Horvath=as.vector(unlist(predictedAges$HorvathClock.output$Horvath.Est)[1]),
                    Hannum=as.vector(unlist(predictedAges$HannumClock.output$Hannum.Clock.Est)[1])
                  )

write(paste(names(predictedAges), as.vector(predictedAges), sep="="), sep="\n", file=stdout())
