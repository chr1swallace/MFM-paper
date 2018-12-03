#!/usr/bin/env Rscript
## IL2RA
## rs61839660 rs62626317 rs2104286
## rs61839660 rs11594656 rs706779
## IL2/IL21 (4q-122973062-123565302)
## rs6837165 rs13122213 rs77516441
source("~/DIRS.txt")
library(annotSnpStats)
library(magrittr)
devtools::load_all("~/RP/snpHaps")
devtools::load_all("~/RP/simGWAS")

library(yaml)
todo <- read_yaml("~/Projects/cvs/sunbeams/figs.yaml")
## todo <- list(list(name="IL2IL21",
##                   d=file.path(CVSROOT,"4q-122973062-123565302"),
##                   snps=c("rs6837165", "rs13122213", "rs77516441")),
##              list(name="IL2RA_AITD",
##                   d=file.path(CVSROOT,"10p-6030000-6220000"),
##                   snps=c( "rs61839660", "rs11594656", "rs706779" )),
##              list(name="IL2RA_MS",
##                   d=file.path(CVSROOT,"10p-6030000-6220000"),
##                   snps=c( "rs61839660", "rs62626317", "rs2104286" ))
## 		list(name="14q",
## d=file.path(CVSROOT,"14q-101290463-101328739"),
## snps=c("")))
x <- todo[[1]]

for(x in todo) {
    of <- file.path("~/share/Projects/cvs/sunbeams",x$data)
    if(file.exists(of))
        next
    f.geno <- file.path(CVSROOT,x$region,"GUESS","all-data.RData")
    (load(f.geno))
    p <- samples(DATA)
    use <- which(p$country=="UK" & p$phenotype=="CONTROL")

    hsnps <- sapply(x$snps, grep,colnames(DATA),value=TRUE)

    XX <- as(DATA[use,hsnps],"SnpMatrix")
    DATA@snps$A1 <- 1
    DATA@snps$A2 <- 2

    f <- tempfile()
    haps <- genhaps(XX,slist=hsnps,
                snps=DATA@snps,
                samples=df,
                A1="A1",A2="A2",
                cols.copy=c("phenotype","PC1","PC2","PC3","PC4"),
                f.in=f,
                redo=TRUE)
    freq <- read.table(paste0(f,".out1"),header=TRUE)
    save(XX,freq,file=of)
}
