
library(annotSnpStats)
library(snpStats)
library(GUESSFM)
library(magrittr)
library(data.table)
snp1 <- "rs706779"
msnps <- c("rs2476491","rs61839660")

regions <- list.files("/rds/project/cew54/rds-cew54-wallace-share/Projects/cvs/",
                      pattern="[0-9][pq]-")
reg <- "10p-6030000-6220000"

setwd("~/Projects/cvs/tagflag")
library(Rcpp)
sourceCpp("tag-pattern-flag.cpp")

devtools::load_all("~/RP/coloc")


## library(parallel)
## options(mc.cores=4)
results <- vector("list",length(regions))
names(results) <- regions
reg <- regions[[1]]
for(reg in sample(regions)) {
	ofile = paste0(reg,".csv")
if(file.exists(ofile))
next
    message(reg)
    load(paste("/rds/project/cew54/rds-cew54-wallace-share/Projects/cvs/",reg,"/GUESS/all-data.RData",sep="")) # international DATA, pcs
    DATA <- DATA[which(samples(DATA)$country == "UK"),] 
    if("aSnpMatrix" %in% class(DATA))
        DATA <- sm(DATA)
    cs <- col.summary(DATA)
    DATA <- DATA[,!is.na(cs$z.HWE) & cs$MAF>0.05]
    Gmat <- as(DATA,"numeric")
    S <- coloc:::cor2(Gmat)
    maf <- cs$MAF[ !is.na(cs$z.HWE) & cs$MAF>0.05]
    ## ret <- lapply(c(0.01,0.05,0.1),function(thr) {
    thr <- 0.05
    ##     use <- which(maf > thr)
    n <- nrow(S)
    if(n<2)
        return(NULL)
    maxtags <- n*(n-1)/2
    atleast <- countatleastone(S)
    patt <- countpatt(S)
    atleast.f <- countatleastonemaf(S,maf)
    patt.f <- countpattmaf(S,maf)
    ret<- list(reg=reg,nsnps=n,thr=thr,maxtags=maxtags,tags=atleast,tags.f=atleast.f,maxpatt=n*(n-1)*(n-2)/2,patt=patt,patt.f=patt.f)  %>%
      as.data.table()
    fwrite(ret,file=ofile)
     results[[reg]] <- ret
}

res <- rbindlist(results)
fwrite(res,file="tagpatt.csv")

