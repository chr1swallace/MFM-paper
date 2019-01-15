#!/usr/bin/env Rscript
#redoing run 13/07/18, wants standard sunbeams (sunbeam_whichtop)
library(devtools)
## install_github("mdfortune/simGWAS")
library(simGWAS)
library(corpcor)
library(mvtnorm)
library(ggplot2)
library(data.table)
library(parallel)
library(cowplot)
library(oce)


GRIDSIZE <- 0.005
PLOTSIZE <- 4

source("~/Projects/cvs/sunbeams/sunbeam.R")
d <- "~/share/Projects/cvs/sunbeams"

## read scenarios
library(yaml)
scen <- read_yaml("~/Projects/cvs/sunbeams/figs.yaml")

nm <- names(scen)[[1]]
i <- 1
ORgrid=GRIDSIZE
## SNPint=Tag

library(parallel)
options(mc.cores=8)
mclapply(names(scen), function(nm) {
## for(nm in names(scen)) {
    message(nm)
    load(file.path(d,scen[[nm]]$data))
    snps<-colnames(freq)[-ncol(freq)]
    LD <- snpStats::ld(XX,XX,stat="R",symmetric=TRUE)
    LD<-as.matrix(make.positive.definite(LD))
    CV1<-scen[[nm]]$snps[[ 1 ]]
    CV2<-scen[[nm]]$snps[[ 2 ]]
    Tag<-scen[[nm]]$snps[[ 3 ]]
    gammatrue<-scen[[nm]]$gammatrue
    gammagt<-scen[[nm]]$gammagt
    gammalt<-scen[[nm]]$gammalt
    for(i in seq_along(scen[[nm]]$N0)) {
        pngprob <- file.path(d,"../figures",paste0("sunbeam_prob_",nm,i,".png"))
        pdfprob <- file.path(d,"../figures",paste0("sunbeam_prob_",nm,i,".pdf"))
	pdfwhich <- file.path(d,"../figures",paste0("sunbeam_which_",nm,i,".pdf"))
	pngwhich <- file.path(d,"../figures",paste0("sunbeam_which_",nm,i,".png"))
	## if(file.exists(pngprob) && file.exists(pdfprob)) 
	## 	next
	N0<-scen[[nm]]$N0[i]
        N1<-scen[[nm]]$N1[i]
        message("\t",N0,"\t",N1) 
        ## system.time(tmp1 <- calc_whichtop(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.01)) ## 22
        ## system.time(tmp2 <- calc_whichtop2(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.01)) ## 16
        ## system.time(tmp3 <- calc_whichtop3(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.01)) ## 2
        ## system.time(calc_whichtop(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.005)) ## 83
        ## system.time(calc_whichtop2(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.005)) ## 
        ## system.time(calc_whichtop3(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.005)) ## 4
        ## system.time(calc_whichtop3(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.001)) ## 24
                       
        ## system.time(tmp1 <- calc_probtop(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.1)) ## 0.48
        ## system.time(tmp3 <- calc_probtop3(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.1)) ## 0.32
        ## system.time(tmp4 <- calc_probtop4(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.1)) ## 0.2

        ## par(mfrow=c(2,2))
        ## image(tmp1); image(tmp3); image(tmp4)
        

        ## system.time(tmp1 <- calc_probtop(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.01)) ## 49
        ## system.time(tmp3 <- calc_probtop3(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.01)) ## 17
        ## system.time(tmp4 <- calc_probtop4(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.01)) ## 7
        ## system.time(tmp3 <- calc_probtop3(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.005)) ## 65
        ## system.time(tmp4 <- calc_probtop4(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=0.005)) ## 26
      
        tmp <- calc_whichtop3(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=GRIDSIZE)
        p <- plot_whichtop(tmp,snps=names(scen[[nm]]$snps),
                           ORmax=1,gammatrue=gammatrue) +
          theme(legend.position=c(0.8,0.5),
                legend.background=element_rect(fill="white",colour="black",size=0.5,linetype=1)) +
          ggtitle(paste(paste(names(scen[[nm]]$snps)[1:2],collapse="+"), "causal"))
        p
        
         ggsave(p,file=pdfwhich,height=PLOTSIZE,width=PLOTSIZE)
        ggsave(p,file=pngwhich, height=PLOTSIZE,width=PLOTSIZE,dpi=300)

        tmp <- calc_probtop4(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORgrid=GRIDSIZE)
        p <- plot_probtop(tmp,snps=names(scen[[nm]]$snps),
                           ORmax=1,gammatrue=gammatrue) +
          theme(legend.position=c(0.8,0.5),
                legend.background=element_rect(fill="white",colour="black",size=0.5,linetype=1)) +
          ## theme(legend.position="right") +
          ggtitle(paste(paste(names(scen[[nm]]$snps)[1:2],collapse="+"), "causal"))
               if(!is.null(gammagt))
                   p=p+annotate("text",x=gammagt[1],y=gammagt[2],label=">",size=6)
        if(!is.null(gammalt))
            p= p+annotate("text",x=gammalt[1],y=gammalt[2],label="<",size=6)
        ggsave(p,file=pdfprob,height=PLOTSIZE,width=PLOTSIZE)
        ggsave(p,file=pngprob, height=PLOTSIZE,width=PLOTSIZE,dpi=300)
    }

    ## if(nm=="DEXI_MS") {
    ##     .probtop <- function(x,y) {
    ##     GenoProbList<-make_GenoProbList(snps,c(CV1,CV2),freq)
    ##     est_z_score<-abs(est_statistic(N0,N1,snps,c(CV1,CV2),gammatrue,freq,GenoProbList))
    ##     if (max(abs(est_z_score))<z_score_sig)
    ##         return(0)
    ##     if(bound_prob_top(est_z_score[which(snps==SNPint)],max(est_z_score),LD[which.max(est_z_score),which(snps==SNPint)])<0.1)
    ##         return(-1)
    ##     prob_snp_top(which(snps==SNPint),est_z_score,LD,num_itt=100)
    ## }
    ## }
    
})
