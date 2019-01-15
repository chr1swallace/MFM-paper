library(gdata)
f <- "~/ASE.xlsx"
## nm <- sheetNames(f,verbose=TRUE)
gdna <- read.xls(f,1)
cm <- read.xls(f,2)
na <- read.xls(f,3)

library(data.table)
library(ggplot2)

################################################################################

## gdna - what is null?
d <- melt(as.data.table(gdna)[,c("Donor",grep("count",names(gdna),value=TRUE)),with=FALSE],
          "Donor",
          list(A=grep("A.count",names(gdna),value=TRUE),
               G=grep("G.count",names(gdna),value=TRUE)))
head(d)
d[,logratio:=log2(G)-log2(A)]
d <- d[!is.na(A),]



ggplot(d,aes(x=Donor,y=logratio)) + geom_point()
m0 <- d #[,.(logratio=mean(logratio)),by="Donor"]
m0 <- cbind("Label"="gDNA","Cell"="gDNA",m0)
mu <- mean(m0$logratio)
mu.se <- sd(m0$logratio)/sqrt(nrow(m0))
mu + c(-1,1)*1.96*mu.se

################################################################################

## cm
d <- melt(as.data.table(cm)[,c("Label","Donor","Direction",grep("count",names(gdna),value=TRUE)),with=FALSE],
          c("Label","Donor","Direction"),
          list(A=grep("A.count",names(gdna),value=TRUE),
               G=grep("G.count",names(gdna),value=TRUE)))
head(d)
d <- d[!is.na(A),]
table(d$Direction)
d[,sw:=grepl("on A read out",Direction)]
## from Dan, 6/11/2018:
## I have added how i have phased the direction of effect for the rare haplotype donors to the end of the ASE methods.  In doing this, i see that the direction is wrong in the figure for donors 1 and 2, please could you swap them. 
d[Label %in% c("Donor 1","Donor 2"),sw:=!sw]
## Q: should I swap the order of chromosomes in the schematic for A+D-het and Donor 4, and swap the direction of imbalance shown for Donor 4 in central memory?  Then we can say that the imbalance shown is top chromosome relative to bottom chromosome
## A: I like the idea of showing the ratio as the top allele to the bottom allele, and so yes, you are right, the alleles of donor 4 and A+D het need changing over. 
d[Label %in% c("Donor 4"),sw:=!sw]
d[,logratio:=ifelse(sw,-1,1) * (log2(G)-log2(A))]
ggplot(d,aes(x=Donor,y=logratio,col=Label)) + geom_point()
m <- d #[,.(logratio=mean(logratio)),by=c("Label","Donor")]
m[,Cell:="Central memory T cells"]

################################################################################

## naive
d <- melt(as.data.table(na)[,c("Label","Donor","Direction",grep("count",names(gdna),value=TRUE)),with=FALSE],
          c("Label","Donor","Direction"),
          list(A=grep("A.count",names(gdna),value=TRUE),
               G=grep("G.count",names(gdna),value=TRUE)))
head(d)
d <- d[!is.na(A),]
table(d$Direction)
d[,sw:=grepl("on A read out",Direction)]
## from Dan, 6/11/2018:
## I have added how i have phased the direction of effect for the rare haplotype donors to the end of the ASE methods.  In doing this, i see that the direction is wrong in the figure for donors 1 and 2, please could you swap them. 
d[Label %in% c("Donor 1","Donor 2"),sw:=!sw]
## Q: should I swap the order of chromosomes in the schematic for A+D-het and Donor 4, and swap the direction of imbalance shown for Donor 4 in central memory?  Then we can say that the imbalance shown is top chromosome relative to bottom chromosome
## A: I like the idea of showing the ratio as the top allele to the bottom allele, and so yes, you are right, the alleles of donor 4 and A+D het need changing over. 
d[Label %in% c("Donor 4"),sw:=!sw]
d[,logratio:=ifelse(sw,-1,1) * (log2(G)-log2(A))]
ggplot(d,aes(x=Donor,y=logratio,col=Label)) + geom_point()
n <- d #[,.(logratio=mean(logratio)),by=c("Label","Donor")]
n[,Cell:="Naive T cells"]

m <- rbind(m[,.(Label,Donor,Cell,logratio)],m0[,.(Label,Donor,Cell,logratio)],n[,.(Label,Donor,Cell,logratio)])
m[,Label:=factor(Label,levels=c("gDNA","A-het","D-het","A+D-het","Donor 1","Donor 2","Donor 3","Donor 4"))]
m[,Cell:=factor(Cell,levels=c("gDNA","Central memory T cells","Naive T cells"))]

library(cowplot)

## ggplot(m[!grepl("Donor",Label),],aes(x=Label,y=2^logratio,col=Label)) +
##   geom_boxplot() + geom_point() +  geom_hline(yintercept=1,lty="dashed") +
##   facet_grid(.~Cell,scales="free_x",space="free_x") +
##   labs(x="Genotype group",y="Allelic Ratio") +
##   scale_y_continuous(breaks=c(0.5,1,1.5,2.0)) +
##   theme(legend.position="none") +
##   background_grid()
## ggsave("~/share/Projects/cvs/figures/ASE-groups.pdf",height=4,width=6)

m1 <- m[grepl("gDNA|Donor",Label),]

m1 <- m1[order(Cell,Donor),]
m1[,NLabel:=as.character(Label)]
m1[,OLabel:=as.character(Label)]
m1[Cell=="gDNA",OLabel:=""]
m1[Cell=="gDNA",NLabel:=as.character(as.numeric(Donor))]
ggplot(m1,aes(x=NLabel,y=2^logratio,col=Label)) +
  geom_hline(yintercept=1,lty="dashed") +
  geom_boxplot() + geom_point(pch=17) +  
  facet_grid(.~Cell,scales="free_x") +
  labs(x="Genotype group",y="Allelic Ratio") +
  scale_y_continuous(breaks=c(0.5,1,1.5,2.0),limits=c(0.5,2)) +
  scale_x_discrete(breaks=unique(m1$NLabel),labels=m1$OLabel[!duplicated(m1$NLabel)]) +
  theme(legend.position="none") +
  background_grid()
ggsave("~/share/Projects/cvs/figures/ASE-donors.pdf",height=5,width=7.5)

################################################################################

### tests
m <- m[grepl("DNA|het",Label),.(logratio=mean(logratio)),by=c("Cell","Label","Donor")]
dna <- m[Cell=="gDNA",]$logratio
res <- m[grep("het",Label), t.test(dna,logratio)[c("p.value")],by=c("Label","Cell")]
mp <- merge(m,res,by=c("Label","Cell"),all.x=TRUE)

myfp <- function(x) {
    x <- sapply(x,format,digits=2)
    wh <- grep("e",x)
    if(length(wh)) 
        x[ wh ]  <- paste0(sub("e-0","*'x'*10^{-", x[wh] ),"}")
    paste0("'p ='~",x)
}
myfp(res$p.value)

ggplot(m[!grepl("Donor",Label),],aes(x=Label,y=2^logratio,col=Label)) +
  geom_hline(yintercept=1,lty="dashed") +
  geom_boxplot() + geom_point() +
  geom_text(aes(x=Label,y=0.15,label=myfp(p.value)),data=res,parse=TRUE,col="black") +
  facet_grid(.~Cell,scales="free_x",space="free_x") +
  labs(x="Genotype group",y="Allelic Ratio") +
  scale_y_continuous(breaks=c(0.5,1,1.5,2.0)) +
  theme(legend.position="none") +
  background_grid()
ggsave("~/share/Projects/cvs/figures/ASE-groups.pdf",height=5,width=7.5)


