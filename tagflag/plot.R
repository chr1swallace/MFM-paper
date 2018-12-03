
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
setwd("~/Projects/cvs/tagflag")
library(data.table)
files <- list.files(".",pattern=".csv")

x <- lapply(files,fread)  %>%
  lapply(., function(z)
      z[,.(reg,nsnps,thr,tags,tags.f,patt,patt.f)])  %>%
  rbindlist()

m <- melt(x,c("reg","nsnps","thr"))
m[grepl("tags",variable),max:=nsnps*(nsnps-1)/2]
m[grepl("patt",variable),max:=nsnps*(nsnps-1)*(nsnps-2)/2]
m[,group:=ifelse(grepl("tags",variable),"at least 1 tag","tagged trios")]
m[,effect:=ifelse(grepl(".f",variable),"equal OR","equal Z")]


library(ggplot2)
library(cowplot)
library(RColorBrewer)
p1 <- ggplot(m[grepl("tags",variable),],aes(x=nsnps,y=100*value/max,col=effect,fill=effect)) +
  geom_smooth() +
  geom_point() +
  ## scale_x_log10() +
  xlim(0,1750) +
  background_grid() +
  labs(x="# SNPs",y="% tagging") +
  scale_colour_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(legend.position="bottom") +
  ggtitle("At least one tag")

p0 <- ggplot(m[grepl("patt",variable),],aes(x=nsnps,y=100*value/max,col=effect,fill=effect)) +
  geom_smooth() +
  geom_point() +
  ## scale_x_log10() +
  xlim(0,1750) +
  background_grid() +
  labs(x="# SNPs",y="% tagging") +
  scale_colour_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(legend.position="bottom") +
  ggtitle("Trios that are tags")

plot_grid(p0,p1,labels=c("d","e"))

ggsave(file="~/share/Projects/cvs/figures/tagging-frequency.pdf",width=8,height=6)

head(m)
m[nsnps<2000,sum(value)/sum(max),by="variable"]
