rm(list= ls())

library(ggstatsplot)
library(data.table)
library(ggplot2)
library(ggpubr)

################################################################
#GSCA 

gsca<-read.csv('GdscIC50AndExprTable.csv')
gsca<-gsca[gsca$symbol=='PF4',]
gsca<-gsca[order(gsca$fdr),]


library(ggplot2)
library(dplyr)


pdf("./STEP4_Drug_sensitivity/PF4_GSCA_cor.pdf",width = 8,height = 6)
ggplot(gsca[1:10,],aes(x=reorder(drug,cor),y=cor))+
  geom_point()+
  geom_segment(aes(x=drug,xend=drug,y=0,yend=cor),size=1.5,color="#C9CACA",linetype="solid")+
  geom_point(size=5,aes(color=-log10(fdr),fill=-log10(fdr)),shape=21)+
  theme_light(base_size = 22)+theme(panel.grid.major.x=element_blank(),
                       panel.border=element_blank(),
                       axis.ticks.x=element_blank())+
  xlab("")+ylab("")+coord_flip()+theme(aspect.ratio=1)
dev.off()

###############################################
#免疫检查点

genelist<-c('PF4','PD-L1','CD274','PD-1','PDCD1','CTLA4')

load('TCGA_LAML_ready.RData')

checkpoint<-exp[row.names(exp) %in% genelist,]
checkpoint<-as.data.frame(t(checkpoint))

p<-list()

p[[1]]<-
ggscatterstats(data = checkpoint, 
               x = PF4,
               y = PDCD1,
               xlab = "PF4 expression (log2 FPKM)",
               ylab = "PDCD1 expression (log2 FPKM)",
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram",
               marginal = TRUE)

p[[2]]<-
ggscatterstats(data = checkpoint, 
               x = PF4,
               y = CD274,
               xlab = "PF4 expression (log2 FPKM)",
               ylab = "CD274 expression (log2 FPKM)",
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram",
               marginal = TRUE)

p[[3]]<-
ggscatterstats(data = checkpoint, 
               x = PF4,
               y = CTLA4,
               xlab = "PF4 expression (log2 FPKM)",
               ylab = "CTLA4 expression (log2 FPKM)",
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram",
               marginal = TRUE)

pdf("./STEP4_Drug_sensitivity/Checkpoint_cor.pdf",width = 18,height = 6)
patchwork::wrap_plots(p,nrow=1) 
dev.off()



