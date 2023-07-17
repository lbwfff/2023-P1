rm(list= ls())

library(survminer)
library(survival)
library(ggpubr)
library(gridExtra)
library(patchwork)

########查看TCGA和GSE37642中PF4的预后影响###############
p<-list()

#TCGA

load('TCGA_LAML_ready.RData')

exp<-2^exp
rpkmTOtpm <- function(rpkm){
  exp(log(rpkm) - log(sum(rpkm)) + log(1e6))
}
exp <- apply(exp, 2, rpkmTOtpm)
exp<-log2(exp)

meta<-data.frame(name=c(colnames(exp)),
                 exp=c(as.numeric(exp[which(rownames(exp)=='PF4'),])),
                 OS=c(NA),OS_TIME=c(NA))

survival<-survival[match(meta$name,survival$sample),]
meta$OS<-survival$OS
meta$OS_TIME<-survival$OS.time

mat<-meta[!is.na(meta$OS)&
           !is.na(meta$OS_TIME),]
matt<-mat
med.exp<-median(matt$exp)
more.med.exp.index<-which(matt$exp>=med.exp)
less.med.exp.index<-which(matt$exp< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS_TIME/30,OS) ~ status, data = matt)
s.diff<-survdiff(Surv(OS_TIME/30,OS) ~ status, data = matt)

sdata.plot3<-ggsurvplot(s.fit, data=matt,
                        palette="Pastel1",
                        pval = TRUE,pval.method = TRUE, conf.int = TRUE,
                        xlab = 'Time (Month)',ggtheme = theme_survminer(),
                        surv.median.line = 'hv',size=2,
                        title=paste0("TCGA survival"))

p[[1]]<-
sdata.plot3$plot+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())


#GSE37642

load('GSE37642_ready.RData')

# genes_expr<-limma::normalizeBetweenArrays(genes_expr, method="quantile")

genes_expr<-mergeexp
phenoDat<-meta

meta<-data.frame(name=c(colnames(genes_expr)),
                 exp=c(as.numeric(genes_expr[which(rownames(genes_expr)=='PF4'),])),
                 OS=c(NA),OS_TIME=c(NA))

phenoDat<-phenoDat[match(meta$name,phenoDat$name),]
meta$OS<-as.numeric(ifelse(phenoDat$os=='dead','1','0'))
meta$OS_TIME<-as.numeric(phenoDat$ostime)
  
mat<-meta[!is.na(meta$OS)&
            !is.na(meta$OS_TIME),]
matt<-mat
med.exp<-median(matt$exp)
more.med.exp.index<-which(matt$exp>=med.exp)
less.med.exp.index<-which(matt$exp< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS_TIME/30,OS) ~ status, data = matt)
s.diff<-survdiff(Surv(OS_TIME/30,OS) ~ status, data = matt)

sdata.plot3<-ggsurvplot(s.fit, data=matt,
                        palette="Pastel1",
                        pval = TRUE,pval.method = TRUE, conf.int = TRUE,
                        xlab = 'Time (Month)',ggtheme = theme_survminer(),
                        surv.median.line = 'hv',size=2,
                        title=paste0("GSE37642 survival"))
p[[2]]<-
sdata.plot3$plot+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  theme(aspect.ratio=1)

############################################
#
# dir.create('./STEP2_Prognosis')

pdf("./STEP2_Prognosis/survplot.pdf",width = 16,height = 6)
wrap_plots(p,nrow=1) 
dev.off()

