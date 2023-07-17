rm(list= ls())

library('ggplot2')
library('MetBrewer')
library('ggpubr')
##################################################
#GSE13159、GSE9476、TCGAplusGTEx，三种数据的PF4差异分析

p<-list()

#

load('GSE13159_ready.RData')

pf4<-as.data.frame(t(genes_expr[which(rownames(genes_expr)=='PF4'),]))
  
phenoDat<-phenoDat[match(rownames(pf4),phenoDat$geo_accession),]
pf4$group<-ifelse(phenoDat$characteristics_ch1.1=='leukemia class: Non-leukemia and healthy bone marrow',
                  'Non-Tumor','Tumor')
pf4$group<-factor(pf4$group,levels = c('Non-Tumor','Tumor'))


p[[1]]<-
  ggplot(pf4, aes(x=group, y=PF4))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(title="GSE13159",y="log2 Expression")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5)

  
##############
#
load('GSE9476_ready.RData')
  
pf4<-as.data.frame(t(genes_expr[which(rownames(genes_expr)=='PF4'),]))
  
phenoDat<-phenoDat[match(rownames(pf4),phenoDat$geo_accession),]
pf4$group<-ifelse(phenoDat$source_name_ch1=='Bone Marrow',
                    'Non-Tumor','Tumor')
pf4$group<-factor(pf4$group,levels = c('Non-Tumor','Tumor'))
  
  
p[[2]]<-
  ggplot(pf4, aes(x=group, y=PF4))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="GSE9476",y="log2 Expression")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)

######################################################
  
#
load('tcgaplusgtex_ready.RData')
rownames(mergearray)<-mergearray$Ensembl_ID
mergearray<-mergearray[,-1]
  
pf4<-as.data.frame(t(mergearray[which(rownames(mergearray)==
                                          id$id[which(id$gene=='PF4')]),]))
  
pf4$group<-ifelse(grepl('GTEX',rownames(pf4)),'Non-Tumor','Tumor')
pf4$group<-factor(pf4$group,levels = c('Non-Tumor','Tumor'))
colnames(pf4)[1]<-c('PF4')
  
p[[3]]<-
  ggplot(pf4, aes(x=group, y=PF4))+
    geom_boxplot(width=0.2,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title="TCGA & GTEx",y="log2 FPKM")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5)

#拼图
# dir.create('./STEP1_DEG')

library(patchwork)
pdf("./STEP1_DEG/DEG_box_THREEDATA.pdf",width = 12,height = 8)
wrap_plots(p,nrow=1, guides="collect") 
dev.off()

#############################################################http://127.0.0.1:38145/graphics/plot_zoom_png?width=1920&height=1017
#TCGA数据库中不同临床特征下PF4的差异表达
rm(list= ls())

load('TCGA_LAML_ready.RData')

exp<-2^exp #吧RPKM改为TPM1
rpkmTOtpm <- function(rpkm){
  exp(log(rpkm) - log(sum(rpkm)) + log(1e6))
}
exp <- apply(exp, 2, rpkmTOtpm)

exp<-log2(exp)
#三个突变，三个细胞数

meta<-data.frame(name=c(colnames(exp)),
                 exp=c(as.numeric(exp[which(rownames(exp)=='PF4'),])),
                 WBC=c(NA),PB_B=c(NA),BM_B=c(NA),
                 NPM1=c(NA),IDHIR132=c(NA),FLT3=c(NA))

for (i in 1:nrow(meta)){
  meta$WBC[i]<-c(pheno$lab_procedure_leukocyte_result_unspecified_value[which(pheno$sampleID==meta$name[i])])
  meta$PB_B[i]<-c(pheno$lab_procedure_blast_cell_outcome_percentage_value[which(pheno$sampleID==meta$name[i])])
  meta$BM_B[i]<-c(pheno$lab_procedure_bone_marrow_blast_cell_outcome_percent_value[which(pheno$sampleID==meta$name[i])])
  
  mutant<-pheno$molecular_analysis_abnormality_testing_result[which(pheno$sampleID==meta$name[i])]
  mutant<-unlist(strsplit(mutant,'[|]'))
  
  # if(length(grep('NPMc R132',mutant))==0) {meta$NPM1[i]<-c('Negative')} else{ 
  # if(grepl('Negative',mutant[grep('NPMc',mutant)])) {meta$NPM1[i]<-c('Negative')} else{meta$NPM1[i]<-c('Positive')}
  # }
  
  if(length(grep('IDH1 R132',mutant))==0) {meta$IDHIR132[i]<-c('Negative')} else{ 
  if(grepl('Negative',mutant[grep('IDH1 R132',mutant)])) {meta$IDHIR132[i]<-c('Negative')} else{meta$IDHIR132[i]<-c('Positive')}
  }
  
  # if(length(grep('FLT3 R132',mutant))==0) {meta$FLT3[i]<-c('Negative')} else{ 
  # if(grepl('Negative',mutant[grep('FLT3',mutant)])) {meta$FLT3[i]<-c('Negative')} else{meta$FLT3[i]<-c('Positive')}
  # }
}


tcgamutant<-fread('tcga/TCGA.LAML.sampleMap_mutation_wustl_gene.gz')
tcgamutant<-as.data.frame(tcgamutant)

for (i in 1:nrow(meta)){
  id<-meta$name[i]
  
  if (length(grep(id,colnames(tcgamutant)))==0) {
    
    meta[i,c(6,8)]<-c(NA)
    
  } else{
  
  mutantinf<-tcgamutant[,colnames(tcgamutant) %in% c('sample',id)]
  
  mutantinf<-mutantinf[mutantinf$sample %in% c('FLT3','NPM1'),]
  
  if(mutantinf[which(mutantinf$sample=='FLT3'),2]==0) {meta$FLT3[i]<-c('Negative')} else{meta$FLT3[i]<-c('Positive')}
  if(mutantinf[which(mutantinf$sample=='NPM1'),2]==0) {meta$NPM1[i]<-c('Negative')} else{meta$NPM1[i]<-c('Positive')}
  }
}

#然后开始分组看PF4的表达
p<-list()
#WBC
wbc<-meta[!is.na(meta$WBC),]
wbc$group<-ifelse(wbc$WBC>median(wbc$WBC),'Hight','Low')
wbc$group<-factor(wbc$group,levels = c('Low','Hight'))

p[[1]]<-
ggplot(wbc, aes(x=group, y=exp))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(title="WBC",y="log2 TPM")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5)

#PB_blast
pb<-meta[!is.na(meta$PB_B),]
pb$group<-ifelse(pb$PB_B>median(pb$PB_B),'Hight','Low')
pb$group<-factor(pb$group,levels = c('Low','Hight'))

p[[2]]<-
ggplot(pb, aes(x=group, y=exp))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(title='PB Blast cell',y="log2 TPM")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5)

#BM_blast

bm<-meta[!is.na(meta$BM_B),]
bm$group<-ifelse(bm$BM_B>median(bm$BM_B),'Hight','Low')
bm$group<-factor(bm$group,levels = c('Low','Hight'))

p[[3]]<-
ggplot(bm, aes(x=group, y=exp))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(title='BM Blast cell',y="log2 TPM")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5)

#NPM1突变
nmp1<-meta[!is.na(meta$NPM1),]
nmp1$group<-ifelse(nmp1$NPM1 == 'Negative' ,'Non-mutant','Mutant')
nmp1$group<-factor(nmp1$group,levels = c('Non-mutant','Mutant'))

p[[4]]<-
  ggplot(nmp1, aes(x=group, y=exp))+
  geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  labs(title='NMP1',y="log2 TPM")+ 
  theme_classic(base_size = 22)+ 
  scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab('')+
  stat_compare_means(label.x = 1,size=5)+
  theme(aspect.ratio=1.5)

#IDH1 R132
nmp1<-meta[!is.na(meta$IDHIR132),]
nmp1$group<-ifelse(nmp1$IDHIR132 == 'Negative' ,'Non-mutant','Mutant')
nmp1$group<-factor(nmp1$group,levels = c('Non-mutant','Mutant'))
  
p[[5]]<-
ggplot(nmp1, aes(x=group, y=exp))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title='IDHIR132',y="log2 TPM")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)

#FLT3
nmp1<-meta[!is.na(meta$FLT3),]
nmp1$group<-ifelse(nmp1$FLT3 == 'Negative' ,'Non-mutant','Mutant')
nmp1$group<-factor(nmp1$group,levels = c('Non-mutant','Mutant'))
  
p[[6]]<-
ggplot(nmp1, aes(x=group, y=exp))+
    geom_boxplot(width=0.1,aes(fill=group),colour='black',alpha = 0.4,linetype="dashed")+
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=group),color="black")+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
    labs(title='FLT3',y="log2 FPKM")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=met.brewer("Hokusai1", 7)[c(6,3)])+
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    xlab('')+
    stat_compare_means(label.x = 1,size=5)+
    theme(aspect.ratio=1.5)

##########################
#拼图
pdf("./STEP1_DEG/DEG_box_WITHfeatures.pdf",width = 12,height = 16)
wrap_plots(p,nrow=2, guides="collect") 
dev.off()
    