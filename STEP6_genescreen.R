rm(list= ls())
# dir.create('./STEP6_genescreen')

library('limma')
###############################################
#筛选基因，取交集，功能富集

#DEG ，GSE13159,GSE9476

load('GSE13159_ready.RData')

# exp<-genes_expr
exp <- normalizeBetweenArrays(genes_expr, method="quantile")
exp<-as.data.frame(exp)

phenoDat<-phenoDat[match(colnames(exp),phenoDat$geo_accession),]

group<-ifelse(phenoDat$characteristics_ch1.1=='leukemia class: Non-leukemia and healthy bone marrow',
              'nontumor','tumor')

design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=phenoDat$geo_accession
design

contrast.matrix<-makeContrasts(tumor-nontumor,levels = design)
contrast.matrix 

fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 

#可视化

table<-nrDEG
table<-table[!is.na(table$adj.P.Val),]
table$group = factor(ifelse(table$P.Value < 0.05 & abs(table$logFC) >= 1, ifelse(table$logFC>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
table(table$group)
table$X<-rownames(table)

library(ggplot2)
pdf("./STEP6_genescreen/GSE13159_volcano.pdf",width = 9,height = 6)
ggplot(table,aes(x=logFC,y=-log10(P.Value),color=group))+
  geom_point(size=2)+
  geom_text_repel(data=table[table$group!='NoSignifi',],aes(label=X),size=3,segment.color='black',show.legend=F)+
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_blank())+
  ylab('-log10 (padj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
  scale_color_manual(values=c("#dd5129","#808fe1","#b9b9b8"),
                     breaks=c("Up", "Down", "NoSignifi"),
                     labels=c(paste0('Up',' (',nrow(table[table$group=='Up',]),')'),
                              paste0('Down',' (',nrow(table[table$group=='Down',]),')'),
                              paste0('NoSignifi',' (',nrow(table[table$group=='NoSignifi',]),')')))
dev.off()

#热图
normalized_counts <- exp
n<-50
pheatmap<-normalized_counts[rownames(normalized_counts) %in% 
                              rownames(nrDEG)[1:n],]

library(ComplexHeatmap)
library(circlize)

mat_scaled = t(scale(t(pheatmap)))
col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4]))

df <- data.frame(group = c(group))
df$group<-factor(df$group,levels = c('nontumor','tumor'))
df$what<-c('whatever')

mat_scaled<-cbind(mat_scaled[,which(df$group!='tumor')],mat_scaled[,which(df$group=='tumor')])
df<-rbind(df[which(df$group!='tumor'),],df[which(df$group=='tumor'),])

split = c(ifelse(df$group=='nontumor',1,2))
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2),labels = c("Non-tumor","Tumor"), 
                                         labels_gp = gpar(col = "white", fontsize = 14)))

pdf('./STEP6_genescreen/GSE13159_heatmap.pdf',width = 40,height = 10)
draw(
  Heatmap(mat_scaled,
          col = col_fun1,name = "Score",
          column_split = split,
          column_title = NULL,
          # column_gap = unit(c(0, 0), "mm"),
          cluster_columns = F,
          cluster_rows =T,
          show_row_dend = T,
          show_column_names =F,
          show_heatmap_legend = T,
          row_names_side = "right",
          width = ncol(mat_scaled)*unit(1, "mm"), 
          height = nrow(mat_scaled)*unit(4, "mm"),
          top_annotation = ha))
dev.off()

GSE13159<-table$X[table$group!='NoSignifi']

###############
#GSE9476
load('GSE9476_ready.RData')

phenoDat<-phenoDat[match(colnames(genes_expr),phenoDat$geo_accession),]

# exp<-genes_expr
exp <- normalizeBetweenArrays(genes_expr, method="quantile")
exp<-as.data.frame(exp)

group<-ifelse(phenoDat$source_name_ch1=='Bone Marrow',
              'nontumor','tumor')

design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))

identical(phenoDat$geo_accession,colnames(exp))
rownames(design)=phenoDat$geo_accession
design

contrast.matrix<-makeContrasts(tumor-nontumor,levels = design)
contrast.matrix 

fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
nrDEG<-nrDEG[order(nrDEG$adj.P.Val),]

#可视化

table<-nrDEG
table<-table[!is.na(table$adj.P.Val),]
table$group = factor(ifelse(table$adj.P.Val < 0.05 & abs(table$logFC) >= 1, ifelse(table$logFC>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
table(table$group)
table$X<-rownames(table)

library(ggplot2)
pdf("./STEP6_genescreen/GSE9476_volcano.pdf",width = 9,height = 6)
ggplot(table,aes(x=logFC,y=-log10(adj.P.Val),color=group))+
  geom_point(size=2)+
  geom_text_repel(data=table[table$group!='NoSignifi',],aes(label=X),size=3,segment.color='black',show.legend=F)+
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_blank())+
  ylab('-log10 (padj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
  scale_color_manual(values=c("#dd5129","#808fe1","#b9b9b8"),
                     breaks=c("Up", "Down", "NoSignifi"),
                     labels=c(paste0('Up',' (',nrow(table[table$group=='Up',]),')'),
                              paste0('Down',' (',nrow(table[table$group=='Down',]),')'),
                              paste0('NoSignifi',' (',nrow(table[table$group=='NoSignifi',]),')')))
dev.off()

#热图
normalized_counts <- exp
n<-50
pheatmap<-normalized_counts[rownames(normalized_counts) %in% 
                              rownames(nrDEG)[1:n],]

library(ComplexHeatmap)
library(circlize)

mat_scaled = t(scale(t(pheatmap)))
col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4]))

df <- data.frame(group = c(group))
df$group<-factor(df$group)

split = c(ifelse(df$group=='nontumor',1,2))
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2),labels = c("Non-tumor", "Tumor"), 
                                         labels_gp = gpar(col = "white", fontsize = 14)))

pdf('./STEP6_genescreen/GSE9476_heatmap.pdf',width = 16,height = 10)
draw(
  Heatmap(mat_scaled,
          col = col_fun1,name = "Score",
          column_split = split,
          column_title = NULL,
          column_gap = unit(c(0.1, 0.1), "mm"),
          cluster_columns = F,
          cluster_rows =T,
          show_row_dend = T,
          show_column_names =F,
          show_heatmap_legend = T,
          row_names_side = "right",
          width = ncol(mat_scaled)*unit(4, "mm"), 
          height = nrow(mat_scaled)*unit(4, "mm"),
          top_annotation = ha))
dev.off()

GSE9476<-table$X[table$group!='NoSignifi']


###################################################
#venn1

library(ggvenn)

venn<-list(GSE13159=c(GSE13159),GSE9476=c(GSE9476))

pdf("./STEP6_genescreen/VENN1.pdf",width = 6,height = 4)
ggvenn(venn,fill_color = c("#0073C2FF", "#EFC000FF"),
          stroke_size = 1, set_name_size = 6)
dev.off()
geolist<-GSE13159[GSE13159%in% GSE9476]

#venn2

deseq<-read.csv('./STEP3_GSEA/TCGA_PF4HIGHT_VS_PF4LOW.csv')
deseq<-deseq[deseq$baseMean>5,]
deseq<-deseq[-1,]
table<-deseq
table<-table[!is.na(table$padj),]
table$group = factor(ifelse(table$pvalue < 0.05 & abs(table$log2FoldChange) >= 1, ifelse(table$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
table(table$group)

venn<-list(TCGA=c(table$X[table$group!='NoSignifi']),GEO=c(geolist))

pdf("./STEP6_genescreen/VENN2.pdf",width = 6,height = 4)
ggvenn(venn,fill_color = c("#0073C2FF", "#EFC000FF"),
          stroke_size = 1, set_name_size = 6)
dev.off()
#
list<-geolist[geolist %in% table$X[table$group!='NoSignifi']]

write.csv(list,file = 'genelist.csv',row.names = F)

#####################################################
#富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGG.db)

list<-read.csv('genelist.csv')
list<-list$x

enrich<-bitr(list, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)[2]

kk.negative <- enrichKEGG(gene  = enrich$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.5,
                          qvalueCutoff =0.5)

edox <- setReadable(kk.negative, 'org.Hs.eg.db', 'ENTREZID')

pdf("./STEP6_genescreen/kegg.pdf",width = 10,height =6)
cnetplot(edox,  circular = TRUE, colorEdge = TRUE) 
dev.off()

go<-enrichGO(gene          = enrich$ENTREZID,
         OrgDb         = org.Hs.eg.db,
         ont           = 'ALL' ,
         pAdjustMethod = "BH",
         pvalueCutoff  = 0.2,
         qvalueCutoff  = 0.5,
         readable      = TRUE)

edox <- setReadable(go, 'org.Hs.eg.db', 'ENTREZID')

pdf("./STEP6_genescreen/go.pdf",width = 10,height = 6)
cnetplot(edox,circular = TRUE, colorEdge = TRUE)
dev.off()





