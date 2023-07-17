rm(list= ls())

library('data.table')
library('DESeq2')
library('ggplot2')
library('MetBrewer')
library('clusterProfiler')
library('tidyr')
library('tibble')
library('corrplot')
library('ggrepel')

######################################################
#使用DESeq2对PF4高低表达的样本进行差异分析，后使用ClusterProfiler进行GSEA富集

load('TCGA_LAML_ready.RData')
exp<-count
group<-data.frame(name=c(colnames(exp)),
                  exp=c(as.numeric(exp[which(rownames(exp)=='PF4'),])),
                  group=c(NA))
group$group<-ifelse(group$exp>median(group$exp),
                    'PF4 Hight','PF4 Low')

exp<-(2^exp-1) #先把矩阵还原到Count
expadj<-apply(exp,2,round)

condition <- factor(group$group)
dds <- DESeqDataSetFromMatrix(expadj, DataFrame(condition), design= ~ condition )
dds

dds$condition<- relevel(dds$condition, ref = "PF4 Low") # 指定哪一组作为对照组
dds <- DESeq(dds)
dds

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)

# dir.create('./STEP3_GSEA')

write.csv(res,file="./STEP3_GSEA/TCGA_PF4HIGHT_VS_PF4LOW.csv")

deseq<-read.csv('./STEP3_GSEA/TCGA_PF4HIGHT_VS_PF4LOW.csv')
deseq<-deseq[deseq$baseMean>5,]

deseq<-deseq[-1,]

##############可视化#####################
#对应8.0的内容，火山图和热图
table<-deseq
table<-table[!is.na(table$padj),]
table$group = factor(ifelse(table$pvalue < 0.05 & abs(table$log2FoldChange) >= 1, ifelse(table$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
table(table$group)

library(ggplot2)
pdf("./STEP3_GSEA/TCGA_volcano.pdf",width = 9,height = 6)
ggplot(table,aes(x=log2FoldChange,y=-log10(pvalue),color=group))+
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
normalized_counts <- counts(dds,normalized=T) 
n<-50
pheatmap<-normalized_counts[rownames(normalized_counts) %in% 
  deseq$X[1:n],]

library(ComplexHeatmap)
library(circlize)

mat_scaled = t(scale(t(pheatmap)))
col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4]))

df <- data.frame(group = c(group$group))
df$group<-factor(df$group)

split = c(ifelse(df$group=='PF4 Hight',1,2))
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4),labels = c("PF4 Hight", "PF4 Low"), 
                                         labels_gp = gpar(col = "white", fontsize = 14)))

pdf('./STEP3_GSEA/TCGA_heatmap.pdf',width = 16,height = 10)
  draw(
  Heatmap(mat_scaled,rect_gp = gpar(col = "white", lwd = 1),
          col = col_fun1,name = "Score",
          column_split = split,
          column_title = NULL,
          column_gap = unit(c(0.1, 0.1), "mm"),
          cluster_columns = T,
          cluster_rows =T,
          show_row_dend = T,
          show_column_names =F,
          show_heatmap_legend = T,
          row_names_side = "right",
          width = ncol(mat_scaled)*unit(2, "mm"), 
          height = nrow(mat_scaled)*unit(4, "mm"),
          top_annotation = ha))
dev.off()
  
  #这个聚类为什么这么奇怪？
  

###############GSEA######################

library(org.Hs.eg.db)
# columns(org.Hs.eg.db)
# k <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
# all_gene<-AnnotationDbi::select(org.Hs.eg.db,keys = k,columns = c("SYMBOL",'ENTREZID','GENENAME','GENETYPE'),keytype="ENSEMBL")

#
# match<-all_gene[match(deseq$X,all_gene$SYMBOL),]
# deseq$ENTREZID<-match$ENTREZID
# 
gsea<-deseq

gsea<-gsea[order(gsea$log2FoldChange,decreasing = T),]

test<-as.numeric(gsea$log2FoldChange)
names(test) = as.character(gsea$X) #记住names这个属性

{ 
  #选择gmt文件
  gmtfile ='c2.cp.v2023.1.Hs.symbols.gmt'
  #GSEA分析
  library(GSEABase) 
  geneset <- clusterProfiler::read.gmt( gmtfile )  
  length(unique(geneset$term))
  egmt <- GSEA(test, TERM2GENE=geneset, 
               minGSSize = 1,
               pvalueCutoff = 0.99,
               verbose=FALSE)
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  rownames(gsea_results_df)
  write.csv(gsea_results_df,file = 'STEP3_GSEA/gsea_results_df.csv')
  
}

#展示前五的基因集
ridgeplot(egmt,showCategory = 20)

library(enrichplot)
gsea_results_df<-gsea_results_df[-grep('COV2',gsea_results_df$Description),]

for (i in 1:5){
geneset_id<-gsea_results_df$ID[i]
p<-gseaplot2(egmt, geneSetID = geneset_id, pvalue_table=F,subplots = 1:3,rel_heights = c(2, 0.8, 1.2),
             title = egmt$Description[which(egmt@result$ID==geneset_id)]) #标题有点太冗余了，要改倒是也容易
p[[1]] <- p[[1]]+ annotate("text", x = 12000, y = 0.4, label = paste0('P.value = ',signif(egmt@result$pvalue[which(egmt@result$ID==geneset_id)],3),
                                                                       '\n NES = ',round(egmt@result$NES[which(egmt@result$ID==geneset_id)],2),
                                                                      '\n FDR = ',signif(egmt@result$p.adjust[which(egmt@result$ID==geneset_id)],3)),size = 5)

PATH<-paste0('./STEP3_GSEA/',geneset_id,'_GSEA.pdf')
pdf(PATH,width = 6,height = 5)
print(p)
dev.off()
} 

###############################################################
#免疫浸润

write.table(exp,file = 'TCGA_exp.txt',sep = '\t',quote = F)

source('Cibersort_source.R')

result <- CIBERSORT('LM22.txt','TCGA_exp.txt', perm = 100, QN = F)

group<-data.frame(name=c(colnames(exp)),
                  exp=c(as.numeric(exp[which(rownames(exp)=='PF4'),])),
                  group=c(NA))
group$group<-ifelse(group$exp>median(group$exp),
                    'PF4 Hight','PF4 Low')

#合并分组信息
result1 = as.data.frame(result)
result1$name<-rownames(result1)
result1 <- merge(group,result1,by="name")

row.names(result1) <- result1[,1]
result1 <- result1[,-1]
write.csv(result1,"./STEP3_GSEA/CIBERSORT_result.csv")

#根据pf4表达水平组排序
re1 <- result1[order(result1[,1],decreasing = T),]
colnames(re1)[2] <- "Type"

re1 = subset(re1,re1$`P-value`<0.05)
re1 = re1[,-c(1,25:27)]
re2<-re1[,-1]

mypalette <- colorRampPalette(brewer.pal(8,"Set3"))
#提取数据，多行变成多列，要多学习‘tidyr’里面的三个函数
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
#提取数据
dat_group = gather(re1,Cell_type,Proportion,-Type )
#合并分组
dat = cbind(dat_cell,dat_group$Type)
colnames(dat)[4] <- "Type"

##2.2柱状图##############
p1 <- ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(28))
p1
#画分组bar
Var1 = re1$Type #high risk & low risk的数量
Var2=c(1:length(row.names(re1))) #high risk & low risk长度，从1开始
annotation_col<- data.frame(Var1 = re1$Type,
                            Var2=c(1:length(row.names(re1))),
                            value=1)
p2 <- ggplot(dat, aes(Var2, Var1),size=0.5)+
  geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
  labs(x = "", y = "")+
  theme_void()+
  theme(legend.position = "top", legend.title = element_blank())+
  scale_x_discrete(name="")+
  scale_fill_manual(labels=c("High PF4 group(n=59)","Low PF4 group(n=47)"),
                    values =c("#DC0000FF","#00A087FF"))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p2
#拼图
pdf('./STEP3_GSEA/CIBERSORT_radio.pdf',height = 6,width = 10)
cowplot::plot_grid(p2,p1,
          rel_heights = c(0.08,1),
          label_x=0,
          label_y=1,
          align = 'v',ncol = 1,greedy = F)
dev.off()

##2.3 免疫细胞热图#############
data <- as.data.frame(t(re1[,-1]))

data_1 <- as.data.frame(lapply(data,as.numeric))
row.names(data_1) <- row.names(data)
data=data_1

library(pheatmap)
pdf(file="./STEP3_GSEA/immunecell_pheatmap.pdf",width=10,height=8,onefile=FALSE)
annotation_col=data.frame(Group=rep(c("Low risk group(n=59)","High risk group(n=47)"),c(59,47)))
annColors <- list(Group = c("Low risk group(n=59)" = "#00A087FF","High risk group(n=47)" ="#DC0000FF"))
rownames(annotation_col)=colnames(data)
pheatmap(data,
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cluster_rows=TRUE,
         show_rownames=TRUE,
         fontsize_row=10,
         show_colnames=FALSE,
         scale="row",
         cluster_cols=FALSE,
         main="")
dev.off()

##2.4组免疫细胞的差异###################
library(ggpubr)
library(ggplot2)
box=dat

#图片美化
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小     
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗 
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.8),#图例在绘图区域的位置
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}

e1 <- ggplot(box,aes(x=reorder(Cell_type,-Proportion),y=Proportion),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=Type),position=position_dodge(0.5),width=0.6)+
  labs(x = "Cell Type", y = "Estimated Proportion")+
  scale_fill_manual(values = c("#DC0000FF","#00A087FF")) +
  theme_zg()
e1 = e1 + stat_compare_means(aes(group = Type),label = "p.signif")
e1
#保存图片
ggsave('./STEP3_GSEA/box_Immune.pdf', plot = e1,width=15,height = 6)

#######################################
#还缺少相关性热图和散点
M<-re1[,-1]
M<-cor(M)

pdf(file="./STEP3_GSEA/immunecell_corrplot.pdf",width=8,height=6,onefile=FALSE)
corrplot(M, type = 'lower', order = 'hclust', tl.col = 'black', cl.ratio = 0.2, tl.srt = 45)
dev.off()

#散点

coran<-result1[,-2]
coran<-coran[coran$`P-value`<0.05,]
coran<-coran[,-c(24:26)]

p<-list()
for (i in 1:22){
  corinf<-coran[,c(1,i+1)]
  colnames(corinf)<-c('exp','cell')
  p[[i]]<-
  ggplot(data = corinf,aes(x=log2(exp),y=cell))+
    geom_point(shape=21,size=4)+#把点再加大了一倍
    stat_smooth(method = lm)+
    stat_cor(data=corinf,method='pearson',size=6)+
    theme_classic(base_size = 22)+
    theme(aspect.ratio=1)+
    labs(y=paste0(colnames(coran)[i+1]), 
         x="PF4 Expression")
}


pdf("./STEP3_GSEA/PF4_COR_immucell.pdf",width = 28,height = 24)
wrap_plots(p,nrow=4) 
dev.off()

#
