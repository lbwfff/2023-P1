rm(list= ls())
# dir.create('./STEP7_prognosis')

library('forestmodel')
library("survival")
library("survminer")
##############################################################
#
load('TCGA_LAML_ready.RData')

list<-read.csv('genelist.csv')
list<-list$x

survival<-survival[match(colnames(exp),survival$sample),]
array<-t(exp[rownames(exp) %in% list ,])

identical(rownames(array),survival$sample)

merge<-as.data.frame(cbind(array,survival[,c(3,4)]))
  
rownames(merge)<-survival$sample
#

vars_for_table<-c(list)
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(OS.time, OS)~', x))) #单次的cox回归分析循环
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = merge)}) #把循环分析的结果整合在一起

pdf('./STEP7_prognosis/allfactor_forest.pdf',width =12,height = 18)
print(forest_model(model_list = univ_models,covariates = vars_for_table,
                   merge_models =T)) #绘图，forest_model做出来这个还行
dev.off()

#整理为表格
table<-data.frame(gene=c(list),
                  p=c(NA),
                  hz=c(NA))
for (i in 1:58){
  table$p[i]<-summary(univ_models[[i]])$coefficients[, 5]
  table$hz[i]<-exp(coef(univ_models[[i]]))
}

filter<-table[table$p<0.05,]
filter<-filter[order(filter$hz,decreasing = T),]

vars_for_table<-c(filter$gene)
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(OS.time, OS)~', x))) #单次的cox回归分析循环
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = merge)}) #把循环分析的结果整合在一起

pdf('./STEP7_prognosis/filter_forest.pdf',width =8,height = 6)
print(
test<-  forest_model(model_list = univ_models,covariates = vars_for_table,
                   merge_models =T)
  ) #绘图，forest_model做出来这个还行
dev.off()


##########################################################
#LASSO
library(glmnet)

lassoarray<-t(exp[rownames(exp) %in% filter$gene,])
net<-survival[,c(4,3)]
colnames(net)<-c('time','status')
net$time<-net$time+1
identical(colnames(exp),survival$sample)

set.seed(23714)

cv<-cv.glmnet(as.matrix(lassoarray),as.matrix(net),nfold=10,family='cox')
fit.train <- cv$glmnet.fit
pdf(file = "./STEP7_prognosis/lasso1.pdf",width =8,height = 6)
plot(cv,las=1)
dev.off()

fit <- glmnet(as.matrix(lassoarray),as.matrix(net),family = "cox")
pdf(file = "./STEP7_prognosis/lasso2.pdf",width =8,height = 6)
plot(fit,xvar = "lambda",label = TRUE, las=1)
dev.off()

##

coef.min = coef(cv, s = "lambda.min")  
active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
lasso_suv_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] 

lasso_coef<-data.frame(gene=c(lasso_suv_geneids),
                       coef=c(coef.min@x))

#计算得分

merge<-as.data.frame(merge)
cox.data = merge[,colnames(merge) %in% c('OS','OS.time',lasso_coef$gene)]


if (ncol(cox.data) < 3) next
cox.data.step <- na.omit(cox.data[, c("OS.time", "OS", colnames(cox.data)[!(colnames(cox.data) %in% c('OS','OS.time'))])])
fml <- as.formula(paste0("Surv(OS.time,OS)~", paste0(colnames(cox.data)[!(colnames(cox.data) %in% c('OS','OS.time'))],
                                                     collapse = "+")))
f <- coxph(fml, data = cox.data.step, id = rownames(cox.data.step))
cox <- f

riskScore = predict(cox, type = "risk", newdata = cox.data.step)
cox.data.plot <- cbind(cox.data.step, riskScore)
cox.data.plot$OS <- as.numeric(as.character(cox.data.plot$OS))
cox.data.plot$riskScore <- as.numeric(cox.data.plot$riskScore)
cox.data.plot <- cox.data.plot[order(cox.data.plot$riskScore), ]

write.csv(cox.data.plot, file = "./STEP7_prognosis/risk_score.csv")
save(cox,file = './STEP7_prognosis/model.RData')
#########################################################
#时间依赖ROC

library(survivalROC)
heatmap_train <- cox.data.plot

nobs<-nrow(heatmap_train)

pdf(file="./STEP7_prognosis/time_ROC.pdf",height = 6,width = 6)
roc=survivalROC(Stime=heatmap_train$OS.time, status=heatmap_train$OS,lambda=0.0001,
                marker = heatmap_train$riskScore,predict.time =1*365)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#FFC1C1',
     xlab="False positive rate", ylab="True positive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=c()
rocCol <- c('#FFC1C1','#63B8FF','#FA8072','#ADFF2F','#FFFF00')
aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
j =0
for (i in c(2,4)){
  roc1=survivalROC(Stime=heatmap_train$OS.time, status=heatmap_train$OS,lambda=0.0001,
                   marker = heatmap_train$riskScore,predict.time =i*365)
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
abline(0,1)
dev.off()

#生存曲线
matt<-cox.data.plot
med.exp<-median(matt$riskScore)
more.med.exp.index<-which(matt$riskScore>=med.exp)
less.med.exp.index<-which(matt$riskScore< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS.time/30,OS) ~ status, data = matt)
s.diff<-survdiff(Surv(OS.time/30,OS) ~ status, data = matt)

sdata.plot3<-ggsurvplot(s.fit, data=matt,
                        palette="Pastel1",
                        pval = TRUE,pval.method = TRUE, conf.int = TRUE,
                        xlab = 'Time (Month)',ggtheme = theme_survminer(),
                        surv.median.line = 'hv',size=2,
                        title=paste0("TCGA survival"))

pdf("./STEP7_prognosis/survplot.pdf",width = 8,height = 6)
sdata.plot3$plot+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  theme(aspect.ratio=1)
dev.off()
#####################################################
#曲线

plot<-matt[order(matt$riskScore),]
plot$patientid<-1:nrow(plot)

library(ggplot2)
p1=ggplot(plot,aes(x=patientid,y=riskScore))+geom_point(aes(color=status))+
  scale_colour_manual(values = c("#990033","#003366"))+
  theme_bw()+labs(x=" ",y="Risk score")+theme(legend.position = "top")+
  geom_hline(yintercept=median(plot$riskScore),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(plot$status==names(table(plot$status))[1]),colour="black", linetype="dotted",size=0.8)+
  scale_x_continuous(expand = c(0,0))+  #把X轴和Y轴的间隙去除
  scale_y_continuous(expand = c(0,0))
p1

#散点
plot$OS<-factor(plot$OS,levels = c(1,0))
p2=ggplot(plot,aes(x=patientid,y=OS.time/365))+geom_point(aes(col=OS))+theme_bw()+
  scale_colour_manual(values = c("#990033","#003366"))+theme(legend.position = "top")+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(plot$status==names(table(plot$status))[1]),colour="black", linetype="dotted",size=0.8)+
  scale_x_continuous(expand = c(0,0))+  #把X轴和Y轴的间隙去除
  scale_y_continuous(expand = c(0,0))
p2

#热图
library(pheatmap)
mycolors <- colorRampPalette(c("white", "#003366", "#990033"), bias = 1.2)(100)

exp_dat<-plot[,colnames(plot) %in% lasso_coef$gene]
tmp=t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1

lasso_coef<-lasso_coef[order(lasso_coef$coef),]
tmp<-tmp[match(lasso_coef$gene,rownames(tmp)),]

test <- tmp %>% as.data.frame()%>% dplyr::mutate(B=row.names(.)) %>% reshape::melt()

test$B<-factor(test$B,levels = rev(lasso_coef$gene))
test$variable<-factor(test$variable,levels = colnames(tmp))

p3<-ggplot(test,aes(x=variable,y=B,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white")+
  labs(x=NULL, y = NULL)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),
        legend.position = "none")


#拼图

pdf("./STEP7_prognosis/riskscore.pdf",width = 8,height = 8)
cowplot::plot_grid(p1,p2,p3,ncol=1,align = "hv")
dev.off()



