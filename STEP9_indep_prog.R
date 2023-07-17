rm(list= ls())
# dir.create('./STEP9_indep_prog')

library('forestmodel')
library("survival")
library("survminer")
##############################################################

load('TCGA_LAML_ready.RData')

list<-c('sampleID','gender','age_at_initial_pathologic_diagnosis',
        'acute_myeloid_leukemia_calgb_cytogenetics_risk_category',
        'history_of_neoadjuvant_treatment')
pheno<-as.data.frame(pheno)
multi<-pheno[,colnames(pheno) %in% list]

risk<-read.csv("./STEP7_prognosis/risk_score.csv")

multi<-multi[match(risk$X,multi$sampleID),]
identical(risk$X,multi$sampleID)

multi$riskScore<-risk$riskScore
multi$OS<-risk$OS
multi$OS.time<-risk$OS.time

colnames(multi)[c(2:5)]<-c('CytogeneticsRisk','Age','Gender','NeoadjuvantTreatment')
multi$CytogeneticsRisk<-factor(multi$CytogeneticsRisk,levels = c('Favorable','Intermediate/Normal','Poor'))
multi$NeoadjuvantTreatment<-factor(multi$NeoadjuvantTreatment,levels = c('No','Yes'))

#COX，单因素
rownames(multi)<-multi$sampleID
multi<-multi[,-1]

vars_for_table<-colnames(multi)[1:5]
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(OS.time, OS)~', x))) #单次的cox回归分析循环
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = multi)}) #把循环分析的结果整合在一起

pdf('./STEP9_indep_prog/sig_COX.pdf',width =10,height = 6)
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)) 
dev.off()

#多因素

res.cox<-coxph(Surv(OS.time,OS)~CytogeneticsRisk+Age+Gender+NeoadjuvantTreatment+riskScore, data=multi)
summary(res.cox)

pdf('./STEP9_indep_prog/mul_COX.pdf',width =10,height = 6)
forest_model(res.cox)
dev.off()

##########################################
#列线图

library(rms)

dd<-datadist(multi)
options(datadist="dd")
options(na.action="na.delete")

coxpbc<-cph(formula = Surv(OS.time,OS) ~  CytogeneticsRisk+Age+Gender+NeoadjuvantTreatment+riskScore,data=multi,
            x=T,y=T,surv = T,na.action=na.delete)

print(coxpbc)
surv<-Survival(coxpbc) 
surv3<-function(x) surv(365,x)
surv4<-function(x) surv(730,x) #如果时间超过数据里最高的生存时间就会报错
surv5<-function(x) surv(1460,x)

x<-nomogram(coxpbc,fun = list(surv3,surv4,surv5),lp=T,
            funlabel = c('1-year survival Probability','2-year survival Probability','4-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("./STEP9_indep_prog/nomogram.pdf",width = 12,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

#校准曲线

f1<-cph(formula = Surv(OS.time,OS) ~  CytogeneticsRisk+Age+Gender+NeoadjuvantTreatment+riskScore,data=multi,
        x=T,y=T,surv = T,na.action=na.delete,time.inc = 365)

cal1<-calibrate(f1, cmethod="KM", method="boot",u=365,m=50,B=1000) 

f3<-cph(formula = Surv(OS.time,OS) ~  CytogeneticsRisk+Age+Gender+NeoadjuvantTreatment+riskScore,data=multi,
        x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 
cal3<-calibrate(f3, cmethod="KM", method="boot",u=730,m=50,B=1000)

f5<-cph(formula = Surv(OS.time,OS) ~  CytogeneticsRisk+Age+Gender+NeoadjuvantTreatment+riskScore,data=multi,
        x=T,y=T,surv = T,na.action=na.delete,time.inc = 1460) 
cal5<-calibrate(f5, cmethod="KM", method="boot",u=1460,m=50,B=1000)


#
pdf("./STEP9_indep_prog/calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 1,errbar.col = c("#f7aa58"),
     xlim = c(0,1),ylim= c(0,1),col = c("#f7aa58"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#f7aa58"), pch = 16)

plot(cal5,lwd = 2,lty = 1,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year",'4-year'), #图例文字
       col =c("#2166AC",'#f7aa58',"#B2182B"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()

#############################################################
#DCA

library(ggDCA)

m1 <- cph(formula = Surv(OS.time,OS) ~  CytogeneticsRisk,data=multi,
          x=T,y=T,surv = T,na.action=na.delete)
m2 <- cph(formula = Surv(OS.time,OS) ~  riskScore,data=multi,
          x=T,y=T,surv = T,na.action=na.delete)
m3 <- cph(formula = Surv(OS.time,OS) ~  CytogeneticsRisk++riskScore,data=multi,
          x=T,y=T,surv = T,na.action=na.delete)

data  <- dca(m1,m2,m3,
             model.names =c('CytogeneticsRisk','riskScore','CytogeneticsRisk+riskScore'))


pdf("./STEP9_indep_prog/DCA.pdf",width = 10,height = 8)
ggplot(data,linetype = T,lwd=1.5,
       color = c(MetBrewer::met.brewer("Manet",n=5)))+
  coord_fixed() +theme_classic(base_size = 18)+theme(aspect.ratio=0.85)
dev.off()




