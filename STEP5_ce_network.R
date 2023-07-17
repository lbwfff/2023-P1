rm(list= ls())

########################################################

targetscan<-fread('TargetScan_8.0_ENST00000296029.3_predicted_targeting_details.csv')
targetscan<-targetscan[targetscan$`context++ score percentile`>95,]

milnc<-fread('human.all6.lnc.mi')
milnc<-milnc[milnc$clipExpNum>1,]
milnc<-milnc[milnc$miRNAname %in% targetscan$miRNA]
milnc<-milnc[!duplicated(milnc$geneName),]

targetscan<-targetscan[targetscan$miRNA %in% milnc$miRNAname,]

network<-as.data.frame(array(NA,c(1,4)))
colnames(network)<-c('name1','name2','inter','type')

mi_m<-data.frame(name1=c(targetscan$miRNA),
                 name2=c('PF4'),
                 inter=c(paste0(targetscan$miRNA,'_PF4')),
                 type=c('mi_m'))

mi_lnc<-data.frame(name1=c(milnc$miRNAname),
                 name2=c(milnc$geneName),
                 inter=c(paste0(milnc$miRNAname,'_',milnc$geneName)),
                 type=c('mi_lnc'))

merge<-rbind(mi_m,mi_lnc)

# dir.create('./STEP5_ce_network')
write.table(merge,file = './STEP5_ce_network/network_ce.txt',sep = '\t',quote = F,row.names = F)

node<-data.frame(node=c('PF4',unique(targetscan$miRNA),
                        unique(milnc$geneName)),
                 attribute=c('mRNA',rep('miRNA',length(unique(targetscan$miRNA))),
                             rep('lncRNA',length(unique(milnc$geneName)))
                   
                 ))

write.table(node,file = './STEP5_ce_network/network_node.txt',sep = '\t',quote = F,row.names = F)



