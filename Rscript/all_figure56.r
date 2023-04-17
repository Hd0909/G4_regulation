
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_genes_program/G4_final/")

source('./functions/functions_evolution_5UTR.r')
load(file="./data/all_meth_exp_signal.rdata")
load(file="./data/predict_G4_up_down.rdata")
all_meth=readRDS(file="./data/aver_meth_multiple_pos_continues2.rds")


cor_meth_dhs <- readRDS("./data/aver_meth_dhs_multiple_pos_cor_all.rds")
all_cg=readRDS(file="./data/aver_cg_multiple_pos_continues2.rds")
all_dhs=readRDS(file="./data/aver_dhs_multiple_pos_continues2.rds")

cor_all<-readRDS(file="./data/aver_meth_multiple_pos_cor_all.rds")

cor_exp<-readRDS(file="./data/aver_exp_multiple_pos_cor_all.rds")
cor_exp_meth<-readRDS(file="./data/aver_exp_meth_multiple_pos_cor_all.rds")
cor_all_dhs<-readRDS(file="./data/aver_dhs_multiple_pos_cor_all.rds")
cor_exp_dhs <- readRDS("./data/aver_exp_dhs_multiple_pos_cor_all.rds")

all_meth$pos1=str_split_fixed(all_meth$pos," ",2)[,1]
temp1=str_split_fixed(all_meth$pos1,"_",2)
neg=ifelse(as.numeric(temp1[,2])>as.numeric(temp1[,1]),1,-1)
all_meth$position<-(as.numeric(temp1[,1])/2+as.numeric(temp1[,2])/2+50)*neg

all_cg$pos1=str_split_fixed(all_cg$pos," ",2)[,1]
temp1=str_split_fixed(all_cg$pos1,"_",2)
neg=ifelse(as.numeric(temp1[,2])>as.numeric(temp1[,1]),1,-1)
all_cg$position<-(as.numeric(temp1[,1])/2+as.numeric(temp1[,2])/2+50)*neg


all_dhs$pos1=str_split_fixed(all_dhs$pos," ",2)[,1]
temp1=str_split_fixed(all_dhs$pos1,"_",2)
neg=ifelse(as.numeric(temp1[,2])>as.numeric(temp1[,1]),1,-1)
all_dhs$position<-(as.numeric(temp1[,1])/2+as.numeric(temp1[,2])/2+50)*neg


load(file="./data/predict_G4_up_down.rdata")

x1=all_meth[all_meth$position<=1000&all_meth$position>=-1500 & all_meth$cancer=="LUAD" ,]
x1$type_2=ifelse(x1$ensembl_gene_id %in% promoter_G4_down_pos$trans_ID[promoter_G4_down_pos[,9]<=100],"+","-")
p_preG4<-ggplot(x1,aes(x=as.factor(position),y=norm_meth,fill=type_2))+geom_boxplot()+theme_hd()+xlab("Distance of methylation regions to TSS")+ylab("DNA methylation")+theme(axis.text.x = element_text(angle = 90))+scale_fill_manual("Predicted G4s in NT strand of TSS downstream 100bp",values =   colorRampPalette(c("white" ,"#F07F16C7"))(6)[c(2,6)])+ggtitle("TCGA-LUAD")

x1_dhs=all_dhs[all_dhs$position<=1000&all_dhs$position>=-1500 & all_dhs$cancer=="LUAD" ,]
x1_dhs$type_2=ifelse(x1_dhs$genes_main_transcript %in% promoter_G4_down_pos$trans_ID[promoter_G4_down_pos[,9]<=100],"+","-")
p_preG4_dhs<-ggplot(x1_dhs,aes(x=as.factor(position),y=signal,fill=type_2))+geom_boxplot()+theme_hd()+xlab("Distance of DHS signal regions to TSS")+ylab("DHS signal")+theme(axis.text.x = element_text(angle = 90))+scale_fill_manual("Predicted G4s in NT strand of TSS downstream 100bp",values =   colorRampPalette(c("white" ,"#F07F16C7"))(6)[c(2,6)])+ggtitle("ENCODE-Lung")


p1=ggplot(all_meth[all_meth$position<=1000&all_meth$position>=-1500 & all_meth$cancer=="LUAD" ,],aes(x=as.factor(position),y=norm_meth))+geom_boxplot()+theme_hd()+xlab("Distance to TSS")+ylab("DNA methylation")+theme(axis.text.x = element_text(angle = 90))

colnames(all_meth)<-gsub("ensembl_gene_id","ensembl_trans_id",colnames(all_meth))


all_cg$pos1=str_split_fixed(all_cg$pos," ",2)[,1]
temp1=str_split_fixed(all_cg$pos1,"_",2)
neg=ifelse(as.numeric(temp1[,2])>as.numeric(temp1[,1]),1,-1)
all_cg$position<-(as.numeric(temp1[,1])/2+as.numeric(temp1[,2])/2+50)*neg

all_cg_melt<-melt(all_cg,id.vars = c("pos1","pos","position","ensembl_gene_id","ensembl_trans_id"))
all_cg_melt$value<-as.numeric(all_cg_melt$value)
p1_cg=ggplot(all_cg_melt[all_cg_melt$position<=1000&all_cg_melt$position>=-1500  ,],aes(x=as.factor(position),y=value,fill=variable))+geom_boxplot()+theme_hd()+xlab("Distance to TSS")+ylab("Frequency of G-runs")+theme(axis.text.x = element_text(angle = 90))+scale_fill_manual(name="Strand",labels=c("Non-template","Template"),values=c( "#F07F16C7","#40A4DE"))
#all_cg_melt$position<-as.factor(all_cg_melt$position)
p1_cg_sig=plotSubSigcg(all_cg_melt[all_cg_melt$position<=1000&all_cg_melt$position>=-1500  ,],"position", "value","variable","", cols1 =  c("#F07F16C7","#40A4DE"),"","Frequency of G-runs") + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1), legend.position = "none")+xlab("Distance to TSS")

p1_dhs=ggplot(all_dhs[all_dhs$position<=1000&all_dhs$position>=-1500 & all_dhs$cancer=="LUAD" ,],aes(x=as.factor(position),y=signal))+geom_boxplot()+theme_hd()+xlab("Distance to TSS")+ylab("DHS signal")+theme(axis.text.x = element_text(angle = 90))+ggtitle("ENCODE-Lung")

p1_meth=ggplot(all_meth[all_meth$position<=1000&all_meth$position>=-1500 & all_meth$cancer=="LUAD" ,],aes(x=as.factor(position),y=norm_meth))+geom_boxplot()+theme_hd()+xlab("Distance to TSS")+ylab("DNA methylation")+theme(axis.text.x = element_text(angle = 90))+ggtitle("TCGA-LUAD")



p1=ggarrange(p1_cg_sig,p1_meth,p1_dhs,nrow=1,ncol=3)

pcor_meth_gg=ggplot(cor_all[cor_all$pos_gg==cor_all$pos_meth& cor_all$position_meth<=1000 &cor_all$position_meth>=-1500,],aes(x=as.factor(position_gg),y=cor))+geom_boxplot(aes(fill=type_x,alpha=I(0.1)),show.legend=F)+geom_jitter(width = 0.25,aes(col=type_x))+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of methylation \nand G-run Frequency")+xlab("Distance of the methylated and Gruns regions to TSS")+labs(col="Strand of the G-runs")+ggtitle(" G-run  VS methylation in the same region")+scale_color_manual(values = c("#F07F16C7","#40A4DE"))

pcor_signal_gg=ggplot(cor_all_dhs[cor_all_dhs$pos_gg==cor_all_dhs$pos_meth& cor_all_dhs$position_meth<=1000 &cor_all_dhs$position_meth>=-1500,],aes(x=as.factor(position_gg),y=cor))+geom_boxplot(aes(fill=type_x,alpha=I(0.1)),show.legend=F)+geom_jitter(width = 0.25,aes(col=type_x))+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of DHS signal \nand G-run Frequency")+xlab("Distance of the DHS signal and Gruns regions to TSS")+labs(col="Strand of the G-runs")+ggtitle(" G-run  VS DHS signal in the same region")+scale_color_manual(values = c("#F07F16C7","#40A4DE"))


pcor_signal_meth=ggplot(cor_meth_dhs[cor_meth_dhs$pos_meth==cor_meth_dhs$pos_signal & cor_meth_dhs$position_other<=1000 &cor_meth_dhs$position_other>=-1500,],aes(x=as.factor(position_other),y=cor))+geom_boxplot()+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+ylab("Correlation of methylation  and DHS signal")+xlab("Distance of the DHS regions to TSS")+ggtitle("Expression  VS DHS signal in different regions")#+geom_hline(yintercept = c(0),linetype="dashed")

p_cor_exp_dhs=ggplot(cor_exp_dhs[cor_exp_dhs$position_other<=1000 &cor_exp_dhs$position_other>=-1500,],aes(x=as.factor(position_other),y=cor))+geom_boxplot()+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of expression and DHS signal")+xlab("Distance of the DHS regions to TSS")+ggtitle("Expression  VS DHS signal in different regions")

p_cor_exp=ggplot(cor_exp[cor_exp$position_gg<=1000 &cor_exp$position_gg>=-1500,],aes(x=as.factor(position_gg),y=cor))+geom_boxplot(aes(fill=type_x,alpha=I(0.1)),show.legend=F)+geom_jitter(width = 0.25,aes(col=type_x))+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of expression \nand G-run Frequency")+xlab("Distance of the Gruns regions to TSS")+ggtitle("Expression  VS Gruns in different regions")+scale_color_manual(values = c("#F07F16C7","#40A4DE"))+labs(col="Strand of the G-runs")

p_cor_exp_meth=ggplot(cor_exp_meth[cor_exp_meth$position_meth<=1000 &cor_exp_meth$position_meth>=-1500,],aes(x=as.factor(position_meth),y=cor))+geom_boxplot()+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+ylab("Correlation of expression and DHS signal")+xlab("Distance of the methylation regions to TSS")+ggtitle("Expression  VS methylation in different regions")


p_cor_exp_dhs=ggplot(cor_exp_dhs[cor_exp_dhs$position_other<=1000 &cor_exp_dhs$position_other>=-1500,],aes(x=as.factor(position_other),y=cor))+geom_boxplot()+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of expression and DHS signal")+xlab("Distance of the DHS regions to TSS")+ggtitle("Expression  VS DHS signal in different regions")
#+geom_hline(yintercept = c(0),linetype="dashed")




pcor_meth_gg100=ggplot(cor_all[cor_all$pos_gg=="0_100"& cor_all$position_meth<=1000 &cor_all$position_meth>=-1500,],aes(x=as.factor(position_meth),y=cor))+geom_boxplot(aes(fill=type_x,alpha=I(0.1)),show.legend=F)+geom_jitter(width = 0.25,aes(col=type_x))+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of methylation and\n G-run Frequency")+xlab("Distance of the methylated regions to TSS")+labs(col="Strand of the G-runs")+ggtitle("G-run in TSS downstream 100bp  VS\n methylation in different regions")+scale_color_manual(values = c("#F07F16C7","#40A4DE"))


pcor_signal_gg100=ggplot(cor_all_dhs[cor_all_dhs$pos_gg=="0_100"& cor_all_dhs$position_meth<=1000 &cor_all_dhs$position_meth>=-1500,],aes(x=as.factor(position_meth),y=cor))+geom_boxplot(aes(fill=type_x,alpha=I(0.1)),show.legend=F)+geom_jitter(width = 0.25,aes(col=type_x))+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of DHS signal\n and G-run Frequency")+xlab("Distance of the DHS signal regions to TSS")+labs(col="Strand of the G-runs")+ggtitle(" G-run in TSS downstream 100bp  VS\n DHS signal in different regions")+scale_color_manual(values = c("#F07F16C7","#40A4DE"))

pcor_signal_meth100=ggplot(cor_meth_dhs[cor_meth_dhs$pos_meth=="0_100" & cor_meth_dhs$position_otherdhs<=1000 &cor_meth_dhs$position_otherdhs>=-1500,],aes(x=as.factor(position_otherdhs),y=cor))+geom_boxplot()+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+ylab("Correlation of methylation\n and DHS signal")+xlab("Distance of the DHS regions to TSS")+ggtitle("Methylation in TSS downstream 100bp VS\n DHS signal in different regions")#+geom_hline(yintercept = c(0),linetype="dashed")

pcor_signal_methsignal100=ggplot(cor_meth_dhs[cor_meth_dhs$pos_signal=="0_100" & cor_meth_dhs$position_other<=1000 &cor_meth_dhs$position_other>=-1500,],aes(x=as.factor(position_other),y=cor))+geom_boxplot()+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+ylab("Correlation of methylation\n and DHS signal")+xlab("Distance of the methylation regions to TSS")+ggtitle("DHS signal in TSS downstream 100bp VS\n Methylation in different regions")#+geom_hline(yintercept = c(0),linetype="dashed")


p0=ggarrange(p1_meth,p1_dhs,nrow=2,ncol=1,labels=c(letters[2:3]),font.label = list(size = 18, color = "black", face = "bold"))

p1=ggarrange(p1_cg_sig,p0,p_cor_exp,nrow=1,ncol=3,labels=c(letters[1],"",letters[4]),font.label = list(size = 18, color = "black", face = "bold"))

p2=ggarrange(pcor_meth_gg,pcor_signal_gg,pcor_meth_gg100,nrow=1,ncol=3,labels=c(letters[5:7]),font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend = "bottom")

p3=ggarrange(pcor_signal_gg100,p_preG4,p_preG4_dhs,nrow=1,ncol=3,labels=c(letters[8:10]),font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend = "bottom")
tiff("./all_figures/all_figure5.tiff",width = 20,height = 23,res=300,units="in",compression = "lzw")

print(ggarrange(p1,p2,p3,nrow=3,ncol=1,heights = c(1,1,1)))
dev.off()


library(raster)
library(ggplot2)




x1_lable=c("DNA methyltransferase","RNA polymerase","Transcription\ninitiation complex","transcription factor","Activator","Enhancer","Template strand","Non-template strand","Pre-mRNA","Methylated CpG","Unmethylated CpG","Repression")

img <- readPNG("./all_figures/model1.png")
g1 <- rasterGrob(img, interpolate=TRUE)
s1=qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g1, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(aes(col=I("white")))+theme_void()
img <- readPNG("./all_figures/model2.png")
g2 <- rasterGrob(img, interpolate=TRUE)
s2=qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(aes(col=I("white")))+theme_void()




s1=s1+annotate("text",x=8,y=c(seq(1.5,7,1.6),7.5,8.8),size=5,label=x1_lable[6:1],fontface="bold",hjust=0)
s2=s2+annotate("text",x=8,y=c(seq(2.2,9,1.55),9.6),size=5,label=x1_lable[12:7],fontface="bold",hjust=0)



p5=ggarrange(s1,s2,nrow=2,ncol=1,font.label = list(size = 18, color = "black", face = "bold"))

tiff("./all_figures/all_figure6.tiff",width = 13,height = 16,res=300,units="in",compression = "lzw")

print(p5)
dev.off()


p1_dhsall=ggplot(unique(all_dhs[all_dhs$position<=1000&all_dhs$position>=-1500&all_dhs$tissue!="lung"  ,c("tissue","position","signal")]),aes(x=as.factor(position),y=signal))+geom_boxplot()+xlab("Distance to TSS")+ylab("DHS signal")+theme(axis.text.x = element_text(angle = 90))+facet_wrap(~tissue)+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90,hjust=1))

p1_methall=ggplot(all_meth[all_meth$position<=1000&all_meth$position>=-1500 &all_meth$cancer!="LUAD" ,],aes(x=as.factor(position),y=norm_meth))+geom_boxplot()+xlab("Distance to TSS")+ylab("DNA methylation")+theme(axis.text.x = element_text(angle = 90))+facet_wrap(~cancer,nrow=4,ncol=3)+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90,hjust=1))

tiff("./all_figures/supply_all_figure6.tiff",width = 18,height = 23,res=300,units="in",compression = "lzw")
print(p1_methall)
dev.off()

tiff("./all_figures/supply_all_figure7.tiff",width = 16,height = 23,res=300,units="in",compression = "lzw")
print(p1_dhsall)
dev.off()






temp=str_split_fixed(cor_all$pos_gg,"_",2)
cor_all$pos_gg_label=paste("Grun in ",ifelse(as.numeric(temp[,2])>as.numeric(temp[,1]),paste("TSS downstream\n",temp[,1],"-",temp[,2],"bp"),paste("TSS upstream\n",temp[,2],"-",temp[,1],"bp")),sep="")
x=unique(cor_all[,c("position_gg","pos_gg_label")])
x=x[order(x[,1]),]
cor_all$pos_gg_label<-factor(cor_all$pos_gg_label,levels = x[,2])
pcor_meth_ggall=ggplot(cor_all[ cor_all$position_meth<=1000 &cor_all$position_meth>=-1500 & cor_all$position_gg <=1000 &cor_all$position_gg>=-1500 ,],aes(x=as.factor(position_meth),y=cor))+geom_boxplot(aes(fill=type_x,alpha=I(0.1)),show.legend=F)+geom_jitter(width = 0.25,aes(col=type_x))+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of methylation and\n G-run Frequency")+xlab("Distance of the methylated regions to TSS")+labs(col="Strand of the G-runs")+scale_color_manual(values = c("#F07F16C7","#40A4DE"))+facet_wrap(~pos_gg_label,nrow=6)

tiff("./all_figures/supply_all_figure8.tiff",width = 20,height = 24,res=300,units="in",compression = "lzw")
print(pcor_meth_ggall)
dev.off()

temp=str_split_fixed(cor_all_dhs$pos_gg,"_",2)
cor_all_dhs$pos_gg_label=paste("Grun in ",ifelse(as.numeric(temp[,2])>as.numeric(temp[,1]),paste("TSS downstream\n",temp[,1],"-",temp[,2],"bp"),paste("TSS upstream\n",temp[,2],"-",temp[,1],"bp")),sep="")
x=unique(cor_all_dhs[,c("position_gg","pos_gg_label")])
x=x[order(x[,1]),]
cor_all_dhs$pos_gg_label<-factor(cor_all_dhs$pos_gg_label,levels = x[,2])

pcor_signal_ggall=ggplot(cor_all_dhs[ cor_all_dhs$position_meth<=1000 &cor_all_dhs$position_meth>=-1500& cor_all_dhs$position_gg <=1000 &cor_all_dhs$position_gg>=-1500,],aes(x=as.factor(position_meth),y=cor))+geom_boxplot(aes(fill=type_x,alpha=I(0.1)),show.legend=F)+geom_jitter(width = 0.25,aes(col=type_x))+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90))+geom_hline(yintercept = c(0),linetype="dashed")+ylab("Correlation of DHS signal\n and G-run Frequency")+xlab("Distance of the DHS signal regions to TSS")+labs(col="Strand of the G-runs")+scale_color_manual(values = c("#F07F16C7","#40A4DE"))+facet_wrap(~pos_gg_label,nrow=6)


tiff("./all_figures/supply_all_figure9.tiff",width = 20,height = 24,res=300,units="in",compression = "lzw")
print(pcor_signal_ggall)
dev.off()

x1=all_meth[all_meth$position<=1000&all_meth$position>=-1500 &all_meth$cancer!="LUAD" ,]
x1$type_2=ifelse(x1$ensembl_gene_id %in% promoter_G4_down_pos$trans_ID[promoter_G4_down_pos[,9]<=100],"+","-")
p_preG4All<-ggplot(x1,aes(x=as.factor(position),y=norm_meth,fill=type_2))+geom_boxplot()+xlab("Distance of methylation regions to TSS")+ylab("DNA methylation")+theme(axis.text.x = element_text(angle = 90))+scale_fill_manual("Predicted G4s in NT strand of TSS downstream 100bp",values =   colorRampPalette(c("white" ,"#F07F16C7"))(6)[c(2,6)])+facet_wrap(~cancer,nrow=4,ncol=3)+theme_hd_minimal2()+theme(axis.text.x = element_text(angle = 90,hjust=1))



x1_dhs=unique(all_dhs[all_dhs$position<=1000&all_dhs$position>=-1500&all_dhs$tissue!="lung" ,c("signal","position","tissue","ensembl_transcript_id")])
x1_dhs$type_2=ifelse(x1_dhs$ensembl_transcript_id %in% promoter_G4_down_pos$trans_ID[promoter_G4_down_pos[,9]<=100],"+","-")
p_preG4_dhsall<-ggplot(x1_dhs,aes(x=as.factor(position),y=signal,fill=type_2))+geom_boxplot()+xlab("Distance of DHS signal regions to TSS")+ylab("DHS signal")+theme(axis.text.x = element_text(angle = 90))+scale_fill_manual("Predicted G4s in NT strand of TSS downstream 100bp",values =   colorRampPalette(c("white" ,"#F07F16C7"))(6)[c(2,6)])+facet_wrap(~tissue)+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 90,hjust=1))

tiff("./supply_all_figure10.tiff",width = 18,height = 23,res=300,units="in",compression = "lzw")
print(p_preG4All)
dev.off()
tiff("./supply_all_figure11.tiff",width = 18,height = 23,res=300,units="in",compression = "lzw")
print(p_preG4_dhsall)
dev.off()



