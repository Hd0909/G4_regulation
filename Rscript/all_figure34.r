
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_genes_program/G4_final/")

source('functions/functions_evolution_5UTR.r')
load(file="./data/all_meth_exp_signal.rdata")


meth_exp<-merge(aver_downstream_meth$aver_downstream_meth,all_data_tsg,by=c("cancer","ensembl_gene_id"))



meth2=PlotMutipleCor(aver_meth_gene,"CC_density","GG_density","Template strand","Non-template strand","norm_meth","Correlation of methylation and frequency\nof G-runs in different strands in TSS downstream","pearson")
#+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes = list(size=5)))

meth_cg<-as.data.frame(rbind(cbind(aver_meth_gene[aver_meth_gene$cancer=="LUAD",c("G4_predict_type","norm_meth")],value=aver_meth_gene[aver_meth_gene$cancer=="LUAD","CC_density"],grun="T"),cbind(aver_meth_gene[aver_meth_gene$cancer=="LUAD",c("G4_predict_type","norm_meth")],value=aver_meth_gene[aver_meth_gene$cancer=="LUAD","GG_density"],grun="NT")))
meth_cg$type_all<-paste("G-run-",meth_cg$grun,sep="")
sp1_meth1=plot_scatter_2(meth_cg,"value","norm_meth","Frequency of G-runs in TSS downstream", "Methylation in TSS downstream",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+theme_hd_minimal()+ggtitle("TCGA-LUAD Normal")+facet_wrap(~type_all,nrow=1,ncol=2,scales = "free_y")+ylim(0,1.25)




cor_age_meth1<-get_cor_agetype(aver_meth_gene,"norm_meth","CC_density","age_type2")
cor_age_meth2<-get_cor_agetype(aver_meth_gene,"norm_meth","GG_density","age_type2")


all_cor_meth_age<-rbind(cbind(strand="G-run-T",cor_age_meth1),cbind(strand="G-run-NT",cor_age_meth2))
p_cor_meth=ggplot(data=all_cor_meth_age,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd_minimal()+xlab("Gene Age Groups (million years)")+theme(axis.text.x =element_text(angle=15,hjust=1))+guides(color=guide_legend(nrow=1,byrow=TRUE))+ylab("Pearson's r(methylation VS Frequency\nof G-runs)")+ggtitle("TCGA Datasets Normal samples")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))+facet_wrap(~strand,nrow=1,scales = "free_y")+ylim(-0.5,0)

cor_age_meth_exp1<-get_cor_agetype(meth_exp,"norm_meth","norm_exp","age_type2")

ggplot(data=cor_age_meth_exp1,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age Groups\n(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+ylab("Pearson's r(Expression VS methylation)")+ggtitle("TCGA-LUAD normal")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))

genes_chromatin$cancer<-gsub(" ","\n",genes_chromatin$cancer)
signal=PlotMutipleCor(genes_chromatin,"CC_density","GG_density","Template","Non-template","signal","Correlation of DHS signal and frequency of\nG-runs in different strand in TSS downstream","pearson")+ylab("ENCODE Datasets normal samples")#+guides(color=guide_legend(nrow=1,byrow=TRUE,override.aes = list(size=5)))+theme(legend.position = "bottom",legend.box ="vertical")

cor_age_signal1<-get_cor_agetype(genes_chromatin,"signal","GG_density","age_type2")
cor_age_signal2<-get_cor_agetype(genes_chromatin,"signal","CC_density","age_type2")

all_cor_signal_age<-rbind(cbind(strand="Non-template",cor_age_signal1),cbind(strand="Template",cor_age_signal2))
p_signal_age<-ggplot(data=all_cor_signal_age,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd_minimal()+xlab("Gene Age Groups (million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(color=guide_legend(nrow=1,byrow=TRUE))+ylab("Pearson's r(DHS signal\nVS Frequency of G-runs)")+ggtitle("ENCODE Datasets normal samples")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))+facet_wrap(~strand,nrow=1,scales = "free_y")+ylim(-0.2,0.4)



signal_cg<-as.data.frame(rbind(cbind(genes_chromatin[genes_chromatin$tissue=="lung",c("G4_predict_type","signal")],value=genes_chromatin[genes_chromatin$tissue=="lung","CC_density"],grun="T"),cbind(genes_chromatin[genes_chromatin$tissue=="lung",c("G4_predict_type","signal")],value=genes_chromatin[genes_chromatin$tissue=="lung","GG_density"],grun="NT")))
signal_cg$type_all<-paste("G-run-",signal_cg$grun,sep="")
sp1_signal1=plot_scatter_2(signal_cg,"value","signal","Frequency of G-runs in TSS downstream", "DHS signal in TSS downstream",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+theme_hd_minimal()+ggtitle("ENCODE-Lung")+facet_wrap(~type_all,nrow=1,ncol=2,scales = "free_y")+ylim(0,120)








sp1signal_cc=plot_scatter_2(genes_chromatin[genes_chromatin$tissue=="lung",],"CC_density","signal","Frequency of G-runs\nin TSS downstream", "DHS signal",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+ggtitle("ENCODE-Lung\nTemplate")+theme_hd()+ylim(0,120)
print(sp1signal_cc)


sp1signal_gg=plot_scatter_2(genes_chromatin[genes_chromatin$tissue=="lung",],"GG_density","signal","Frequency of G-runs\nin TSS downstream", "DHS signal",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+theme_hd()+ylim(0,120)+ggtitle("ENCODE-Lung\nNon-template")
print(sp1signal_gg)

pcor_meth_exp<-plotForest(meth_exp,"norm_exp","norm_meth","Correlation of  gene expression and\n methylation in the TSS downstream",0.5,1.3,"0.5cm")


pcor_signal_exp<-plotForest(all_data_tsg[!is.na(all_data_tsg$DHS_signal),],"norm_exp","DHS_signal","Correlation of  gene expression and\n DHS_signal in the TSS downstream",0.5,1.7,"0.5cm")

pcor_signal_meth<-plotForest(aver_meth_gene[!is.na(aver_meth_gene$DHS_signal),],"norm_meth","DHS_signal","Correlation of  DNA methylation and\n DHS_signal in the TSS downstream",0.5,1.7,"0.5cm")




pgg_1=ggarrange(meth2,sp1_meth1,nrow=1,ncol=2,labels = c(letters[1:2])  ,font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend = "none",widths = c(1,2))
pgg_2=ggarrange(signal,sp1_signal1,nrow=1,ncol=2,labels = c(letters[3:4])  ,font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,2),common.legend = T,legend = "bottom")
pgg_12<-ggarrange(pgg_1,pgg_2,nrow=2,ncol=1,common.legend = T,legend = "bottom")
pgg_3=ggarrange(pcor_meth_exp,pcor_signal_exp,pcor_signal_meth,nrow=1,ncol=3,labels = c(letters[5:7])  ,font.label = list(size = 18, color = "black", face = "bold"))

tiff("./all_figures/all_figure3.tiff",width = 22,height = 22,res=300,units="in",compression = "lzw")
print(ggarrange(pgg_12,pgg_3,nrow=2,ncol=1,heights = c(2,1)))
dev.off()



p1_g4=ggplot(genes_infor,aes(x=G4_down_predict_neg,y=CC_density))  + geom_signif(comparisons = getComp(c("+","-")),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),outlier.colour = NA,aes(fill=G4_down_predict_neg))+ylab("Frequency of G-runs")+theme_hd_minimal() + guides(fill="none")+scale_fill_manual(values = colorRampPalette(c("white" ,"#40A4DE"))(6)[c(2,6)])+ggtitle("Templated strand of \nTSS downstream")+xlab("Predicted G4")

p2_g4=ggplot(genes_infor,aes(x=G4_down_predict_pos,y=GG_density))  + geom_signif(comparisons = getComp(c("+","-")),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),outlier.colour = NA,aes(fill=G4_down_predict_pos))+ylab("Frequency of G-runs")+theme_hd_minimal() + guides(fill="none")+scale_fill_manual(values = colorRampPalette(c("white" ,"#F07F16C7"))(6)[c(2,6)])+xlab("Predicted G4")+ggtitle("Non-templated strand\nof TSS downstream")


ggarrange(p1_g4,p2_g4)



x=list(Predicted_G4_T=genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_neg=="+"],`G4-seq_validated_G4`=genes_infor$ensembl_gene_id[genes_infor$G4_down_exp=="+"])
title_x="Overlap of G4s in TSS downstream"

p=venn.diagram(x,fill = c( "#00B2EE","#104E8B"),
               alpha = c(0.5, 0.5), cex = 1,cat.fontface = "bold",lty =2, fontfamily =3,cat.cex=1,sub.fontface="bold",margin=0.1,sub="", sub.pos = c(0.3,-0.05),filename = NULL,main=title_x,main.fontface ="bold",  main.cex = 1.2,col=NA,	ext.text = TRUE,ext.line.lwd = 2, ext.dist = -0.15,ext.length = 0.9,ext.pos = -14,rotation.degree = 185,cat.pos=c(120,240))
pVNN2<-ggarrange(p)
pVNN2

x=list(Predicted_G4_NT=genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_pos=="+"],`G4-seq_validated_G4`=genes_infor$ensembl_gene_id[genes_infor$G4_down_exp=="+"])
title_x="Overlap of G4s in TSS downstream"

p=venn.diagram(x,fill = c( "#FFC125","#104E8B"),
               alpha = c(0.5, 0.5), cex = 1,cat.fontface = "bold",lty =2, fontfamily =3,cat.cex=1,sub.fontface="bold",margin=0.1,sub="", sub.pos = c(0.3,-0.05),filename = NULL,main=title_x,main.fontface ="bold",  main.cex = 1.2,col=NA,	ext.text = TRUE,ext.line.lwd = 2, ext.dist = -0.15,ext.length = 0.9,ext.pos = -14,rotation.degree = 5)
pVNN<-ggarrange(p)
pVNN







df_all<-as.data.frame(table(genes_infor$G4_predict_type,genes_infor$Gene.Type,genes_infor$age_type2))
p_freq=ggplot(df_all, aes(x = Var3, y = Freq, fill = Var1)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))+ylab("Frequency")+theme_hd_minimal2()+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "bottom",axis.title.x =element_blank() )+guides(fill="none")+ggtitle("Genes with predicted G4s in different strand")

types=sort(as.character(unique(genes_infor$G4_predict_type)))
p2_exp=ggplot(data = all_data_tsg[all_data_tsg$cancer=="LUAD"& !is.na(all_data_tsg$Gene.Type),],aes(y=norm_exp,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Expression ")+theme_hd_minimal()+ggtitle("TCGA-LUAD normal")+theme(axis.text.x =element_blank(),legend.position = "none",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))


p3_signal=ggplot(data = genes_chromatin[genes_chromatin$cancer=="lung"& !is.na(genes_chromatin$Gene.Type),],aes(y=signal,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("DHS signal in TSS downstream  ")+theme_hd_minimal()+ggtitle("ENCODE-Lung")+theme(axis.text.x = element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))

p4_meth=ggplot(data = aver_meth_gene[aver_meth_gene$cancer=="LUAD",],aes(y=norm_meth,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Methylation in TSS downstream ")+ggtitle("TCGA-LUAD normal")+theme_hd_minimal()+theme(axis.text.x = element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))



pall0_strand=ggarrange(pVNN,pVNN2,nrow=2,ncol=1,labels = c(letters[7:8])  ,font.label = list(size = 18, color = "black", face = "bold"))
pall1_strand=ggarrange(p1_g4,p2_g4,p_freq,nrow=1,ncol=3,labels=c(letters[1:3]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1,2))
pall12_strand=ggarrange(p2_exp,p3_signal,nrow=1,ncol=2,labels=c(letters[4:5]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1))
pall2_strand=ggarrange(p4_meth,pall0_strand,nrow=1,ncol=2,labels=c(letters[6],""),font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend = "bottom")

pp<-ggarrange(pall1_strand,pall12_strand,pall2_strand,nrow=3,ncol=1,font.label = list(size = 18, color = "black", face = "bold"))
tiff("./all_figures/all_figure4.tiff",width = 15,height = 18,res=300,units="in",compression = "lzw")
print(pp)
dev.off()



tiff("./all_figures/supply_all_figure3.tiff",width = 16,height = 20,res=300,units="in",compression = "lzw")
types=sort(as.character(unique(genes_infor$G4_predict_type)))
p2_exp_all=ggplot(data = all_data_tsg[all_data_tsg$cancer!="LUAD"& !is.na(all_data_tsg$Gene.Type),],aes(y=norm_exp,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Expression ")+theme_hd_minimal()+theme(axis.text.x =element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))+facet_wrap(~cancer)
print(p2_exp_all)
dev.off()

tiff("./all_figures/supply_all_figure4.tiff",width = 16,height = 20,res=300,units="in",compression = "lzw")
p3_signal_all=ggplot(data = genes_chromatin[genes_chromatin$cancer!="lung"& !is.na(genes_chromatin$Gene.Type),],aes(y=signal,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("DHS signal in TSS downstream  ")+theme_hd_minimal()+theme(axis.text.x = element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))+facet_wrap(~cancer)
print(p3_signal_all)
dev.off()
tiff("./all_figures/supply_all_figure5.tiff",width = 16,height = 20,res=300,units="in",compression = "lzw")
p_meth_all<-ggplot(data = aver_meth_gene[aver_meth_gene$cancer!="LUAD",],aes(y=norm_meth,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Methylation in TSS downstream ")+theme_hd_minimal()+theme(axis.text.x = element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))+facet_wrap(~cancer)
print(p_meth_all)
dev.off()

print("over")
