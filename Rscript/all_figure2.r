
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_genes_program/G4_final/")

source('functions/functions_evolution_5UTR.r')
load(file="./data/all_meth_exp_signal.rdata")

GTEXmedian_tpm_infor<-merge(GTEXmedian_tpm,genes_infor,by.x="Gene.ID",by.y="ensembl_gene_id")

p1=plotStrandCor(all_data_tsg,"norm_exp","GG_density","CC_density","age_type2",method="pearson")+ggtitle("TCGA Cohort")
p2=plotStrandCor(GTEXmedian_tpm_infor,"norm_exp","GG_density","CC_density","age_type2",method="pearson")+ggtitle("GTEX Cohort")

exp2=PlotMutipleCor(GTEXmedian_tpm_infor,"CC_density","GG_density","Template","Non-template","norm_exp","Correlation of expression and frequency\nof G-runs in different strand in TSS downstream","pearson")+theme(legend.position = "bottom",legend.box ="vertical")+ylab("GTEX Cohort")+labs(fill="Strand")
exp1=PlotMutipleCor(all_data_tsg,"CC_density","GG_density","Template","Non-template","norm_exp","Correlation of expression and frequency\nof G-runs in different strand in TSS downstream","pearson")+theme(legend.position = "bottom",legend.box ="vertical")

sp1_cpd_exp1=plot_scatter_2(all_data_tsg[  all_data_tsg$cancer=="LUAD",],"CC_density","norm_exp","Frequency of G-runs\nin template strand", "Expression ",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+ggtitle("TCGA-LUAD normal")+theme_hd()+ylim(0,16)
sp1_cpd_exp2=plot_scatter_2(all_data_tsg[which(all_data_tsg$cancer=="LUAD"),],"GG_density","norm_exp","Frequency of G-runs\nin non-template strand", "Expression ",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+ggtitle("TCGA-LUAD normal")+theme_hd()+ylim(0,16)
sp1_cpd_exp1GTEX=plot_scatter_2(GTEXmedian_tpm_infor[ GTEXmedian_tpm_infor$cancer=="Lung",],"CC_density","norm_exp","Frequency of G-runs\nin template strand", "Expression ",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+ggtitle("GTEX-Lung normal")+theme_hd()+ylim(0,16)


sp1_cpd_exp2GTEX=plot_scatter_2(GTEXmedian_tpm_infor[GTEXmedian_tpm_infor$cancer=="Lung",],"GG_density","norm_exp","Frequency of G-runs\nin non-template strand", "Expression ",0.01,15,"black")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+ggtitle("GTEX-Lung normal")+theme_hd()+ylim(0,16)


ggarrange(sp1_cpd_exp1,sp1_cpd_exp2,nrow=2,ncol=1)



pall1=ggarrange(sp1_cpd_exp2GTEX,sp1_cpd_exp1GTEX,nrow=1,ncol=2,labels=c(letters[2:3]),font.label = list(size = 18, color = "black", face = "bold"))
pall3=ggarrange(pall1,p2+ylim(0,0.37),nrow=2,ncol=1,labels=c("",letters[4]),font.label = list(size = 18, color = "black", face = "bold"))

pall4<-ggarrange(exp2,pall3,nrow=1,ncol=2,labels=c(letters[1],""),font.label = list(size = 18, color = "black", face = "bold"))
tiff("./all_figures/all_figure2.tiff",width = 17,height = 17,res=300,units="in",compression = "lzw")
print(pall4)
dev.off()


pallS1=ggarrange(sp1_cpd_exp2,sp1_cpd_exp1,nrow=2,ncol=1,labels=c(letters[2:3]),font.label = list(size = 18, color = "black", face = "bold"))

pallS2<-ggarrange(exp1,pallS1,nrow=1,ncol=2,labels=c(letters[1],""),font.label = list(size = 18, color = "black", face = "bold"))

pallS=ggarrange(pallS2,p1+ylim(0,0.36),nrow=2,ncol=1,labels=c("",letters[4]),font.label = list(size = 18, color = "black", face = "bold"))

tiff("./all_figures/supply_all_figure2.tiff",width = 12,height = 14,res=300,units="in",compression = "lzw")
print(pallS)
dev.off()


