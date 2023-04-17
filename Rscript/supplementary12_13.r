
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_genes_program/G4_final/")

source('functions/functions_evolution_5UTR.r')
load(file="./data/all_meth_exp_signal.rdata")

load(file="./data/translation_Data.rdata")

pcor_gg<-plotForest(all_TE_translate,"norm_exp_log2","utrRNA_GG_density","Correlation of The Frequency of\nG-run in 5UTR and log2(TE+1)",0.5,1.5,"5.5cm")
ptr2gr=plotForest(all_TR_translate,"norm_exp_log2","utrRNA_GG_density","Correlation of Frequency of G-runs in 5UTR\nand log2(TR+1) (Translation Ratio)",0.5,1.5,"5.5cm")
pevi2gr=plotForest(all_EVI_translate,"norm_exp_log2","utrRNA_GG_density","Correlation of Frequency of G-runs in 5UTR\nand log2(EVI+1) (Elongation Velocity Index)",0.5,1.5,"5.5cm")

tiff("./all_figures/supply_all_figure12.tiff",width = 12,height = 13,res=300,units="in",compression = "lzw")
print(ggarrange(pcor_gg,ptr2gr,pevi2gr,nrow=3,ncol=1,labels=c(letters[1:3]),font.label = list(size = 18, color = "black", face = "bold")))
dev.off()
sig_G4_TE<-get_diff(all_TE_translate,"rG4_5utr_pre","norm_exp_log2","+")
pff_G4_TE<-plotForesttwogroup(sig_G4_TE,"Translational Efficiency(log2(TE+1))","4.2mm","5UTR_G4s+","5UTR_G4s-")
ggarrange(pff_G4_TE$p1)
sig_G4_TR<-get_diff(all_TR_translate,"rG4_5utr_pre","norm_exp_log2","+")
pff_G4_TR<-plotForesttwogroup(sig_G4_TR,"log2(TR+1) (Translation Ratio)","4.2mm","5UTR_RNA_G4s+","5UTR_RNA_G4s-")
ggarrange(pff_G4_TR$p1)

sig_G4_EVI<-get_diff(all_EVI_translate,"rG4_5utr_pre","norm_exp_log2","+")
pff_G4_EVI<-plotForesttwogroup(sig_G4_EVI,"all_EVI_translate","4.2mm","5UTR_RNA_G4s+","5UTR_RNA_G4s-")
ggarrange(pff_G4_EVI$p1)
pall_tr_evi=ggarrange( pff_G4_TE$p1,pff_G4_TR$p1,pff_G4_EVI$p1,nrow=3,ncol=1,labels=c(letters[1:3]),font.label = list(size = 18, color = "black", face = "bold"))
tiff("./supply_all_figure13.tiff",width = 12,height = 12,res=300,units="in",compression = "lzw")
print(pall_tr_evi)
dev.off()
