
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_genes_program/G4_final/")
genes_infor<-readRDS("./data/genes_infor_downloaded_from_ensembl105_gene_infor.rds")
all_data_genes<-readRDS("./data/genes_infor_downloaded_from_ensembl105_all_data_genes.rds")
load(file="./data/other_species_all_data.rdata")
source('functions/functions_evolution_5UTR.r')
load(file="./data/result_promoter_window_ACTG_1000_1000.rdata")
load(file="./data/all_meth_exp_signal.rdata")

gene_pos<-genes_infor[,c("genes_main_transcript","chromosome_name","transcript_start","transcript_end","strand","ensembl_gene_id")]
colnames(gene_pos)<-c("trans_ID","seqnames","start","end","strand","gene_id")
gene_pos$strand[gene_pos$strand=="1"]="+"
gene_pos$strand[gene_pos$strand=="-1"]="-"
gene_pos$seqnames<-paste("chr",gene_pos$seqnames,sep="")
gene_pos_data<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=T)

pos_all_up<-getFreq(gene_pos_data,-1000,0,Hsapiens,"+")
pos_all_down<-getFreq(gene_pos_data,0,1000,Hsapiens,"+")


x=rbind(cbind(type="TSS_upstream",melt(pos_all_up[1:13])),cbind(type="TSS_downstream" ,melt(pos_all_down[1:13])))
x$variable<-gsub("_density","",as.character(x$variable))
x$type<-factor(x$type,levels =c("TSS_upstream","TSS_downstream") )
x$Gene.Type=genes_infor[x$Group.1,"Gene.Type"]

pGC_content=ggplot(x[x$variable=="GC_content",],aes(x=type,y=value))  + geom_signif(comparisons = getComp(unique(as.character(x$type))),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),outlier.colour = NA,aes(fill=type))+ylab("GC content")+theme_hd_minimal() + guides(fill="none")+scale_fill_manual(values = c("#7FFF00", "#458B00"))+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 30,hjust = 1))
pG_C=ggplot(x[x$type=="TSS_downstream" & x$variable %in% c("G","C"),],aes(x=variable,y=value))  + geom_signif(comparisons = getComp(c("G","C")),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),outlier.colour = NA,aes(fill=variable))+ylab("Nucleotide Frequency")+theme_hd_minimal() + guides(fill="none")+scale_fill_manual(values = c("#7FFF00", "#458B00"))+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 30,hjust = 1))+ggtitle("TSS downstream")


getPairedPlot(genes_infor,"CC_density","GG_density","Template ","Non-template")
p_age=plotSubSigParied(genes_infor[!is.na(genes_infor$age_type2),], "CC_density","GG_density","age_type2","",cols1=c("#40A4DE", "#F07F16C7"),"Template ","Non-template","Age groups(Million years)","Frequency of G-runs\nin TSS downstream",cols_x="black")+ylim(0,0.18)+labs(fill="Strand")

all_data<-all_data[all_data$species != "Drosophila.melanogaster",]  
all_data$species<-gsub("\\.","\n",all_data$species)
all_data$species=factor(all_data$species,levels=c("Human","Mouse","Pig","Cow","Zebrafish","Caenorhabditis\nelegans")[6:1])

p_spe=plotSubSigParied(all_data, "CC_density","GG_density","species","",cols1=c("#40A4DE", "#F07F16C7"),"Template ","Non-template","species","Frequency of G-runs\nin TSS downstream",cols_x="black")+xlab("Species") +ylim(0,0.2)+labs(fill="Strand")+theme(axis.text.x = element_text(angle = 60))



a_site=SlidingWindow(median,1:result_window_1000_1000$length_x,window=50,step=10)
bm_4=plot_density_TSS7(result_window_1000_1000,a_site,"density",result_window_1000_1000$genes_promoter_df$gene_id %in%genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000$genes_promoter_df)



gmt<-read_gmt("./data/h.all.v2022.1.Hs.symbols.gmt",tidy = T)
term2gene=as.data.frame(gmt)
term2gene$gene_set<-gsub("HALLMARK_","",term2gene$gene_set)
x=merge(term2gene,genes_infor,by.x="gene",by.y="hgnc_symbol")
temp_x=aggregate(x[,c("CC_density","GG_density")],list(gene_set=x$gene_set),mean)
x$gene_set<-factor(x$gene_set,levels=temp_x$gene_set[order(temp_x$GG_density-temp_x$CC_density)])
terms_x<-plotSubSigParied(x[!is.na(x$age_type2),], "GG_density","CC_density","gene_set","",cols1=c("#F07F16C7","#40A4DE" ),"Non-template","Template ","HallMark Geneset","Frequency of G-runs\nin TSS downstream",cols_x="black")+labs(fill="Strand")+coord_flip()


pall1=ggarrange(bm_4+theme_hd()+scale_color_manual(labels=c("A","C","G","T"),values = c("#104E8B", "#B23AEE", "#CD2626", "#458B00"))+ylab("Nucleotide Frequency")+xlab("Distance from TSS(bp)")+theme(axis.text.x = element_text(angle = 90)),pGC_content,nrow=1,ncol=2,labels=c(letters[1:2]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1))
pall2=ggarrange(pall1,p_spe+guides(fill="none"),p_age+guides(fill="none"),nrow=3,ncol=1,labels=c("",letters[3:4]),font.label = list(size = 18, color = "black", face = "bold"))
pall3=ggarrange(pall2,terms_x, nrow=1,ncol=2,labels=c("",letters[5]),font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend ="bottom" ,widths = c(1,1.5))

tiff("./all_figures/all_figure1.tiff",width = 15,height = 18,res=300,units="in",compression = "lzw")
print(pall3)
dev.off()




pt=ggplot(genes_infor[!is.na(genes_infor$age_type2),],aes(x=age_type2,y=CC_density,fill=age_type2))+geom_boxplot(show.legend = FALSE)+geom_signif(comparisons=getComp(levels(genes_infor$age_type2[!is.na(genes_infor$age_type2)])), test='wilcox.test', map_signif_level=TRUE,col="black",step_increase=0.1) +theme_hd()+scale_fill_manual(values=colorRampPalette(c("white" ,"#1C84C9"))(6)[2:6])+theme(axis.text.x = element_text(angle = 30,hjust = 1))+ylab("Frequency of G-runs\nin template strand")+xlab("Age groups(Million years)")+ggtitle("All genes")#+stat_summary(fun.data =  countFunction, geom="text", size = 4,col="black")
pnt=ggplot(genes_infor[!is.na(genes_infor$age_type3),],aes(x=age_type2,y=GG_density,fill=age_type2))+geom_boxplot(show.legend = FALSE)+geom_signif(comparisons=getComp(levels(genes_infor$age_type2[!is.na(genes_infor$age_type2)])), test='wilcox.test', map_signif_level=TRUE,col="black",step_increase=0.1) +theme_hd()+scale_fill_manual(values=colorRampPalette(c("white" ,"#CF8310"))(6)[2:6])+theme(axis.text.x = element_text(angle = 30,hjust = 1))+ylab("Frequency of G-runs\nin non-template strand")+xlab("Age groups(Million years)")+ggtitle("All genes")#+stat_summary(fun.data =  countFunction, geom="text", size = 4,col="black")
p_age_all1<-ggarrange(pG_C,pt,pnt,nrow=1,ncol=3,labels=c(letters[1:3]),font.label = list(size = 18, color = "black", face = "bold"))

tiff("./all_figures/supply_all_figure1.tiff",width =7,height = 7,res=300,units="in",compression = "lzw")
print(p_age_all1)
dev.off()








