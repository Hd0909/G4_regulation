options(stringsAsFactors = FALSE)
devtools::install_github(repo = "tiagochst/ELMER.data")
library("ELMER.data")
library("GenomicRanges")
library(rstatix)
library(ggplotify)
library("ggExtra")
library(ggsignif)
library("cowplot")
library("gridExtra")
library("RColorBrewer")
library("ggpubr")
#devtools::install_github("jhrcook/readgmt")
library(readgmt)
library(clusterProfiler)
library(ggplot2)
library(VennDiagram)
library(stringr)
library("ggvenn")
require(plyr)
#library(metacor)
library(foreach)
library(doParallel)
library(meta)
library(parallel)
library(GenomicRanges)
library(reshape2)
library(Repitools)
library(EnrichedHeatmap)
library(evobiR)
library(pqsfinder)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(proto)
library(GGally)
require(plyr)
library(png)
library(grid)
library(cocor)
library(meta)
library(GSEABase)
library("ggVennDiagram")
#library(BioVenn)
#install.packages("extrafont")
library(extrafont)
library(pqsfinder)
#font_import()
### get the expression wildth
pvalueTosig<-function(pvalue_fish){
  Sig <- ifelse(pvalue_fish>0.05,"NS",ifelse(pvalue_fish>0.01, "*" ,  ifelse(pvalue_fish >0.001, "**" , "***")))
  return(Sig) 
  }

get_cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
#get_cv_log<-function(x){sqrt(exp(sd(x,na.rm =T)^2)-1)}
get_cv_log<-function(x){
  
  x=2^x-1##get the original TPM value
  return(sd(x,na.rm=T)/mean(x,na.rm=T))
}
get_TSEI_log2<-function(data_x){
  ##transform the log(TPM+1) back to TPM
  data_x=2^data_x-1
  sum(1-data_x/max(data_x))/(length(data_x)-1)
}


get_TSEI_log2_array<-function(data_x){
  ##transform the log(TPM+1) back to TPM
  data_x=2^data_x
  sum(1-data_x/max(data_x))/(length(data_x)-1)
}


get_TSEI_pre<-function(data_x){
  sum(1-2^(data_x-max(data_x)))/(length(data_x)-1)
}


get_TSEI<-function(data_x){
  sum(1-data_x/max(data_x))/(length(data_x)-1)
}

countFunction <- function(x){
  return(data.frame(y=0,label=round(length(x),2)))}

get_substring<-function(x,sp,n){
  x_1=strsplit(x,sp,fixed = T)
  x_2=do.call(rbind.data.frame, x_1)
  return(as.character(x_2[,n]))
}
theme_hd<-function(base_size = 14, base_family = "Helvetica"){
  theme_classic(base_size = base_size, base_family = base_family) %+replace%  theme(plot.title = element_text(lineheight=.8, size=14, face="bold"),
                                                                                    axis.line = element_line(colour = "black",size = 1),
                                                                                    legend.background = element_rect(colour = "white"),
                                                                                    legend.key = element_rect(colour = "white"),legend.position = "bottom",
                                                                                    axis.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    axis.title = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    strip.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.text =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.title =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    plot.margin = unit(c(1,1,1,1), "cm"))}


theme_hd_minimal<-function(base_size = 14, base_family = "Helvetica"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%  theme(plot.title = element_text(lineheight=.8, size=14, face="bold"),
                                                                                    axis.line = element_line(colour = "black",size = 1),
                                                                                    legend.background = element_rect(colour = "white"),
                                                                                    legend.key = element_rect(colour = "white"),legend.position = "bottom",
                                                                                    axis.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    axis.title = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    strip.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.text =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.title =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    plot.margin = unit(c(1,1,1,1), "cm"))}

theme_hd_minimal2<-function(base_size = 14, base_family = "Helvetica"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%  theme(plot.title = element_text(lineheight=.8, size=14, face="bold"),
                                                                                    legend.background = element_rect(colour = "white"),
                                                                                    legend.key = element_rect(colour = "white"),legend.position = "bottom",
                                                                                    axis.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    axis.title = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.text =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.title =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    plot.margin = unit(c(1,1,1,1), "cm"))}



get_main_transcript<-function(get_gene_pos_temp){
  # get the transcript which is principal1
  main_trans<-subset( get_gene_pos_temp,APPRIS.annotation=="principal1")
  genes=unique(get_gene_pos_temp$ensembl_gene_id)
  genes=genes[! genes %in% main_trans$ensembl_gene_id] ## genes which have no principal transcript
  get_gene_pos_temp=get_gene_pos_temp[!get_gene_pos_temp$ensembl_gene_id  %in% main_trans$ensembl_gene_id,]
  get_gene_pos_temp<-rbind(main_trans,get_gene_pos_temp) # genes with  principal1 and genes that do not have principal1
  max_trans_length<-aggregate(get_gene_pos_temp$trans_length,list(get_gene_pos_temp$ensembl_gene_id),function(x){max(x,na.rm = T)}) # get the max transcript length
  f_select=paste(get_gene_pos_temp$ensembl_gene_id,get_gene_pos_temp$trans_length,sep="_") %in% paste(max_trans_length[,1],max_trans_length[,2],sep="_")
  select_gene_trans<-get_gene_pos_temp[f_select,]
  # for genes with transcript which have the same length randomly select one
  return(select_gene_trans[!duplicated(select_gene_trans$ensembl_gene_id),])
}


get_main_transcriptENSEMBL<-function(get_gene_pos_temp){
  # get the transcript which is principal1
  main_trans<-subset( get_gene_pos_temp,Ensembl.Canonical==1)
  genes=unique(get_gene_pos_temp$Gene.stable.ID)
  genes=genes[! genes %in% main_trans$Gene.stable.ID] ## genes which have no principal transcript
  get_gene_pos_temp=get_gene_pos_temp[!get_gene_pos_temp$Gene.stable.ID  %in% main_trans$Gene.stable.ID,]
  get_gene_pos_temp<-rbind(main_trans,get_gene_pos_temp) # genes with  principal1 and genes that do not have principal1
  max_trans_length<-aggregate(get_gene_pos_temp$Transcript.length..including.UTRs.and.CDS.,list(get_gene_pos_temp$Gene.stable.ID),function(x){max(x,na.rm = T)}) # get the max transcript length
  f_select=paste(get_gene_pos_temp$Gene.stable.ID,get_gene_pos_temp$Transcript.length..including.UTRs.and.CDS.,sep="_") %in% paste(max_trans_length[,1],max_trans_length[,2],sep="_")
  select_gene_trans<-get_gene_pos_temp[f_select,]
  # for genes with transcript which have the same length randomly select one
  return(select_gene_trans[!duplicated(select_gene_trans$Gene.stable.ID),])
}




plot_age_trans_count<-function(all_data_age_x,types_x){
  sts <- boxplot.stats(all_data_age_x$transcript_count)$stats
  p=ggplot(all_data_age_x, aes(x=factor(Age_class), y=transcript_count, fill=Gene.Type)) + 
    geom_boxplot(colour = "#E5E5E5", outlier.color = "#B3B3B3",outlier.shape = NA) + 
    coord_cartesian(ylim = c(sts*1.3,sts/1.3))+
    scale_fill_manual(values=c("#1874CD","#EE2C2C"), name = "Gene Group") +ggtitle(types_x)+
    theme_bw()+
    theme(plot.title = element_text(lineheight=.8, size=14, face="bold.italic"),
          axis.line.x = element_line(colour = "black",size = 0.7,lineend = 2),
          axis.line.y = element_line(colour = "black",size = 0.7),
          axis.title.x = element_text( size=17, face="bold"),
          axis.title.y = element_text( size=17, face="bold"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x=element_text( size=12,face="bold",vjust = 0.3),
          axis.text.y=element_text( size=12,face="bold",vjust = 0.3,angle = 90),
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(colour = "white"),
          plot.margin = unit(c(0.5,1,0.5,0.3), "cm")) + labs(x = "Origin Time", y = "transcript count")
  print(p)
}



#####################33 Only calculated one time if two TSG or non-cancer genes are paralogs    ; calculated the mean transcript count for the paralogs 
calculate_paralog<-function(data_x,all_genes,trans_count){
  paralog_transcount<-c()
  genes=unique(data_x$ensembl_gene_id)
  gene_select=c()
  gene_dup=c()
  paralog_unique=c()
  for(gene in genes){
    if(gene %in% gene_dup){next}
    para_gene=data_x[data_x[,1]==gene,]
    types_gene=all_genes[gene,"Gene.Type"]
    temp_cancer=c(gene,para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[all_genes$Gene.Type==types_gene]])
    temp_non_cancer=para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[all_genes$Gene.Type=="Non-Cancer"]]
    temp_other_cancer=para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[!all_genes$Gene.Type %in% c("Non-Cancer",types_gene)]]
    temp_cancer_count=c(mean(trans_count$all_trans[temp_cancer,"transcript_count"],rm.na=T),
                        mean(trans_count$pro_trans[temp_cancer,"transcript_count"],rm.na=T),
                        mean(trans_count$non_pro_trans[temp_cancer,"transcript_count"],rm.na=T),
                        mean(trans_count$non_pro_trans[temp_cancer,"Age_class"],rm.na=T))
    temp_non_cancer_count=c(mean(trans_count$all_trans[temp_non_cancer,"transcript_count"],rm.na=T),
                            mean(trans_count$pro_trans[temp_non_cancer,"transcript_count"],rm.na=T),
                            mean(trans_count$non_pro_trans[temp_non_cancer,"transcript_count"],rm.na=T),
                            mean(trans_count$non_pro_trans[temp_non_cancer,"Age_class"],rm.na=T))
    temp_other_cancer_count=c(mean(trans_count$all_trans[temp_other_cancer,"transcript_count"],rm.na=T),
                              mean(trans_count$pro_trans[temp_other_cancer,"transcript_count"],rm.na=T),
                              mean(trans_count$non_pro_trans[temp_other_cancer,"transcript_count"],rm.na=T),
                              mean(trans_count$non_pro_trans[temp_other_cancer,"Age_class"],rm.na=T))
    paralog_transcount<-rbind(paralog_transcount,c(length(temp_cancer),temp_cancer_count,
                                                   length(temp_non_cancer),temp_non_cancer_count,
                                                   length(temp_other_cancer),temp_other_cancer_count))
    gene_dup=c(gene_dup,temp_cancer)
    gene_select=c(gene_select,gene)
    paralog_unique<-rbind(paralog_unique,para_gene)
  }
  rownames(paralog_transcount)<-gene_select
  colnames(paralog_transcount)<-paste(rep(c("cancer","non_cancer_paralog","other_patalog"),each=5),c("length","all","protein","non_protein","Age_class"),sep="_")
  return(list(a=as.data.frame(paralog_transcount),b=paralog_unique))
}

calculate_paralog_2<-function(data_x,all_genes){
  paralog_transcount<-c()
  genes=unique(data_x$ensembl_gene_id)
  gene_select=c()
  gene_dup=c()
  paralog_unique=c()
  for(gene in genes){
    if(gene %in% gene_dup){next}
    para_gene=data_x[data_x[,1]==gene,]
    types_gene=all_genes[gene,"Gene.Type"]
    temp_cancer=c(gene,para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[all_genes$Gene.Type==types_gene]])
    gene_dup=c(gene_dup,temp_cancer)
    paralog_unique<-rbind(paralog_unique,para_gene)
  }
  return(paralog_unique)
}


get_paralog_same_age_type<-function(x_paralog,all_gene_trans_infor){
  genes<-unique(x_paralog[,1])
  paralog_pairs<-c()
  for(gene in genes){
    
    temp_p<-x_paralog[x_paralog[,1]==gene,]
    temp_infor<-all_gene_trans_infor[unique(c(gene,temp_p[,2])),]
    temp_age=temp_infor[which(temp_infor$age_type==temp_infor[gene,"age_type"]),]
    if(length(temp_age[,1])==0){next}
    trans_count=aggregate(temp_age$all_trans_count,list(temp_age$Gene.Type),median)
    gene_count=aggregate(temp_age$all_trans_count,list(temp_age$Gene.Type),length)
    pro_count=aggregate(temp_age$pro_count,list(temp_age$Gene.Type),median)
    non_pro_count=aggregate(temp_age$non_pro_trans_count,list(temp_age$Gene.Type),median)
    promoter_cpg_oe=aggregate(temp_age$promoter_cpg_OE,list(temp_age$Gene.Type),median)
    temp_result=cbind(gene,pro_count[,1],temp_infor[gene,"age_type"],gene_count[,2],trans_count[,2], pro_count[,2],non_pro_count[,2],promoter_cpg_oe[,2])
    paralog_pairs<-rbind(paralog_pairs,temp_result)
  }
  rownames(paralog_pairs)<-paralog_pairs[,1]
  colnames(paralog_pairs)<-c("gene","Gene.Type","age_type","Gene.Type.count","all_trans_count","pro_trans_count","non_pro_trans_count","promoter_cpg_OE")
  return(as.data.frame(paralog_pairs))
}





plot_promoter_data<-function(genes_infor,title_x,y_name){
  TSG_non_cancer<-genes_infor[genes_infor$Gene.Type!="Oncogene",]
  TSG_non_cancer$HK_types=factor(paste(TSG_non_cancer$HK,TSG_non_cancer$Gene.Type,sep=":"),levels = c("HK_gene:Non-Cancer","non_HK_gene:TSG","HK_gene:TSG","non_HK_gene:Non-Cancer"))
  comp3=list(c("1_13", "13_20"),c("1_13", "20_23"),c("1_13", "23+"),c("13_20", "20_23"),c("13_20", "23+"),c("20_23", "23+"))
  x=quantile(TSG_non_cancer[,y_name], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  s1=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=as.factor(age_type) ,y=get(y_name),color=HK))  + 
    facet_wrap(~Gene.Type+HK,nrow=1) + geom_signif(comparisons = comp3,map_signif_level=TRUE,
                                                   textsize=3,tip_length = 2,col="black", test='wilcox.test', y_position=seq(x[2]+steps*3,x[2]+steps*9,steps))+
    geom_point(size=I(0.5),position=position_jitterdodge(dodge.width=0.9)) +
    geom_boxplot(alpha=I(0.7),outlier.colour = NA, 
                 position = position_dodge(width=0.9))+scale_y_continuous(limits =x+c(0,steps*9))+ggtitle(title_x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_type")
  
  s2=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=as.factor(Age_class) ,y=get(y_name),color=I("white"),fill=HK))  + 
    facet_wrap(~Gene.Type+HK,nrow=1) +
    geom_boxplot(alpha=I(0.7),position = position_dodge(width=0.9))+scale_y_continuous(limits =x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_class")+ggtitle("")
  p_all<-ggarrange(s1, s2, nrow = 2, labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )
  print(p_all)
}


plot_promoter_data_2<-function(genes_infor,title_x,y_name){
  TSG_non_cancer<-genes_infor[genes_infor$Gene.Type!="Oncogene",]
  comp3=list(c("Non-Cancer", "TSG"))
  x=quantile(TSG_non_cancer[,y_name], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  s1=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=Gene.Type ,y=get(y_name),col=Gene.Type))  + 
    facet_wrap(~age_type,nrow=1) + geom_signif(comparisons = comp3,map_signif_level=TRUE,
                                               textsize=3,tip_length = 2,col="black", test='wilcox.test', y_position=seq(x[2]+steps*6,x[2]+steps*9,steps))+
    geom_point(size=I(0.5),position=position_jitterdodge(dodge.width=0.9)) +
    geom_boxplot(alpha=I(0.7),outlier.colour = NA, 
                 position = position_dodge(width=0.9))+scale_y_continuous(limits =x+c(0,steps*7))+ggtitle(title_x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_type")
  
  s2=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=Gene.Type ,y=get(y_name),col=Gene.Type))  + 
    facet_wrap(~HK+age_type,nrow=1) + geom_signif(comparisons = comp3,map_signif_level=TRUE,
                                                  textsize=3,tip_length = 2,col="black", test='wilcox.test', y_position=seq(x[2]+steps*3,x[2]+steps*9,steps))+
    geom_point(size=I(0.5),position=position_jitterdodge(dodge.width=0.9)) +
    geom_boxplot(alpha=I(0.7),outlier.colour = NA, 
                 position = position_dodge(width=0.9))+scale_y_continuous(limits =x+c(0,steps*7))+ggtitle(title_x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90),strip.text = element_text(size=6))+xlab("Age_type")
  s3=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=as.factor(Age_class) ,y=get(y_name),color=I("white"),fill=Gene.Type))  + 
    facet_wrap(~Gene.Type,nrow=2) +
    geom_boxplot(alpha=I(0.7),position = position_dodge(width=0.9))+scale_y_continuous(limits =x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_class")+ggtitle("")
    p_all<-ggarrange(s1, s2, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "bottom"  )
    p_all1<-ggarrange(p_all, s3, ncol  = 2, common.legend = TRUE, legend = "bottom" )
  print(p_all1)
}




firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}



plot_boxplot<-function( all_data_tsg,y_lab,title_x){
  comp=list(c("low","median"),c("low","high"),c("median","high"))
  x=quantile(all_data_tsg[,y_lab], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  p=ggplot(data = all_data_tsg,aes(y=get(y_lab),x=promoter_cpg_OE_category,fill=promoter_cpg_OE_category,col=I("grey")))+
    geom_boxplot(outlier.colour = NA)+facet_wrap(~cancer,nrow=1)
  s3=p+theme(axis.text.x = element_text(angle = 90),strip.text = element_text(size=25),axis.title =element_text(size=20))+ 
    scale_y_continuous(limits =x+c(0,steps*5))+
    geom_signif(comparisons = comp,map_signif_level=TRUE,col="black",tip_length = 4, test='wilcox.test', y_position=seq(x[2]+steps*2,x[2]+steps*9,steps))+ggtitle(title_x)
  
  return(s3+ylab(y_lab))
}

plot_boxplot_2<-function( all_data_tsg,y_lab,title_x){
  comp=list(c("low","median"),c("low","high"),c("median","high"))
  x=quantile(all_data_tsg[,y_lab], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  p=ggplot(data = all_data_tsg,aes(y=get(y_lab),x=promoter_cpg_OE_category,fill=Gene.Type,col=I("grey")))+
    geom_boxplot(outlier.colour = NA)+facet_wrap(~cancer,nrow=1)
  s3=p+theme(axis.text.x = element_text(angle = 90),strip.text = element_text(size=25),axis.title =element_text(size=20))+ 
    scale_y_continuous(limits =x+c(0,steps*5))+
    geom_signif(comparisons = comp,map_signif_level=TRUE,col="black",tip_length = 4, test='wilcox.test', y_position=seq(x[2]+steps*2,x[2]+steps*9,steps))+ggtitle(title_x)
  
  return(s3+ylab(y_lab))
}






### delete low expression
get_cor_for_each_cancer<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    print(cancer)
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff_n=quantile(temp_data$norm_exp,0.2,na.rm=T)
    quantile_cutoff_t=quantile(temp_data$tumor_exp,0.2,na.rm=T)
    if (cancer %in% c("CRC_GSE137327_RNA","Lung_IMA_GSE86958_RNA")){
      quantile_cutoff_n=1
      quantile_cutoff_t=1   
    }
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff_n & temp_data$tumor_exp< quantile_cutoff_t )),]
     temp_data$diff_exp=abs(temp_data$diff_exp)
        temp_cor_non_cancer=cor.test(temp_data[temp_data$Gene.Type=="Non-Cancer",name1],temp_data[temp_data$Gene.Type=="Non-Cancer",name2],method = method_x)
    temp_cor_TSG=cor.test(temp_data[temp_data$Gene.Type=="TSG",name1],temp_data[temp_data$Gene.Type=="TSG",name2],method = method_x)
    len_non_cancer=length(temp_data[temp_data$Gene.Type=="Non-Cancer",name1])
    len_TSG=length(temp_data[temp_data$Gene.Type=="TSG",name1])
    if(method_x=="kendall"){
      r1=sin(pi*0.5*temp_cor_non_cancer$estimate)
      r2=sin(pi*0.5*temp_cor_TSG$estimate)      
    }else{
      r1=temp_cor_non_cancer$estimate
      r2=temp_cor_TSG$estimate    
    }

    if(r1>r2){temp="greater"
    temp_x="<"}else{temp="less"
    temp_x=">"}
    cocor_x=cocor.indep.groups(r1.jk=r1, r2.hm=r2, n1=len_non_cancer, n2=len_TSG, alternative=temp, alpha=0.05, conf.level=0.95, null.value=0)
    result_cor<-rbind(result_cor,c(cancer,"Non_Cancer",round(temp_cor_non_cancer$estimate,3),temp_cor_non_cancer$p.value,cocor_x@fisher1925$p.value))
    result_cor<-rbind(result_cor,c(cancer,"TSG",round(temp_cor_TSG$estimate,3),temp_cor_TSG$p.value,cocor_x@fisher1925$p.value))
  }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue","Diff_pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  result_cor$Pvalue=-log10(as.numeric(result_cor$Pvalue))
  result_cor$Diff_pvalue=as.numeric(result_cor$Diff_pvalue)
  return(as.data.frame(result_cor))
}

get_cor_for_each_cancer_tcga<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    #print(cancer)
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=1 
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff & temp_data$tumor_exp< quantile_cutoff )),]
    temp_data$diff_exp=abs(temp_data$diff_exp)
    temp_cor_non_cancer=cor.test(temp_data[temp_data$Gene.Type=="Non-Cancer",name1],temp_data[temp_data$Gene.Type=="Non-Cancer",name2],method = method_x)
    temp_cor_TSG=cor.test(temp_data[temp_data$Gene.Type=="TSG",name1],temp_data[temp_data$Gene.Type=="TSG",name2],method = method_x)
    len_non_cancer=length(temp_data[temp_data$Gene.Type=="Non-Cancer",name1])
    len_TSG=length(temp_data[temp_data$Gene.Type=="TSG",name1])
    if(method_x=="kendall"){
      r1=sin(pi*0.5*temp_cor_non_cancer$estimate)
      r2=sin(pi*0.5*temp_cor_TSG$estimate)      
    }else{
      r1=temp_cor_non_cancer$estimate
      r2=temp_cor_TSG$estimate    
    }
    
    if(r1>r2){temp="greater"
    temp_x="<"}else{temp="less"
    temp_x=">"}
    cocor_x=cocor.indep.groups(r1.jk=r1, r2.hm=r2, n1=len_non_cancer, n2=len_TSG, alternative=temp, alpha=0.05, conf.level=0.95, null.value=0)
    result_cor<-rbind(result_cor,c(cancer,"Non_Cancer",round(temp_cor_non_cancer$estimate,3),temp_cor_non_cancer$p.value,cocor_x@fisher1925$p.value))
    result_cor<-rbind(result_cor,c(cancer,"TSG",round(temp_cor_TSG$estimate,3),temp_cor_TSG$p.value,cocor_x@fisher1925$p.value))
  }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue","Diff_pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  result_cor$Pvalue=-log10(as.numeric(result_cor$Pvalue))
  result_cor$Diff_pvalue=as.numeric(result_cor$Diff_pvalue)
  return(as.data.frame(result_cor))
}

get_filter_low_expression_DW_tcga<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=1
    print(quantile_cutoff)   
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff & temp_data$tumor_exp< quantile_cutoff )),]
    temp_data$diff_exp=abs(temp_data$diff_exp)
    result=rbind(result,temp_data)
  }
  return(result[result$Gene.Type=="TSG",])
}


get_filter_low_expression_DW<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff_n=quantile(temp_data$norm_exp,0.2,na.rm=T)
    quantile_cutoff_t=quantile(temp_data$tumor_exp,0.2,na.rm=T)
    if (cancer %in% c("CRC_GSE137327_RNA","Lung_IMA_GSE86958_RNA")){
      quantile_cutoff_n=1
      quantile_cutoff_t=1   
    }
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff_n & temp_data$tumor_exp< quantile_cutoff_t )),]
    temp_data$diff_exp=abs(temp_data$diff_exp)
    result=rbind(result,temp_data)
     }
  return(result[result$Gene.Type=="TSG",])
}

get_cor_for_each_cancer0<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    temp_cor_non_cancer=cor.test(temp_data[temp_data$Gene.Type=="Non-Cancer",name1],temp_data[temp_data$Gene.Type=="Non-Cancer",name2],method = method_x)
    temp_cor_TSG=cor.test(temp_data[temp_data$Gene.Type=="TSG",name1],temp_data[temp_data$Gene.Type=="TSG",name2],method = method_x)
    len_non_cancer=length(temp_data[temp_data$Gene.Type=="Non-Cancer",name1])
    len_TSG=length(temp_data[temp_data$Gene.Type=="TSG",name1])
    if(method_x=="kendall"){
      r1=sin(pi*0.5*temp_cor_non_cancer$estimate)
      r2=sin(pi*0.5*temp_cor_TSG$estimate)      
    }else{
      r1=temp_cor_non_cancer$estimate
      r2=temp_cor_TSG$estimate    
    }

    if(r1>r2){temp="greater"
    temp_x="<"}else{temp="less"
    temp_x=">"}
    cocor_x=cocor.indep.groups(r1.jk=r1, r2.hm=r2, n1=len_non_cancer, n2=len_TSG, alternative=temp, alpha=0.05, conf.level=0.95, null.value=0)
    result_cor<-rbind(result_cor,c(cancer,"Non_Cancer",round(temp_cor_non_cancer$estimate,3),temp_cor_non_cancer$p.value,cocor_x@fisher1925$p.value))
    result_cor<-rbind(result_cor,c(cancer,"TSG",round(temp_cor_TSG$estimate,3),temp_cor_TSG$p.value,cocor_x@fisher1925$p.value))
  }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue","Diff_pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[result_cor$Pvalue==0]=2.2e-16
  result_cor$Pvalue=-log10(as.numeric(result_cor$Pvalue))
  result_cor$Diff_pvalue=as.numeric(result_cor$Diff_pvalue)
  return(as.data.frame(result_cor))
}


get_cor_for_each_cancer_all<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    temp_cor_all=cor.test(temp_data[,name1],temp_data[,name2],method = method_x)
    result_cor<-rbind(result_cor,c(cancer,"All_genes",round(temp_cor_all$estimate,3),formatC(temp_cor_all$p.value, format = "e", digits = 2)))
 }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  
  if(result_cor$cor[1]<0){order_cancer=result_cor[order(result_cor$cor,decreasing = T),1]}else{order_cancer=result_cor[order(result_cor$cor,decreasing = F),1]}
  
  result_cor$Cancer=factor(result_cor$Cancer,levels=order_cancer)
  return(as.data.frame(result_cor))
}

get_cor_for_each_age<-function(all_data_tsg,name1,name2,ylab_x,method_x){
  cancers=unique(all_data_tsg$cancer)
  ages=unique(all_data_tsg[!is.na(all_data_tsg[,ylab_x]),ylab_x])
  result_cor=c()
  for(cancer in cancers){
    for(age in ages){
      temp_data=all_data_tsg[all_data_tsg[,ylab_x]==age &all_data_tsg$cancer==cancer ,]
      temp_cor_all=cor.test(temp_data[,name1],temp_data[,name2],method = method_x)
      result_cor<-rbind(result_cor,c(cancer,age,"All_genes",round(temp_cor_all$estimate,3),formatC(temp_cor_all$p.value, format = "e", digits = 2)))
      
    }
  }
  colnames(result_cor)<-c("Cancer","age","Gene.Type","cor","Pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  
  return(as.data.frame(result_cor))
}



plot_cor<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer0(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  #print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="TSG",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                 ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                             yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))
  
  cor_cpg_exp$Significance<- ifelse(round(10^-cor_cpg_exp$Pvalue,5)>0.05,"NS", "Pvalue<0.05")
  
   cor_plot1=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue,shape=Significance))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))+
     scale_shape_manual("Significance",values=c(1,19),label=c("No","yes"))
  
   #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot1=cor_plot1,cor_plot_matrix=cor_cpg_exp))
}



plot_cor_all<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer_all(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="All_genes",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                            ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))
  
  cor_cpg_exp$Significance<- ifelse(round(10^-cor_cpg_exp$Pvalue,5)>0.05,"NS", "Pvalue<0.05")
  
  cor_plot1=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                     yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue,shape=Significance))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))+
    scale_shape_manual("Significance",values=c(1,19),label=c("No","yes"))
  
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot1=cor_plot1,cor_plot_matrix=cor_cpg_exp))
}


plot_cor_1<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer0(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="TSG",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                            ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("-log10(0.05)","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1.3,5,10,15),guide="legend")+xlab(firstup(method_x))
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot_matrix=cor_cpg_exp))
}


plot_cor_no_order<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer0(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  order_cancer=levels(all_data_tsg$cancer)
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="TSG",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                            ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("-log10(0.05)","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1.3,5,10,15),guide="legend")+xlab(firstup(method_x))
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot_matrix=cor_cpg_exp))
}

plot_cor_TSG<-function(all_data_tsg,name1,name2,title_x){
  cor_cpg_exp<-get_cor_for_each_cancer(all_data_tsg,name1,name2)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>30]=29
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=Pearson_Cor))+geom_segment(data=cor_matrix,aes(x=TSG,xend=Non_Cancer, y=factor(Cancer), yend=factor(Cancer)),colour="#BFEFFF", size = 1)+
    scale_fill_manual(values=c("#1874CD","#EE2C2C"))+ geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label="***")+
    geom_point(data=cor_cpg_exp,aes(y=Cancer,x=Pearson_Cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(cor_plot)
}




plot_non_paired_paralog<-function(TSG_paralog_non_cancer,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=-0.5,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non_Cancer_paralogs"))
  p1=ggplot(data = TSG_paralog_non_cancer[!is.na(TSG_paralog_non_cancer[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+
    ggtitle("Paralogs")+ theme(plot.title = element_text(face = "bold"),axis.text.x =element_blank(),axis.title.x = element_blank())
  p1=p1#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")+ylab(y_name)
  return(p1+ylab(y_lab))
}

plot_class<-function(genes_infor,x_name,y_name,x_lab,y_lab){
  countFunction <- function(x){
    return(data.frame(y=-0.5,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non-Cancer"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=get(x_name),y=get(y_name),fill=get(x_name)))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = getComp(unique(genes_infor[,x_name])), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#B2D0F7", "#0A67BF"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"))
  p1=p1#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  return(p1+ylab(y_lab)+xlab(x_lab))
}

plot_class_legend<-function(genes_infor,x_name,y_name,x_lab,y_lab){
  countFunction <- function(x){
    return(data.frame(y=-0.5,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non-Cancer"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=get(x_name),y=get(y_name),fill=get(x_name)))+geom_boxplot(width=0.5,size=1)+
    geom_signif(comparisons = getComp(unique(genes_infor[,x_name])), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#B2D0F7", "#0A67BF"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"))
  p1=p1#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  return(p1+ylab(y_lab)+xlab(x_lab))
}

plot_TSG_non_cancer<-function(genes_infor,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=-0.5,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non-Cancer"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"),axis.text.x =element_blank(),axis.title.x = element_blank())
  p1=p1#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  return(p1+ylab(y_lab))
}

plot_TSG_non_cancer_withlegend<-function(genes_infor,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=0,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non-Cancer"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = TRUE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=function(p)sprintf("p = %.2g", p),col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"),axis.text.x =element_blank(),axis.title.x = element_blank())
  p1=p1+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")+ylab(y_name)
  return(p1+ylab(y_lab))
}

plot_all_genes<-function(genes_infor,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=0,label=round(length(x),2)))}
  genes_infor$Gene.Type=factor(genes_infor$Gene.Type,levels = c("Non-Cancer","TSG","Oncogene"))
  comp=list(c("TSG", "Non-Cancer"),c("TSG","Oncogene"),c("Non-Cancer","Oncogene"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=function(p)sprintf("p = %.2g", p),col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#2E8B57", "#EE2C2C"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"))
  p1=p1+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")+ylab(y_name)
  return(p1+ylab(y_lab))
}


plotPairedParalog<-function(TSG_paralog_age,col.x,lab.x){
  
  tsg_genes=TSG_paralog_age[TSG_paralog_age[,2]=="TSG",]
  non_cancer_genes=TSG_paralog_age[TSG_paralog_age[,2]=="Non-Cancer",]
  overlap_gene=as.character(tsg_genes[as.character(tsg_genes[,1]) %in% as.character(non_cancer_genes[,1]),1])
  
  data_x=data.frame(TSG=as.numeric(as.character(tsg_genes[overlap_gene,col.x])),Non_Cancer_Paralogs=as.numeric(as.character(non_cancer_genes[overlap_gene,col.x])),age=as.numeric(as.character(non_cancer_genes[overlap_gene,]$Age_class)))
  
  p2=ggpaired(data_x[!(is.na(data_x$TSG) | is.na(data_x$Non_Cancer_Paralogs)),], cond2 = "TSG", cond1 = "Non_Cancer_Paralogs",
              color  = "condition", line.color = "gray78", line.size = 0.1,palette = "npg",size=1,width = 0.5)+ylab(lab.x)+theme_hd()+ theme(plot.title = element_text(hjust = 0.5,face = "plain"),axis.text.x  = element_blank(),axis.title.x = element_blank())+
    stat_compare_means(method="wilcox.test", paired=TRUE, aes(label = paste0("p = ", ..p.format..)),size=4.8,fontface="bold")+ggtitle("Paralogs with the same age")+scale_color_manual("Condition",values=c("#1874CD","#EE2C2C"))
  p2=p2+theme(axis.title.y = element_text(size = 14),plot.title = element_text(size=14,face="bold"))
  return(p2)
}

pairs_list<-function(x){
  re_all<-list()
  
  for(i in 1:(length(x)-1)){
    for(j in (i+1):length(x)){
    re_all[[paste(i,j)]]<-c(x[i],x[j])      
    }
  }
  return(re_all)
}


plotAgelength<-function(genes_infor,y.name,y_lab,pos,y_lab_2){
  genes_infor=genes_infor[!is.na(genes_infor$Age_class),]
  genes_infor$age_type=genes_infor[,y_lab_2]
  countFunction <- function(x){
    return(data.frame(y=pos,label=round(length(x),2)))}
  in.x1<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.7
  in.x<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.5
  label_x<-unique(genes_infor$age_type)
  comp3=pairs_list(label_x)
  s3=ggplot(genes_infor[!is.na(genes_infor$Age_class),],aes(x=as.factor(age_type) ,y=get(y.name)))  + geom_signif(comparisons = comp3,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(aes(fill=Gene.Type),width=0.4)+ylab(y_lab)
  label.df <- data.frame(age_type=label_x,y.n = rep(in.x1,length(label_x)))
  colnames(label.df)<-c("age_type",y.name)
  r <- .1
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- r*2*in.x1/3 * sin(t)
  arc.df <- data.frame(age_type = x, x= y)
  colnames(arc.df)<-c("age_type",y.name)
  t.p<-c()
  for(i in label.df$age_type){
    t1<-wilcox.test(genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$age_type==i,y.name],genes_infor[genes_infor$Gene.Type=="Non-Cancer" & genes_infor$age_type==i,y.name])
    t.p<-c(t.p,t1$p.value)
  }
  
  Sig <- ifelse(t.p>0.05,"NS",ifelse(t.p >0.01, "*" , ifelse(t.p >0.001, "**" , "***")))
  
  s3=s3 + geom_text(data = label.df, label = Sig,col=I("blue"))+ geom_line(data = arc.df, aes(age_type+1, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+2, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+3, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+4, get(y.name)+in.x))+xlab("Gene Age (million years)")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+theme(axis.text.x =   element_text(angle = 30,hjust = 1))
  return(s3)
         #+stat_summary(fun.data =  countFunction, col="black",geom="text", size = 4.8,fontface="bold"))
}



plotAgelength5utrratio<-function(genes_infor,y.name,y_lab,pos){
  countFunction <- function(x){
    return(data.frame(y=pos,label=round(length(x),2)))}
  in.x1<-quantile(genes_infor[,y.name],0.996,na.rm = T)*1.7
  in.x<-quantile(genes_infor[,y.name],0.996,na.rm = T)*1.5
  comp3=list(c("1_13", "13_20"),c("1_13", "20_23"),c("1_13", "23+"),c("13_20", "20_23"),c("13_20", "23+"),c("20_23", "23+"))
  s3=ggplot(genes_infor[!is.na(genes_infor$Age_class),],aes(x=as.factor(age_type) ,y=get(y.name)))  + geom_signif(comparisons = comp3,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),aes(fill=Gene.Type),width=0.4)+ylab(y_lab)
  label.df <- data.frame(age_type= c("1_13", "13_20","20_23", "23+"),y.n = rep(in.x1,4))
  colnames(label.df)<-c("age_type",y.name)
  r <- .1
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- r*2*in.x1/3 * sin(t)
  arc.df <- data.frame(age_type = x, x= y)
  colnames(arc.df)<-c("age_type",y.name)
  t.p<-c()
  for(i in label.df$age_type){
    t1<-wilcox.test(genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$age_type==i,y.name],genes_infor[genes_infor$Gene.Type=="Non-Cancer" & genes_infor$age_type==i,y.name])
    t.p<-c(t.p,t1$p.value)
  }
  
  Sig <- ifelse(t.p>0.05,"NS",ifelse(t.p >0.01, "*" , ifelse(t.p >0.001, "**" , "***")))
  
  s3=s3 + geom_text(data = label.df, label = Sig,col=I("blue"))+ geom_line(data = arc.df, aes(age_type+1, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+2, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+3, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+4, get(y.name)+in.x))+xlab("Gene Age class")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()
  return(s3+stat_summary(fun.data =  countFunction, col="black",geom="text", size = 4.8,fontface="bold"))
}


plotForest<-function(all_data,x.name,y.name,title.x,squaresize,spacing,colgap,method_x="pearson"){
  cor_plot1<-get_cor_for_each_cancer_all(all_data,x.name,y.name,method_x)
  x<-table(all_data$cancer)
  cols.x<-rep("black",length(cor_plot1$Cancer))
  cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
  m1<<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
  p.ff<-grid.grabExpr(print(forest(m1,fontsize = 11,squaresize=squaresize,plotwidth="5cm",colgap = colgap,spacing = spacing,layout ="JAMA" ,ff.fixed="bold",fs.fixed = 13,col.square=cols.x,xlab=title.x,ff.study="bold",fs.heading  = 13,ff.test.overall="bold",fs.xlab = 12,ff.xlab = "bold")))
  print(forest(m1,fontsize = 11,squaresize=0.9,plotwidth="5cm",colgap = colgap,spacing = spacing,layout ="JAMA" ,ff.fixed="bold",fs.fixed = 13,col.square=cols.x,xlab=title.x,ff.study="bold",fs.heading  = 13,ff.test.overall="bold",col.diamond="red",fs.xlab = 12,ff.xlab = "bold"))
  return(p.ff)
}



get_cor_age<-function(all_data,x.name,y.name){
  age<-unique(all_data[!is.na(all_data$Age_class),]$Age_class)
  all_cor<-c()
  for( age_x in age){
      all_data_sub=all_data[which(all_data$Age_class==age_x),]
      cor_plot1<-get_cor_for_each_cancer_all(all_data_sub,x.name,y.name,"pearson")
      x<-table(all_data_sub$cancer)
      cols.x<-rep("black",length(cor_plot1$Cancer))
      cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
      m1<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
      all_cor<-rbind(all_cor,cbind(cancer=m1$studlab,age_number=m1$n, age=age_x,fix=m1$TE.fixed,random=m1$TE.random,cor= m1$cor))
  }
  all_cor<-as.data.frame(all_cor)
  all_cor$age=factor(all_cor$age,levels = 1:26)
  return(all_cor)
}

get_cor_agetypeprevious<-function(all_data,x.name,y.name,y_lab,method="pearson"){
  all_data$age_type=all_data[,y_lab]
  age<-unique(all_data[!is.na(all_data$age_type),]$age_type)
  all_cor<-c()
  for( age_x in age){
    all_data_sub=all_data[which(all_data$age_type==age_x),]
    cor_plot1<-get_cor_for_each_cancer_all(all_data_sub,x.name,y.name,"pearson")
    x<-table(all_data_sub$cancer)
    cols.x<-rep("black",length(cor_plot1$Cancer))
    cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
    m1<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
    print(cor_plot1)
    all_cor<-rbind(all_cor,cbind(cancer=m1$studlab,age_number=m1$n, age=age_x,fix=m1$TE.fixed,random=m1$TE.random,cor= m1$cor,p_value=cor_plot1$Pvalue))
  }
  all_cor<-as.data.frame(all_cor)
  all_cor[,"cor"]<-as.numeric(all_cor[,"cor"])
  all_cor[,"p_value"]<-as.numeric(all_cor[,"p_value"])
  all_cor$Significance <- ifelse(all_cor[,"p_value"]>0.05,"NS",ifelse(all_cor[,"p_value"] >0.01, "*" , ifelse(all_cor[,"p_value"] >0.001, "**" , "***")))
  all_cor$Significance <- factor(all_cor$Significance ,levels = c("NS","*","**","***"))
  
  all_cor$age=factor(all_cor$age,levels = levels(all_data$age_type))
  return(all_cor)
}


get_cor_agetype<-function(all_data,x.name,y.name,y_lab,method="pearson"){
  all_data$age_type=as.factor(all_data[,y_lab])
  age<-unique(all_data[!is.na(all_data$age_type),]$age_type)
  all_cor<-c()
  print(age)
  for( age_x in age){
    all_data_sub=all_data[which(all_data$age_type==age_x),]
    cor_plot1<-get_cor_for_each_cancer_all(all_data_sub,x.name,y.name,method)
    x<-table(all_data_sub$cancer)
    cols.x<-rep("black",length(cor_plot1$Cancer))
    cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
    m1<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
    print(cor_plot1)
    print(age_x)
    all_cor<-rbind(all_cor,cbind(cancer=m1$studlab,age_number=m1$n, age=age_x,fix=m1$TE.fixed,random=m1$TE.random,cor= m1$cor,p_value=cor_plot1$Pvalue))
  }
  all_cor<-as.data.frame(all_cor)
  all_cor[,"cor"]<-as.numeric(all_cor[,"cor"])
  all_cor[,"p_value"]<-as.numeric(all_cor[,"p_value"])
  all_cor$Significance <- ifelse(all_cor[,"p_value"]>0.05,"NS",ifelse(all_cor[,"p_value"] >0.01, "*" , ifelse(all_cor[,"p_value"] >0.001, "**" , "***")))
  all_cor$Significance <- factor(all_cor$Significance ,levels = c("NS","*","**","***"))
  
    all_cor$age=factor(all_cor$age,levels = levels(all_data$age_type))
  return(all_cor)
}


get_cor_agetype_single<-function(all_data,x.name,y.name,y_lab){
  all_data$age_type=all_data[,y_lab]
  age<-unique(all_data[!is.na(all_data$age_type),]$age_type)
  all_cor<-c()
  for( age_x in age){
    all_data_sub=all_data[which(all_data$age_type==age_x),]
    cor_x<-cor.test(all_data_sub[,x.name],all_data_sub[,y.name])
    all_cor<-rbind(all_cor,cbind(age=age_x,cor= cor_x$estimate,p_value=cor_x$p.value))
  }

all_cor<-as.data.frame(all_cor)
all_cor[,"cor"]<-as.numeric(all_cor[,"cor"])
all_cor[,"p_value"]<-as.numeric(all_cor[,"p_value"])
all_cor$Significance <- ifelse(all_cor[,"p_value"]>0.05,"NS",ifelse(all_cor[,"p_value"] >0.01, "*" , ifelse(all_cor[,"p_value"] >0.001, "**" , "***")))
all_cor$Significance <- factor(all_cor$Significance ,levels = c("NS","*","**","***"))

all_cor$age=factor(all_cor$age,levels = levels(all_data$age_type))
  return(all_cor)
}


get_diff<-function(data_x,x_lab,y_lab,pos_x){
  re.result<-c()
  cancers=unique(data_x$cancer)
  data_x=data_x[which(data_x$Gene.Type!="Oncogene"),]
  for(cancer in cancers){
    temp=data_x[which(data_x$cancer==cancer),c(y_lab,x_lab)]
    colnames(temp)<-c("yvalue","x_lab")
    x.p<-wilcox.test(temp$yvalue ~ temp$x_lab)
    temp$x_lab=ifelse(temp$x_lab==pos_x,"pos","neg")
    sum_temp<-temp %>% dplyr::group_by(x_lab) %>% dplyr::summarise(dplyr::across("yvalue",list(mean=~mean(.x, na.rm = TRUE),length=length,sd=~sd(.x, na.rm = TRUE))))
    temp.sum<-t(data.frame(as.numeric(c(sum_temp[1,-1],sum_temp[2,-1]))))
    colnames(temp.sum)<-paste(rep(as.character(unlist(sum_temp[,1])),each=3),colnames(sum_temp)[-1],sep = "_")
    re.result<-rbind(re.result,data.frame(cancer=cancer,wilcoxP=x.p$p.value,temp.sum))
  }
  return(re.result)
}


getsigm6a_diff<-function(data_x,m6a_utr,n,y_lab){
  re.result<-c()
  x=table(m6a_utr[m6a_utr$m6a=="+",]$ensembl_gene_id)
  data_x$m6a=0
  data_x$m6a=x[data_x$ensembl_gene_id]
  data_x$m6a[is.na(data_x$m6a)]=0
  data_x$m6a_in<-ifelse(data_x$m6a>=n,"+","no")
  data_x$m6a_in[data_x$m6a==0]<-"-"
  data_x=data_x[which(data_x$m6a_in!="no"),]
  cancers=unique(data_x$cancer)
  for(cancer in cancers){
    temp=data_x[data_x$cancer==cancer,c(y_lab,"m6a_in")]
    colnames(temp)<-c("yvalue","m6a_in")
    temp$m6a_in=ifelse(temp$m6a_in=="+","pos","neg")
    x.p<-wilcox.test(temp$yvalue ~ temp$m6a_in)
    
    sum_temp<-temp %>% dplyr::group_by(m6a_in) %>% dplyr::summarise(dplyr::across("yvalue",list(mean=~mean(.x, na.rm = TRUE),length=length,sd=~sd(.x, na.rm = TRUE))))
    temp.sum<-t(data.frame(as.numeric(c(sum_temp[1,-1],sum_temp[2,-1]))))
    colnames(temp.sum)<-paste(rep(as.character(unlist(sum_temp[,1])),each=3),colnames(sum_temp)[-1],sep = "_")
    re.result<-rbind(re.result,data.frame(cancer=cancer,wilcoxP=x.p$p.value,temp.sum))
  }
  return(re.result)
}


getsigm6a_difflen<-function(data_x,n,y_lab,ylab2){
  re.result<-c()
  cancers=unique(data_x$cancer)
  for(cancer in cancers){
    temp=data_x[data_x$cancer==cancer,c(y_lab,ylab2)]
    colnames(temp)<-c("yvalue","m6a_in")
    temp$m6a_in=ifelse(temp$m6a_in=="+","pos","neg")
    x.p<-wilcox.test(temp$yvalue ~ temp$m6a_in)
    print(table(temp$m6a_in))
    sum_temp<-temp %>% dplyr::group_by(m6a_in) %>% dplyr::summarise(dplyr::across("yvalue",list(mean=~mean(.x, na.rm = TRUE),length=length,sd=~sd(.x, na.rm = TRUE))))
    temp.sum<-t(data.frame(as.numeric(c(sum_temp[1,-1],sum_temp[2,-1]))))
    colnames(temp.sum)<-paste(rep(as.character(unlist(sum_temp[,1])),each=3),colnames(sum_temp)[-1],sep = "_")
    re.result<-rbind(re.result,data.frame(cancer=cancer,wilcoxP=x.p$p.value,temp.sum))
  }
  return(re.result)
}




#Potential 3UTR m6A+

plotForesttwogroup<-function(sif_m6a,title_x, colgap,pos,control,spacingx =0.9){
  sif_m6a[,"pos_yvalue_mean"]=round(sif_m6a[,"pos_yvalue_mean"],2)
  sif_m6a[,"neg_yvalue_mean"]=round(sif_m6a[,"neg_yvalue_mean"],2)
  m1 <<- metacont(pos_yvalue_length,  pos_yvalue_mean,pos_yvalue_sd ,neg_yvalue_length,  neg_yvalue_mean,neg_yvalue_sd,studlab=cancer,data =sif_m6a, sm = "SMD",label.e=pos,label.c=control)
  sig1<-sif_m6a[,"wilcoxP"]
  cols.x<-rep("black",length(sif_m6a[,"wilcoxP"]))
  cols.x[as.numeric(sif_m6a[,"wilcoxP"])<0.05]="red"
  cols.x<<-cols.x ### forest only use global varaint
  title_x<<-title_x
  colgap <<- colgap
  spacingx <<-spacingx
  SMD=m1$TE
  names(SMD)<-m1$data$cancer
  p1_3utr_diff <- as.grob(~forest(m1,xlab=title_x,fontsize = 11,squaresize=0.9,spacing = spacingx, col.square =cols.x,colgap = colgap,fs.xlab = 12,ff.xlab = "bold"))
  return(list(p1=p1_3utr_diff,SMD=SMD))
}

getSMDforestcutoff<-function(data_x,m6a_3utr,y_lab,pos,control){
  num<-min(16,length(unique(m6a_3utr[,"cancer"])))
  re_all<-c()
  for(n in 1:num){
    temp_sig<-getsigm6a_diff(data_x,m6a_3utr,n,y_lab)
    print(n)
   m1 <- metacont(pos_yvalue_length,  pos_yvalue_mean,pos_yvalue_sd ,neg_yvalue_length,  neg_yvalue_mean,neg_yvalue_sd,studlab=cancer,data =temp_sig, sm = "SMD",label.e=pos,label.c=control)
   re_all<-rbind(re_all,cbind(cutoff=paste(">=",n),smd=m1$TE,cancer=m1$studlab))
     }

  re_all<-as.data.frame(re_all)
  re_all[,"cutoff"]<-factor(re_all[,"cutoff"],levels =paste(">=",1:num) )
  re_all[,"smd"]<-as.numeric(re_all[,"smd"])
  return(as.data.frame(re_all))
}




plotCorMultiple<-function(all_data_tsg_all_genes,y.name,x1,x1.name,title.x,cols){
  all_norm_plot<-c()
  for(i in 1:length(x1)){
    cor_norm_exp<-get_cor_for_each_cancer_all(all_data_tsg_all_genes,y.name,x1[i],"pearson")
    all_norm_plot<-rbind(all_norm_plot,cbind(type_x=x1.name[i],cor_norm_exp))
  }
  #col=c("#FFC125","#00B2EE","#EE6363", "#458B00", "purple")
  all_norm_plot$Pvalue<-as.numeric(all_norm_plot$Pvalue)
  cor_matrix1= aggregate(.~type_x,all_norm_plot[,-c(2:3)],median)
  print(cor_matrix1)
  cor_matrix2= aggregate(.~type_x,all_norm_plot[,-c(2:3)],mad)
  print(cor_matrix2)
  all_norm_plot$type_x<-factor(all_norm_plot$type_x,levels = x1.name)
  all_norm_plot$Pvalue<--log10(as.numeric(all_norm_plot$Pvalue))
  all_norm_plot$Significance<- ifelse(round(10^-all_norm_plot$Pvalue,5)>0.05,"NS", "Pvalue<0.05")
  if(length(unique(all_norm_plot$Significance))>1){
    cor_plot2= ggplot(all_norm_plot,aes(y=Cancer,x=cor,col=type_x))+  geom_point(aes(size=Pvalue,shape=Significance))+scale_shape_manual("Significance",values=c(1,19),label=c("No","yes"))
  }else{
    cor_plot2= ggplot(all_norm_plot,aes(y=Cancer,x=cor,col=type_x))+  geom_point(aes(size=Pvalue))
    
  }
  cor_plot2=cor_plot2+labs(size = "-Log10(P-value)")+scale_color_manual("",values=cols,label=x1.name)+ggtitle(title.x)+  scale_size_area(labels=c("1","2",">5"),max_size=7,limits=c(0,16),breaks=c(1,2,5))+guides(color=guide_legend(nrow=2,byrow=TRUE),size=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+xlab("Pearson's r")+geom_hline(yintercept = 0.4,size=1.5)+geom_vline(xintercept = 0,size=1)+theme_hd_minimal2()+ylab("Dataset")+guides(size = guide_legend(order = 1),col = guide_legend(order = 3))
  return(cor_plot2)
}


plotcancercount<-function(genes_infor,y.name,y_lab,pos){
  countFunction <- function(x){
    return(data.frame(y=pos,label=round(length(x),2)))}
  in.x1<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.7
  in.x<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.5
  comp3=list(c("0", "1"),c("1", "2"),c("2", ">2"))
  s3=ggplot(genes_infor[!is.na(genes_infor$cancer_count2) & genes_infor$Gene.Type!="Oncogene" ,],aes(x=as.factor(cancer_count2) ,y=get(y.name)))  + geom_signif(comparisons = comp3,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),outlier.colour = NA,aes(fill=Gene.Type),width=0.4)+ylab(y_lab)
  label.df <- data.frame(cancer_count2= c("0", "1","2", ">2"),y.n = rep(in.x1,4))
  colnames(label.df)<-c("cancer_count2",y.name)
  r <- .1
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- r*2*in.x1/3 * sin(t)
  arc.df <- data.frame(cancer_count2 = x, x= y)
  colnames(arc.df)<-c("cancer_count2",y.name)
  t.p<-c()
  for(i in label.df$cancer_count2){
    t1<-wilcox.test(genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$cancer_count2==i,y.name],genes_infor[genes_infor$Gene.Type=="Non-Cancer" & genes_infor$cancer_count2==i,y.name])
    t.p<-c(t.p,t1$p.value)
  }
  
  Sig <- ifelse(t.p>0.05,"NS",ifelse(t.p >0.01, "*" , ifelse(t.p >0.001, "**" , "***")))
  
  s3=s3 + geom_text(data = label.df, label = Sig,col=I("blue"))+ geom_line(data = arc.df, aes(cancer_count2+1, get(y.name)+in.x))+ geom_line(data = arc.df, aes(cancer_count2+2, get(y.name)+in.x))+ geom_line(data = arc.df, aes(cancer_count2+3, get(y.name)+in.x))+ geom_line(data = arc.df, aes(cancer_count2+4, get(y.name)+in.x))+xlab("Gene Age class")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()
  return(s3+stat_summary(fun.data =  countFunction, col="black",geom="text", size = 4.8,fontface="bold"))
}



plotSubSig<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black"){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group,col=I(cols_x)))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(length(cols1)>1){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"]="NS"
  if(size_x==0.5){
     stat.test=stat.test %>% filter(p.adj.signif!="ns")
     #print(stat.test)
  }
  if(dim(stat.test)[1]>0){
     p1=p1 +  geom_signif(tip_length = 0.01,xmin=stat.test$xmin,size=size_x, xmax=stat.test$xmax, annotations=stat.test$p.adj.signif, y_position=max(stat.test$y.position),fontface="bold") 
 
  }
  return(p1)
}

plotSubSigcg<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black"){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  x.data$x<-as.factor(x.data$x)
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(col=group))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(length(cols1)>1){
    p1=p1+scale_color_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"]="NS"
  if(size_x==0.5){
    stat.test=stat.test %>% filter(p.adj.signif!="ns")
    #print(stat.test)
  }
  if(dim(stat.test)[1]>0){
    p1=p1 +  geom_signif(tip_length = 0.01,xmin=stat.test$xmin,size=size_x, xmax=stat.test$xmax, annotations=stat.test$p.adj.signif, y_position=max(stat.test$y.position),fontface="bold",textsize = 2.5,angle=30) 
    
  }
  return(p1)
}


plotSubSigOneSide<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black",side){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group,col=I(cols_x)))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(length(cols1)>1){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group,alternative=side)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"]="NS"
  if(size_x==0.5){
    stat.test=stat.test %>% filter(p.adj.signif!="ns")
    #print(stat.test)
  }
  if(dim(stat.test)[1]>0){
    p1=p1 +  geom_signif(tip_length = 0.01,xmin=stat.test$xmin,size=size_x, xmax=stat.test$xmax, annotations=stat.test$p.adj.signif, y_position=max(stat.test$y.position),step_increase = 0.1,fontface="bold") 
    
  }
  return(p1)
}



plotSubSigall4<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black",add_x=1){
  p<-list()
  cancer=unique(x.data$cancer)
  print(cancer)
  i=1
  max_x<-max(x.data[,ylab1])+add_x
  if(length(cancer) %% 4!=0){
    nr<-length(cancer) %/% 4 
  }else{nr<-length(cancer) %/% 4 -1 }
  j=1
  while(i <length(cancer)){
    print(length(cancer))
    if(i%/% 4 !=nr){
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])
      p1=p1+theme(axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")+ scale_y_continuous(limits=c(0,max_x+add_x),expand = c(0, 0))
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+add_x),expand = c(0, 0))
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+add_x),expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+add_x),expand = c(0, 0))
     p_1<-ggarrange(p1,p2,p3,p4,nrow=1,ncol=4,widths = c(1.2,1,1,1))
      
    }else{
      
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])+ scale_y_continuous( limits=c(0,max_x+add_x),expand = c(0, 0))
      p1=p1+xlab("")
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+add_x),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+add_x),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+add_x),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p_1<-ggarrange(p1,p2,p3,p4,nrow=1,ncol=4,widths = c(1.2,1,1,1),common.legend = T,legend = "bottom")
      
    }
    print(i)
   print(p_1)
    p[[j]]=p_1
    i=i+4
    j=j+1
    print(i)
  }
  p_all<-ggarrange(p[[1]],p[[2]],p[[3]],nrow=3,ncol=1,common.legend = T,legend="bottom",heights = c(1,1,1.8))
  return(p_all)  
}

plotSubSigall5<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black"){
  p<-list()
  cancer=unique(x.data$cancer)
  i=1
  max_x<-max(x.data[,ylab1])
  if(length(cancer) %% 5!=0){
    nr<-length(cancer) %/% 5 
  }else{nr<-length(cancer) %/% 5 -1 }
  j=1
  while(i <length(cancer)){
    print(length(cancer))
    if(i%/% 5 !=nr){
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])
      p1=p1+theme(axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p5<-plotSubSig(x.data[x.data$cancer==cancer[i+4],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+4])
      p5=p5+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p_1<-ggarrange(p1,p2,p3,p4,p5,nrow=1,ncol=5,widths = c(1.2,1,1,1,1))
      
    }else{
      
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])
      p1=p1+xlab("")
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+add_x),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+add_x),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+add_x),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p5<-plotSubSig(x.data[x.data$cancer==cancer[i+4],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+4])
      p5=p5+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous(limits=c(0,max_x+add_x),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p_1<-ggarrange(p1,p2,p3,p4,p5,nrow=1,ncol=5,widths = c(1.2,1,1,1,1),common.legend = T,legend = "bottom")
      
         }
    print(i)
   print(p_1)
    p[[j]]=p_1
    i=i+5
    j=j+1
    print(i)
  }
  p_all<-ggarrange(p[[1]],p[[2]],p[[3]],nrow=3,ncol=1,common.legend = T,legend="bottom",heights = c(1,1,1.8))
  return(p_all)  
}

plotSubSig_flip<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,cols_x="black"){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group,col=I(cols_x)))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(cols1!=""){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"]="NS"
  if(dim(stat.test)[1]>0){
    p1=p1 + coord_flip()+annotate("text",1:length(stat.test$p.adj.signif),max(stat.test$y.position),label=stat.test$p.adj.signif,col="red")
    
  }
  return(p1)
}


plotSubSig_pvalue<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(cols1!=""){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test=stat.test %>% filter(p.adj.signif!="ns")
  if(dim(stat.test)[1]>0){
    p1=p1 +  geom_signif(tip_length = 0.01,xmin=stat.test$xmin, xmax=stat.test$xmax, annotations=stat.test$p.adj, y_position=stat.test$y.position,col="black",textsize = 4.8,fontface="bold") 
    
  }
  return(p1)
}

get_main_transcript<-function(get_gene_pos_temp){
  # get the transcript which is principal1
  main_trans<-subset( get_gene_pos_temp,APPRIS.annotation=="principal1")
  genes=unique(get_gene_pos_temp$ensembl_gene_id)
  genes=genes[! genes %in% main_trans$ensembl_gene_id] ## genes which have no principal transcript
  get_gene_pos_temp=get_gene_pos_temp[!get_gene_pos_temp$ensembl_gene_id  %in% main_trans$ensembl_gene_id,]
  get_gene_pos_temp<-rbind(main_trans,get_gene_pos_temp) # genes with  principal1 and genes that do not have principal1
  max_trans_length<-aggregate(get_gene_pos_temp$trans_length,list(get_gene_pos_temp$ensembl_gene_id),function(x){max(x,na.rm = T)}) # get the max transcript length
  f_select=paste(get_gene_pos_temp$ensembl_gene_id,get_gene_pos_temp$trans_length,sep="_") %in% paste(max_trans_length[,1],max_trans_length[,2],sep="_")
  select_gene_trans<-get_gene_pos_temp[f_select,]
  # for genes with transcript which have the same length randomly select one
  return(select_gene_trans[!duplicated(select_gene_trans$ensembl_gene_id),])
}


plot3utrFreq<-function(genes_infor,temp_gene){
  
  
  temp_gene_trans<-aggregate(.~ensembl_gene_id+ensembl_transcript_id,data=temp_gene,max)
  
  temp_trans<-unique(unique(temp_gene[,1:2]))
  trans_count=table(temp_trans$ensembl_gene_id)
  
  temp_gene_uniq<-unique(temp_gene_trans[,-2])
  utr3_count<-table(temp_gene_uniq$ensembl_gene_id)
  genes_infor$utr3_count<-as.numeric(utr3_count[genes_infor$ensembl_gene_id])
  genes_infor$trans_count<-as.numeric(trans_count[genes_infor$ensembl_gene_id])
  genes_infor$utr3_freq<-genes_infor$utr3_count/genes_infor$trans_count
  
  
  temp_data=c()
  for(i in 0:10){
    temp_data<-rbind(temp_data,cbind(transript_count=paste(">",i),genes_infor[genes_infor$trans_count>i,]))
  }  
  temp_data[,"transript_count"]<-factor(temp_data[,"transript_count"],levels = paste(">",0:10))
  return(temp_data)
}



plot_histoligical<-function(genes_infor,col_x,comps,ylab){
  overlap.data=rbind(cbind(cancer="Carcinoma",genes_infor[which(genes_infor$Carcinoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Leukemia",genes_infor[which(genes_infor$Leukemia!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Lymphoma",genes_infor[which(genes_infor$Lymphoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Sarcoma",genes_infor[which(genes_infor$Sarcoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Myeloma",genes_infor[which(genes_infor$Myeloma!="" & genes_infor$Gene.Type=="TSG"),col_x]))
  overlap.data<-as.data.frame(overlap.data)
  overlap.data$cancer<-factor(overlap.data$cancer,levels=c("Myeloma","Carcinoma","Sarcoma","Leukemia","Lymphoma"))
  colnames(overlap.data)<-c("cancer_types",col_x)
  overlap.data[,col_x]<-as.numeric(overlap.data[,col_x])
  
  
  comp<-comps
    #list(c("Carcinoma","Leukemia"),c("Carcinoma","Lymphoma"),c("Carcinoma","Sarcoma"),c("Leukemia","Lymphoma"),c("Leukemia","Sarcoma"),c("Lymphoma","Sarcoma"),c("Myeloma","Leukemia"),c("Myeloma","Lymphoma"),c("Myeloma","Sarcoma"),c("Carcinoma","Myeloma"))
  #,test.args = c(alternative = "greater")
  px2_dif=ggplot(data = overlap.data,aes(x=cancer_types,y=get(col_x),fill=cancer_types))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold") +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"),axis.title.x =element_blank(),axis.text.x = element_text(angle =90))# 30,vjust=0.8
  px2_dif=px2_dif#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold",angle=90)
  print(px2_dif)
  genes_infor$cancer_count2=as.character(genes_infor$cancer_count2)
  genes_infor$cancer_count2[genes_infor$cancer_count2 %in% c("2",">2")]=">1"
  genes_infor$cancer_count2<-factor(genes_infor$cancer_count2,levels = c("1",">1"))
  px3_dif=ggplot(data = genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$cancer_count2 %in% c("1",">1"), ],aes(x=cancer_count2,y=get(col_x),fill=cancer_count2))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = list(c("1",">1")), test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold",test.args = c(alternative = "less")) +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"))
  px3_dif=px3_dif+xlab("Number of cancer type")+scale_fill_manual(values = c("#8EE5EE", "#20B2AA"))
  print(px3_dif)
  #+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  
  pall5<-list(p1=px2_dif,p2=px3_dif)
  return(pall5)
}



CutFile <- function(Input.name, Output.name, Annotation.name, Exception.symbols = "!Ser" ){
  Exception.symbols=toupper(Exception.symbols)
  if (file.exists(Output.name)){
    cat("Cut-File is already exited!\n")
    return(Output.name)
  }
  Input <- file(Input.name, "r")
  cat("",file=Output.name,append = FALSE)
  cat("",file=Annotation.name,append = FALSE)
  N.Char <- nchar(Exception.symbols)
  while (TRUE) {
    Line <- readLines(Input, n = 1)
    if (length(Line) == 0 ) {
      break
    }
    if ( toupper(substr(Line, start = 1, stop = N.Char)) == Exception.symbols ){
      cat(Line,"\n",file=Annotation.name,append = TRUE)
      next
    }
    if (!nzchar(Line)){
      next
    }
    cat(Line,"\n",file=Output.name,append = TRUE)
  }
  close(Input)
  return(Output.name)
}

getComp<-function(x){
  re_list<-list()
  for(i in 1:(length(x)-1)){
    for(j in i:(length(x)-1)){
      re_list[[paste(i,j)]]=c(x[i],x[j+1])
    }
    
  }
  return(re_list)
}


plot_histoligical<-function(genes_infor,col_x,comps,ylab){
  overlap.data=rbind(cbind(cancer="Carcinoma",genes_infor[which(genes_infor$Carcinoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Leukemia",genes_infor[which(genes_infor$Leukemia!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Lymphoma",genes_infor[which(genes_infor$Lymphoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Sarcoma",genes_infor[which(genes_infor$Sarcoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Myeloma",genes_infor[which(genes_infor$Myeloma!="" & genes_infor$Gene.Type=="TSG"),col_x]))
  overlap.data<-as.data.frame(overlap.data)
  overlap.data$cancer<-factor(overlap.data$cancer,levels=c("Myeloma","Carcinoma","Sarcoma","Leukemia","Lymphoma"))
  colnames(overlap.data)<-c("cancer_types",col_x)
  overlap.data[,col_x]<-as.numeric(overlap.data[,col_x])
  
  
  comp<-comps
  #list(c("Carcinoma","Leukemia"),c("Carcinoma","Lymphoma"),c("Carcinoma","Sarcoma"),c("Leukemia","Lymphoma"),c("Leukemia","Sarcoma"),c("Lymphoma","Sarcoma"),c("Myeloma","Leukemia"),c("Myeloma","Lymphoma"),c("Myeloma","Sarcoma"),c("Carcinoma","Myeloma"))
  #,test.args = c(alternative = "greater")
  px2_dif=ggplot(data = overlap.data,aes(x=cancer_types,y=get(col_x),fill=cancer_types))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold") +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"),axis.title.x =element_blank(),axis.text.x = element_text(angle =90))# 30,vjust=0.8
  px2_dif=px2_dif#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold",angle=90)
  print(px2_dif)
  genes_infor$cancer_count2=as.character(genes_infor$cancer_count2)
  genes_infor$cancer_count2[genes_infor$cancer_count2 %in% c("2",">2")]=">1"
  genes_infor$cancer_count2<-factor(genes_infor$cancer_count2,levels = c("1",">1"))
  px3_dif=ggplot(data = genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$cancer_count2 %in% c("1",">1"), ],aes(x=cancer_count2,y=get(col_x),fill=cancer_count2))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = list(c("1",">1")), test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold",test.args = c(alternative = "less")) +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"))
  px3_dif=px3_dif+xlab("Number of cancer type")+scale_fill_manual(values = c("#8EE5EE", "#20B2AA"))
  print(px3_dif)
  #+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  
  pall5<-list(p1=px2_dif,p2=px3_dif)
  return(pall5)
}

plot_classified_hist<-function(draws,cut_off,lables,x_lab){
  dens <- density(draws)
  dd <- with(dens,data.frame(x,y))
  cols_x=brewer.pal(n = 8, name = "Set2")
  # cols_x=lables
  s1=ggplot(data=dd,aes(x,y))+geom_line()+
    geom_ribbon(data=subset(dd,x<=cut_off[2]),aes(ymax=y,fill=I(cols_x[1])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd, x>cut_off[2] & x<=cut_off[3]),aes(ymax=y,fill=I(cols_x[2])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd, x>cut_off[3] & x<=cut_off[4]),aes(ymax=y,fill=I(cols_x[3])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd, x>cut_off[4] & x<=cut_off[5]),aes(ymax=y,fill=I(cols_x[4])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd,x>cut_off[5]),aes(ymax=y,fill=I(cols_x[5])),ymin=0,
                colour=NA) + labs(fill = "Groups")+ylab("Density")+xlab(x_lab)
  
  print(s1)
}


plot_scatter<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,file_name){
  jpeg(file_name,width = 1300,height = 1200,res = 120,quality = 100)
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.1),color="cancer",size=0.6,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,label.y =label_y_pos)+facet_wrap(~cancer,nrow=3)  + theme(legend.position = "none",strip.text = element_text(size=25),axis.title =element_text(size=20))
  print(sp)
  dev.off()
}

plot_scatter_3<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=1,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,aes(color = get(col_x)),fontface="bold",output.type="text")  + theme(strip.text = element_text(size=25),axis.title =element_text(size=15))
  sc1<-cor.test(all_data_tsg[,y_lab],all_data_tsg[,x_lab])
  pvalue=sc1$p.value
  if(pvalue==0){pvalue="2.2e-16"}
  data_text_1= paste("All genes: R=",round(sc1$estimate,3),"   p < ", formatC(pvalue, format = "e", digits = 2),sep="")
  sp=sp+annotate(geom="text", x=label_x_pos-0.05, y=label_y_pos,size=4, label=data_text_1, color="black")
  return(sp+scale_color_manual(values=c("#1874CD","#EE2C2C")))
}

plot_scatter_4<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=1,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,aes(color = get(col_x)),fontface="bold",output.type="text")  + theme(strip.text = element_text(size=25),axis.title =element_text(size=15))
  sc1<-cor.test(all_data_tsg[,y_lab],all_data_tsg[,x_lab])
  pvalue=sc1$p.value
  if(pvalue==0){pvalue="2.2e-16"}
  data_text_1= paste("All genes: R=",round(sc1$estimate,3),"   p < ", pvalue,sep="")
  sp=sp+annotate(geom="text", x=label_x_pos+0.1, y=label_y_pos,size=4, label=data_text_1, color="black")
  
  return(sp+scale_color_manual(values=c("#1874CD","#EE2C2C")))
}


plot_scatter_5<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  all_data_tsg$Gene.Type=factor(all_data_tsg$Gene.Type,levels=c("Non-Cancer","TSG"))
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=1,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,aes(color = get(col_x)),fontface="bold",output.type="text")  + theme(strip.text = element_text(size=25),axis.title =element_text(size=15))
  sc1<-cor.test(all_data_tsg[,y_lab],all_data_tsg[,x_lab])
  pvalue=sc1$p.value
  if(pvalue==0){pvalue="2.2e-16"}
  data_text_1= paste("All genes: R=",round(sc1$estimate,3),"   p < ", pvalue,sep="")
  sp=sp+facet_wrap(~cancer,nrow=5)+theme_hd()
  
  return(sp+scale_color_manual(values=c("#1874CD","#EE2C2C")))
}



plot_scatter_2<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.8),color=col_x,size=1.5,shape=16,
                  add = "reg.line",  add.params = list(color = "red"),  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}


plot_scatter_2kendall<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.6),color=col_x,size=1,shape=16,
                  add = "reg.line",  add.params = list(color = "blue"),  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "kendall",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}






plot_scatter_6<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=2,
                  add = "reg.line", add.params = list(color = "red"), # Add regressin line#size=1.5,shape=16,
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}


plot_scatter_7<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.1),color=col_x,size=0.7,shape=16,
                  add = "reg.line",  # Add regressin line#
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}



plot_scatter_2speaman<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.8),color=col_x,size=1.5,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "spearman",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}
getSigRatio<-function(genes_infor,ylab_x,xlab_x,class_x,title_y,title_x,side="two.sided"){
  genes_infor<-genes_infor[!is.na(genes_infor[,xlab_x]),]
  genes_infor[,xlab_x]=as.factor(genes_infor[,xlab_x])
  motif_count=as.data.frame(table(genes_infor[which(genes_infor[,ylab_x]==class_x),xlab_x],genes_infor[which(genes_infor[,ylab_x]==class_x),]$Gene.Type))
  motif_all=as.data.frame(table(genes_infor[,xlab_x],genes_infor$Gene.Type))
  merge_motif=merge(motif_count,motif_all,by=c("Var1","Var2"))
  genes_infor[,xlab_x]<-factor(genes_infor[,xlab_x])
  pvalue_fish<-c()
  print(ylab_x)
  for(x in levels(genes_infor[,xlab_x])){
    temp=genes_infor[as.character(genes_infor[,xlab_x])==x& genes_infor$Gene.Type!="Oncogene",]
    temp_x<-fisher.test( table(temp$Gene.Type,temp[,ylab_x]),alternative = side)
    pvalue_fish<-c(pvalue_fish,temp_x$p.value)
    print(x)
    print(table(temp$Gene.Type,temp[,ylab_x]))
    print(temp_x$p.value)
  }
  Sig <- ifelse(pvalue_fish>0.05,"NS",ifelse(pvalue_fish>0.01, "*" ,  ifelse(pvalue_fish >0.001, "**" , "***")))
  n=length(unique(merge_motif$Var1))  
  merge_motif$ratio=merge_motif$Freq.x/merge_motif$Freq.y
  p_motif<-ggplot(data=merge_motif[merge_motif$Var2!="Oncogene",],aes(x=Var1,y=ratio,fill=Var2))+geom_bar(stat = "identity", position=position_dodge())+scale_fill_manual("Gene.Type",values=c("#1874CD","#EE2C2C"))+xlab(title_x)+ylab(title_y)+theme_hd()+theme(axis.text.x = element_text(angle = 90))+annotate("text",x=1:n,y=max(merge_motif$ratio)+0.1,label=Sig,size=5)+annotate("segment",x = c(1:n)-0.25, xend = c(1:n)+0.25, y = max(merge_motif$ratio)+0.05, yend = max(merge_motif$ratio)+0.05)
  
  return(p_motif)
}




GroupbyMinMaxBed<-function(pos_data){
  x=unique(pos_data[,4])
  re_pos<-c()
  for(i in x){
    temp=pos_data[pos_data[,4]==i,]
    if(dim(temp)[1]==1){
      temp_x=temp
    }else{
      temp_x=c(temp[1,1],min(temp[,2]),max(temp[,3]),temp[1,4])
    }
    names(temp_x)<-colnames(pos_data)
    re_pos<-rbind(re_pos,temp_x)
  }
  re_pos<-as.data.frame(re_pos)
  colnames(re_pos)<-colnames(pos_data)
  return(re_pos)
}



getUtrSeqDNA<-function(t1,t2,name1,pos_UTR,genes_infor){
  pos_UTR$strand[pos_UTR$strand=="1"]="+"
  pos_UTR$strand[pos_UTR$strand=="-1"]="-"
  pos_UTR$chromosome_name=paste("chr",pos_UTR$chromosome_name,sep="")
  seq_data_utr<-getSeq(Hsapiens, pos_UTR$chromosome_name,
                       start=pos_UTR[,t1],end=pos_UTR[,t2],strand=pos_UTR$strand)
  seq_data_df_utr=as.data.frame(seq_data_utr)
  all_data_seq_utr=cbind(pos_UTR,seq_data_df_utr)
  
  select_seq=all_data_seq_utr[all_data_seq_utr$ensembl_transcript_id %in% genes_infor[genes_infor$Gene.Type==name1,]$genes_main_transcript,]
  uniq_trans=unique(select_seq$ensembl_transcript_id)
  fasta_seq=c()
  for(i in uniq_trans){
    temp_seq=select_seq[select_seq$ensembl_transcript_id ==i,]
    temp_p=c()
    for(temp in temp_seq$x){temp_p=paste(temp_p,temp,sep="")}
    if(nchar(temp_p)<3){next}
    fasta_seq<-rbind(fasta_seq,paste(">",i,"_",temp_seq[1,1],sep=""))
    fasta_seq<-rbind(fasta_seq,temp_p)
  }
  write.table(fasta_seq,file=paste("./result/",name1,"_",t1,"_DNA.fasta",sep=""),quote=F,col.names=F,row.names = F)
}


get_utr_seqRaw<-function(t1,t2,all_data_genes,genes_infor){
  all_trans_pos=all_data_genes$get_gene_pos[all_data_genes$get_gene_pos$ensembl_transcript_id %in% genes_infor$genes_main_transcript,c("ensembl_gene_id","ensembl_transcript_id",t1,t2,"strand","chromosome_name")]
  pos_UTR<-unique(all_trans_pos[!is.na(all_trans_pos[,t1]),])
  pos_UTR<-checkPosition(pos_UTR)
  pos_UTR$strand[pos_UTR$strand=="1"]="+"
  pos_UTR$strand[pos_UTR$strand=="-1"]="-"
  pos_UTR$chromosome_name=paste("chr",pos_UTR$chromosome_name,sep="")
  seq_data_utr<-getSeq(Hsapiens, pos_UTR$chromosome_name,
                       start=pos_UTR[,t1],end=pos_UTR[,t2],strand=pos_UTR$strand)
  seq_data_df_utr=as.data.frame(seq_data_utr)
  all_data_seq_utr=cbind(pos_UTR,seq_data_df_utr)
  select_seq=all_data_seq_utr[all_data_seq_utr$ensembl_transcript_id %in% genes_infor$genes_main_transcript,]
  uniq_trans=unique(select_seq$ensembl_transcript_id)
  fasta_seq=c()
  len_x=c()
  for(i in uniq_trans){
    temp_seq=select_seq[select_seq$ensembl_transcript_id ==i,]
    temp_p=c()
    for(temp in temp_seq$x){temp_p=paste(temp_p,temp,sep="")}
    if(nchar(temp_p)<1){next}
    fasta_seq<-rbind(fasta_seq,temp_p)
    len_x<-c(len_x,nchar(temp_p))
  }
  return(data.frame(trans=uniq_trans,len=len_x,seq=fasta_seq))
}





GroupbyMinMax<-function(UTR5_pos){
  x=unique(UTR5_pos$ensembl_transcript_id)
  re_pos<-c()
  for(i in x){
    temp=UTR5_pos[UTR5_pos$ensembl_transcript_id==i,]
    if(dim(temp)[1]==1){
      temp_x=temp
    }else{
      temp_x=c(temp[1,1:2],min(temp[,3]),max(temp[,4]),temp[1,5:6])
    }
    names(temp_x)<-colnames(UTR5_pos)
    re_pos<-rbind(re_pos,temp_x)
  }
  re_pos<-as.data.frame(re_pos)
  colnames(re_pos)<-colnames(UTR5_pos)
  return(re_pos)
}



checkPosition<-function(UTR5_pos){
  x=unique(UTR5_pos$ensembl_transcript_id)
  re_pos<-c()
  num=0
  for(i in x){
    temp=UTR5_pos[UTR5_pos$ensembl_transcript_id==i,]
    if(dim(temp)[1]==1){
      temp_x=temp
    }else{
      if(temp[1,5]=="1"){
        temp_order<-order(as.numeric(temp[,3]),decreasing = F)
        temp_x<-temp[order(as.numeric(temp[,3]),decreasing = F),]
        if(sum(temp_order ==1:length(temp_order))!=length(temp_order)){
          num=num+1
        }
      }else{
        temp_order<-order(as.numeric(temp[,3]),decreasing = T)
        temp_x<-temp[order(as.numeric(temp[,3]),decreasing = T),]
        if(sum(temp_order ==1:length(temp_order))!=length(temp_order)){
          num=num+1
          print(num)
        }
      }
      
    }
    names(temp_x)<-colnames(UTR5_pos)
    re_pos<-rbind(re_pos,temp_x)
  }
  re_pos<-as.data.frame(re_pos)
  colnames(re_pos)<-colnames(UTR5_pos)
  print(num)
  return(re_pos)
}


getCompositionAll<-function(pos_UTR,seq_data_5utr,genes_infor,n){
  genes_infor<-genes_infor[genes_infor$Gene.Type!="Oncogene",]
  bas_freq=oligonucleotideFrequency(seq_data_5utr,as.prob = T,width = n,step = n)
  bas_has<-aggregate(bas_freq,list(pos_UTR$ensembl_gene_id),sum)
  rownames(bas_has)<-bas_has[,1]
  no_5utr_genes<-genes_infor[!genes_infor$ensembl_gene_id %in% rownames(bas_has),"ensembl_gene_id"]
  bas_has=bas_has[,-1]
  no_5utr_Data<-matrix(0,nrow = length(no_5utr_genes),ncol=length(bas_has))
  rownames(no_5utr_Data)<-no_5utr_genes
  colnames(no_5utr_Data)<-colnames(bas_has)
  bas_has<-rbind(bas_has,no_5utr_Data)
  bas_has[bas_has==0]=0
  bas_has[bas_has>0]=1
  x_has<-melt(as.matrix(bas_has))
  x_has$Gene.Type=genes_infor[as.character(x_has$Var1),]$Gene.Type
  x_has$age_type2<-genes_infor[as.character(x_has$Var1),]$age_type2
  x_t<-as.data.frame(table(x_has$Gene.Type,x_has$value,x_has$Var2))
  print(x_t)
  all_re<-c()
  for(i in unique(x_t$Var3)){
    all_re<-rbind(all_re,c(i,funFisher(x_t[x_t$Var3==i,"Freq"])))
  }
  all_re<-as.data.frame(all_re)
  all_re[,2]<-as.numeric(all_re[,2])
  all_re[,3]<-as.numeric(all_re[,3])
  return(all_re)
}


getComposition<-function(pos_UTR,seq_data_5utr,genes_infor,n){
  genes_infor<-genes_infor[genes_infor$Gene.Type!="Oncogene",]
  bas_freq=oligonucleotideFrequency(seq_data_5utr,as.prob = T,width = n,step = n)
  bas_has<-aggregate(bas_freq,list(pos_UTR$ensembl_gene_id),sum)
  rownames(bas_has)<-bas_has[,1]
  no_5utr_genes<-genes_infor[!genes_infor$ensembl_gene_id %in% rownames(bas_has),"ensembl_gene_id"]
  bas_has=bas_has[,-1]
  bas_has[bas_has==0]=0
  bas_has[bas_has>0]=1
  x_has<-melt(as.matrix(bas_has))
  x_has$Gene.Type=genes_infor[as.character(x_has$Var1),]$Gene.Type
  x_has$age_type2<-genes_infor[as.character(x_has$Var1),]$age_type2
  x_t<-as.data.frame(table(x_has$Gene.Type,x_has$value,x_has$Var2))
  print(x_t)
  all_re<-c()
  for(i in unique(x_t$Var3)){
    all_re<-rbind(all_re,c(i,funFisher(x_t[x_t$Var3==i,"Freq"])))
  }
  all_re<-as.data.frame(all_re)
  print(head(all_re))
  all_re[,2]<-as.numeric(all_re[,2])
  all_re[,3]<-as.numeric(all_re[,3])
  return(all_re)
}

funFisher<-function(x){
  x_p<-fisher.test(matrix(x,2,byrow = T))
  return(c(x_p$p.value, x_p$estimate))
}
  
  
  
getSeqMotif<-function(motif_x,pos_UTR,seq_data_5utr,genes_infor,all_data_tsg,all_data_tsg_filtered){
  TATA_pos=vmatchPattern(motif_x,seq_data_5utr)
  TATA_count_all=as.numeric(sapply(TATA_pos@ends, length))
  sum_tata<-aggregate(TATA_count_all,list(pos_UTR$ensembl_gene_id),sum)
  rownames(sum_tata)<-sum_tata[,1]
  genes_infor$tata_count<-sum_tata[genes_infor$ensembl_gene_id,2]
  genes_infor$tata="-"
  genes_infor$tata[genes_infor$tata_count>0]="+"
  print(table(genes_infor$tata))
  pMotif1<-getSigRatio(genes_infor,"tata","age_type2","+","Ratio of genes with  Motif","Gene Age Groups(million years)")
  print(pMotif1)
  all_data_tsg$tata=genes_infor[all_data_tsg$ensembl_gene_id,"tata"]
  genes_infor1=genes_infor[genes_infor$Gene.Type!="Oncogene",]
  print(table(genes_infor1$Gene.Type,genes_infor1$tata))
  print(fisher.test( table(genes_infor1$Gene.Type,genes_infor1$tata)))
  sig_rbp<-get_diff(all_data_tsg,"tata","norm_exp","+")
  pff_rbp_all0<-plotForesttwogroup(sig_rbp,"normal expression in TCGA","4.2mm","+","-")
  ggarrange(pff_rbp_all0$p1)
  
  all_data_tsg_filtered$tata="-"
  all_data_tsg_filtered$tata=genes_infor[all_data_tsg_filtered$ensembl_gene_id,"tata"]
  sig_rbp<-get_diff(all_data_tsg_filtered,"tata","diff_exp","+")
  pff_rbp_all<-plotForesttwogroup(sig_rbp,"Abs(Log2FC) of Down-regulated TSG in TCGA","4.2mm","+","-")
  ggarrange(pff_rbp_all$p1)
  return(list(p1=pff_rbp_all0$p1,p2=pff_rbp_all$p1))
}  

getSeqMotifMulti<-function(motif_xs,pos_UTR,seq_data_5utr,genes_infor,all_data_tsg,all_data_tsg_filtered){
  motif_genes<-c()
  for( motif_x in motif_xs){
       
  TATA_pos=vmatchPattern(motif_x,seq_data_5utr)
  TATA_count_all=as.numeric(sapply(TATA_pos@ends, length))
  sum_tata<-aggregate(TATA_count_all,list(pos_UTR$ensembl_gene_id),sum)
  rownames(sum_tata)<-sum_tata[,1]
  motif_genes<-rbind(motif_genes,sum_tata[sum_tata[,2]>0,])
   }


  genes_infor$tata=ifelse(genes_infor$ensembl_gene_id %in% motif_genes[,1],"+","-")
  print(table(genes_infor$tata))
  pMotif1<-getSigRatio(genes_infor,"tata","age_type2","+","Ratio of genes with  Motif","Gene Age Groups(million years)")
  print(pMotif1)
  all_data_tsg$tata=genes_infor[all_data_tsg$ensembl_gene_id,"tata"]
  genes_infor1=genes_infor[genes_infor$Gene.Type!="Oncogene",]
  print(table(genes_infor1$Gene.Type,genes_infor1$tata))
  print(fisher.test( table(genes_infor1$Gene.Type,genes_infor1$tata)))
  sig_rbp<-get_diff(all_data_tsg,"tata","norm_exp","+")
  pff_rbp_all0<-plotForesttwogroup(sig_rbp,"normal expression in TCGA","4.2mm","+","-")
  ggarrange(pff_rbp_all0$p1)
  
  all_data_tsg_filtered$tata="-"
  all_data_tsg_filtered$tata=genes_infor[all_data_tsg_filtered$ensembl_gene_id,"tata"]
  sig_rbp<-get_diff(all_data_tsg_filtered,"tata","diff_exp","+")
  pff_rbp_all<-plotForesttwogroup(sig_rbp,"Abs(Log2FC) of Down-regulated TSG in TCGA","4.2mm","+","-")
  ggarrange(pff_rbp_all$p1)
  return(list(p1=pff_rbp_all0$p1,p2=pff_rbp_all$p1))
}  


getG4finder<-function(pos_UTR,t1,t2,extend_bp){
  pos_UTR$strand[pos_UTR$strand=="1"]="+"
  pos_UTR$strand[pos_UTR$strand=="-1"]="-"
  pos_UTR$chromosome_name=paste("chr",pos_UTR$chromosome_name,sep="")
  seq_data_utr<-getSeq(Hsapiens, pos_UTR$chromosome_name,
                       start=pos_UTR[,t1]-extend_bp,end=pos_UTR[,t2]+extend_bp,strand=pos_UTR$strand)
  seq_data_df_utr=as.data.frame(seq_data_utr)
  all_data_seq_utr=cbind(pos_UTR,seq_data_df_utr)
  select_seq=all_data_seq_utr[all_data_seq_utr$ensembl_transcript_id %in%      genes_infor$genes_main_transcript,]
  all_re<-c()
  for( i in 1:dim(select_seq)[1]){
    pv <- pqsfinder(DNAString(select_seq$x[i]))
    if(length(pv)>0){
      temp<-as.data.frame(pv)
      temp2<-elementMetadata(pv)
      all_re<-rbind(all_re,cbind(select_seq[i,1:5],temp,temp2))
    }
  }
  return(all_re)
  
}

getG4finder2<-function(pos_UTR,t1,t2,extend_bp){
  pos_UTR$strand[pos_UTR$strand=="1"]="+"
  pos_UTR$strand[pos_UTR$strand=="-1"]="-"
  pos_UTR$chromosome_name=paste("chr",pos_UTR$chromosome_name,sep="")
  seq_data_utr<-getSeq(Hsapiens, pos_UTR$chromosome_name,
                       start=pos_UTR[,t1]-extend_bp,end=pos_UTR[,t2]+extend_bp,strand=pos_UTR$strand)
  seq_data_df_utr=as.data.frame(seq_data_utr)
  all_data_seq_utr=cbind(pos_UTR,seq_data_df_utr)
  select_seq=all_data_seq_utr[all_data_seq_utr$ensembl_transcript_id %in%      genes_infor$genes_main_transcript,]
  all_re<-c()
  for( i in 1:dim(select_seq)[1]){
    pv <- pqsfinder(DNAString(select_seq$x[i]),min_score=0)
    if(length(pv)>0){
      temp<-as.data.frame(pv)
      temp2<-elementMetadata(pv)
      all_re<-rbind(all_re,cbind(select_seq[i,1:5],temp,temp2))
    }
  }
  return(all_re)
  
}

getG4finderRNA<-function(utr5_seq){
  all_re<-c()
  for( i in 1:length(utr5_seq$len)){
      pv <- pqsfinder(DNAString(utr5_seq$seq[i]),strand = "+")
      if(length(pv)>0){
        temp<-as.data.frame(pv)
        temp2<-elementMetadata(pv)
        all_re<-rbind(all_re,cbind(utr5_seq[i,1:2],temp,temp2))
      }
    }
  return(all_re)
  
}


get_CpG_promoter_CpG_OE<-function(gene_pos_data,length_up,length_down,species_x){
  genes_promoter=promoters(gene_pos_data,length_up,length_down)
  genes_promoter_df=as.data.frame(genes_promoter)
  print("delete these rows whose promoter start is less than 0")
  print(genes_promoter_df[genes_promoter_df$start<0,])
  genes_promoter_df<-genes_promoter_df[genes_promoter_df$start>0,]
  # delete genes with no promoter
  #seqnames start  end width strand           trans_ID            gene_id
  # 1               chrM  -999  100  1100      + ENSMUST00000082387 ENSMUSG00000064336
  length_x=length_up+length_down
  print("Here is the promoter length")
  print(length_x)
  seq_data<-getSeq(species_x, genes_promoter_df$seqnames,
                   start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  CG_pos=vmatchPattern("CG",seq_data)
  GC_pos=vmatchPattern("GC",seq_data)
  G_pos=vmatchPattern("G",seq_data)
  C_pos=vmatchPattern("C",seq_data)
  CG_count_all=sapply(CG_pos@ends, length)
  GC_count_all=sapply(GC_pos@ends, length)
  G_count_all=sapply(G_pos@ends, length)
  C_count_all=sapply(C_pos@ends, length)
  ### the length will be delete (C_count/length-G_count/length)/(C_count/length+G_count/length)
  CpG_OE_all=CG_count_all*length_x/(C_count_all*G_count_all)
  GpC_OE_all=GC_count_all*length_x/(C_count_all*G_count_all)
  GC_content_all=(C_count_all+G_count_all)/length_x
  GC_density_all=GC_count_all/length_x
  CG_density_all=CG_count_all/length_x
  GC_skew_all=(G_count_all-C_count_all)/(C_count_all+G_count_all)
  CG_GC_ratio_all=CG_count_all/GC_count_all
  all_data<-cbind(CpG_OE_all,GpC_OE_all,GC_content_all,GC_density_all,CG_density_all,GC_skew_all,CG_GC_ratio_all)
  rownames(all_data)<-genes_promoter_df$trans_ID
  data_temp<-data.frame(gene=genes_promoter_df$gene_id,trans=genes_promoter_df$trans_ID,all_data[genes_promoter_df$trans_ID,])
  #print(dim(data_temp))
  #print(head(data_temp))
  promoter_stat=aggregate(data_temp[,-c(1:2)],list(data_temp$gene),mean)
  #print(dim(promoter_stat))
  #print(head(promoter_stat))
  rownames(promoter_stat)<-promoter_stat$Group.1
  return(promoter_stat)
}
filter_low_data<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=1
    #print(quantile_cutoff)
    temp_data=temp_data[which(temp_data$Gene.Type=="TSG" & temp_data$diff_exp<0  & !(temp_data$norm_exp < quantile_cutoff & temp_data$tumor_exp < quantile_cutoff )),]  
    result=rbind(result,temp_data)
  }
  
  return(result)
}

get_slide_w_count<-function(x,window=200,step=50,length_a=length_x){
  a=rep(0,length_a)
  a[x]=1
  a_count=SlidingWindow(sum,a,window,step)
  return(a_count)
}


get_slide_w_count2<-function(x,window=50,step=10,length_a=length_x){
  a=rep(0,length_a)
  a[x]=1
  a_count=SlidingWindow(sum,a,window,step)
  return(a_count)
}

get_slide_w_count3<-function(x,window=50,step=10,length_a=length_x){
  a=rep(0,length_a)
  a[x]=1
  x2=x-c(0,x[1:(length(x)-1)])
  a[x[x2==1]]=0
  a_count=SlidingWindow(sum,a,window,step)
  return(a_count)
}

lengthMultiple<-function(x){
  x2=x-c(0,x[1:(length(x)-1)]);
  dup_len=length(x2[x2==1]);
  return(length(x)-dup_len)}





plot_density_TSS3<-function(x,a_site,types_x,f_select,type_genes="",get_gene_pos_genes,genes_promoter_df){
  print(length(a_site))
  print(dim(x))
  colnames(x)<-a_site
  genes_promoter_df=genes_promoter_df[f_select,]
  data_temp<-data.frame(gene=genes_promoter_df$gene_id,trans=genes_promoter_df$trans_ID,type=get_gene_pos_genes[genes_promoter_df$gene_id,]$age_type,Gene.Type=get_gene_pos_genes[genes_promoter_df$gene_id,]$Gene.Type, x[f_select,])
  data_temp=data_temp[!is.na(data_temp$type),]
  data_x=melt(data_temp,id=c("gene","type","trans","Gene.Type"))
  data_x$variable=as.numeric(gsub("X","",data_x$variable))
  length_x= (round(max(a_site)/1000)*1000)/2 #max(a_site)=3850.5 
  levCat.plot <- ggplot(data_x,aes(variable-length_x,value))+
    scale_color_brewer(palette = "Set1")+
    theme_classic()
  p=levCat.plot+
    stat_summary(fun.data = mean_cl_boot,geom = "ribbon",size = 1,aes(fill =  type),alpha = 0.3)+
    guides(fill = "none")+
    stat_summary(fun.y = mean,geom = "line",size = 1,aes(colour =  type))+
    labs(x = "distance from TSS (bp)",y = types_x,colour = "")+
    geom_vline(xintercept = 0,linetype = "dashed" )#+ylim(y1,y2)
  print(p+facet_wrap(~Gene.Type,nrow=1))
}
plot_density_TSS6<-function(x,a_site,types_x,f_select,type_genes="",get_gene_pos_genes,genes_promoter_df){
  print(length(a_site))
  print(dim(x))
  colnames(x)<-a_site
  genes_promoter_df=genes_promoter_df[f_select,]
  data_temp<-data.frame(gene=genes_promoter_df$gene_id,trans=genes_promoter_df$trans_ID,type=get_gene_pos_genes[genes_promoter_df$gene_id,]$age_type3,Gene.Type=get_gene_pos_genes[genes_promoter_df$gene_id,]$Gene.Type, x[f_select,])
  data_temp=data_temp[!is.na(data_temp$type),]
  data_x=melt(data_temp,id=c("gene","type","trans","Gene.Type"))
  data_x$variable=as.numeric(gsub("X","",data_x$variable))
  length_x= (round(max(a_site)/1000)*1000)/2 #max(a_site)=3850.5 
  levCat.plot <- ggplot(data_x,aes(variable-length_x,value))+
    scale_color_brewer(palette = "Set1")+
    theme_classic()
  p=levCat.plot+
    stat_summary(fun.data = mean_cl_boot,geom = "ribbon",size = 1,aes(fill =  Gene.Type),alpha = 0.3)+
    guides(fill = "none")+
    stat_summary(fun.y = mean,geom = "line",size = 1,aes(colour =  Gene.Type))+
    labs(x = "distance from TSS (bp)",y = types_x,colour = "")+
    geom_vline(xintercept = 0,linetype = "dashed" )#+ylim(y1,y2)
  print(p+facet_wrap(~type,nrow=1)+scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_fill_manual(values=c("#1874CD","#EE2C2C")))
}

plot_density_TSS4<-function(x,a_site,types_x,f_select,type_genes="",get_gene_pos_genes,genes_promoter_df){
  print(length(a_site))
  print(dim(x))
  colnames(x)<-a_site
  genes_promoter_df=genes_promoter_df[f_select,]
  data_temp<-data.frame(gene=genes_promoter_df$gene_id,trans=genes_promoter_df$trans_ID,type=get_gene_pos_genes[genes_promoter_df$gene_id,]$age_type,Gene.Type=get_gene_pos_genes[genes_promoter_df$gene_id,]$Gene.Type, x[f_select,])
  data_temp=data_temp[!is.na(data_temp$type),]
  data_x=melt(data_temp,id=c("gene","type","trans","Gene.Type"))
  data_x$variable=as.numeric(gsub("X","",data_x$variable))
  length_x=  (round(max(a_site)/1000)*1000)/2#max(a_site)=3850.5 
  levCat.plot <- ggplot(data_x,aes(variable-length_x,value))+
    scale_color_brewer(palette = "Set1")+
    theme_classic()
  p=levCat.plot+
    stat_summary(fun.data = mean_cl_boot,geom = "ribbon",size = 1,aes(fill =  Gene.Type),alpha = 0.3)+
    guides(fill = "none")+
    stat_summary(fun.y = mean,geom = "line",size = 1,aes(colour =  Gene.Type))+
    labs(x = "distance from TSS (bp)",y = types_x,colour = "")+
    geom_vline(xintercept = 0,linetype = "dashed" )#+ylim(y1,y2)
  return(p+scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_fill_manual(values=c("#1874CD","#EE2C2C")))
}

plot_density_TSS5<-function(x_all,a_site,types_x,f_select,type_genes="",get_gene_pos_genes,genes_promoter_df){
  print(length(a_site))
  all_result<-c()
  genes_promoter_df=genes_promoter_df[f_select,]
  cc=names(x_all)
  cc=cc[grepl("density",cc)]
  for(c1 in cc){
    
    x=x_all[[c1]]
    print(dim(x))
  colnames(x)<-a_site
  data_temp<-data.frame(gene=genes_promoter_df$gene_id,trans=genes_promoter_df$trans_ID,type=get_gene_pos_genes[genes_promoter_df$gene_id,]$age_type,Gene.Type=get_gene_pos_genes[genes_promoter_df$gene_id,]$Gene.Type, x[f_select,],type_actg=c1)
  data_temp=data_temp[!is.na(data_temp$type),]
  all_result<-rbind(all_result,data_temp)
  
  }
  data_x=melt(all_result,id=c("gene","type","trans","Gene.Type","type_actg"))
  
  data_x$variable=as.numeric(gsub("X","",data_x$variable))
  length_x= (round(max(a_site)/1000)*1000)/2 #max(a_site)=3850.5 
  levCat.plot <- ggplot(data_x,aes(variable-length_x,value))+
    scale_color_brewer(palette = "Set1")+
    theme_classic()
  p=levCat.plot+
    stat_summary(fun.data = mean_cl_boot,geom = "ribbon",size = 1,aes(fill =  type_actg),alpha = 0.3)+
    guides(fill = "none")+
    stat_summary(fun.y = mean,geom = "line",size = 1,aes(colour =  type_actg))+
    labs(x = "Distance from TSS (bp)",y = "Density",colour = "")+
    geom_vline(xintercept = 0,linetype = "dashed" )#+ylim(y1,y2)
  print(p+facet_wrap(~Gene.Type,nrow=1))
}



plot_density_TSS7<-function(x_all,a_site,types_x,f_select,type_genes="",get_gene_pos_genes,genes_promoter_df){
  print(length(a_site))
  all_result<-c()
  genes_promoter_df=genes_promoter_df[f_select,]
  cc=names(x_all)
  cc=cc[grepl("density",cc)]
  for(c1 in cc){
    
    x=x_all[[c1]]
    print(dim(x))
    colnames(x)<-a_site
    data_temp<-data.frame(gene=genes_promoter_df$gene_id,trans=genes_promoter_df$trans_ID,type=get_gene_pos_genes[genes_promoter_df$gene_id,]$age_type,Gene.Type=get_gene_pos_genes[genes_promoter_df$gene_id,]$Gene.Type, x[f_select,],type_actg=c1)
    data_temp=data_temp[!is.na(data_temp$type),]
    all_result<-rbind(all_result,data_temp)
    
  }
  data_x=melt(all_result,id=c("gene","type","trans","Gene.Type","type_actg"))
  
  data_x$variable=as.numeric(gsub("X","",data_x$variable))
  length_x= (round(max(a_site)/1000)*1000)/2 #max(a_site)=3850.5 
  levCat.plot <- ggplot(data_x,aes(variable-length_x,value))+
    scale_color_brewer(palette = "Set1")+
    theme_classic()
  p=levCat.plot+
    stat_summary(fun.data = mean_cl_boot,geom = "ribbon",size = 1,aes(fill =  type_actg),alpha = 0.3)+
    guides(fill = "none")+
    stat_summary(fun.y = mean,geom = "line",size = 1,aes(colour =  type_actg))+
    labs(x = "distance from TSS (bp)",y = "density",colour = "")+
    geom_vline(xintercept = 0,linetype = "dashed" )#+ylim(y1,y2)
  print(p)
}



get_CpG_promoter_window<-function(gene_pos_data,length_up,length_down,Hsapiens){
  genes_promoter=promoters(gene_pos_data,length_up,length_down)
  genes_promoter_df=as.data.frame(genes_promoter)
  length_w=200
  length_x<<-length_down+length_up
  seq_data<-getSeq(Hsapiens, genes_promoter_df$seqnames,
                   start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  
  
  # get_count<-function(type_x,seq_data){
  #   c1 <- makeCluster(detectCores()-3)
  #   on.exit(stopCluster(c1))
  #   clusterExport(c1,list("SlidingWindow","get_slide_w_count"))
  #   x_pos=vmatchPattern(type_x,seq_data)
  #   count_x<-parSapply(c1,x_pos@ends, get_slide_w_count)
  #   return(count_x)
  # }
  
  library(parallel)
  c1 <- makeCluster(detectCores()-3)
  on.exit(stopCluster(c1))
  clusterExport(c1,list("SlidingWindow","get_slide_w_count","length_x"))
  
  GC_pos=vmatchPattern("GC",seq_data)
  GC_count<-parSapply(c1,GC_pos@ends, FUN=get_slide_w_count)
  print("gc")
  
  CG_pos=vmatchPattern("CG",seq_data)
  CG_count<-parSapply(c1,CG_pos@ends, FUN=get_slide_w_count)    
  print(str(CG_count))
  
  print("cG")
  
  G_pos=vmatchPattern("G",seq_data)
  G_count<-parSapply(c1,G_pos@ends, FUN=get_slide_w_count)
  
  print("G")
  C_pos=vmatchPattern("C",seq_data)
  C_count<-parSapply(c1,C_pos@ends, FUN=get_slide_w_count)
  print("C")
  print(length_w)
  print(str(CG_count))
  CpG_OE=t(CG_count*length_w/(C_count*G_count))
  GpC_OE=t(GC_count*length_w/(C_count*G_count))
  GC_content=t((C_count+G_count)/length_w)
  GC_density=t(GC_count/length_w)
  CG_density=t(CG_count/length_w)
  GC_skew=t((G_count-C_count)/(C_count+G_count))
  CG_GC_ratio=t(CG_count/GC_count)
  result_windows=list(CpG_OE=CpG_OE,GpC_OE=GpC_OE,GC_content=GC_content,GC_density=GC_density,CG_density=CG_density,GC_skew=GC_skew,CG_GC_ratio=CG_GC_ratio,length_x=length_x,genes_promoter_df=genes_promoter_df)
  return(result_windows)
}



get_sequence_window<-function(gene_pos_data,length_up,length_down,Hsapiens){
  genes_promoter=promoters(gene_pos_data,length_up,length_down)
  genes_promoter_df=as.data.frame(genes_promoter)
  length_w=50
  length_x<<-length_down+length_up
  seq_data<-getSeq(Hsapiens, genes_promoter_df$seqnames,
                   start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  
  
  # get_count<-function(type_x,seq_data){
  #   c1 <- makeCluster(detectCores()-3)
  #   on.exit(stopCluster(c1))
  #   clusterExport(c1,list("SlidingWindow","get_slide_w_count"))
  #   x_pos=vmatchPattern(type_x,seq_data)
  #   count_x<-parSapply(c1,x_pos@ends, get_slide_w_count)
  #   return(count_x)
  # }
  
  library(parallel)
  c1 <- makeCluster(detectCores()-3)
  on.exit(stopCluster(c1))
  clusterExport(c1,list("SlidingWindow","get_slide_w_count2","length_x"))
  
  A_pos=vmatchPattern("A",seq_data)
  A_count<-parSapply(c1,A_pos@ends, FUN=get_slide_w_count2)
  print("gc")
  
  T_pos=vmatchPattern("T",seq_data)
  T_count<-parSapply(c1,T_pos@ends, FUN=get_slide_w_count2)    
  
  
  print("cG")
  
  G_pos=vmatchPattern("G",seq_data)
  G_count<-parSapply(c1,G_pos@ends, FUN=get_slide_w_count2)
  
  print("G")
  C_pos=vmatchPattern("C",seq_data)
  C_count<-parSapply(c1,C_pos@ends, FUN=get_slide_w_count2)
  print("C")
  print(length_w)
  
  C_density=t(C_count/length_w)
  G_density=t(G_count/length_w)
  A_density=t(A_count/length_w)
  T_density=t(T_count/length_w)
  result_windows=list(C_density=C_density,G_density=G_density,A_density=A_density,T_density=T_density,length_x=length_x,genes_promoter_df=genes_promoter_df)
  return(result_windows)
}


get_sequence_window_Mutiple<-function(gene_pos_data,length_up,length_down,Hsapiens){
  genes_promoter=promoters(gene_pos_data,length_up,length_down)
  genes_promoter_df=as.data.frame(genes_promoter)
  length_w=50
  length_x<<-length_down+length_up
  seq_data<-getSeq(Hsapiens, genes_promoter_df$seqnames,
                   start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  
  
  # get_count<-function(type_x,seq_data){
  #   c1 <- makeCluster(detectCores()-3)
  #   on.exit(stopCluster(c1))
  #   clusterExport(c1,list("SlidingWindow","get_slide_w_count"))
  #   x_pos=vmatchPattern(type_x,seq_data)
  #   count_x<-parSapply(c1,x_pos@ends, get_slide_w_count)
  #   return(count_x)
  # }
  
  library(parallel)
  c1 <- makeCluster(detectCores()-3)
  on.exit(stopCluster(c1))
  clusterExport(c1,list("SlidingWindow","get_slide_w_count3","length_x"))
  
  GG_pos=vmatchPattern("GG",seq_data)
  GG_count<-parSapply(c1,GG_pos@ends, FUN=get_slide_w_count3)
  print("gc")
  
  GGG_pos=vmatchPattern("GGG",seq_data)
  GGG_count<-parSapply(c1,GGG_pos@ends, FUN=get_slide_w_count3)    
  
  
  print("cG")
  
  CC_pos=vmatchPattern("CC",seq_data)
  CC_count<-parSapply(c1,CC_pos@ends, FUN=get_slide_w_count3)
  
  print("G")
  CCC_pos=vmatchPattern("CCC",seq_data)
  CCC_count<-parSapply(c1,CCC_pos@ends, FUN=get_slide_w_count3)
  print("C")
  print(length_w)
  
  CC_density=t(CC_count/length_w)
  GG_density=t(GG_count/length_w)
  GGG_density=t(GGG_count/length_w)
  CCC_density=t(CCC_count/length_w)
  result_windows=list(CC_density=CC_density,GG_density=GG_density,CCC_density=CCC_density,GGG_density=GGG_density,length_x=length_x,genes_promoter_df=genes_promoter_df)
  return(result_windows)
}



get_gene_region = function(x, upstream, downstream) {
  x=as.data.frame(x)
  t1=x[x$strand=="+","start"]
  x[x$strand=="+","start"]=t1+ upstream
  x[x$strand=="+","end"]=t1+ downstream-1
  t2=x[x$strand=="-","end"]
  x[x$strand=="-","start"]=t2- downstream+1
  x[x$strand=="-","end"]=t2- upstream
  
  
  x<-makeGRangesFromDataFrame(x,keep.extra.columns=T)
  
  return(x)
}

getFreq<-function(gene_pos_data,length_up,length_down,species_x,strand_x){
  length_x=abs(length_up)+abs(length_down)
  genes_promoter=get_gene_region(gene_pos_data,length_up,length_down)
  genes_promoter_df=as.data.frame(genes_promoter)
  print("delete these rows whose promoter start is less than 0")
  print(genes_promoter_df[genes_promoter_df$start<0,])
  genes_promoter_df<-genes_promoter_df[genes_promoter_df$start>0,]
  # delete genes with no promoter
  #seqnames start  end width strand           trans_ID            gene_id
  # 1               chrM  -999  100  1100      + ENSMUST00000082387 ENSMUSG00000064336
  seq_data_utr<-getSeq(species_x, genes_promoter_df$seqnames,
                       start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  
  G_pos=vmatchPattern("G",seq_data_utr)
  GG_pos=vmatchPattern("GG",seq_data_utr)
  GGG_pos=vmatchPattern("GGG",seq_data_utr)
  C_pos=vmatchPattern("C",seq_data_utr)
  CC_pos=vmatchPattern("CC",seq_data_utr)
  CCC_pos=vmatchPattern("CCC",seq_data_utr)
  A_pos=vmatchPattern("A",seq_data_utr)
  T_pos=vmatchPattern("T",seq_data_utr)
  G_count_all=as.numeric(sapply(G_pos@ends, length))
  C_count_all=as.numeric(sapply(C_pos@ends, length))
  GG_count_all=as.numeric(sapply(GG_pos@ends, lengthMultiple))
  CC_count_all=as.numeric(sapply(CC_pos@ends, lengthMultiple))
  GGG_count_all=as.numeric(sapply(GGG_pos@ends, lengthMultiple))
  CCC_count_all=as.numeric(sapply(CCC_pos@ends, lengthMultiple))
  genes_promoter_df$G_density<-G_count_all/length_x
  genes_promoter_df$C_density<-C_count_all/length_x
  genes_promoter_df$GG_density<-GG_count_all/length_x
  genes_promoter_df$CC_density<-CC_count_all/length_x
  genes_promoter_df$GGG_density<-GGG_count_all/length_x
  genes_promoter_df$CCC_density<-CCC_count_all/length_x
  genes_promoter_df$CCC_GGG<-(CCC_count_all+GGG_count_all)/length_x
  genes_promoter_df$CC_GG<-(CC_count_all+GG_count_all)/length_x
  genes_promoter_df$GCbias<-(G_count_all-C_count_all)/length_x
  genes_promoter_df$GC_content=(C_count_all+G_count_all)/length_x
  T_count_all=as.numeric(sapply(T_pos@ends, length))
  A_count_all=as.numeric(sapply(A_pos@ends, length))
  genes_promoter_df$T_density<-T_count_all/length_x
  genes_promoter_df$A_density<-A_count_all/length_x
  genes_promoter_df$GTbias<-(G_count_all-T_count_all)/length_x
  pos_all<-aggregate(genes_promoter_df[,c("G_density","C_density","A_density","T_density","GGG_density","GG_density","CCC_density","CC_density","CC_GG","GC_content","CCC_GGG","GCbias","GTbias")],list(genes_promoter_df$gene_id),sum)
  rownames(pos_all)<-pos_all$Group.1
  
  
  
  
  return(pos_all)
  
}

getFreq2<-function(gene_pos_data,species_x,strand_x){
  genes_promoter_df=as.data.frame(gene_pos_data)
  length_x=genes_promoter_df$end-genes_promoter_df$start+1
  print("delete these rows whose promoter start is less than 0")
  print(genes_promoter_df[genes_promoter_df$start<0,])
  genes_promoter_df<-genes_promoter_df[genes_promoter_df$start>0,]
  # delete genes with no promoter
  #seqnames start  end width strand           trans_ID            gene_id
  # 1               chrM  -999  100  1100      + ENSMUST00000082387 ENSMUSG00000064336
  seq_data_utr<-getSeq(species_x, genes_promoter_df$seqnames,
                       start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  
  G_pos=vmatchPattern("G",seq_data_utr)
  GG_pos=vmatchPattern("GG",seq_data_utr)
  GGG_pos=vmatchPattern("GGG",seq_data_utr)
  C_pos=vmatchPattern("C",seq_data_utr)
  CC_pos=vmatchPattern("CC",seq_data_utr)
  CCC_pos=vmatchPattern("CCC",seq_data_utr)
  A_pos=vmatchPattern("A",seq_data_utr)
  T_pos=vmatchPattern("T",seq_data_utr)
  G_count_all=as.numeric(sapply(G_pos@ends, length))
  C_count_all=as.numeric(sapply(C_pos@ends, length))
  GG_count_all=as.numeric(sapply(GG_pos@ends, lengthMultiple))
  CC_count_all=as.numeric(sapply(CC_pos@ends, lengthMultiple))
  GGG_count_all=as.numeric(sapply(GGG_pos@ends, lengthMultiple))
  CCC_count_all=as.numeric(sapply(CCC_pos@ends, lengthMultiple))
  genes_promoter_df$G_density<-G_count_all/length_x
  genes_promoter_df$C_density<-C_count_all/length_x
  genes_promoter_df$GG_density<-GG_count_all/length_x
  genes_promoter_df$CC_density<-CC_count_all/length_x
  genes_promoter_df$GGG_density<-GGG_count_all/length_x
  genes_promoter_df$CCC_density<-CCC_count_all/length_x
  genes_promoter_df$CCC_count<-CCC_count_all
  T_count_all=as.numeric(sapply(T_pos@ends, length))
  A_count_all=as.numeric(sapply(A_pos@ends, length))
  genes_promoter_df$T_density<-T_count_all/length_x
  genes_promoter_df$A_density<-A_count_all/length_x
  pos_all<-aggregate(genes_promoter_df[,c("G_density","C_density","A_density","T_density","GGG_density","GG_density","CCC_density","CC_density","CCC_count")],list(genes_promoter_df$gene_id),sum)
  rownames(pos_all)<-pos_all$Group.1
  return(pos_all)
}

getFreqRaw<-function(utr5_seq){
  genes_promoter_df=utr5_seq
  length_x=utr5_seq$len
  seq_data_utr<-DNAStringSet(utr5_seq$seq)
  G_pos=vmatchPattern("G",seq_data_utr)
  GG_pos=vmatchPattern("GG",seq_data_utr)
  GGG_pos=vmatchPattern("GGG",seq_data_utr)
  C_pos=vmatchPattern("C",seq_data_utr)
  CC_pos=vmatchPattern("CC",seq_data_utr)
  CCC_pos=vmatchPattern("CCC",seq_data_utr)
  A_pos=vmatchPattern("A",seq_data_utr)
  T_pos=vmatchPattern("T",seq_data_utr)
  G_count_all=as.numeric(sapply(G_pos@ends, length))
  C_count_all=as.numeric(sapply(C_pos@ends, length))
  GG_count_all=as.numeric(sapply(GG_pos@ends, lengthMultiple))
  CC_count_all=as.numeric(sapply(CC_pos@ends, lengthMultiple))
  GGG_count_all=as.numeric(sapply(GGG_pos@ends, lengthMultiple))
  CCC_count_all=as.numeric(sapply(CCC_pos@ends, lengthMultiple))
  genes_promoter_df$G_density<-G_count_all/length_x
  genes_promoter_df$C_density<-C_count_all/length_x
  genes_promoter_df$GG_density<-GG_count_all/length_x
  genes_promoter_df$CC_density<-CC_count_all/length_x
  genes_promoter_df$GGG_density<-GGG_count_all/length_x
  genes_promoter_df$CCC_density<-CCC_count_all/length_x
  genes_promoter_df$CCC_count<-CCC_count_all
  T_count_all=as.numeric(sapply(T_pos@ends, length))
  A_count_all=as.numeric(sapply(A_pos@ends, length))
  genes_promoter_df$T_density<-T_count_all/length_x
  genes_promoter_df$A_density<-A_count_all/length_x

  return(genes_promoter_df)
}


getHist<-function(col_x,genes_infor,title_x){
  
  x1=table(genes_infor[genes_infor$cancer_types3!=0 & genes_infor$Gene.Type=="TSG",col_x],genes_infor[genes_infor$cancer_types3!=0& genes_infor$Gene.Type=="TSG",]$cancer_types3)
  
  temp=fisher.test(x1,alternative = "greater")
  sig=pvalueTosig(temp$p.value)
  x1[2,]/apply(x1,2,sum)
  x2=x1[2,]/apply(x1,2,sum)
  cancer_ratio<-data.frame(ratio=x2,cancer_type=factor(c("1",">1"),levels =c("1",">1") ))
  y_pos=max(cancer_ratio$ratio)+0.05
  p_cancer1<-ggplot(cancer_ratio,aes(x=cancer_type,y=ratio,fill=cancer_type))+geom_bar(stat = "identity",width = 0.6)+ylab(title_x)+theme_hd()+theme(legend.position = "none")+annotate("text",x=1.5,y=y_pos+0.05,label=sig,size=8)+scale_fill_manual(values=c("#ADD8E6", "#1874CD"))+annotate("segment",x = 1, xend = 2, y = y_pos, yend = y_pos)+scale_fill_manual(values=c("#ADD8E6", "#1874CD"))
  
  
  overlap.data=rbind(cbind(cancer="Carcinoma",genes_infor[which(genes_infor$Carcinoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Leukemia",genes_infor[which(genes_infor$Leukemia!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Lymphoma",genes_infor[which(genes_infor$Lymphoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Sarcoma",genes_infor[which(genes_infor$Sarcoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Myeloma",genes_infor[which(genes_infor$Myeloma!="" & genes_infor$Gene.Type=="TSG"),col_x]))
  overlap.data<-as.data.frame(overlap.data)
  overlap.data$cancer<-factor(overlap.data$cancer,levels=c("Myeloma","Carcinoma","Sarcoma","Leukemia","Lymphoma"))
  colnames(overlap.data)<-c("cancer_types",col_x)
  
  x1=table(overlap.data[,col_x],overlap.data$cancer_types)
  temp=fisher.test(x1)
  pvalue<-c()
  for(tt in getComp(colnames(x1))){
    print(tt)
    temp=fisher.test(x1[,tt],alternative = "greater")
    temp2=fisher.test(x1[,tt],alternative = "less")
    temp1=fisher.test(x1[,tt],alternative = "two.side")
    print(fisher.test(x1[,tt]))
    pvalue<-rbind(pvalue,cbind(temp1$p.value,temp$p.value,temp2$p.value,paste(tt,collapse = "")))
  }
  
  print(pvalue)
  sig=pvalueTosig(temp$p.value)
  x1[2,]/apply(x1,2,sum)
  x2=x1[2,]/apply(x1,2,sum)
  cancer_ratio_overlap<-data.frame(ratio=x2,cancer_type=colnames(x1))
  
  
  
  return(list(p1=p_cancer1,data=cancer_ratio_overlap,pvalue=pvalue))
  
}


getSigPlot<-function(col_x,genes_infor,title_x){
  
  x1=table(genes_infor[,col_x],genes_infor$Gene.Type)
  
  temp=fisher.test(x1)
  sig=pvalueTosig(temp$p.value)
  x1[2,]/apply(x1,2,sum)
  x2=x1[2,]/apply(x1,2,sum)
  cancer_ratio<-data.frame(ratio=x2,type=names(x2))
  y_pos=max(cancer_ratio$ratio)+0.05
  p_cancer1<-ggplot(cancer_ratio,aes(x=type,y=ratio,fill=type))+geom_bar(stat = "identity",width = 0.6)+ylab("Ratio of genes")+theme_hd()+theme(legend.position = "none")+annotate("text",x=1.5,y=y_pos+0.05,label=sig,size=8)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+annotate("segment",x = 1, xend = 2, y = y_pos, yend = y_pos)+xlab(title_x)
  return(p_cancer1)
}


plotGoTerm<-function(kegg_gene_infor,col_x,x,title_x,side){
  term_ratio=as.data.frame(kegg_gene_infor[kegg_gene_infor$Gene.Type!="Oncogene" ,] %>% dplyr::count(get(col_x),term,Gene.Type) %>% dplyr::group_by(term,Gene.Type) %>%dplyr::mutate(freq = n / sum(n)))
  
  p_term<-ggplot(data=term_ratio[term_ratio[,1]==x,],aes(x=term,y=freq,fill=Gene.Type))+geom_bar(stat = "identity", position=position_dodge())+scale_fill_manual("Gene.Type",values=c("#1874CD","#EE2C2C"))+xlab("Gene Age Groups(million years)")+ylab(title_x)+theme_hd()+theme(axis.text.x = element_text(size=I(8),angle = 90,hjust=0.95,vjust=0.2))
  
  temp_x=term_ratio[term_ratio[,1]==x & term_ratio$Gene.Type=="Non-Cancer",]
  levels<-temp_x[order(temp_x$freq),"term"]
  kegg_gene_infor$term<-factor(kegg_gene_infor$term,levels = levels)
  pMotif_term<-getSigRatio(kegg_gene_infor[kegg_gene_infor$Gene.Type!="Oncogene" ,],col_x,"term",x,title_x,"GO TERMS",side)+theme(axis.text.x = element_text(angle = 55,hjust = 1,vjust = 1))
  
  
  pMotif_term_paralog<-getSigRatio(kegg_gene_infor[kegg_gene_infor$Gene.Type!="Oncogene"  & kegg_gene_infor$Gene_paralogs=="paralog",],col_x,"term",x,title_x,"GO TERMS",side)+ scale_fill_manual(name = "Gene.Type",labels = c("Non-Cancer-Paralogs", "TSG"),values=c("#1874CD","#EE2C2C"))+theme(axis.text.x = element_text(angle = 55,hjust = 1,vjust = 1))
  
  return(list(p1=pMotif_term,p2=pMotif_term_paralog))
}


plotStatPromoter2<-function(genes_infor,all_data_tsg,TSG_paralog_non_cancer,promoter_stat_up,promoter_stat_down,type_x){
  genes_infor[,paste(type_x,"_pos",sep="")]=promoter_stat_up[rownames(genes_infor),type_x]
  genes_infor[,paste(type_x,"_neg",sep="")]=promoter_stat_down[rownames(genes_infor),type_x]
  genes_infor[,paste(type_x,"_pos",sep="")][is.na(genes_infor[,paste(type_x,"_pos",sep="")])]=0
  genes_infor[,paste(type_x,"_neg",sep="")][is.na(genes_infor[,paste(type_x,"_neg",sep="")])]=0
  genes_infor[,paste(type_x,"_negPosratio",sep="")]=log2(genes_infor[,paste(type_x,"_neg",sep="")]+1)-log2(genes_infor[,paste(type_x,"_pos",sep="")]+1)
  min_x=min(c(genes_infor[,paste(type_x,"_pos",sep="")],genes_infor[,paste(type_x,"_neg",sep="")]),na.rm = T)
  max_x=max(c(genes_infor[,paste(type_x,"_pos",sep="")],genes_infor[,paste(type_x,"_neg",sep="")]),na.rm = T)*1.2
  TSG_paralog_non_cancer[,paste(type_x,"_pos",sep="")]=promoter_stat_up[rownames(TSG_paralog_non_cancer),type_x]
  TSG_paralog_non_cancer[,paste(type_x,"_neg",sep="")]=promoter_stat_down[rownames(TSG_paralog_non_cancer),type_x]
  TSG_paralog_non_cancer[,paste(type_x,"_negPosratio",sep="")]=genes_infor[rownames(TSG_paralog_non_cancer),paste(type_x,"_negPosratio",sep="")]
  
  all_data_tsg[,paste(type_x,"_pos",sep="")]<-genes_infor[all_data_tsg$ensembl_gene_id,paste(type_x,"_pos",sep="")]
  all_data_tsg[,paste(type_x,"_neg",sep="")]<-genes_infor[all_data_tsg$ensembl_gene_id,paste(type_x,"_neg",sep="")]
  all_data_tsg[,paste(type_x,"_negPosratio",sep="")]<-genes_infor[all_data_tsg$ensembl_gene_id,paste(type_x,"_negPosratio",sep="")]
  
  
  p1<-plot_non_paired_paralog(TSG_paralog_non_cancer,paste(type_x,"_negPosratio",sep=""),paste(type_x,"_negPosratio",sep=""))
  p0<-plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],paste(type_x,"_negPosratio",sep=""),paste(type_x,"_negPosratio",sep=""))
  
  
  p2<-plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],paste(type_x,"_pos",sep=""),paste(type_x,"_pos",sep=""))+ylim(min_x,max_x)
  p3<-plot_non_paired_paralog(TSG_paralog_non_cancer,paste(type_x,"_pos",sep=""),paste(type_x,"_pos",sep=""))+ylim(min_x,max_x)
  
  p4<-plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],paste(type_x,"_neg",sep=""),paste(type_x,"_neg",sep=""))+ylim(min_x,max_x)
  p5<-plot_non_paired_paralog(TSG_paralog_non_cancer,paste(type_x,"_neg",sep=""),paste(type_x,"_neg",sep=""))+ylim(min_x,max_x)
  
  
  print(ggarrange(p0,p2,p4,p1,p3,p5,nrow=2,ncol=3))
  
  
  p1=plotSubSig(genes_infor[!is.na(genes_infor$age_type3),],"age_type3", paste(type_x,"_negPosratio",sep=""),"Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age class\n(million years)",paste(type_x,"_negPosratio",sep="")) + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")
  
  
  p2=plotSubSig(genes_infor[!is.na(genes_infor$age_type3),],"age_type3", paste(type_x,"_pos",sep=""),"Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age class\n(million years)",paste(type_x,"_pos",sep="")) + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+ylim(min_x,max_x)
  
  p3=plotSubSig(genes_infor[!is.na(genes_infor$age_type3),],"age_type3", paste(type_x,"_neg",sep=""),"Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age class\n(million years)",paste(type_x,"_neg",sep="")) + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+ylim(min_x,max_x)
  print(ggarrange(p1,p2,p3,nrow=1))
  
  
  
  all_data_tsg_filtered=filter_low_data(all_data_tsg )
  all_data_tsg_filtered$diff_exp=abs(all_data_tsg_filtered$diff_exp)
  
  
  
  p.exp<-plotForest(all_data_tsg,"norm_exp",paste(type_x,"_pos",sep=""),paste("Pearson's r(",type_x,"_pos  VS gene expression in TCGA",sep=""),0.9,0.9,"1cm")
  print(ggarrange(p.exp))
  p.ff<-plotForest(all_data_tsg_filtered,"diff_exp",paste(type_x,"_pos",sep=""),paste("Pearson's r(",type_x,"_pos  VS abs(log2FC))\n Down-regulated TSG in TCGA",sep=""),0.9,0.9,"1cm")
  
  p.exp_downupratio<-plotForest(all_data_tsg,"norm_exp",paste(type_x,"_negPosratio",sep=""),paste("Pearson's r(",type_x,"_negPosratio  VS gene expression in TCGA",sep=""),0.9,0.9,"1cm")
  p.ff_downupratio<-plotForest(all_data_tsg_filtered,"diff_exp",paste(type_x,"_negPosratio",sep=""),paste("Pearson's r(",type_x,"_negPosratio  VS abs(log2FC))\n Down-regulated TSG in TCGA",sep=""),0.9,0.9,"1cm")
  
  
  p.exp_down<-plotForest(all_data_tsg,"norm_exp",paste(type_x,"_neg",sep=""),paste("Pearson's r(",type_x,"_neg  VS gene expression in TCGA",sep=""),0.9,0.9,"1cm")
  p.ff_down<-plotForest(all_data_tsg_filtered,"diff_exp",paste(type_x,"_neg",sep=""),paste("Pearson's r(",type_x,"_neg  VS abs(log2FC))\n Down-regulated TSG in TCGA",sep=""),0.9,0.9,"1cm")
  print(ggarrange(p.ff_down,p.exp_down),nrow=2,ncol=1)
  print(ggarrange(p.ff_downupratio,p.exp_downupratio),nrow=2,ncol=1)
  print(ggarrange(p.ff,p.exp),nrow=2,ncol=1)
  
  
  p=ggpaired(genes_infor[!is.na(genes_infor[,paste(type_x,"_neg",sep="")]*genes_infor[paste(type_x,"_pos",sep="")]),], cond1 = paste(type_x,"_pos",sep=""), cond2 = paste(type_x,"_neg",sep=""),fill = "condition", palette = "jco",line.color =  "gray90")+stat_compare_means(method = "t.test", paired = TRUE)
  p$layers <- p$layers[-2]
  print(p+facet_wrap(~Gene.Type))
  
  all_re<-list(genes_infor=genes_infor,all_data_tsg=all_data_tsg)
  return(all_re)
}


plotStatPromoter3<-function(genes_infor,all_data_tsg,TSG_paralog_non_cancer,promoter_stat_up,type_x){
  genes_infor[,type_x]=promoter_stat_up[rownames(genes_infor),type_x]
 genes_infor[,type_x][is.na(genes_infor[,type_x])]=0
  min_x=min(c(genes_infor[,type_x],genes_infor[,type_x]),na.rm = T)
  max_x=max(c(genes_infor[,type_x],genes_infor[,type_x]),na.rm = T)*1.2
  TSG_paralog_non_cancer[,type_x]=promoter_stat_up[rownames(TSG_paralog_non_cancer),type_x]
 
  all_data_tsg[,type_x]<-genes_infor[all_data_tsg$ensembl_gene_id,type_x]
 
  
  p1<-plot_non_paired_paralog(TSG_paralog_non_cancer,type_x,type_x)
  p2<-plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],type_x,type_x)

  p3=plotSubSig(genes_infor[!is.na(genes_infor$age_type3),],"age_type3",type_x,"Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age class\n(million years)",type_x) + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+ylim(min_x,max_x)
 
  
  
  
  all_data_tsg_filtered=filter_low_data(all_data_tsg )
  all_data_tsg_filtered$diff_exp=abs(all_data_tsg_filtered$diff_exp)
  
  
  
  p.exp<-plotForest(all_data_tsg,"norm_exp",type_x,paste("Pearson's r(",type_x,"  VS gene expression in TCGA",sep=""),0.9,0.9,"1cm")
  print(ggarrange(p.exp))
  p.ff<-plotForest(all_data_tsg_filtered,"diff_exp",type_x,paste("Pearson's r(",type_x," VS abs(log2FC))\n Down-regulated TSG in TCGA",sep=""),0.9,0.9,"1cm")
  pall1=ggarrange(p1,p2,p3,nrow=1,ncol=3)
  pall2=ggarrange(p.ff,p.exp,nrow=1,ncol=2)
  print(ggarrange(pall1,pall2,nrow=2,ncol=1))
  return(list(p1=p1,p2=p2,p3=p3,p.ff=p.ff,p.exp=p.exp))

}
#delete NA
plotStatPromoterNoNA<-function(genes_infor,all_data_tsg,TSG_paralog_non_cancer,promoter_stat_up,type_x){
  genes_infor[,type_x]=promoter_stat_up[rownames(genes_infor),type_x]
  genes_infor=genes_infor[!is.na(genes_infor[,type_x]),]
  min_x=min(c(genes_infor[,type_x],genes_infor[,type_x]),na.rm = T)
  max_x=max(c(genes_infor[,type_x],genes_infor[,type_x]),na.rm = T)*1.2
  TSG_paralog_non_cancer[,type_x]=promoter_stat_up[rownames(TSG_paralog_non_cancer),type_x]
  
  all_data_tsg[,type_x]<-genes_infor[all_data_tsg$ensembl_gene_id,type_x]
  all_data_tsg=all_data_tsg[!is.na(all_data_tsg[,type_x]),]
  
  p1<-plot_non_paired_paralog(TSG_paralog_non_cancer,type_x,type_x)
  p2<-plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],type_x,type_x)
  
  p3=plotSubSig(genes_infor[!is.na(genes_infor$age_type3),],"age_type3",type_x,"Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age class\n(million years)",type_x) + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+ylim(min_x,max_x)
  
  
  
  
  all_data_tsg_filtered=filter_low_data(all_data_tsg )
  all_data_tsg_filtered$diff_exp=abs(all_data_tsg_filtered$diff_exp)
  
  
  
  p.exp<-plotForest(all_data_tsg,"norm_exp",type_x,paste("Pearson's r(",type_x,"  VS gene expression in TCGA",sep=""),0.9,0.9,"1cm")
  print(ggarrange(p.exp))
  p.ff<-plotForest(all_data_tsg_filtered,"diff_exp",type_x,paste("Pearson's r(",type_x," VS abs(log2FC))\n Down-regulated TSG in TCGA",sep=""),0.9,0.9,"1cm")
  pall1=ggarrange(p1,p2,p3,nrow=1,ncol=3)
  pall2=ggarrange(p.ff,p.exp,nrow=1,ncol=2)
  print(ggarrange(pall1,pall2,nrow=2,ncol=1))
  return(list(p1=p1,p2=p2,p3=p3,p.ff=p.ff,p.exp=p.exp))
  
}

GetExpressionResultCor<-function(all_data_tsg,col_x){
  all_data_tsg_filtered=filter_low_data(all_data_tsg )
  all_data_tsg_filtered$diff_exp=abs(all_data_tsg_filtered$diff_exp)

  p.exp<-plotForest(all_data_tsg,"norm_exp",col_x,paste("Pearson's r(",col_x,"  VS gene expression in TCGA",sep=""),0.9,0.9,"1cm")
  print(ggarrange(p.exp))
  p.ff<-plotForest(all_data_tsg_filtered,"diff_exp",col_x,paste("Pearson's r(",col_x,"  VS abs(log2FC))\n Down-regulated TSG in TCGA",sep=""),0.9,0.9,"1cm")
  

return(ggarrange(p.ff,p.exp,nrow=1,ncol=2))

}



de_factor<-function(data_x){
  data_x<-as.data.frame(data_x)
  f_column<-is.na(as.numeric(as.character(data_x[1,])))# check the first row to see if this coloum is numeric or character
  if(length(f_column[!f_column])>1){
    data_1<-apply(data_x[,!f_column],2,function(x) as.numeric(as.character(x)))
  }else{data_1<- as.numeric(as.character(data_x[,!f_column]))}
  
  if(length(f_column[f_column])>1){
    data_2<-apply(data_x[,f_column],2,function(x) as.character(x))
  }else{data_2<- as.character(data_x[,f_column])}
  
  data_x[,f_column]<-data_2
  data_x[,!f_column]<-data_1
  return(data_x)
}


getG4finderNew<-function(gene_pos_data,length_up,length_down,species_x,strand_x){
  genes_promoter=promoters(gene_pos_data,length_up,length_down)
  genes_promoter_df=as.data.frame(genes_promoter)
  print("delete these rows whose promoter start is less than 0")
  print(genes_promoter_df[genes_promoter_df$start<0,])
  genes_promoter_df<-genes_promoter_df[genes_promoter_df$start>0,]
  # delete genes with no promoter
  #seqnames start  end width strand           trans_ID            gene_id
  # 1               chrM  -999  100  1100      + ENSMUST00000082387 ENSMUSG00000064336
  seq_data_utr<-getSeq(species_x, genes_promoter_df$seqnames,
                       start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  seq_data_df_utr=as.data.frame(seq_data_utr)
  all_data_seq_utr=cbind(genes_promoter_df,seq_data_df_utr)
  select_seq=all_data_seq_utr
  all_re<-c()
  for( i in 1:dim(select_seq)[1]){
    pv <- pqsfinder(DNAString(select_seq$x[i]),strand = strand_x)
    if(length(pv)>0){
      temp<-as.data.frame(pv)
      temp2<-elementMetadata(pv)
      all_re<-rbind(all_re,cbind(select_seq[i,1:7],temp,temp2))
    }
  }
  return(all_re)
  
}



getG4finderNew2<-function(gene_pos_data,length_up,length_down,species_x,strand_x){
  genes_promoter=promoters(gene_pos_data,length_up,length_down)
  genes_promoter_df=as.data.frame(genes_promoter)
  print("delete these rows whose promoter start is less than 0")
  print(genes_promoter_df[genes_promoter_df$start<0,])
  genes_promoter_df<-genes_promoter_df[genes_promoter_df$start>0,]
  # delete genes with no promoter
  #seqnames start  end width strand           trans_ID            gene_id
  # 1               chrM  -999  100  1100      + ENSMUST00000082387 ENSMUSG00000064336
  seq_data_utr<-getSeq(species_x, genes_promoter_df$seqnames,
                       start=genes_promoter_df$start, end=genes_promoter_df$end,strand=genes_promoter_df$strand)
  seq_data_df_utr=as.data.frame(seq_data_utr)
  all_data_seq_utr=cbind(genes_promoter_df,seq_data_df_utr)
  select_seq=all_data_seq_utr
  all_re<-c()
  for( i in 1:dim(select_seq)[1]){
    pv <- pqsfinder(DNAString(select_seq$x[i]),strand = strand_x,min_score=1)
    if(length(pv)>0){
      temp<-as.data.frame(pv)
      temp2<-elementMetadata(pv)
      all_re<-rbind(all_re,cbind(select_seq[i,1:7],temp,temp2))
    }
  }
  return(all_re)
  
}


getGenePosData<-function(temp_data,column_x,addCHR=T){
  gene_pos<-temp_data[,column_x]
  colnames(gene_pos)<-c("trans_ID","seqnames","start","end","strand","gene_id")
  gene_pos$strand[gene_pos$strand=="1"]="+"
  gene_pos$strand[gene_pos$strand=="-1"]="-"
  if(addCHR){
    gene_pos$seqnames=paste("chr",gene_pos[,2],sep="")
  }
  
  gene_pos<-gene_pos[!(is.na(gene_pos$start) | is.na(gene_pos$end)),]
  gene_pos_data<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=T)
  return(gene_pos_data)
}



getIntron<-function(exon_data,genes_infor,n){
  temp_exon<-exon_data[exon_data$ensembl_transcript_id %in%  genes_infor$ensembl_transcript_id & exon_data$rank==n, ]
  temp_exon$chromosome_name=genes_infor[temp_exon$ensembl_gene_id,]$chromosome_name
  temp_exon$strand=genes_infor[temp_exon$ensembl_gene_id,]$strand
  rownames(temp_exon)<-temp_exon$ensembl_gene_id
  print(paste(" exon ",n))
  
  temp_exon2<-exon_data[exon_data$ensembl_transcript_id %in%  genes_infor$ensembl_transcript_id & exon_data$rank==n+1, ]
  rownames(temp_exon2)<-temp_exon2$ensembl_gene_id
  print(paste(" exon ",n+1))
  
  temp_exon$exon2_start<-temp_exon2[temp_exon$ensembl_gene_id,"exon_chrom_start"]
  temp_exon$exon2_end<-temp_exon2[temp_exon$ensembl_gene_id,"exon_chrom_end"]
  temp_exon$intron_start<-temp_exon$exon_chrom_end
  temp_exon$intron_end<-temp_exon$exon2_start
  temp_exon$intron_start[temp_exon$strand=="-1"]=temp_exon$exon2_end[temp_exon$strand=="-1" ]
  temp_exon$intron_end[temp_exon$strand=="-1" ]=temp_exon$exon_chrom_start[temp_exon$strand=="-1" ]
  temp_exon$intron_start=temp_exon$intron_start+1
  temp_exon$intron_end=temp_exon$intron_end-1
  print(head(temp_exon[,c(2,4:5,12:13)]))
  print(head(temp_exon2[,c(2,4:6)]))
  return(temp_exon)
}


PlotDensity<-function(result_window_1000_1000_multiple,genes_infor,title_x){
  a_site=SlidingWindow(median,1:result_window_1000_1000_multiple$length_x,window=50,step=10)
  pdf(paste("./result/",title_x,"_mutiple_density_plot.pdf",sep=""),width = 15,height = 10)
  sm_1=plot_density_TSS6(result_window_1000_1000_multiple[["GG_density"]],a_site,"Non-template GG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  
  
  sm_2=plot_density_TSS6(result_window_1000_1000_multiple[["GGG_density"]],a_site,"Non-template GGG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  
  sm_3=plot_density_TSS6(result_window_1000_1000_multiple[["CC_density"]],a_site,"Template GG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  
  sm_4=plot_density_TSS6(result_window_1000_1000_multiple[["CCC_density"]],a_site,"Template GGG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  
  bm_1=plot_density_TSS4(result_window_1000_1000_multiple[["GG_density"]],a_site,"Non-template GG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  
  
  bm_2=plot_density_TSS4(result_window_1000_1000_multiple[["GGG_density"]],a_site,"Non-template GGG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  
  bm_3=plot_density_TSS4(result_window_1000_1000_multiple[["CC_density"]],a_site,"Template GG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  
  bm_4=plot_density_TSS4(result_window_1000_1000_multiple[["CCC_density"]],a_site,"Template GGG_density",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)
  all_list=list(gg1=sm_1,ggg1=sm_2,cc1=sm_3,cc1=sm_4,gg2=bm_1,ggg2=bm_2,cc2=bm_3,cc2=bm_4)
  dev.off()
  return(all_list)
}

PlotDensityIntron<-function(result_window_1000_1000_multiple,genes_infor,title_x,num){
  a_site=SlidingWindow(median,1:result_window_1000_1000_multiple$length_x,window=50,step=10)
 xlab_x=paste("Distance from exon",num,"-intron",num," boundary(bp)",sep="")
  bm_1=plot_density_TSS4(result_window_1000_1000_multiple[["GG_density"]],a_site,"Frequency of G-runs in non-template strand",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)+xlab(xlab_x)+theme_hd()
  bm_3=plot_density_TSS4(result_window_1000_1000_multiple[["CC_density"]],a_site,"Frequency of G-runs in template strand",result_window_1000_1000_multiple$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type!="Oncogene"],"",genes_infor,result_window_1000_1000_multiple$genes_promoter_df)+xlab(xlab_x)+theme_hd()
  tiff(paste("./result/",title_x,"_mutiple_density_plot.tiff",sep=""),width = 15,height = 7,res=300,units="in",compression = "lzw")
   print(ggarrange(bm_1,bm_3,nrow=1,ncol=2,labels = c(letters[1:2])  ,font.label = list(size = 18, color = "black", face = "bold")))
  dev.off()

}



GetAllresult<-function(genes_infor,kegg_data_select,col_x,main_x){
  
  prG41<-getSigRatio(genes_infor,col_x,"age_type3","+",paste("Ratio of genes with",main_x),"Gene Age Groups(million years)")
  prG42<-getSigRatio(genes_infor[genes_infor$Gene_paralogs=="paralog",],col_x,"age_type3","+",paste("Ratio of genes with",main_x),"Gene Age Groups(million years)")+ scale_fill_manual(name = "Gene.Type",labels = c("Non-Cancer-Paralogs", "TSG"),values=c("#1874CD","#EE2C2C"))
  p1=ggarrange(prG41,prG42,nrow=1,ncol=2)
  
  
  kegg_gene_infor=merge(kegg_data_select,genes_infor)
  kegg_gene_infor$term<-gsub("GO_","",kegg_gene_infor$term)
  temp_kegg<-aggregate(.~term+Gene.Type,data=kegg_gene_infor[,c("term","Gene.Type","hg38_5UTR_length_log")],median)
  temp_kegg_tsg=temp_kegg[temp_kegg$Gene.Type=="TSG",]
  temp_kegg_non_cancer=temp_kegg[temp_kegg$Gene.Type=="Non-Cancer",]
  rownames(temp_kegg_non_cancer)<-temp_kegg_non_cancer[,1]
  diff<-temp_kegg_tsg[,3]-temp_kegg_non_cancer[temp_kegg_tsg[,1],3]
  names(diff)<-temp_kegg_tsg[,1]
  kegg_gene_infor$term<-factor(kegg_gene_infor$term,levels = names(sort(diff)))
  
  pG4_term<-plotGoTerm(kegg_gene_infor,col_x,"+",paste("Ratio of genes with",main_x),"greater")
  p3=ggarrange(pG4_term$p1,pG4_term$p2)
  
  p_G4<-getHist(col_x,genes_infor,paste("Ratio of genes with",main_x))
  p_cancerG42<-ggplot(p_G4$data,aes(x=cancer_type,y=ratio,fill=cancer_type))+geom_bar(stat = "identity")+ylab(paste("Ratio of genes with",main_x))+theme_hd()+theme(legend.position = "none",axis.text.x = element_text(angle = 90))#+annotate("text",x=c(2.5),y=c(0.97),label="**",size=8)+annotate("segment",x = c(1), xend = c(5), y = c(0.95), yend =  c(0.95))
  
  p2=ggarrange(p_G4$p1,p_cancerG42,nrow=1,ncol=2)
  
  return(list(p1=p1,p2=p_G4$p1,p3=p_cancerG42,p4=p3))
}



GetExpResult<-function(all_data_tsg,col_x){
  all_data_tsg_filtered=filter_low_data(all_data_tsg )
  all_data_tsg_filtered$diff_exp=abs(all_data_tsg_filtered$diff_exp)
  
  sig_exp<-get_diff(all_data_tsg,col_x,"norm_exp","+")
  pff_exp<-plotForesttwogroup(sig_exp,"expression in TCGA","4.2mm",paste(col_x,"+"),paste(col_x,"-"))
  
  sig_diff<-get_diff(all_data_tsg_filtered,col_x,"diff_exp","+")
  pff_diff<-plotForesttwogroup(sig_diff,"abs(log2FC) of  DOWN-TSG in TCGA","4.2mm",paste(col_x,"+"),paste(col_x,"-"))
  
  ggarrange(pff_exp$p1,pff_diff$p1,nrow=2,ncol=1)
}


get_evolution_chromatin<-function(genes_infor_data_all){
  tissue=c("kidney","lung","prostate","thyroid gland","colon","liver","breast epithelium")
  tissue=paste("class_",tissue,sep="")
  result_chro=c()
  for(i in 1:length(tissue) ){
    x1_non_cancer=as.data.frame(group_by(genes_infor_data_all[genes_infor_data_all$Gene.Type=="Non-Cancer",],age_type3,get(tissue[i])) %>% 
                                  dplyr::summarise(count_n = dplyr::n()) %>% mutate(ratio=count_n/sum(count_n)))
    x1_tsg=as.data.frame(group_by(genes_infor_data_all[genes_infor_data_all$Gene.Type=="TSG",],age_type3,get(tissue[i])) %>% 
                           dplyr::summarise(count_n = dplyr::n()) %>% mutate(ratio=count_n/sum(count_n)))
    result_chro=rbind(result_chro,cbind(tissue[i],type="Non-Cancer",x1_non_cancer[x1_non_cancer$`get(tissue[i])`=="More open",]))
    result_chro=rbind(result_chro,cbind(tissue[i],type="TSG",x1_tsg[x1_tsg$`get(tissue[i])`=="More open",]))
  }
  colnames(result_chro)<-c("tissue","Gene.Type","age_type3","type","count","Ratio")
  result_chro=as.data.frame(result_chro)
  return(as.data.frame(result_chro))
}


get_evolution_chromatin2<-function(genes_infor_data_all){
  tissue=c("kidney","lung","prostate","thyroid gland","colon","liver","breast epithelium")
  tissue=paste("class_",tissue,sep="")
  result_chro=c()
  for(i in 1:length(tissue) ){
    x1_non_cancer=as.data.frame(group_by(genes_infor_data_all[genes_infor_data_all$Gene.Type=="Non-Cancer",],get(tissue[i])) %>% 
                                  dplyr::summarise(count_n = dplyr::n()) %>% mutate(ratio=count_n/sum(count_n)))
    x1_tsg=as.data.frame(group_by(genes_infor_data_all[genes_infor_data_all$Gene.Type=="TSG",],get(tissue[i])) %>% 
                           dplyr::summarise(count_n = dplyr::n()) %>% mutate(ratio=count_n/sum(count_n)))
    result_chro=rbind(result_chro,cbind(tissue[i],type="Non-Cancer",x1_non_cancer[x1_non_cancer$`get(tissue[i])`=="More open",]))
    result_chro=rbind(result_chro,cbind(tissue[i],type="TSG",x1_tsg[x1_tsg$`get(tissue[i])`=="More open",]))
  }
  colnames(result_chro)<-c("tissue","Gene.Type","type","count","Ratio")
  result_chro=as.data.frame(result_chro)
  return(as.data.frame(result_chro))
}


plot_chromatin_evo<-function(chromatin_data){
  chromatin_data$tissue=sub(" DHS signal","", chromatin_data$tissue)
  chromatin_2= ggplot(data = chromatin_data,aes(x=factor(age_type3),y=Ratio,col=Gene.Type,group=Gene.Type))+geom_point(size=4)+geom_line()+facet_wrap(~tissue,ncol=3)+theme_hd()+scale_color_manual(name = "Gene.Type",values=c("#1874CD","#EE2C2C"))+xlab("Gene age (Million years)")+guides(col=guide_legend(nrow=2,byrow=TRUE))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90),strip.text = element_text(face = "bold"),
          strip.background = element_blank())+ylab("Ratio of genes with more open chromatin status")
  return(chromatin_2)
}

PlotSpecies<-function(all_data,col_x){
  comp=list(c("TSG", "Non-Cancer"))
  max_x=max(all_data[,col_x])*1.2
  p=ggplot(data=all_data[all_data$species!="Cow",],aes(x=Gene.Type,y=get(col_x),fill=Gene.Type))+geom_boxplot(show.legend = FALSE,width=0.4)+
    geom_signif(comparisons = comp,map_signif_level=TRUE,col="black", test='wilcox.test',color="red", textsize = 6.88,size=0, tip_length = 0.00)+
    facet_wrap(~species,ncol = 1)+coord_flip()+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd_minimal()+ylab(col_x)+ylim(0,max_x)+theme(axis.text.x=element_text(angle = 90) )#+stat_summary(fun.data =  countFunction, geom="text", size = 3,col="black")
  return(p)
}

PlotMutipleCor<-function(all_data_tsg,t1,t2,t1_name,t2_name,x,title_x,cor_x){
  cor_cc<-get_cor_for_each_cancer_all(all_data_tsg,t1,x,cor_x)
  cor_gg<-get_cor_for_each_cancer_all(all_data_tsg,t2,x,cor_x)
  all_plot=rbind(cbind(type_x=t2_name,cor_gg),cbind(type_x=t1_name,cor_cc))
  cor_matrix1= aggregate(.~type_x,all_plot[,-c(2:3)],median)
  print("median")
  print(cor_matrix1)
  print(all_plot)
 
  #cor_matrix2= aggregate(.~type_x,all_plot[,-c(2:3)],BiocGenerics::mad)
  print("MAD")
  #print(cor_matrix2)
  all_plot$Pvalue<- -log10(as.numeric(all_plot$Pvalue))
  
  cor_plot2=ggplot(all_plot,aes(y=Cancer,x=cor,col=type_x))+  geom_point(aes(size=Pvalue))+labs(size = "-Log10(P-value)")+scale_color_manual("Strand",values=c( "#F07F16C7","#40A4DE"))+ggtitle(title_x)+  scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15))+guides(color=guide_legend(nrow=1,byrow=TRUE,override.aes = list(size=5)))+xlab("Pearson's r")+geom_hline(yintercept = 0.4,size=1.5)+geom_vline(xintercept = 0,size=1)+theme_hd_minimal2()+theme(legend.position = "bottom")+ylab("Dataset")
  return(cor_plot2)
}




GetMutipleCorMethAll<-function(all_data_tsg,t1,t2,t1_name,t2_name,x,cor_x){
  cor_cc<-get_cor_for_each_cancer_all(all_data_tsg,t1,x,cor_x)
  cor_gg<-get_cor_for_each_cancer_all(all_data_tsg,t2,x,cor_x)
  all_plot=rbind(cbind(type_x=t2_name,cor_gg),cbind(type_x=t1_name,cor_cc))
  all_plot$Pvalue<- -log10(as.numeric(all_plot$Pvalue))
  
  return(all_plot)
}




GetMutipleCorMeth<-function(all_data_tsg,t1,t2,t1_name,t2_name,x,cor_x,type_x){
  s_x<-unique(all_data_tsg[,type_x])
  all_plot<-c()
  for(s in s_x){
    cor_cc<-get_cor_for_each_cancer_all(all_data_tsg[all_data_tsg[,type_x]==s,],t1,x,cor_x)
    cor_gg<-get_cor_for_each_cancer_all(all_data_tsg[all_data_tsg[,type_x]==s,],t2,x,cor_x)
    all_plot=rbind(all_plot,cbind(type_x=paste(s,"_","G-run-",t2_name,sep=""),cor_gg),cbind(type_x=paste(s,"_","G-run-",t1_name,sep=""),cor_cc))
  }
  all_plot$type_x=gsub("Template","T", gsub("Non-template","NT",all_plot$type_x))


  all_plot$Pvalue<- -log10(as.numeric(all_plot$Pvalue))
  
  return(all_plot)
}


PlotMutipleCorMeth<-function(all_data_tsg,t1,t2,t1_name,t2_name,x,title_x,cor_x,type_x){
  s_x<-unique(all_data_tsg[,type_x])
  all_plot<-c()
  for(s in s_x){
    cor_cc<-get_cor_for_each_cancer_all(all_data_tsg[all_data_tsg[,type_x]==s,],t1,x,cor_x)
    cor_gg<-get_cor_for_each_cancer_all(all_data_tsg[all_data_tsg[,type_x]==s,],t2,x,cor_x)
    all_plot=rbind(all_plot,cbind(type_x=paste(s,"_","G-run-",t2_name,sep=""),cor_gg),cbind(type_x=paste(s,"_","G-run-",t1_name,sep=""),cor_cc))
  }
  all_plot$type_x=gsub("Template","T", gsub("Non-template","NT",all_plot$type_x))
  cor_matrix1= aggregate(.~type_x,all_plot[,-c(2:3)],median)
  print("median")
  print(cor_matrix1)
  
  #cor_matrix2= aggregate(.~type_x,all_plot[,-c(2:3)],BiocGenerics::mad)
  print("MAD")
  #print(cor_matrix2)
  all_plot$Pvalue<- -log10(as.numeric(all_plot$Pvalue))
  
  cor_plot2=ggplot(all_plot,aes(y=Cancer,x=cor,col=type_x))+  
    geom_point(aes(size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+  scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15))+guides(color=guide_legend(nrow=1,byrow=TRUE,override.aes = list(size=15)))+xlab("Pearson's r")+geom_hline(yintercept = 0.4,size=1.5)+geom_vline(xintercept = 0,size=1)+theme_hd_minimal2()+theme(legend.position = "bottom")+ylab("Dataset")
  return(cor_plot2)
}

PlotMutipleCorMethExp<-function(all_data_tsg,t1,t1_name,x,title_x,cor_x,type_x){
  s_x<-unique(all_data_tsg[,type_x])
  all_plot<-c()
  for(s in s_x){
    cor_cc<-get_cor_for_each_cancer_all(all_data_tsg[all_data_tsg[,type_x]==s,],t1,x,cor_x)
    all_plot=rbind(all_plot,cbind(type_x=paste("Meth-",s,"_",t1_name,sep=""),cor_cc))
  }
  cor_matrix1= aggregate(.~type_x,all_plot[,-c(2:3)],median)
  print("median")
  print(cor_matrix1)
  all_plot$type_x=gsub("Template","T", gsub("Non_template","NT",all_plot$type_x))
  
  #cor_matrix2= aggregate(.~type_x,all_plot[,-c(2:3)],BiocGenerics::mad)
  print("MAD")
  #print(cor_matrix2)
  all_plot$Pvalue<- -log10(as.numeric(all_plot$Pvalue))
  
  cor_plot2=ggplot(all_plot,aes(y=Cancer,x=cor,col=type_x))+  
    geom_point(aes(size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+  scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15))+guides(color=guide_legend(nrow=1,byrow=TRUE,override.aes = list(size=5)))+xlab("Pearson's r")+geom_hline(yintercept = 0.4,size=1.5)+geom_vline(xintercept = 0,size=1)+theme_hd_minimal2()+theme(legend.position = "bottom")+ylab("Dataset")
  return(cor_plot2)
}

get_CpG_promoter_other_species<-function(promoter_cpg_OE_x,unique_orthologs,title_x){
  if(title_x=="Human"){unique_orthologs_x=unique(unique_orthologs[,c("Gene.stable.ID","Gene.Type")])}else{
    unique_orthologs_x=unique(unique_orthologs[,c(paste(title_x,".gene.stable.ID",sep=""),"Gene.Type")])
  }
  TSG_non_cancer<-unique_orthologs_x[unique_orthologs_x$Gene.Type!="Oncogene",]
  non_cancer_data=unique_orthologs_x[unique_orthologs_x$Gene.Type=="Non-Cancer",]
  ## gene-noncancer pairs where gene-TSG exist at the same time
  non_cancer_dup=non_cancer_data[non_cancer_data[,1] %in% unique_orthologs_x[unique_orthologs_x$Gene.Type=="TSG",1],]
  print(dim(non_cancer_dup))
  print("There are some genes that have both TSG and non cancer orthologs\n but the gene is regarded as TSG as long as it is a orthologs of human TSG ")
  if(title_x!="Human"){
    print( unique(unique_orthologs[unique_orthologs[,paste(title_x,".gene.stable.ID",sep="")]==non_cancer_dup[5,1],c("Gene.stable.ID", paste(title_x,".gene.stable.ID",sep=""),"Gene.Type")]) )
  }
  #### delete the gene-noncancer pairs if gene-TSG exist
  unique_orthologs_select=TSG_non_cancer[! paste(TSG_non_cancer[,1],TSG_non_cancer[,2],sep="_") %in% paste(non_cancer_dup[,1],non_cancer_dup[,2],sep="_"),]
  print("thus I keeep those genes and regard them as TSG")
  rownames(unique_orthologs_select)<-as.character(unique_orthologs_select[,1])
  promoter_cpg_OE_x<-promoter_cpg_OE_x[promoter_cpg_OE_x$Group.1 %in% unique_orthologs_select[,1],]
  promoter_cpg_OE_x$Gene.Type=unique_orthologs_select[promoter_cpg_OE_x$Group.1,2]
  comp=list(c("TSG", "Non-Cancer"))
  promoter_cpg_OE_x$species<-title_x
  return(promoter_cpg_OE_x)
}




GetAlldata<-function(promoter_elegan,promoter_Fruitfly,promoter_zebra,promoter_cow,promoter_pig,promoter_mouse,promoter_human,unique_orthologs){
  data1=get_CpG_promoter_other_species(promoter_elegan,unique_orthologs,"Caenorhabditis.elegans")
  data2=get_CpG_promoter_other_species(promoter_Fruitfly,unique_orthologs,"Drosophila.melanogaster")
  data3=get_CpG_promoter_other_species(promoter_zebra,unique_orthologs,"Zebrafish")
  data4=get_CpG_promoter_other_species(promoter_cow,unique_orthologs,"Cow")
  data5=get_CpG_promoter_other_species(promoter_pig,unique_orthologs,"Pig")
  data6=get_CpG_promoter_other_species(promoter_mouse,unique_orthologs,"Mouse")
  data7=get_CpG_promoter_other_species(promoter_human,unique_orthologs,"Human")
  all_data=rbind(data1,data2,data3,data4,data5,data6,data7)
  level_x=c("Caenorhabditis.elegans","Drosophila.melanogaster","Zebrafish","Cow","Pig","Mouse","Human")
  all_data$species<-factor(all_data$species,levels = level_x[7:1])
  all_data_sd=apply(all_data[,2:5],1,sd)
  all_data$sd=all_data_sd
  return(all_data)
}



plotSubSigParied<- function(x.data,xlab,ylab,sub1,title.x,cols1="",xlab.name,ylab.name,xlab1,ylab2,cols_x="black"){
  x.data<-x.data[,c(xlab,ylab,sub1)]
  colnames(x.data)<-c("x","y","group")
  x_df<-reshape2::melt(x.data)
  print(max(x_df$value))
  x_df$variable<-factor(x_df$variable,levels = c("x","y"))
  p1=ggplot(data = x_df,aes(x=group,y=value))+geom_boxplot(aes(fill=variable,col=I(cols_x)))+xlab(xlab1)+ylab(ylab2)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle =30,hjust = 1))
  if(length(cols1)>1){
    p1=p1+scale_fill_manual(values=cols1,labels=c(xlab.name,ylab.name))}
  print(sub1)
  stat.test<-c()
  for(g in levels(x.data$group)){
    temp_x=x.data[x.data$group==g,]
    wt<-wilcox.test(temp_x$x,temp_x$y,paired=T)
    stat.test<-c(stat.test,wt$p.value)
  }
  stat.test<-pvalueTosig(stat.test)
  print(unique(x.data$group))
  print(stat.test)
    p1=p1 +  geom_signif(tip_length = 0.01,xmin=1:length(stat.test)-0.1 ,xmax=1:length(stat.test)+0.1 ,annotations=stat.test, y_position=max(x_df$value,na.rm=T),fontface="bold") 
  return(p1)
}



plot_scatter_density<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggplot(data = all_data_tsg,aes(x=get(x_lab),y=get(y_lab)))+geom_bin2d(bins=70)+
    scale_fill_continuous(type="viridis")+stat_cor(method = "pearson",fontface="bold",output.type="text")+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",fontface="bold",output.type="text") +geom_smooth(method = "lm")
  return(sp)
}


getPairedPlotMultiHK<-function(data_x1,con1,con2,con1_lab,con2_lab,class_x){
  data_x=rbind(data.frame(value=data_x1[,con1],type=con1_lab,age=data_x1[,class_x]),data.frame(value=data_x1[,con2],type=con2_lab,age=data_x1[,class_x]))

  data_x$type_age<-paste(data_x$age,data_x$type)
  print(aggregate(data_x$value,list(data_x$type_age),median))
  compares_x=getComp(unique(data_x$type_age))[c(2,5,1,6)]
  data_x$type_age<-factor(data_x$type_age,levels=names(table(data_x$type_age))[c(4,3,2,1)])
  p2=ggplot(data_x[!is.na(data_x$age),], aes(x=type_age,y=value,fill=type_age),)+
    geom_boxplot()+theme_hd_minimal()+geom_signif(comparisons=compares_x, test='wilcox.test', map_signif_level=TRUE,col="black",step_increase=0.1) +
    scale_fill_manual(values=c("#BDE1FF", "#FFCDA6","#40A4DE", "#F07F16C7"))+theme(axis.text.x = element_text(angle = 30,hjust = 1))+
    ylab("Frequency of G-runs")+xlab("Genes")
    p2=p2+theme(axis.title.x = element_blank())+guides(fill="none")
  return(p2)
} 


getPairedPlot<-function(data_x1,con1,con2,con1_lab,con2_lab){
  data_x=as.data.frame(rbind(data.frame(value=data_x1[,con1],type=con1_lab),
                             data.frame(value=data_x1[,con2],type=con2_lab)))
  data_x$type<-factor(data_x$type,levels=c(con1_lab,con2_lab))
  p2=ggplot(data_x, aes(x=type,y=value,fill=type,color=I("gray78")),)+geom_boxplot()+theme_hd()+stat_compare_means(method="wilcox.test", paired=TRUE, aes(label =  ..p.signif..),label.x.npc="center",size=4.8,fontface="bold")+scale_fill_manual("Strand",values=c("#40A4DE", "#F07F16C7"))+ylab("Frequency of G-runs")
  p2=p2+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
  return(p2)
}  

getPairedPlotSpe<-function(data_x1,con1,con2,con1_lab,con2_lab,class_x){
  data_x1=data_x1[data_x1$species!="Drosophila.melanogaster",]
  data_x=rbind(data.frame(value=data_x1[,con1],type=con1_lab,species=data_x1[,class_x]),data.frame(value=data_x1[,con2],type=con2_lab,species=data_x1[,class_x]))
  data_x$species<-gsub("\\.","\n",data_x$species)
  data_x$species=factor(data_x$species,levels=c("Human","Mouse","Pig","Cow","Zebrafish","Caenorhabditis\nelegans"))
  data_x$type=factor(data_x$type,levels=c(con1_lab,con2_lab))
  p2=ggplot(data_x, aes(x=type,y=value,fill=type,color=I("gray78")),)+geom_boxplot()+theme_hd_minimal()+stat_compare_means(method="wilcox.test", paired=TRUE, aes(label =  ..p.signif..),label.x.npc="center",size=4.8,fontface="bold")+scale_fill_manual("Strand",values=c("#40A4DE", "#F07F16C7"))+ylab("Frequency of G-runs")
  p2=p2+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
  
  return(p2+facet_wrap(~species,nrow=1))
}  

getPairedPlotage<-function(data_x1,con1,con2,con1_lab,con2_lab,class_x){
  data_x=rbind(data.frame(value=data_x1[,con1],type=con1_lab,age=data_x1[,class_x]),data.frame(value=data_x1[,con2],type=con2_lab,age=data_x1[,class_x]))
  data_x$age=factor(data_x$age,levels=names(table(data_x$age)))
  data_x$type=factor(data_x$type,levels=c(con1_lab,con2_lab))
  p2=ggplot(data_x[!is.na(data_x$age),], aes(x=type,y=value,fill=type,color=I("gray78")),)+geom_boxplot()+theme_hd_minimal()+stat_compare_means(method="wilcox.test", paired=TRUE, aes(label =  ..p.signif..),label.x.npc="center",size=4.8,fontface="bold")+scale_fill_manual("Strand",values=c("#40A4DE", "#F07F16C7"))+ylab("Frequency of G-runs")
  p2=p2+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
  
  return(p2+facet_wrap(~age,nrow=1))
}  

getPairedPlotage2<-function(data_x1,con1,con2,con1_lab,con2_lab,class_x){
  data_x=rbind(data.frame(value=data_x1[,con1],type=con1_lab,age=data_x1[,class_x]),data.frame(value=data_x1[,con2],type=con2_lab,age=data_x1[,class_x]))
  data_x$age=factor(data_x$age,levels=rev(names(table(data_x$age))))
  data_x$type=factor(data_x$type,levels=c(con1_lab,con2_lab))
  p2=ggplot(data_x[!is.na(data_x$age),], aes(x=age,y=value,fill=age,color=I("gray78")),)+geom_boxplot()+theme_hd_minimal()+stat_compare_means(method="wilcox.test", aes(label =  ..p.signif..),label.x.npc="center",size=4.8,fontface="bold")+scale_fill_manual("Strand",values=c( "#D0ECFA","#3D55B3"))+ylab("Frequency of G-runs")
  p2=p2+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
  
  return(p2+facet_wrap(~type,nrow=1))
}  



plotStrandCor<-function(all_data,x.name,y.name,y1.name,y_lab,method="pearson"){
  cor_age_tsg<-get_cor_agetype(all_data,x.name,y.name,y_lab)
  cor_age_tsg1<-get_cor_agetype(all_data,x.name,y1.name,y_lab)
  strand_x=ifelse(y.name=="CC_density","Template","Non-template")
  strand_x1=ifelse(y1.name=="CC_density","Template","Non-template")
  all_cor_exp_age<-rbind(cbind(strand=strand_x,cor_age_tsg),cbind(strand=strand_x1,cor_age_tsg1))
  p1=ggplot(data=all_cor_exp_age,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd_minimal()+xlab("Gene Age Groups\n(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+ylab("Pearson's r(Expression\nVS Frequency of G-runs)")+ggtitle("TCGA cohort  normal")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))+facet_wrap(~strand,scales = "free_y")
  return(p1)
}

plotStrandCorBarplot<-function(all_data,x.name,y.name,y1.name,y_lab,method="pearson"){
  cor_age_tsg<-get_cor_agetype_single(all_data,x.name,y.name,y_lab)
  cor_age_tsg1<-get_cor_agetype_single(all_data,x.name,y1.name,y_lab)
  strand_x=ifelse(y.name=="CC_density","Template","Non-template")
  strand_x1=ifelse(y1.name=="CC_density","Template","Non-template")
  all_cor_exp_age<-rbind(cbind(strand=strand_x,cor_age_tsg),cbind(strand=strand_x1,cor_age_tsg1))
  p1=ggplot(data=all_cor_exp_age,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age Groups\n(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(color=guide_legend(nrow=2,byrow=TRUE))+ylab("Pearson's r(Expression\nVS Frequency of G-runs)")+ggtitle("TCGA cohort  normal")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))+facet_wrap(~strand)
  return(p1)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=TRUE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac) <- gsub("mean",measurevar,colnames(datac))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



getGeneRegionInterval<-function(genes_infor){
  
  re_all<-c()
  up=2500
  down=2500
  for(i in 1:length(genes_infor$ensembl_gene_id)){
    if(i%%100==0)
    {print(i)}
    temp=genes_infor[i,]
    start=temp$transcript_start
    end=temp$transcript_end
    inter1<-levels(cut(start:end,50,dig.lab=40))
    inter2=t(apply(str_split_fixed(gsub("\\(|]","", inter1),",",4)[,1:2],1,as.numeric))
    inter2[1,1]=start
    inter2[length(inter2[,1]),2]=end
    
    start=temp$transcript_end
    end=temp$transcript_end+2500
    inter_down1<-levels(cut(start:end,50,dig.lab=40))
    inter_down2=t(apply(str_split_fixed(gsub("\\(|]","", inter_down1),",",4)[,1:2],1,as.numeric))
    inter_down2[1,1]=start
    inter_down2[length(inter_down2[,1]),2]=end
    
    start=temp$transcript_start-2500
    end=temp$transcript_start
    inter_up1<-levels(cut(start:end,50,dig.lab=40))
    inter_up2=t(apply(str_split_fixed(gsub("\\(|]","", inter_up1),",",4)[,1:2],1,as.numeric))
    inter_up2[1,1]=start
    inter_up2[length(inter_up2[,1]),2]=end
    
    if(temp$strand=="1"){seq_id_up=-50:-1
    seq_id_body=1:50
    seq_id_down=51:100
    }else{seq_id_up=100:51
    seq_id_body=50:1
    seq_id_down=-1:-50
    }
    re_all=rbind(re_all,cbind(chr=paste("chr",temp$chromosome_name,sep=""),rbind(inter_up2,inter2,inter_down2),temp$ensembl_gene_id,c(seq_id_up,seq_id_body,seq_id_down),temp$strand))
  }
  re_all[,6]=ifelse(re_all[,6]==1,"+","-")
  return(re_all)
}




plot_Meth<-function(data_x,title_x,y_name){
  meth_se=summarySE_median(data_x, measurevar="norm_meth", groupvars=c("type","pos"))
  temp=str_split_fixed(meth_se$pos,"_",2)
  neg=ifelse(as.numeric(temp[,2])>as.numeric(temp[,1]),1,-1)
  meth_se$position<-(as.numeric(temp[,1])/2+as.numeric(temp[,2])/2+50)*neg
  pd<-position_dodge(0.1)
  p=ggplot(meth_se,aes(x=position,y=norm_meth,col=type,group=type))  +
    geom_errorbar(aes(ymin=norm_meth-1.5*se, ymax=norm_meth+1.5*se),colour="gray78", position=pd,width=150) +ylab("Methylation")+ #+facet_wrap(~Gene.Type,nrow=2)+
    geom_point(position=pd,size=1, shape=21, fill="white")+geom_line(position=pd,size=1.5)+xlab("distance to TSS")+
    scale_color_manual(values=c("#40A4DE", "#F07F16C7"))+theme_hd()
  
}
# no strand specific metjylation
get_aver_meth_previous_wrong<-function(result_meth,start,end){
  downstream<-read.table(paste("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/data/hm450.hg38.manifest_",start,"_",end,".bed",sep=""))
  result_meth$result_gene_probe=paste(result_meth_all$transcript_id, result_meth$probe)
  downstream$result_gene_probe=paste(downstream[,4], downstream[,10])
  downstream_NT<-downstream[paste(downstream$V6,downstream$V12) %in% c("+ +","- -"),]
  downstream_T<-downstream[!paste(downstream$V6,downstream$V12) %in% c("+ +","- -"),]
  downstream_meth<-result_meth[result_meth$result_gene_probe %in% c(downstream$result_gene_probe),]
  downstream_meth$type=ifelse(downstream_meth$result_gene_probe%in% downstream_NT$result_gene_probe ,"Meth-NT","Meth-T")
  aver_downstream_meth_strand<-aggregate(downstream_meth[,c("diff_meth","tumor_meth","norm_meth")],list(cancer=downstream_meth$cancer,ensembl_trans_id=downstream_meth$transcript_id,type=downstream_meth$type),function(x) mean(x,na.rm=T))
  print(dim(downstream_meth))
  aver_downstream_meth_strand$type=factor(aver_downstream_meth_strand$type,levels=c("Meth-T","Meth-NT"))
  return(aver_downstream_meth_strand)
}
get_aver_meth<-function(result_meth,start,end){
  downstream<-read.table(paste("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/data/hm450.hg38.manifest_",start,"_",end,".bed",sep=""))
  result_meth$result_gene_probe=paste(result_meth_all$transcript_id, result_meth$probe)
  
   downstream$result_gene_probe=paste(downstream[,4], downstream[,10])
   downstream_meth<-result_meth[result_meth$result_gene_probe %in% c(downstream$result_gene_probe),]
   
  aver_downstream_meth<-aggregate(downstream_meth[,c("diff_meth","tumor_meth","norm_meth")],list(cancer=downstream_meth$cancer,ensembl_gene_id=downstream_meth$transcript_id),function(x) mean(x,na.rm=T))
  return(aver_downstream_meth)
}


get_aver_dhs<-function(result_meth,start,end){
  downstream<-read.table(paste("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/data/hm450.hg38.manifest_",start,"_",end,".bed",sep=""))
  result_meth$result_gene_probe=paste(result_meth_all$transcript_id, result_meth$probe)
  
  downstream$result_gene_probe=paste(downstream[,4], downstream[,10])
  downstream_meth<-result_meth[result_meth$result_gene_probe %in% c(downstream$result_gene_probe),]
  
  aver_downstream_meth<-aggregate(downstream_meth[,c("diff_meth","tumor_meth","norm_meth")],list(cancer=downstream_meth$cancer,ensembl_gene_id=downstream_meth$transcript_id),function(x) mean(x,na.rm=T))
  return(aver_downstream_meth)
}





ggCorPlotNew<-function(all_meth,all_cg,ylab_x="norm_meth"){
  pos<-unique(all_meth$pos)
  result_all<-c()
  for(p in pos){
    print(p)
    temp_gg<-all_cg[all_cg$pos==p,c("ensembl_gene_id","ensembl_trans_id","CC_density","GG_density","pos")]
    for(pp in pos){
      print(paste(p,pp))
      temp_gg2<-all_meth[all_meth$pos==pp,c("cancer",ylab_x,"ensembl_trans_id")]
      temp_gg_other<-merge(temp_gg,temp_gg2,by=c("ensembl_trans_id"))
      
      meth=GetMutipleCorMethAll(temp_gg_other,"CC_density","GG_density","Template strand","Non-template strand",ylab_x,"pearson")
      
      result_all<-rbind(result_all,cbind(pos_gg=p,pos_meth=pp,meth))
      
    }
    
  }
  
  
  return(result_all)
  
}

ggCorPlotNewExp<-function(all_data_tsg,all_cg){
  pos<-unique(all_cg$pos)
  result_all<-c()
  for(p in pos){
    print(p)
    temp_gg<-all_cg[all_cg$pos==p,c("ensembl_gene_id","ensembl_trans_id","CC_density","GG_density","pos")]
    temp_gg_other<-merge(temp_gg,all_data_tsg,by.x=c("ensembl_trans_id"),by.y="ensembl_transcript_id")
    exp=GetMutipleCorMethAll(temp_gg_other,"CC_density.x","GG_density.x","Template strand","Non-template strand","norm_exp","pearson")
      result_all<-rbind(result_all,cbind(pos_gg=p,exp))
  }
  return(result_all)
}

MethCorPlotNewExp<-function(all_data_tsg,all_meth,ylab_x="norm_exp",ylab_y="norm_meth"){
  pos<-unique(all_meth$pos)
  result_all<-c()
  colnames(all_data_tsg)<-gsub("ensembl_transcript_id","ensembl_trans_id",colnames(all_data_tsg))
  for(p in pos){
    print(p)
    temp_gg<-all_meth[all_meth$pos==p,c("cancer",ylab_y,"ensembl_trans_id")]
    temp_gg_other<-merge(temp_gg,all_data_tsg,by=c("cancer","ensembl_trans_id"))
    exp=get_cor_for_each_cancer_all(temp_gg_other,ylab_y,ylab_x,"pearson")
    result_all<-rbind(result_all,cbind(pos_meth=p,exp))
  }
  return(result_all)
}

MethCorPlotNewSignal<-function(all_dhs,all_meth){
  pos<-unique(all_dhs$pos)
  pos1<-unique(all_dhs$pos)
  result_all<-c()
  for(p1 in pos1){
    temp= all_dhs[all_dhs$pos==p1,c("cancer","signal","ensembl_trans_id")]
    for(p in pos){
    print(p)
    temp_gg<-all_meth[all_meth$pos==p,c("cancer","norm_meth","ensembl_trans_id")]
    temp_gg_other<-merge(temp_gg,temp,by=c("cancer","ensembl_trans_id"))
    exp=get_cor_for_each_cancer_all(temp_gg_other,"norm_meth","signal","pearson")
    result_all<-rbind(result_all,cbind(pos_meth=p,pos_signal=p1,exp))
  }
  }
  
  return(result_all)
}




ggCorPlotNewSub<-function(all_meth_gene,all_cg){
  pos<-unique(all_meth_gene$pos)
  result_all<-c()
  for(p in pos){
    print(p)
    temp_gg<-all_cg[all_cg$pos==p,c("ensembl_gene_id","ensembl_trans_id","CC_density","GG_density","pos")]
    for(pp in pos){
      print(paste(p,pp))
      temp_gg2<-all_meth_gene[all_meth_gene$pos==pp,c("cancer","norm_meth","ensembl_trans_id","G4_predict_type","G4_down_predict_pos")]
      temp_gg_other<-merge(temp_gg,temp_gg2,by=c("ensembl_trans_id"))
      
      meth=GetMutipleCorMeth(temp_gg_other,"CC_density","GG_density","Template strand","Non-template strand","norm_meth","pearson","type")
      
      result_all<-rbind(result_all,cbind(pos_gg=p,pos_meth=pp,meth))
      
    }
    
  }
  
  
  return(result_all)
  
}



ggCorPlot<-function(all_meth_strand2){
  pos<-unique(all_meth_strand2$pos)
  result_all<-c()
  for(p in pos){
    print(p)
    temp_gg<-all_meth_strand2[all_meth_strand2$pos==p,c("pos_CC_density","pos_GG_density","cancer","type","ensembl_trans_id")]
    for(pp in pos){
      print(paste(p,pp))
      temp_gg2<-all_meth_strand2[all_meth_strand2$pos==pp,c("pos_CC_density","cancer","pos_GG_density","type","norm_meth","ensembl_trans_id")]
      temp_gg_other<-merge(temp_gg,temp_gg2,by=c("cancer", "type","ensembl_trans_id"))
      
      meth=GetMutipleCorMeth(temp_gg_other,"pos_CC_density.x","pos_GG_density.x","Template strand","Non-template strand","norm_meth","pearson","type")
      
      result_all<-rbind(result_all,cbind(pos_gg=p,pos_meth=pp,meth))
      
    }
    
  }
  
  
  return(result_all)
  
}

getMultipleRegionFragment<-function(genes_infor,start,end){
  gene_pos<-genes_infor[,c("genes_main_transcript","chromosome_name","transcript_start","transcript_end","strand","ensembl_gene_id")]
  colnames(gene_pos)<-c("trans_ID","seqnames","start","end","strand","gene_id")
  gene_pos$strand[gene_pos$strand=="1"]="+"
  gene_pos$strand[gene_pos$strand=="-1"]="-"
  gene_pos_data<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=T)
  genes_promoter_up=get_gene_region(gene_pos_data,start,end)
  temp_promoter_up<-as.data.frame(genes_promoter_up)
  temp_promoter_up$seqnames<-paste("chr",temp_promoter_up$seqnames,sep="")
  write.table(temp_promoter_up[,c("seqnames","start","end","trans_ID","width","strand")],file=paste("./data/hg38_gene_promoter_",start,"_",end,"_chr.txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
}

plotMethPos<-function(data_x,cancer,genes=""){
  if(length(genes)<2){genes=unique(data_x$ensembl_gene_id)}
  temp=data_x[data_x$cancer==cancer & data_x$ensembl_gene_id %in% genes,c("ensembl_gene_id","pos","type","norm_meth")]
  temp$x=paste(temp$pos,temp$ensembl_gene_id)
  temp_df<-reshape2::dcast(temp[,c("x","type","norm_meth")],x~type)
  temp_df$pos=str_split_fixed(temp_df$x," ",2)[,1]
  temp1=str_split_fixed(temp_df$pos,"_",2)
  neg=ifelse(as.numeric(temp1[,2])>as.numeric(temp1[,1]),1,-1)
  temp_df$position<-(as.numeric(temp1[,1])/2+as.numeric(temp1[,2])/2+50)*neg
 print(dim(temp))
  temp_x=aggregate(cbind(T=temp_df$`Meth-T`,NT=temp_df$`Meth-NT`,dif=temp_df$`Meth-NT`-temp_df$`Meth-T`),list(Pos=temp_df$position),function(x)median(x,na.rm=T))
  
  x=melt(temp_x,id.vars = "Pos")
  x1=melt(temp_df,id.vars = c("x","pos","position"))
  #ggplot(x1[abs(x1$position)<1000,],aes(x=factor(position),y=value,col=variable))+geom_boxplot()
  
  pd<-position_dodge(0.1)
  p=ggplot(x,aes(x=Pos,y=value,col=variable))+geom_point()+geom_line(position=pd,size=1.5)+theme_hd()
  print(p)
  return(x)
}


