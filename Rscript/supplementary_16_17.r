# Set the working directory to the specified path
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new")

# Load necessary functions and data
source('functions/functions_evolution_5UTR.r')  # Load custom functions for 5' UTR evolution analysis
load(file="./data/all_meth_exp_signal_promoter.rdata")  # Load methylation and expression signal data

# Load necessary libraries
library(readr)

# Load data for various structural features
MR <- read_delim("data/687605ea2bc26_zip/687605ea2bc26_MR.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
DR <- read_delim("data/6875c9fe1c68e_zip/6875c9fe1c68e_DR.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Z <- read_delim("data/687610f726fe0_zip/687610f726fe0_Z.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
STR <- read_delim("data/68761a097a212_zip/68761a097a212_STR.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Add structural feature information to genes_infor
genes_infor$STR <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(STR$Sequence_name, "ENSG", 2)[, 1], "+", "-")
genes_infor$MR <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(MR$Sequence_name, "ENSG", 2)[, 1], "+", "-")
genes_infor$DR <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(DR$Sequence_name, "ENSG", 2)[, 1], "+", "-")
genes_infor$Z <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(Z$Sequence_name, "ENSG", 2)[, 1], "+", "-")

# Load positional data for structural features
MR_pos <- read_delim("data/6878a80c3a7db_zip/6878a80c3a7db_MR_pos.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
DR_pos <- read_delim("data/6878e76edc2f5_zip/6878e76edc2f5_DR_pos.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Z_pos <- read_delim("data/6878befbbdf58_zip/6878befbbdf58_Z_pos.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
STR_pos <- read_delim("data/6878d108b3c8e_zip/6878d108b3c8e_STR_pos.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Load G4 data
G4_pos <- read_delim("data/687c6d13951f3_zip/687c6d13951f3_GQ_pos.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
G4 <- read_delim("data/687c7b979deb3_zip/687c7b979deb3_GQ.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
predict_G4_cds_downstream <- readRDS(file = "./data/predict_G4_cds_downstream.rds")

# Add G4 information to genes_infor
genes_infor$G4_pos <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(G4_pos$Sequence_name, "ENSG", 2)[, 1], "+", "-")
genes_infor$G4 <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(G4$Sequence_name, "ENSG", 2)[, 1], "+", "-")

# Convert positional data to numeric
temp_pos <- predict_G4_cds_downstream$pos
DR <- as.numeric(DR_pos$Start)
Z <- as.numeric(Z_pos$Start)
STR <- as.numeric(STR_pos$Start)
MR <- as.numeric(MR_pos$Start)
G4 <- as.numeric(temp_pos[, 8])  # Confirm that the 8th column is Start

# Create a data frame for density plot
df <- tibble(
  value = c(DR, Z, STR, MR, G4),
  type = rep(c("DR", "Z-DNA", "STR", "MR", "G4"),
             times = c(length(DR), length(Z), length(STR), length(MR), length(G4)))
)

# Plot density of structural feature start sites
p_density <- ggplot(df, aes(value, fill = type)) +
  geom_density(alpha = 0.45, size = 0.8) +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")) +
  labs(x = "Genomic position (Start)",
       y = "Density",
       title = "Density of structural feature start sites") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  ) + facet_wrap(~type, nrow = 1)

# Add positional information to genes_infor
genes_infor$STR_pos <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(STR_pos$Sequence_name, "ENSG", 2)[, 1], "+", "-")
genes_infor$MR_pos <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(MR_pos$Sequence_name, "ENSG", 2)[, 1], "+", "-")
genes_infor$DR_pos <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(DR_pos$Sequence_name, "ENSG", 2)[, 1], "+", "-")
genes_infor$Z_pos <- ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(Z_pos$Sequence_name, "ENSG", 2)[, 1], "+", "-")

# Load necessary libraries for Venn diagram
library(ggVennDiagram)
library(ggplot2)

# Construct a list of gene sets for Venn diagram
x <- list(
  G4 = genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_pos == "+"],
  Short_tandem_repeats = genes_infor$ensembl_gene_id[genes_infor$STR_pos == "+"],
  Direct_repeats = genes_infor$ensembl_gene_id[genes_infor$DR_pos == "+"],
  Mirror_repeats_and_triplex = genes_infor$ensembl_gene_id[genes_infor$MR_pos == "+"],
  Z_DNA = genes_infor$ensembl_gene_id[genes_infor$Z_pos == "+"]
)

other<-unique(c(x$Z_DNA,x$Direct_repeats,x$Mirror_repeats_and_triplex,x$Short_tandem_repeats))
G4_alone<-x$G4[!x$G4%in% other]
G4_cor<-PlotMutipleCor(GTEXmedian_tpm_infor[GTEXmedian_tpm_infor$Gene.ID %in%G4_alone, ],"CC_density","GG_density","Template","Non-template","norm_exp","Correlation of expression and frequency\nof G-runs in different strand in TSS downstream","spearman")+theme(legend.position = "bottom",legend.box ="vertical")+ylab("GTEX Cohort")+labs(fill="Strand")+xlab("spearman's rho (expression vs G-run frequency)")+ggtitle("Genes only contains G4")

# Plot Venn diagram
p_venn <- ggVennDiagram(x,
                        label_alpha = 0,  # Remove grey fill for numbers
                        label_color = "black", label_size = 4, edge_size = 0.8,
                        set_color = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF")
) +
  scale_fill_distiller(palette = "Set3", direction = 1) +  # Soft color palette
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  labs(title = "Overlap of non-B DNA structures")

# Load necessary libraries for boxplot
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

# Prepare data for boxplot
plot_dat <- genes_infor %>%
  dplyr::select(
    ensembl_gene_id,  # Already named ensembl_gene_id
    GG_density,
    Z_pos, MR_pos, DR_pos, STR_pos, G4_down_predict_pos
  ) %>%
  pivot_longer(
    cols = -c(ensembl_gene_id, GG_density),
    names_to = "Feature",
    values_to = "Pos_flag"
  ) %>%
  mutate(
    Feature = recode_factor(
      Feature,
      Z_pos = "Z-DNA",
      MR_pos = "Mirror_repeats_and_triplex",
      DR_pos = "Direct_repeats",
      STR_pos = "Short_tandem_repeats",
      G4_down_predict_pos = "G4"
    ),
    Pos_flag = factor(Pos_flag, levels = c("-", "+"))
  )

# Define colors for boxplot
my_cols <- colorRampPalette(c("white", "#40A4DE"))(6)[c(2, 6)]

# Plot boxplot for structural features
p_box <- ggplot(plot_dat,
                aes(x = Pos_flag, y = GG_density, fill = Pos_flag)) +
  geom_boxplot(alpha = 0.7, outlier.colour = NA, width = 0.6) +
  geom_signif(
    comparisons = list(c("+", "-")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    col = "black",
    step_increase = 0.05,
    tip_length = 0.01
  ) +
  facet_wrap(~ Feature, ncol = 1) +
  scale_fill_manual(values = my_cols) +
  labs(
    x = "Non-B DNA structure",
    y = "Frequency of G-runs",
    title = "Non-templated strand of TSS downstream"
  ) +
  theme_hd_minimal() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1.5, "lines")
  ) + coord_flip()


g4_ids <- genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_pos == "+"]
zdna_ids <- genes_infor$ensembl_gene_id[genes_infor$Z_pos == "+"]
dr_ids <- genes_infor$ensembl_gene_id[genes_infor$DR_pos == "+"]
mr_ids <- genes_infor$ensembl_gene_id[genes_infor$MR_pos == "+"]
str_ids <- genes_infor$ensembl_gene_id[genes_infor$STR_pos == "+"]

other_ids <- unique(c(zdna_ids, dr_ids, mr_ids, str_ids))


G4_alone <- setdiff(g4_ids, other_ids)

ZDNA_alone <- setdiff(zdna_ids, c(g4_ids, dr_ids, mr_ids, str_ids))

DR_alone <- setdiff(dr_ids, c(g4_ids, zdna_ids, mr_ids, str_ids))

MR_alone <- setdiff(mr_ids, c(g4_ids, zdna_ids, dr_ids, str_ids))

STR_alone <- setdiff(str_ids, c(g4_ids, zdna_ids, dr_ids, mr_ids))

# 2. 创建 gene_tbl 数据框
gene_tbl <- data.frame(
  ensembl_gene_id = c(G4_alone, ZDNA_alone, DR_alone, MR_alone, STR_alone),
  category = factor(rep(c("G4", "Z-DNA", "Direct Repeats", "Mirror Repeats", "Short Tandem Repeats"),
                        times = c(length(G4_alone), length(ZDNA_alone), length(DR_alone), length(MR_alone), length(STR_alone))),
                    levels = c( "Direct Repeats", "G4","Mirror Repeats", "Short Tandem Repeats", "Z-DNA"))
)


expr <- merge(all_data_tsg[, c("ensembl_gene_id", "cancer", "norm_exp")],
              gene_tbl, by = "ensembl_gene_id")


p_exp_non_B <- ggplot(expr, aes(x = category, y = norm_exp, fill = cancer)) +
  geom_boxplot() +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(title = "Expression of genes containing only one type of non-B DNA structure",
       x = "Non-B DNA structure",
       y = "Expression level",
       fill = "Cancer Type") +
  theme_hd()+
  theme(legend.position = "bottom")+coord_flip()

print(p_exp_non_B)





# Arrange plots for Supplementary Figure 16
pall_G41 <- ggarrange(p_box + theme(legend.position = "none"), p_venn, nrow = 1, ncol = 2, 
                      labels = c(letters[1:2]), 
                      font.label = list(size = 18, color = "black", face = "bold"), 
                      widths = c(1, 1))

pall_G42 <- ggarrange(p_exp_non_B, G4_cor, nrow = 1, ncol = 2, 
                      labels = c(letters[4:5]), 
                      font.label = list(size = 18, color = "black", face = "bold"), 
                      widths = c(1, 1))

pall_G43 <- ggarrange(pall_G41, p_density, pall_G42, nrow = 3, ncol = 1, 
                      labels = c("", letters[3], ""), 
                      font.label = list(size = 18, color = "black", face = "bold"), 
                      heights = c(1.5, 1, 2.7))

# Save Supplementary Figure 16 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure16.tiff", 
     width = 18, height = 20, res = 300, units = "in", compression = "lzw")
print(pall_G43)
dev.off()






temp_intron1<-getIntron(all_data_genes$get_gene_pos_exon,genes_infor,1)
gene_pos_data_exon1<-unique(getGenePosData(temp_intron1,c("ensembl_transcript_id","chromosome_name","exon_chrom_start","exon_chrom_end","strand","ensembl_gene_id")))
gene_pos_data_intron1<-unique(getGenePosData(temp_intron1,c("ensembl_transcript_id","chromosome_name","intron_start","intron_end","strand","ensembl_gene_id")))
gene_pos_data_5utr<-unique(getGenePosData(all_data_genes$get_gene_pos[all_data_genes$get_gene_pos$ensembl_transcript_id %in% genes_infor$ensembl_transcript_id,],c("ensembl_transcript_id","chromosome_name","5_utr_start","5_utr_end","strand","ensembl_gene_id")))

getRegionInfor<-function(gene_pos_data_exon1,gene_pos_data_5utr,gene_pos_data_intron1, down_gr,labs_x){
  hits <- findOverlaps(gene_pos_data_intron1, down_gr)
  intron_in_1kb <- pintersect(gene_pos_data_intron1[queryHits(hits)],
                              down_gr[subjectHits(hits)])
  
  hits <- findOverlaps(gene_pos_data_exon1, down_gr)
  exon_in_1kb <- pintersect(gene_pos_data_exon1[queryHits(hits)],
                            down_gr[subjectHits(hits)])
  
  hits <- findOverlaps(gene_pos_data_5utr, down_gr)
  utr5_in_1kb <- pintersect(gene_pos_data_5utr[queryHits(hits)],
                            down_gr[subjectHits(hits)])
  pos_intron1<- getFreqDefinedRegion(intron_in_1kb,Hsapiens)
  pos_exon1<- getFreqDefinedRegion(exon_in_1kb,Hsapiens)
  pos_utr5<- getFreqDefinedRegion(utr5_in_1kb,Hsapiens)
  pos_down<- getFreqDefinedRegion(down_gr,Hsapiens)
  
  
  
  title_x=""
  if(labs_x=="TSS_down_1kb"){title_x="Regions in TSS downstream 1kb"}
  if(labs_x=="TSS_down_300bp"){title_x="Regions in TSS downstream 300bp"}
 
  df_all <- bind_rows(
    tibble(GG_density = pos_utr5$GG_density,   Region = "5'UTR"),
    tibble(GG_density = pos_exon1$GG_density,  Region = "First exon"),
    tibble(GG_density = pos_intron1$GG_density,Region = "First intron"),
    tibble(GG_density = pos_down$GG_density,Region = labs_x)
  ) %>% mutate(Region = factor(Region,levels = c(labs_x,"First intron","5'UTR",
                                                 "First exon")))
  pair_list <- combn(levels(df_all$Region), 2, simplify = FALSE)
  p=ggplot(df_all, aes(x = Region, y = GG_density, fill = Region)) +
    geom_boxplot(outlier.colour = "grey90",outlier.shape  = 19,
                 outlier.size   = 1.2,col="grey90") +
    geom_signif(comparisons = pair_list, test        = "wilcox.test",
                map_signif_level = TRUE,col         = "black",
                step_increase = 0.1,y_position = 0.2) +
    scale_fill_manual(values = c("#1B9E77","#D95F02","#7570B3","#E6AB02")) +
    theme_classic(base_size = 14) +theme_hd()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x =title_x, y = "G-runs Frequency")
  return(list(p=p,intron=intron_in_1kb,exon=exon_in_1kb,utr5=utr5_in_1kb,df_all=df_all))
}

gene_pos<-genes_infor[,c("genes_main_transcript","chromosome_name","transcript_start","transcript_end","strand","ensembl_gene_id")]
colnames(gene_pos)<-c("trans_ID","seqnames","start","end","strand","gene_id")
gene_pos$strand[gene_pos$strand=="1"]="+"
gene_pos$strand[gene_pos$strand=="-1"]="-"
gene_pos$seqnames<-paste("chr",gene_pos$seqnames,sep="")
gene_pos_data<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=T)

down1000_gr <- get_gene_region(gene_pos_data,0,1000)
down300_gr <- get_gene_region(gene_pos_data,0,300)
down300_1000_gr <- get_gene_region(gene_pos_data,300,1000)
pdown1000<-getRegionInfor(gene_pos_data_exon1,gene_pos_data_5utr,gene_pos_data_intron1, down1000_gr,"TSS_down_1kb")
pdown300<-getRegionInfor(gene_pos_data_exon1,gene_pos_data_5utr,gene_pos_data_intron1, down300_gr,"TSS_down_300bp")
pdown300_1000<-getRegionInfor(gene_pos_data_exon1,gene_pos_data_5utr,gene_pos_data_intron1, down300_1000_gr,"TSS_down_300bp_1kb")

all_data_down<-rbind(cbind(type="Down 300bp",pdown300$df_all),cbind(type="Down 300bp_1kb",pdown300_1000$df_all))
all_data_down$Region<-as.character(all_data_down$Region)
all_data_down$Region[all_data_down$Region %in% c("TSS_down_300bp","TSS_down_300bp_1kb")]="All"

library(ggpubr)   # 已含 ggplot2

# 1. 计算每组中位数
median_df <- aggregate(GG_density ~ Region + type,
                       data = all_data_down,
                       FUN  = median)   # 也可加 IQR 同理
all_data_down$Region<-factor(all_data_down$Region,levels=c(unique(all_data_down$Region)))
# 2. 绘图 + 统计检验
intron_bar1<- ggbarplot(all_data_down,
                        x = "type", y = "GG_density",
                        fill = "type",
                        facet.by = "Region",
                        nrow = 1,
                        add = "mean_sd") +            # 如要 IQR 可改为 "median_iqr"
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     step.increase = 0.05) +
  scale_fill_manual(values = c("Down 300bp"      = "#1B9E77",
                               "Down 300bp_1kb" = "#D95F02")) +
  theme_hd_minimal() +theme(axis.text.x = element_text(angle=45,hjust=1))+
  labs(x = NULL, y = "G-run frequency in non-template strand")+ylim(0,0.13)


# 2. 生成所有两两比较对
predict_G4_cds_downstream<-readRDS(file="./data/predict_G4_cds_downstream.rds")
temp_pos=predict_G4_cds_downstream$pos
df <- as.data.frame(temp_pos)

names(temp_pos)[c(2,3,8,9)] <- c("winStart","winEnd","g4_txStart","g4_txEnd")
temp_pos300<-temp_pos[temp_pos$g4_txStart<=300,]
temp_pos300_1000<-temp_pos[temp_pos$g4_txStart>300,]




getG4location<-function(pdown,temp_pos){
  ## 1. 计算染色体绝对坐标
  isPos <- temp_pos$strand == "+"
  temp_pos$chromStart <- ifelse(isPos,
                                temp_pos$winStart + temp_pos$g4_txStart - 1,
                                temp_pos$winEnd   - temp_pos$g4_txEnd   + 1)
  
  temp_pos$chromEnd   <- ifelse(isPos,
                                temp_pos$winStart + temp_pos$g4_txEnd - 1,
                                temp_pos$winEnd   - temp_pos$g4_txStart + 1)
  gr <- GRanges(
    seqnames = temp_pos$seqnames,
    ranges   = IRanges(start = temp_pos$chromStart,
                       end   = temp_pos$chromEnd),
    strand   = temp_pos$strand,
    trans_ID = temp_pos$trans_ID,
    gene_id  = temp_pos$gene_id,
    score    = temp_pos$score
  )
  pdown$utr5$type  <- "UTR5"
  pdown$intron$type <- "First Intron"
  pdown$exon$type   <- "First Exon"
  feats <- c(pdown$utr5, pdown$intron, pdown$exon)
  
  ## 3. 现在可以正常 overlap
  hits <- findOverlaps(gr, feats, ignore.strand = TRUE)
  
  ## 4. 统计
  feat_type <- feats$type[subjectHits(hits)]
  tbl  <- table(feat_type)
  prop <- round(tbl / length(gr) * 100, 2)
  pie_df <- data.frame(
    Feature = names(tbl),
    Count   = as.numeric(tbl),
    Percent = as.numeric(prop)
  )
  return(pie_df)
}

g4_down300<-getG4location(pdown300,temp_pos300)
g4_down300_1000<-getG4location(pdown300_1000,temp_pos300_1000)
g4_down1000<-getG4location(pdown1000,temp_pos)

## 准备数据
g4_down1000$group <- "TSS Downstream\n1kb"
g4_down300_1000$group <- "TSS Downstream\n300bp-1kb"
g4_down300$group <- "TSS Downstream\n300bp"

df <- rbind(g4_down1000, g4_down300_1000, g4_down300)
df$Feature<-factor(df$Feature,levels = c("UTR5","First Exon","First Intron"))
## 作图
library(ggplot2)
intron_bar<-ggplot(df, aes(x = Feature, y = Percent, fill = Feature)) +
  geom_col(width = 0.7, color = "white", size = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "black") +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = NULL, y = "Percentage of predicted G4\nin non-template strand") +
  theme_hd_minimal() +
  theme(legend.position = "none")+theme(axis.text.x = element_text(angle=45,hjust=1))





all_genes_reverse_complementary_imotif <- read_csv("data/all_genes_reverse_complementary_imotif.csv")

genes_infor$imotif<-ifelse(genes_infor$ensembl_transcript_id %in% str_split_fixed(all_genes_reverse_complementary_imotif$chr,"_",2)[,1],"+","-")

x=list(Predicted_G4_NT=genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_pos=="+"],imotif=genes_infor$ensembl_gene_id[genes_infor$imotif=="+"])
G4_aloneimotif<-x$Predicted_G4_NT[!x$Predicted_G4_NT %in% x$imotif]
p_cor_imotif<-PlotMutipleCor(GTEXmedian_tpm_infor[GTEXmedian_tpm_infor$Gene.ID %in% G4_aloneimotif,],"CC_density","GG_density","Template","Non-template","norm_exp","Correlation of expression and \n G-run frequency in different strand in TSS downstream","spearman")+theme(legend.position = "bottom",legend.box ="vertical")+ylab("GTEX Cohort")+labs(fill="Strand")+xlab("Spearman's rho (expression vs G-run frequency")+ggtitle("Genes only contains G4(without imotif)")



title_x="Overlap in TSS downstream"

p=venn.diagram(x,fill = c( "#00B2EE","#104E8B"),
               alpha = c(0.5, 0.5), cex = 1,cat.fontface = "bold",lty =2, fontfamily =3,cat.cex=1,sub.fontface="bold",margin=0.1,sub="", sub.pos = c(0.3,-0.05),filename = NULL,main=title_x,main.fontface ="bold",  main.cex = 1.2,col=NA,	ext.text = TRUE,ext.line.lwd = 2, ext.dist = -0.15,ext.length = 0.9,ext.pos = -14,rotation.degree = 185,cat.pos=c(120,240))
pVNN_imotif<-ggarrange(p)
pVNN_imotif

palls0<-ggarrange(pVNN_imotif,intron_bar1,intron_bar,nrow=3,ncol=1,labels=c(letters[c(1,3,4)]),font.label = list(size = 18, color = "black", face = "bold"))

palls1<-ggarrange(palls0,p_cor_imotif, nrow=1,ncol=2,labels=c("",letters[c(2)]),font.label = list(size = 18, color = "black", face = "bold"))

tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure17.tiff",width = 17,height = 17,res=300,units="in",compression = "lzw")
print(palls1)
dev.off()