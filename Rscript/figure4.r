# Set the working directory to the specified path
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new")

# Load necessary data and functions
genes_infor <- readRDS("./data/genes_infor_downloaded_from_ensembl105_gene_infor.rds")  # Load gene information
all_data_genes <- readRDS("./data/genes_infor_downloaded_from_ensembl105_all_data_genes.rds")  # Load all gene data
load(file="./data/other_species_all_data.rdata")  # Load data for other species
source('functions/functions_evolution_5UTR.r')  # Load custom functions for 5' UTR evolution analysis
load(file="./data/result_promoter_window_ACTG_1000_1000.rdata")  # Load promoter window data
load(file="./data/all_meth_exp_signal_promoter.rdata")  # Load methylation and expression signal data
load(file="./data/predict_G4_up_down.rdata")  # Load G4 prediction data

# Filter out oncogenes from the gene information
genes_infor <- genes_infor[genes_infor$Gene.Type != "Oncogene", ]

# Create a new column for age type with combined categories
genes_infor$age_type3 <- as.character(genes_infor$age_type2)
genes_infor$age_type3[genes_infor$age_type3 %in% c("(67.6 - 355.7]", "(0 - 67.6]")] <- "(0 - 355.7]"

# Extract relevant columns and rename them for consistency
gene_pos <- genes_infor[, c("genes_main_transcript", "chromosome_name", "transcript_start", "transcript_end", "strand", "ensembl_gene_id")]
colnames(gene_pos) <- c("trans_ID", "seqnames", "start", "end", "strand", "gene_id")

# Convert strand information to standard format
gene_pos$strand[gene_pos$strand == "1"] <- "+"
gene_pos$strand[gene_pos$strand == "-1"] <- "-"

# Prefix chromosome names with 'chr'
gene_pos$seqnames <- paste("chr", gene_pos$seqnames, sep = "")

# Create a GRanges object from the gene position data
gene_pos_data <- makeGRangesFromDataFrame(gene_pos, keep.extra.columns = TRUE)

# Calculate G-run frequencies for upstream and downstream regions
pos_all_up <- getFreq(gene_pos_data, -1000, 0, Hsapiens, "+")  # Upstream region
pos_all_down <- getFreq(gene_pos_data, 0, 1000, Hsapiens, "+")  # Downstream region

# Add downstream region frequencies to the gene information data frame
genes_infor[, colnames(pos_all_down[, 6:11])] <- pos_all_down[genes_infor$ensembl_gene_id, 6:11]

# Calculate G-run frequencies for a shorter downstream region (300 bp)
pos_all_down <- getFreq(gene_pos_data, 0, 300, Hsapiens, "+")
genes_infor[, paste("down300", colnames(pos_all_down[, c(2, 6:11)]), sep = "_")] <- pos_all_down[genes_infor$ensembl_gene_id, c(2, 6:11)]

# Add specific G-run density and CpG OE values for the 300 bp downstream region
genes_infor$down300_AA_density <- pos_all_down[, "AA_density"]
genes_infor$down300_GG_density <- pos_all_down[, "GG_density"]
genes_infor$down300_CpG_OE <- pos_all_down[, "CpG_OE"]

# Convert Gene.Type to a factor for plotting
genes_infor$Gene.Type2 <- as.factor(genes_infor$Gene.Type)

# Identify paralog genes for TSGs
gene_paralog <- all_data_genes$gene_paralog
gene_paralog <- gene_paralog[gene_paralog$hsapiens_paralog_ensembl_gene != "", ]
TSG_genes <- genes_infor[genes_infor$Gene.Type == "TSG", ]$ensembl_gene_id
TSG_paralog_genes <- gene_paralog[gene_paralog$ensembl_gene_id %in% TSG_genes & 
                                    gene_paralog$hsapiens_paralog_ensembl_gene %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type == "Non-Cancer"], ]
TSG_paralog_non_cancer <- genes_infor[genes_infor$ensembl_gene_id %in% unique(as.character(c(TSG_genes, TSG_paralog_genes$hsapiens_paralog_ensembl_gene))), ]

# Update Gene.Type for non-cancer paralogs
TSG_paralog_non_cancer$Gene.Type[TSG_paralog_non_cancer$Gene.Type == "Non-Cancer"] <- "Non_Cancer_paralogs"

# Add a column to indicate paralog status
genes_infor$Gene_paralogs <- rep("Non_paralog", length(genes_infor$ensembl_gene_id))
genes_infor$Gene_paralogs[genes_infor$ensembl_gene_id %in% unique(as.character(c(TSG_genes, TSG_paralog_genes$hsapiens_paralog_ensembl_gene)))] <- "paralog"

# Plot G-run density for TSGs and non-cancer genes
p11 <- plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type != "Oncogene", ], "GG_density", "G-run frequency\nin NT strand of TSS downstream 1kb")

# Plot G-run density for paralog genes
p21 <- plot_non_paired_paralog(TSG_paralog_non_cancer, "GG_density", "G-run frequency\nin NT strand of TSS downstream 1kb")

# Plot G-run density for the 300 bp downstream region
p11down300 <- plot_TSG_non_cancer_withlegend(genes_infor[genes_infor$Gene.Type != "Oncogene", ], "down300_GG_density", "G-run frequency\nin NT strand of TSS downstream 300bp")
p21down300 <- plot_non_paired_paralog(TSG_paralog_non_cancer, "down300_GG_density", "G-run frequency\nin NT strand of TSS downstream 300bp") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Identify G4 motifs in the 300 bp downstream region
genes_infor$down300_pre_pos <- ifelse(genes_infor$ensembl_transcript_id %in% promoter_G4_down_pos$trans_ID[promoter_G4_down_pos[, 9] <= 300], "+", "-")
genes_infor$down300_pre_neg <- ifelse(genes_infor$ensembl_transcript_id %in% promoter_G4_down_neg$trans_ID[promoter_G4_down_neg[, 9] <= 300], "+", "-")

# Load G4 prediction results for the 300 bp downstream region
temp_down300 <- read.delim("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/result/human_5UTR300plus.hits.max.PDS.w50.35.bed", header = FALSE, stringsAsFactors = FALSE)

# Identify G4 motifs in the 300 bp downstream region
genes_infor$down300G4_down_exp <- ifelse(genes_infor$genes_main_transcript %in% temp_down300[, 4], "+", "-")
genes_infor$pre_type <- ifelse(genes_infor$ensembl_transcript_id %in% promoter_G4_down_pos$trans_ID[promoter_G4_down_pos[, 9] <= 300], "down300bp", ifelse(genes_infor$G4_down_predict_pos == "+", "down300-1kb", "-"))

# Merge GTEX expression data with gene information
GTEXmedian_tpm_infor <- merge(GTEXmedian_tpm, genes_infor, by.x = "Gene.ID", by.y = "ensembl_gene_id")

# Add pre_type and G-run density information to chromatin data
genes_chromatin$pre_type <- genes_infor[genes_chromatin$ensembl_gene_id, "pre_type"]
genes_chromatin$down300_pre_pos <- genes_infor[genes_chromatin$ensembl_gene_id, "down300_pre_pos"]
genes_chromatin$down300_GG_density <- genes_infor[genes_chromatin$ensembl_gene_id, "down300_GG_density"]
genes_chromatin$down300_CC_density <- genes_infor[genes_chromatin$ensembl_gene_id, "down300_CC_density"]

# Add pre_type and G-run density information to methylation data
aver_meth_gene$pre_type <- genes_infor[aver_meth_gene$ensembl_gene_id, "pre_type"]
aver_meth_gene$down300_pre_pos <- genes_infor[aver_meth_gene$ensembl_gene_id, "down300_pre_pos"]
aver_meth_gene$down300_GG_density <- genes_infor[aver_meth_gene$ensembl_gene_id, "down300_GG_density"]
aver_meth_gene$down300_CC_density <- genes_infor[aver_meth_gene$ensembl_gene_id, "down300_CC_density"]

# Plot promoter DHS signal for TSGs in normal tissue
p1_signal_g4all <- plot_classmulti(genes_chromatin[genes_chromatin$cancer == "colon" & genes_chromatin$Gene.Type == "TSG", ], "pre_type", "signal", "", "Promoter DHS signal of TSGs in normal tissue") +
  ggtitle("ENCODE colon") +
  scale_fill_manual(name = "Position of Predicted G4\n in TSS downstream", values = c("#D6D6D6", "#EEB422", "#CD6600")) +
  theme(axis.text.x = element_blank())

# Plot expression of TSGs in normal tissue
p1_exp_g4all <- plot_classmulti(GTEXmedian_tpm_infor[GTEXmedian_tpm_infor$cancer == "Colon.Sigmoid" & GTEXmedian_tpm_infor$Gene.Type.x == "TSG", ], "pre_type", "norm_exp", "", "Expression of TSGs in normal tissue") +
  ggtitle("GTEX-sigmoid colon") +
  scale_fill_manual(name = "Position of Predicted G4\n in TSS downstream", values = c("#D6D6D6", "#EEB422", "#CD6600")) +
  theme(axis.text.x = element_blank())

# Plot promoter methylation of TSGs in normal tissue
p1_meth_g4all <- plot_classmulti(aver_meth_gene[aver_meth_gene$cancer == "COAD" & aver_meth_gene$Gene.Type == "TSG", ], "pre_type", "norm_meth", "", "Promoter Methylation of TSGs in normal tissue") +
  ggtitle("TCGA-COAD") +
  scale_fill_manual(name = "Position of Predicted G4\n in TSS downstream", values = c("#D6D6D6", "#EEB422", "#CD6600")) +
  theme(axis.text.x = element_blank())

# Calculate and plot correlations for methylation, expression, and DHS signal
aver_meth_gene$down300_GG_density <- genes_infor[aver_meth_gene$ensembl_gene_id, "down300_GG_density"]
aver_meth_gene$down300_CC_density <- genes_infor[aver_meth_gene$ensembl_gene_id, "down300_CC_density"]

all_data_tsg$down300_GG_density <- genes_infor[all_data_tsg$ensembl_gene_id, "down300_GG_density"]
all_data_tsg$down300_CC_density <- genes_infor[all_data_tsg$ensembl_gene_id, "down300_CC_density"]

p_meth <- plotCorbox(aver_meth_gene[aver_meth_gene$Gene.Type == "TSG", ], "down300_CC_density", "down300_GG_density", "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", "Template strand", "Non-template strand", "norm_meth", "Correlation of promoter methylation \nand G-run frequency in TSGs(TCGA cohort)", "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

p_exp <- plotCorbox(GTEXmedian_tpm_infor[GTEXmedian_tpm_infor$Gene.Type.x == "TSG", ], "down300_CC_density", "down300_GG_density", "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", "Template strand", "Non-template strand", "norm_exp", "Correlation of expression and\nG-run frequency in TSGs(GTEX cohort)", "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

p_exp2 <- plotCorbox(all_data_tsg[all_data_tsg$Gene.Type == "TSG", ], "down300_CC_density", "down300_GG_density", "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", "Template strand", "Non-template strand", "norm_exp", "Correlation of expression and\nG-run frequency in TSGs(GTEX cohort)", "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

p_signal <- plotCorbox(genes_chromatin[genes_chromatin$Gene.Type == "TSG", ], "down300_CC_density", "down300_GG_density", "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", "Template strand", "Non-template strand", "signal", "Correlation of promoter DHS signal and\nG-run frequency in TSGs(ENCODE cohort)", "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

# Aggregate correlation data
temp_exp <- aggregate(p_exp$data$cor, list(type = p_exp$data$type), mean)
temp_meth <- aggregate(p_meth$data$cor, list(type = p_meth$data$type), mean)
temp_signal <- aggregate(p_signal$data$cor, list(type = p_signal$data$type), mean)

# Combine all correlation data into a single data frame
all_cor <- rbind(
  cbind(type_x = "Expression vs G-runs freqency in TSGs", temp_exp),
  cbind(type_x = "Methylation vs G-runs freqency in TSGs", temp_meth),
  cbind(type_x = "DHS signal vs G-runs freqency in TSGs", temp_signal)
)


heat_long <- all_cor %>%
  mutate(
    strand = str_extract(type, "(Template|Non-template)"),
    region = str_extract(type, "TSS.*")
  ) %>%
  mutate(
    strand = factor(strand, levels = c("Template", "Non-template")),
    region = factor(region, levels = c("TSS down_300bp", "TSS down_1kb")),
    type_x = factor(type_x, levels = c("Expression vs G-runs freqency in TSGs", "Methylation vs G-runs freqency in TSGs", "DHS signal vs G-runs freqency in TSGs"))
  ) %>%
  mutate(
    alpha = ifelse(strand == "Non-template" & region == "TSS down_300bp", 1, 0.7)
  )
p_heat <- ggplot(heat_long,
            aes(region, strand, fill = x, alpha = alpha)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "#3B4992",
                       mid = "white",
                       high = "#E64B35",
                       midpoint = 0,
                       limits = c(-0.25, 0.25),
                       name = "Spearman ρ") +
  scale_alpha_identity() +  # 使用自定义的 alpha 值
  geom_text(aes(label = sprintf("%.2f", x)),
            size = 4) +
  facet_wrap(~ type_x, ncol = 1) +
  labs(x = "Distance from TSS", y = "Strand")+ggtitle("Mean correlations: Molecular feartures vs G-run frequency") +
  theme_hd_minimal2()+theme(strip.text = element_text(size=14,face = "bold"),
                            legend.text = element_text(angle = 90, vjust = 0.5)  )# 图例文字垂直)

print(p_heat )


p11paralog=getPairedPlotMultiHK(TSG_paralog_non_cancer,"CC_density","GG_density","Template ","Non-template","Gene.Type")+ylab("G-run frequency\nin TSS downstream 1kb")

pall51<-ggarrange(p11paralog,p21down300,nrow=1,ncol=2,labels=c(letters[1:2]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1))
pall52<-ggarrange(p1_exp_g4all,p1_meth_g4all,p1_signal_g4all,nrow=1,ncol=3,labels=c(letters[4:6]),font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend = "bottom")
pall53<-ggarrange(pall51,pall52,nrow=2,ncol=1,widths = c(1.5,1))
pall54<-ggarrange(pall53,p_heat ,nrow=1,ncol=2)

tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/all_figure4.tiff",width = 15,height = 15,res=300,units="in",compression = "lzw")
print(pall54)
dev.off()







# Filter out Drosophila melanogaster from the dataset
all_data <- all_data[!all_data$species %in% c("Drosophila.melanogaster"), ]

# Convert species to a factor, excluding the 9th level (Drosophila melanogaster)
all_data$species <- factor(all_data$species, levels = levels(all_data$species)[-9])

# Plot species-specific data for CC_GG
p3 <- PlotSpecies(all_data, "CC_GG")

# Plot species-specific data for GG_density
p3species <- PlotSpecies(all_data, "GG_density") + ylab("G-run frequency in \nNT strand of TSS downstream 1kb")

# Plot gene age class data for GG_density
p4age <- plotSubSigOneSide(genes_infor[!is.na(genes_infor$age_type3) & genes_infor$Gene.Type != "Oncogene", ], 
                           "age_type2", "GG_density", "Gene.Type", "", 
                           cols1 = c("#1874CD", "#EE2C2C"), 
                           "Gene Age class\n(million years)", 
                           "G-run frequency\nin NT strand of TSS downstream 1kb", 
                           side = "two.sided") + 
  labs(fill = "Gene.Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Arrange plots for Supplementary Figure 8
pall5_1 <- ggarrange(p11, p4age, nrow = 1, ncol = 2, 
                     labels = c(letters[1:2]), 
                     font.label = list(size = 18, color = "black", face = "bold"))

pall5_2 <- ggarrange(p3species, p11down300, nrow = 1, ncol = 2, 
                     labels = c(letters[3:4]), 
                     font.label = list(size = 18, color = "black", face = "bold"), 
                     common.legend = TRUE, legend = "bottom", widths = c(2, 1, 1))

# Save Supplementary Figure 8 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure8.tiff", 
     width = 13, height = 18, res = 300, units = "in", compression = "lzw")
print(ggarrange(pall5_1, pall5_2, nrow = 2, common.legend = TRUE, legend = "bottom"))
dev.off()

p1_exp_g4allcancer=plot_classmulti(GTEXmedian_tpm_infor[GTEXmedian_tpm_infor$cancer!="Colon.Sigmoid" &GTEXmedian_tpm_infor$Gene.Type.x=="TSG",],"pre_type","norm_exp","Position of predicted G4 in NT strand of TSS downstream","Expression of TSGs")+scale_fill_manual(values = c("#D6D6D6", "#EEB422", "#CD6600"))+facet_wrap(~cancer)+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 30,hjust = 1),legend.position = "none")

tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure9.tiff",width = 25,height = 27,res=300,units="in",compression = "lzw")
ggarrange(p1_exp_g4allcancer)
dev.off()


p1_meth_g4allcancer=plot_classmulti(aver_meth_gene[aver_meth_gene$cancer!="COAD"& aver_meth_gene$Gene.Type=="TSG",],"pre_type","norm_meth","Position of predicted G4\nin NT strand of TSS downstream","Promoter Methylation of TSGs")+scale_fill_manual(values = c("#D6D6D6", "#EEB422", "#CD6600"))+facet_wrap(~cancer)+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 30,hjust = 1),legend.position = "none")

tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure10.tiff",width = 10,height = 10,res=300,units="in",compression = "lzw")
ggarrange(p1_meth_g4allcancer)
dev.off()

p1_signal_g4allcancer=plot_classmulti(genes_chromatin[genes_chromatin$cancer!="colon" &genes_chromatin$Gene.Type=="TSG",],"pre_type","signal","Position of predicted G4\nin NT strand of TSS downstream","Promoter DHS signal of TSGs")+scale_fill_manual(values = c("#D6D6D6", "#EEB422", "#CD6600"))+facet_wrap(~cancer)+theme_hd_minimal()+theme(axis.text.x = element_text(angle = 30,hjust = 1),legend.position = "none")
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure11.tiff",width = 10,height = 10,res=300,units="in",compression = "lzw")
print(p1_signal_g4allcancer)
dev.off()

# 定义函数
process_and_export_data <- function(data, output_file = "supplementary_table2.txt") {
  data$Pvalue<-round(data$Pvalue,3)
  
  cor_data <- dcast(data, Cancer ~ type, value.var = "cor")
  pvalue_data <- dcast(data, Cancer ~ type, value.var = "Pvalue")
  names(cor_data)[-1] <- gsub("\n", " ", names(cor_data)[-1])
  names(cor_data)[-1] <- paste(names(cor_data)[-1], "Correlation", sep = " ")
  names(pvalue_data)[-1] <- gsub("\n", " ", names(pvalue_data)[-1])
  names(pvalue_data)[-1] <- paste(names(pvalue_data)[-1], "P-value", sep = " ")
  wide_data <- merge(cor_data, pvalue_data, by = "Cancer")
  
  names(wide_data)[1]<-"Tissue"
  print(wide_data[1:3,1:4])
  write.csv(wide_data, output_file,quote = F,row.names = T)
  
}

# 使用函数处理数据

process_and_export_data(p_exp$data, "supplementary_table2_exp.csv")
process_and_export_data(p_exp2$data, "supplementary_table2_exp2.csv")
process_and_export_data(p_meth$data, "supplementary_table2_meth.csv")
process_and_export_data(p_signal$data, "supplementary_table2_signal.csv")




