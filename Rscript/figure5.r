# Set the working directory to the specified path
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new")

# Load necessary functions and data
source('functions/functions_evolution_5UTR.r')  # Load custom functions for 5' UTR evolution analysis
load(file="./data/all_meth_exp_signal_promoter.rdata")  # Load methylation and expression signal data

# Prepare gene position data
gene_pos <- genes_infor[, c("genes_main_transcript", "chromosome_name", "transcript_start", "transcript_end", "strand", "ensembl_gene_id")]
colnames(gene_pos) <- c("trans_ID", "seqnames", "start", "end", "strand", "gene_id")

# Convert strand information to standard format
gene_pos$strand[gene_pos$strand == "1"] <- "+"
gene_pos$strand[gene_pos$strand == "-1"] <- "-"

# Prefix chromosome names with 'chr'
gene_pos$seqnames <- paste("chr", gene_pos$seqnames, sep = "")

# Create a GRanges object from the gene position data
gene_pos_data <- makeGRangesFromDataFrame(gene_pos, keep.extra.columns = TRUE)

# Calculate G-run frequencies for the downstream region (300 bp)
pos_all_down <- getFreq(gene_pos_data, 0, 300, Hsapiens, "+")

# Add downstream region frequencies to the gene information data frame
genes_infor[, paste("down300", colnames(pos_all_down[, c(2, 6:11)]), sep = "_")] <- pos_all_down[genes_infor$ensembl_gene_id, c(2, 6:11)]

# Add specific G-run density values for the 300 bp downstream region
all_data_tsg$down300_GG_density <- genes_infor[all_data_tsg$ensembl_gene_id, "down300_GG_density"]
all_data_tsg$down300_CC_density <- genes_infor[all_data_tsg$ensembl_gene_id, "down300_CC_density"]

# Filter low data for TSGs
all_data_tsg_filtered <- filter_low_data(all_data_tsg)
all_data_tsg_filtered$diff_exp <- abs(all_data_tsg_filtered$diff_exp)

# Add methylation data for the 300 bp downstream region
aver_meth_gene$down300_GG_density <- genes_infor[aver_meth_gene$ensembl_gene_id, "down300_GG_density"]
aver_meth_gene$down300_CC_density <- genes_infor[aver_meth_gene$ensembl_gene_id, "down300_CC_density"]

# Load rG4-seq data for the 300 bp downstream region
temp_down300 <- read.delim("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/result/human_down300_rG4_h38.bed", header = FALSE, stringsAsFactors = FALSE)
genes_infor$rG4_down300_exp <- ifelse(genes_infor$genes_main_transcript %in% temp_down300[, 4], "+", "-")
all_data_tsg_filtered$rG4_down300_exp <- genes_infor[all_data_tsg_filtered$ensembl_gene_id, "rG4_down300_exp"]

# Plot correlation of expression downregulation extent and G-run frequency
p_diff <- plotCorbox(all_data_tsg_filtered, "down300_CC_density", "down300_GG_density", "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", "Template strand", "Non-template strand", "diff_exp", "Correlation of expression downregulation extent\nand G-run frequency in Down-TSGs", "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

# Plot correlation of methylation downregulation extent and G-run frequency
p_diff_meth <- plotCorbox(aver_meth_gene[which(aver_meth_gene$ensembl_gene_id %in% all_data_tsg_filtered$ensembl_gene_id & aver_meth_gene$diff_meth > 0), ], "down300_CC_density", "down300_GG_density", "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", "Template strand", "Non-template strand", "diff_meth", "Correlation of promoter methylation downregulation extent\nand G-run frequency in Down-TSGs", "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

# Plot multiple correlations for methylation
meth_diff <- PlotMutipleCor(aver_meth_gene[which(aver_meth_gene$ensembl_gene_id %in% all_data_tsg_filtered$ensembl_gene_id & aver_meth_gene$diff_meth > 0), ], "down300_CC_density", "down300_GG_density", "Template strand", "Non-template strand", "diff_meth", "Correlation of promoter methylation downregulation extent\nwith G-run Frequency in TSS down 300bp (Down-TSG )", "spearman") +
  theme(legend.position = "bottom", legend.box = "vertical") +
  xlab("Spearman's Rho")

# Plot multiple correlations for expression
exp_diff <- PlotMutipleCor(all_data_tsg_filtered, "down300_CC_density", "down300_GG_density", "Template strand", "Non-template strand", "diff_exp", "Correlation of expression downregulation extent\nwith G-run Frequency in TSS down 300bp (Down-TSG )", "spearman") +
  theme(legend.position = "bottom", legend.box = "vertical") +
  xlab("Spearman's Rho")

# Identify significant G4 motifs in the promoter region
sig_G4_promoter <- get_diff(all_data_tsg_filtered, "rG4_down_exp", "diff_exp", "+")
pff_G4_promoter_all <- plotForesttwogroup(sig_G4_promoter, "Abs(log2FC) of Down-TSG in TCGA(TSS down_1kb)", "4.2mm", "rG4-seq_G4s+", "rG4-seq_G4s-", spacingx = 1.5)
ggarrange(pff_G4_promoter_all$p1)

# Identify significant G4 motifs in the 300 bp downstream region
sig_G4_promoterdown300 <- get_diff(all_data_tsg_filtered, "rG4_down300_exp", "diff_exp", "+")
pff_G4_promoter_alldown300 <- plotForesttwogroup(sig_G4_promoterdown300, "Expression downregulation extent of Down-TSG in TCGA cohort(TSS down_300bp)", "4.2mm", "rG4-seq_G4s+", "rG4-seq_G4s-", spacingx = 1)
ggarrange(pff_G4_promoter_alldown300$p1)

# Load validation datasets
load("./data/processing_validation_datasets_new.RData")

# Filter low data for other cohorts
all_data_other_cohort_filtered <- filter_low_data(all_data_other_cohort)
all_data_other_cohort_filtered$diff_exp <- abs(all_data_other_cohort_filtered$diff_exp)

# Add G-run density values for the 300 bp downstream region
all_data_other_cohort_filtered$down300_CC_density <- genes_infor[all_data_other_cohort_filtered$ensembl_gene_id, "down300_CC_density"]
all_data_other_cohort_filtered$down300_GG_density <- genes_infor[all_data_other_cohort_filtered$ensembl_gene_id, "down300_GG_density"]

# Plot correlation of expression downregulation extent and G-run frequency in validation cohort
p_diff_validate <- plotCorbox(all_data_other_cohort_filtered, "down300_CC_density", "down300_GG_density", "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", "Template strand", "Non-template strand", "diff_exp", "Correlation of expression downregulation extent\nand G-run frequency in Down-TSGs(Validate cohort)", "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

# Identify significant G4 motifs in the validation cohort
all_data_other_cohort_filtered$rG4_down300_exp <- genes_infor[all_data_other_cohort_filtered$ensembl_gene_id, "rG4_down300_exp"]
sig_G4_promoterValidate <- get_diff(all_data_other_cohort_filtered, "rG4_down_exp", "diff_exp", "+")
pff_G4_promoter_allValidate <- plotForesttwogroup(sig_G4_promoterValidate, "Abs(log2FC) of expression in Down-TSG \nin validation cohort(TSS down_1kb)", "4.2mm", "rG4-seq_G4s+", "rG4-seq_G4s-", spacingx = 1.5)
ggarrange(pff_G4_promoter_allValidate$p1)

# Identify significant G4 motifs in the 300 bp downstream region in validation cohort
sig_G4_promoterdown300Validate <- get_diff(all_data_other_cohort_filtered, "rG4_down300_exp", "diff_exp", "+")
pff_G4_promoter_alldown300Validate <- plotForesttwogroup(sig_G4_promoterdown300Validate, "Expression downregulation extent in Down-TSG \nin validation cohort(TSS down_300bp)", "4.2mm", "rG4-seq_G4s+", "rG4-seq_G4s-", spacingx = 1.5)
ggarrange(pff_G4_promoter_alldown300Validate$p1)

# Plot G-run density for RNA-G4 validated by rG4-seq
p_rg4 <- plot_class(genes_infor, "rG4_down_exp", "GG_density", "RNA-G4 validated by rG4-seq", "G-run frequency") +
  ggtitle("TSS downstream") +
  scale_fill_manual(values = c("#99C5FF", "#104E8B")) +
  ylim(0, 0.25)

p_rg4down300 <- plot_class(genes_infor[genes_infor$Gene.Type == "TSG", ], "rG4_down300_exp", "down300_GG_density", "RNA-G4 validated by rG4-seq", "G-run frequency in TSGs") +
  ggtitle("TSS downstream 300bp") +
  scale_fill_manual(values = c("#99C5FF", "#104E8B")) +
  ylim(0, 0.2)

# Arrange plots for Figure 7
pall7_1 <- ggarrange(p_diff, exp_diff, nrow = 1, ncol = 2, labels = c(letters[1:2]), font.label = list(size = 18, color = "black", face = "bold"), widths = c(1, 1.2))
pall7_2 <- ggarrange(p_diff_meth, meth_diff, nrow = 1, ncol = 2, labels = c(letters[3:4]), font.label = list(size = 18, color = "black", face = "bold"), widths = c(1, 1.2))
pall7_3 <- ggarrange(p_rg4down300, pff_G4_promoter_alldown300$p1, nrow = 1, ncol = 2, labels = c(letters[5:6]), font.label = list(size = 18, color = "black", face = "bold"), widths = c(1, 4.6))

# Save Figure 7 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/all_figure5.tiff", width = 15, height = 17, res = 300, units = "in", compression = "lzw")
print(ggarrange(pall7_1, pall7_2, pall7_3, nrow = 3, ncol = 1))
dev.off()


# Plot scatter plots for Supplementary Figures 12 and 13
sp3all <- plot_scatter_2speaman(all_data_tsg_filtered, "down300_GG_density", "diff_exp", "G-run frequencies in\nTSS NT-down 300bp", "Expression downregulation extent in Down-TSG", 0.01, 15, "black") +
  theme(strip.background = element_blank()) +
  guides(colour = "none", fill = "none") +
  scale_color_manual(values = c("black")) +
  theme_hd() +
  facet_wrap(~cancer)

# Save Supplementary Figure 12 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure12.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")
print(sp3all + theme_hd_minimal())
dev.off()

sp3methall <- plot_scatter_2speaman(aver_meth_gene[which(aver_meth_gene$ensembl_gene_id %in% all_data_tsg_filtered$ensembl_gene_id & aver_meth_gene$diff_meth > 0), ], "down300_GG_density", "diff_meth", "G-run frequencies in\nTSS NT-down 300bp", "Promoter methylation downregulation extent in Down-TSG", 0.01, 15, "black") +
  theme(strip.background = element_blank()) +
  guides(colour = "none", fill = "none") +
  scale_color_manual(values = c("black")) +
  theme_hd() +
  facet_wrap(~cancer)

# Save Supplementary Figure 13 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure13.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")
print(sp3methall + theme_hd_minimal())
dev.off()

# Plot scatter plots for Supplementary Figure 14
pall_other <- ggarrange(p_diff_validate, pff_G4_promoter_alldown300Validate$p1, nrow = 2, ncol = 1, labels = c(letters[1:2]), font.label = list(size = 18, color = "black", face = "bold"), heights = c(1.2, 1))

# Save Supplementary Figure 14 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure14.tiff", width = 12, height = 12, res = 300, units = "in", compression = "lzw")
print(pall_other)
dev.off()

# Define a function to filter low data for Oncogenes
filter_low_data_Oncogene <- function(all_data_tsg) {
  # Extract unique cancer types
  cancers <- unique(all_data_tsg$cancer)
  result <- c()  # Initialize an empty result vector
  
  # Loop through each cancer type
  for (cancer in cancers) {
    # Subset data for the current cancer type
    temp_data <- all_data_tsg[all_data_tsg$cancer == cancer, ]
    quantile_cutoff <- 1  # Set a quantile cutoff value
    
    # Filter data for Oncogenes with positive differential expression
    # and exclude those with both normal and tumor expression below the cutoff
    temp_data <- temp_data[which(temp_data$Gene.Type == "Oncogene" & 
                                   temp_data$diff_exp > 0 & 
                                   !(temp_data$norm_exp < quantile_cutoff & 
                                       temp_data$tumor_exp < quantile_cutoff)), ]
    
    # Bind the filtered data to the result
    result <- rbind(result, temp_data)
  }
  
  return(result)  # Return the filtered data
}

# Apply the filtering function to the dataset
all_data_tsg_filtered_Onco <- filter_low_data_Oncogene(all_data_tsg)

# Plot forest plot for expression differences in Oncogenes
p.ff_diff_onco <- plotForest(all_data_tsg_filtered_Onco, "diff_exp", "down300_GG_density", 
                             "(G-run frequencies in TSS NT-down 300bp VS log2FC\nof expression in up-regulated Oncogenes) in TCGA cohort", 
                             0.9, 1.2, "0.8cm", method_x = "spearman")

# Calculate absolute methylation differences
aver_meth_gene$diff_methabs <- abs(aver_meth_gene$diff_meth)

# Plot forest plot for methylation differences in Oncogenes
p.ff_diff_oncometh <- plotForest(aver_meth_gene[which(aver_meth_gene$ensembl_gene_id %in% all_data_tsg_filtered_Onco$ensembl_gene_id & 
                                                        aver_meth_gene$diff_meth < 0), ], 
                                 "diff_methabs", "down300_GG_density", 
                                 "(G-run frequencies in TSS NT-down 300bp VS Abs(log2FC) of\npromoter methylation in up-regulated Oncogenes) in TCGA cohort", 
                                 0.9, 1.2, "0.8cm", method_x = "spearman")

# Plot correlation boxplot for expression differences in Oncogenes
p_diff_Onco <- plotCorbox(all_data_tsg_filtered_Onco, "down300_CC_density", "down300_GG_density", 
                          "TSS down_300bp", "CC_density", "GG_density", "TSS down_1kb", 
                          "Template strand", "Non-template strand", "diff_exp", 
                          "Spearman's Rho(G-run frequencies in TSS NT-down 300bp\nVS log2FC of expression in up-regulated Oncogenes)", 
                          "spearman") +
  ylab("Spearman's Rho") +
  xlab("Regions in different strand of TSS downstream")

# Add rG4-seq data for Oncogenes
all_data_tsg_filtered_Onco$rG4_down300_exp <- genes_infor[all_data_tsg_filtered_Onco$ensembl_gene_id, "rG4_down300_exp"]

# Identify significant G4 motifs in the 300 bp downstream region for Oncogenes
sig_G4_promoterdown300Onco <- get_diff(all_data_tsg_filtered_Onco, "rG4_down300_exp", "diff_exp", "+")
pff_G4_promoter_alldown300Onco <- plotForesttwogroup(sig_G4_promoterdown300Onco, 
                                                     "Abs(log2FC) of expression in up-regulated Oncogenes in TCGA(TSS down_300bp)", 
                                                     "4.2mm", "rG4-seq_G4s+", "rG4-seq_G4s-", spacingx = 1)
ggarrange(pff_G4_promoter_alldown300Onco$p1)

# Define a function to filter low data for Oncogenes in other cohorts
filter_low_dataOncoArray <- function(all_data_tsg) {
  # Extract unique cancer types
  cancers <- unique(all_data_tsg$cancer)
  result <- c()  # Initialize an empty result vector
  
  # Loop through each cancer type
  for (cancer in cancers) {
    # Subset data for the current cancer type
    temp_data <- all_data_tsg[all_data_tsg$cancer == cancer, ]
    quantile_cutoff <- (quantile(temp_data$norm_exp, 0.2, na.rm = TRUE) + 
                          quantile(temp_data$tumor_exp, 0.2, na.rm = TRUE)) / 2
    print(quantile_cutoff)  # Print the quantile cutoff for reference
    
    # Filter data for Oncogenes with positive differential expression
    # and exclude those with both normal and tumor expression below the cutoff
    temp_data <- temp_data[which(temp_data$Gene.Type == "Oncogene" & 
                                   temp_data$diff_exp > 0 & 
                                   !(temp_data$norm_exp < quantile_cutoff & 
                                       temp_data$tumor_exp < quantile_cutoff)), ]
    
    # Bind the filtered data to the result
    result <- rbind(result, temp_data)
  }
  
  return(result)  # Return the filtered data
}

# Apply the filtering function to the dataset
all_data_other_cohort_filteredOnCo <- filter_low_dataOncoArray(all_data_other_cohort)

# Add G-run density values for the 300 bp downstream region
all_data_other_cohort_filteredOnCo$down300_CC_density <- genes_infor[all_data_other_cohort_filteredOnCo$ensembl_gene_id, "down300_CC_density"]
all_data_other_cohort_filteredOnCo$down300_GG_density <- genes_infor[all_data_other_cohort_filteredOnCo$ensembl_gene_id, "down300_GG_density"]

# Add rG4-seq data for Oncogenes in other cohorts
all_data_other_cohort_filteredOnCo$rG4_down300_exp <- genes_infor[all_data_other_cohort_filteredOnCo$ensembl_gene_id, "rG4_down300_exp"]

# Plot forest plot for expression differences in Oncogenes in validation datasets
p.ff_diffveriOnco <- plotForest(all_data_other_cohort_filteredOnCo, "diff_exp", "down300_GG_density", 
                                "Spearman's Rho(G-run frequencies in TSS NT-down 300bp VS log2FC\nof expression in up-regulated Oncogenes) in validation datasets", 
                                0.9, 1.2, "0.8cm", method_x = "spearman")

# Identify significant G4 motifs in the 300 bp downstream region in validation datasets
sig_G4_promoterdown300Onco2 <- get_diff(all_data_other_cohort_filteredOnCo, "rG4_down300_exp", "diff_exp", "+")
pff_G4_promoter_alldown300Onco2 <- plotForesttwogroup(sig_G4_promoterdown300Onco2, 
                                                      "Abs(log2FC) of expression in up-regulated Oncogenes in validation datasets(TSS down_300bp)", 
                                                      "4.2mm", "rG4-seq_G4s+", "rG4-seq_G4s-", spacingx = 1)
ggarrange(pff_G4_promoter_alldown300Onco2$p1)

# Plot G-run density for all genes
downNT_all <- plot_all_genesSig(genes_infor, "down300_GG_density", "G-run frequencies in\n TSS NT-down 300bp")

# Arrange plots for Supplementary Figure 19
p0 <- ggarrange(downNT_all, p_diff_Onco, p.ff_diff_onco, nrow = 1, ncol = 3, 
                labels = c(letters[1:3]), 
                font.label = list(size = 18, color = "black", face = "bold"), 
                widths = c(1.5, 1.5, 2.7))

p1 <- ggarrange(pff_G4_promoter_alldown300Onco$p1, p.ff_diff_oncometh, 
                pff_G4_promoter_alldown300Onco2$p1, p.ff_diffveriOnco, 
                nrow = 2, ncol = 2, labels = c(letters[4:7]), 
                font.label = list(size = 18, color = "black", face = "bold"), 
                widths = c(3, 2.7))

# Save Supplementary Figure 19 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure19.tiff", 
     width = 21, height = 16, res = 300, units = "in", compression = "lzw")
print(ggarrange(p0, p1, nrow = 2, ncol = 1, heights = c(1, 2)))
dev.off()


save.image("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/rdata/figure5.rdata.RData")

