# Set the working directory to the project folder
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new")

# Load custom functions for evolution analysis
source('functions/functions_evolution_5UTR.r')

# Load promoter region data
load(file="./data/all_meth_exp_signal_promoter.rdata")
load(file="./data/other_species_all_data.rdata")
gene_pos<-genes_infor[,c("genes_main_transcript","chromosome_name","transcript_start","transcript_end","strand","ensembl_gene_id")]
colnames(gene_pos)<-c("trans_ID","seqnames","start","end","strand","gene_id")
gene_pos$strand[gene_pos$strand=="1"]="+"
gene_pos$strand[gene_pos$strand=="-1"]="-"
gene_pos$seqnames<-paste("chr",gene_pos$seqnames,sep="")
gene_pos_data<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=T)

# Calculate nucleotide density upstream and downstream of TSS
pos_all_up <- getFreq(gene_pos_data, -1000, 0, Hsapiens, "+")
pos_all_down <- getFreq(gene_pos_data, 0, 1000, Hsapiens, "+")

# Add nucleotide density data to genes_infor
genes_infor$CC_density <- pos_all_down[genes_infor$ensembl_gene_id, "CC_density"]
genes_infor$GG_density <- pos_all_down[genes_infor$ensembl_gene_id, "GG_density"]
genes_infor$AA_density <- pos_all_down[genes_infor$ensembl_gene_id, "AA_density"]
genes_infor$CpG_OE <- pos_all_down[genes_infor$ensembl_gene_id, "CpG_OE"]

# Calculate GG_density and categorize into treatment or control groups
median_GG <- quantile(genes_infor$GG_density, 0.75)
genes_infor$GG_group <- ifelse(genes_infor$GG_density > median_GG, "treatment", "control")
genes_infor$GG_group <- as.factor(genes_infor$GG_group)

# Merge with GTEX median TPM data
GTEXmedian_tpm_infor <- merge(GTEXmedian_tpm, genes_infor, by.x = "Gene.ID", by.y = "ensembl_gene_id")

# Plot multiple correlations for expression vs. G-run frequency
exp2 <- PlotMutipleCor(GTEXmedian_tpm_infor, "CC_density", "GG_density", "Template", "Non-template", "norm_exp", "Expression vs. G-run frequency", "spearman") +
  theme(legend.position = "bottom") +
  ylab("GTEX Cohort") +
  labs(fill = "Strand") +
  xlab("Spearman's rho (TSS downstream 1kb)")

# Plot multiple correlations for promoter methylation vs. G-run frequency
meth1 <- PlotMutipleCor(aver_meth_gene, "CC_density", "GG_density", "Template strand", "Non-template strand", "norm_meth", "Promoter methylation vs. G-run frequency", "spearman") +
  xlab("Spearman's rho (TSS downstream 1kb)") +
  theme(legend.position = "none") +
  ylab("TCGA Cohort")

# Plot multiple correlations for promoter DHS signal vs. G-run frequency
genes_chromatin$cancer <- gsub(" ", "\n", genes_chromatin$cancer)
signal <- PlotMutipleCor(genes_chromatin, "CC_density", "GG_density", "Template", "Non-template", "signal", "Promoter DHS signal vs. G-run frequency", "spearman") +
  ylab("ENCODE Cohort") +
  xlab("Spearman's rho (TSS downstream 1kb)") +
  theme(legend.position = "none")

types=sort(as.character(unique(genes_infor$G4_predict_type)))
# Plot boxplot for G4 density in different strands
p1_g4 <- ggplot(genes_infor, aes(x = G4_down_predict_neg, y = CC_density)) +
  geom_signif(comparisons = getComp(c("+", "-")), map_signif_level = TRUE, col = "black", test = 'wilcox.test', step_increase = 0.1, tip_length = 0.01, y_position = 0.15) +
  geom_boxplot(alpha = I(0.7), outlier.colour = NA, aes(fill = G4_down_predict_neg)) +
  ylab("G-run frequency") +
  theme_hd_minimal() +
  guides(fill = "none") +
  scale_fill_manual(values = colorRampPalette(c("white", "#40A4DE"))(6)[c(2, 6)]) +
  ggtitle("Template strand") +
  xlab("Predicted G4\nin TSS downstream 1kb")

p2_g4 <- ggplot(genes_infor, aes(x = G4_down_predict_pos, y = GG_density)) +
  geom_signif(comparisons = getComp(c("+", "-")), map_signif_level = TRUE, col = "black", test = 'wilcox.test', step_increase = 0.1, tip_length = 0.01, y_position = 0.15) +
  geom_boxplot(alpha = I(0.7), outlier.colour = NA, aes(fill = G4_down_predict_pos)) +
  ylab("G-run frequency") +
  theme_hd_minimal() +
  guides(fill = "none") +
  scale_fill_manual(values = colorRampPalette(c("white", "#F07F16C7"))(6)[c(2, 6)]) +
  xlab("Predicted G4\nin TSS downstream 1kb") +
  ggtitle("Non-template strand")

# Plot boxplot for expression in different G4 prediction types
p2_exp <- ggplot(data = all_data_tsg[all_data_tsg$cancer == "COAD" & !is.na(all_data_tsg$Gene.Type), ], aes(y = norm_exp, x = G4_predict_type, fill = G4_predict_type)) +
  geom_boxplot() +
  geom_signif(comparisons = getComp(types), map_signif_level = TRUE, col = "black", test = 'wilcox.test', step_increase = 0.1) +
  ylab("Expression ") +
  theme_hd_minimal() +
  ggtitle("TCGA-COAD normal") +
  theme(axis.text.x = element_blank(), legend.position = "none", axis.title.x = element_blank()) +
  scale_fill_manual("Predicted G4s in different strand", values = c("#CDCDC1", "#00B2EE", "#FFC125", "#8FBC8F"))

# Plot boxplot for promoter DHS signal in different G4 prediction types
p3_signal <- ggplot(data = genes_chromatin[genes_chromatin$cancer == "colon" & !is.na(genes_chromatin$Gene.Type), ], aes(y = signal, x = G4_predict_type, fill = G4_predict_type)) +
  geom_boxplot() +
  geom_signif(comparisons = getComp(types), map_signif_level = TRUE, col = "black", test = 'wilcox.test', step_increase = 0.1) +
  ylab("Promoter DHS signal") +
  theme_hd_minimal() +
  ggtitle("ENCODE-Colon") +
  theme(axis.text.x = element_blank(), legend.position = "none", axis.title.x = element_blank()) +
  scale_fill_manual("Predicted G4s in different strand", values = c("#CDCDC1", "#00B2EE", "#FFC125", "#8FBC8F"))

# Plot boxplot for promoter methylation in different G4 prediction types
p4_meth <- ggplot(data = aver_meth_gene[aver_meth_gene$cancer == "COAD", ], aes(y = norm_meth, x = G4_predict_type, fill = G4_predict_type)) +
  geom_boxplot() +
  geom_signif(comparisons = getComp(types), map_signif_level = TRUE, col = "black", test = 'wilcox.test', step_increase = 0.1) +
  ylab("Promoter methylation") +
  ggtitle("TCGA-COAD normal") +
  theme_hd_minimal() +
  theme(axis.text.x = element_blank(), legend.position = "bottom", axis.title.x = element_blank()) +
  scale_fill_manual("Predicted G4s in different strand", values = c("#CDCDC1", "#00B2EE", "#FFC125", "#8FBC8F"))

# Arrange plots
pall1 <- ggarrange(meth1, signal, nrow = 2, ncol = 1, labels = c(letters[2:3]), font.label = list(size = 18, color = "black", face = "bold"))
pall2 <- ggarrange(p1_g4, p2_g4, nrow = 1, ncol = 2, labels = c(letters[4:5]), font.label = list(size = 18, color = "black", face = "bold"), heights = c(2, 3))
pall3 <- ggarrange(pall1, pall2, nrow = 2, ncol = 1, font.label = list(size = 18, color = "black", face = "bold"), heights = c(2, 1))
pall4 <- ggarrange(exp2, pall3, nrow = 1, ncol = 2, labels = c(letters[1], ""), font.label = list(size = 18, color = "black", face = "bold"), widths = c(2, 1.3), common.legend = TRUE, legend = "bottom")
pall5 <- ggarrange(p2_exp, p3_signal, p4_meth, nrow = 1, ncol = 3, labels = c(letters[6:8]), font.label = list(size = 18, color = "black", face = "bold"), common.legend = TRUE, legend = "bottom")
pall6 <- ggarrange(pall4, pall5, nrow = 2, ncol = 1, font.label = list(size = 18, color = "black", face = "bold"), heights = c(3, 1))

# Save figure
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/all_figure2.tiff", width = 17, height = 17, res = 300, units = "in", compression = "lzw")
print(pall6)
dev.off()



df_all <- as.data.frame(table(genes_infor$G4_predict_type, genes_infor$Gene.Type, genes_infor$age_type2))
p_freq <- ggplot(df_all, aes(x = Var3, y = Freq, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual("Predicted G4s", values = c("#CDCDC1", "#00B2EE", "#FFC125", "#8FBC8F")) +
  ylab("Frequency") +
  theme_hd_minimal2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom", axis.title.x = element_blank()) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Genes with predicted G4s\nin different strand")

# Plot correlation between G-run frequency and G4hunter score
p_cor <- ggplot(genes_infor, aes(GG_density, g4hunterPosscore)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  stat_cor(method = "spearman", size = 4, color = "black") +
  labs(x = "G-run frequency", y = "G4hunter score") +
  theme_hd_minimal()

# Create Venn diagrams for overlap of G4s
x <- list(Predicted_G4_T = genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_neg == "+"], `G4-seq_validated_G4` = genes_infor$ensembl_gene_id[genes_infor$G4_down_exp == "+"])
title_x <- "\n\n\n\n\n\n\n\n\n\n\nOverlap of G4s"
p <- venn.diagram(x, fill = c("#00B2EE", "#104E8B"), alpha = c(0.5, 0.5), cex = 1, cat.fontface = "bold", lty = 2, fontfamily = 3, cat.cex = 1, sub.fontface = "bold", margin = 0.1, sub = "", sub.pos = c(0.3, -0.05), filename = NULL, main = title_x, main.fontface = "bold", main.cex = 1.2, col = NA, ext.text = TRUE, ext.line.lwd = 2, ext.dist = -0.15, ext.length = 0.9, ext.pos = -14, rotation.degree = 185, cat.pos = c(120, 240))
pVNN2 <- ggarrange(p)

x <- list(Predicted_G4_NT = genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_pos == "+"], `G4-seq_validated_G4` = genes_infor$ensembl_gene_id[genes_infor$G4_down_exp == "+"])
title_x <- "\n\n\n\n\n\n\n\n\n\n\nOverlap of G4s"
p <- venn.diagram(x, fill = c("#FFC125", "#104E8B"), alpha = c(0.5, 0.5), cex = 1, cat.fontface = "bold", lty = 2, fontfamily = 3, cat.cex = 1, sub.fontface = "bold", margin = 0.1, sub = "", sub.pos = c(0.3, -0.05), filename = NULL, main = title_x, main.fontface = "bold", main.cex = 1.2, col = NA, ext.text = TRUE, ext.line.lwd = 2, ext.dist = -0.15, ext.length = 0.9, ext.pos = -14, rotation.degree = 5)
pVNN <- ggarrange(p)

# Load G4hunter results and process
library(readr)
g4hunter <- read_delim("data/g4hunter/Results_all_gene_0_1000_seq/all_gene_0_1000_seq-Merged_add_gene_names.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
g4hunter <- as.data.frame(g4hunter[g4hunter$Start != "Start", ])
g4hunter$Score <- as.numeric(g4hunter$Score)
g4hunter_max <- aggregate(g4hunter$Score, list(trans_id = str_split_fixed(g4hunter[, 1], "_", 2)[, 1]), max)
g4hunter_max <- g4hunter_max[g4hunter_max$x > 1, ]
rownames(g4hunter_max) <- g4hunter_max[, 1]
genes_infor$g4hunterPos <- ifelse(genes_infor$ensembl_transcript_id %in% g4hunter_max$trans_id[g4hunter_max$x > 1], "+", "-")
genes_infor$g4hunterPosscore <- g4hunter_max[genes_infor$ensembl_transcript_id, 2]

# Create Venn diagram for overlap of G4s
x <- list(pqsfinder = genes_infor$ensembl_gene_id[genes_infor$G4_down_predict_pos == "+"], G4Hunter = genes_infor$ensembl_gene_id[genes_infor$g4hunterPos == "+"])
title_x <- "\n\n\n\n\n\n\n\n\n\n\nOverlap of G4s"
p <- venn.diagram(x, fill = c("#BFEFFF", "#00B2EE"), alpha = c(0.5, 0.5), cex = 1, cat.fontface = "bold", lty = 2, fontfamily = 3, cat.cex = 1, sub.fontface = "bold", margin = 0.1, sub = "", sub.pos = c(0.3, -0.05), filename = NULL, main = title_x, main.fontface = "bold", main.cex = 1.2, col = NA, ext.text = TRUE, ext.line.lwd = 2, ext.dist = -0.15, ext.length = 0.9, ext.pos = -14, rotation.degree = 185, cat.pos = c(100, 260))
pVNN3 <- ggarrange(p)

# Plot regression results for expression vs. G-run frequency
p1 <- plotStrandCor(all_data_tsg, "norm_exp", "GG_density", "CC_density", "age_type2", method = "spearman") + ggtitle("TCGA Cohort(Expression vs. G-run frequency)") + ylab("Spearman's rho") + ylim(-0.02, 0.36)

p2 <- plotStrandCor(GTEXmedian_tpm_infor, "norm_exp", "GG_density", "CC_density", "age_type2", method = "spearman") + ggtitle("GTEX Cohort(Expression vs. G-run frequency)") + ylab("Spearman's rho") + ylim(-0.01, 0.36)

exp1 <- PlotMutipleCor(all_data_tsg, "CC_density", "GG_density", "Template", "Non-template", "norm_exp", "Expression vs. G-run frequency", "spearman") + theme(legend.position = "bottom", legend.box = "vertical") + labs(fill = "Strand") + xlab("Spearman's rho (TSS downstream 1kb)") + ylab("TCGA Cohort")

# Arrange plots
palls1 <- ggarrange(exp1, p2, nrow = 1, ncol = 2, labels = c(letters[1:2]), font.label = list(size = 18, color = "black", face = "bold"), widths = c(1, 2))
palls2 <- ggarrange(p1, p_freq, nrow = 1, ncol = 2, labels = c(letters[3:4]), font.label = list(size = 18, color = "black", face = "bold"), widths = c(2, 1))
palls3 <- ggarrange(pVNN, pVNN2, pVNN3, p_cor, nrow = 1, ncol = 4, labels = c(letters[5:8]), font.label = list(size = 18, color = "black", face = "bold"))

palls4 <- ggarrange(palls1, palls2, palls3, nrow = 3, ncol = 1, heights = c(1.2, 1.2, 0.8))

# Save figure
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure2.tiff", width = 15, height = 15, res = 300, units = "in", compression = "lzw")
print(palls4)
dev.off()


# Identify alternative transcripts
alternative_pos <- all_data_genes$get_gene_pos[
  !paste(all_data_genes$get_gene_pos$ensembl_gene_id, all_data_genes$get_gene_pos$transcript_start) %in% 
    paste(genes_infor$ensembl_gene_id, genes_infor$transcript_start) & 
    all_data_genes$get_gene_pos$transcript_biotype == "protein_coding", 
]

# Calculate the difference in transcript start positions
alternative_pos$diff <- abs(alternative_pos$transcript_start - genes_infor[alternative_pos$ensembl_gene_id, ]$transcript_start)

# Select the top transcript for each gene based on the maximum difference
library(dplyr)
top_tx <- alternative_pos %>%
  group_by(ensembl_gene_id) %>%
  slice_max(diff, n = 1) %>%   # Select the transcript with the maximum difference
  ungroup()

# Create a data frame with unique alternative gene positions
alternative_gene_pos <- unique(as.data.frame(top_tx[, c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "strand", "ensembl_gene_id")]))
colnames(alternative_gene_pos) <- c("trans_ID", "seqnames", "start", "end", "strand", "gene_id")

# Convert strand information to standard format
alternative_gene_pos$strand[alternative_gene_pos$strand == "1"] <- "+"
alternative_gene_pos$strand[alternative_gene_pos$strand == "-1"] <- "-"

# Prefix chromosome names with 'chr'
alternative_gene_pos$seqnames <- paste("chr", alternative_gene_pos$seqnames, sep = "")

# Create a GRanges object from the alternative gene position data
alternative_gene_pos_data <- makeGRangesFromDataFrame(alternative_gene_pos, keep.extra.columns = TRUE)

# Calculate G-run frequencies for the downstream region (1000 bp)
alternative_pos_all_down <- getFreq(alternative_gene_pos_data, 0, 1000, Hsapiens, "+")

# Add G-run density values to the alternative gene positions
alternative_gene_pos$CC_density <- alternative_pos_all_down[alternative_gene_pos$gene_id, "CC_density"]
alternative_gene_pos$GG_density <- alternative_pos_all_down[alternative_gene_pos$gene_id, "GG_density"]
alternative_gene_pos$GG_density_select <- genes_infor[alternative_gene_pos$gene_id, "GG_density"]

# Prepare data for plotting
plot_df <- data.frame(
  gene_id = alternative_gene_pos$gene_id,
  Type = rep(c("CC_density", "GG_density"), each = nrow(alternative_gene_pos)),
  Density = c(alternative_gene_pos$CC_density, alternative_gene_pos$GG_density)
)

plot_df2 <- data.frame(
  gene_id = genes_infor$ensembl_gene_id,
  Type = rep(c("CC_density", "GG_density"), each = nrow(genes_infor)),
  Density = c(genes_infor$CC_density, genes_infor$GG_density)
)

plot_all <- rbind(cbind(classes = "Canonical isoforms", plot_df2), cbind(classes = "Isoforms with \nalternative TSS", plot_df))

# Convert Type to factor for plotting
plot_all$Type <- factor(plot_all$Type, levels = c("CC_density", "GG_density"), labels = c("Template", "Non_template"))

# Plot G-run frequency for canonical and alternative isoforms
p_alternative <- ggplot(plot_all, aes(x = Type, y = Density, fill = Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.9) +
  scale_fill_manual(values = c("#40A4DE", "#F07F16C7")) +
  labs(x = NULL, y = "G-runs frequency\nin TSS downstream") +
  theme_hd_minimal() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  ggsignif::geom_signif(
    comparisons = list(c("Non_template", "Template")),
    map_signif_level = TRUE,
    textsize = 4, vjust = 0.2
  ) +
  scale_y_log10() +
  facet_wrap(~classes, nrow = 1)

# Save alternative gene positions to a CSV file
alternative_gene_pos_select <- alternative_gene_pos[, 1:6]
names(alternative_gene_pos_select)[1] <- "ensembl_transcript_id"
write.table(alternative_gene_pos_select, file = './data/alternative_gene_pos_select.csv', quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

# Function to get GTEX data
get_GTEX_data <- function(file_name, gene_trans) {
  rownames(gene_trans) <- genes_infor$genes_main_transcript
  GTEX_data <- read.csv(file_name)
  GTEX_data <- GTEX_data[GTEX_data$X0 != 0 & GTEX_data$X0 %in% genes_infor$genes_main_transcript, ]
  GTEX_data <- GTEX_data[!(duplicated(GTEX_data$X0) | apply(GTEX_data[, -1], 1, max) == "0.0"), ]
  
  tissue_tpm <- as.data.frame(apply(GTEX_data[, -1], 2, as.numeric))
  rownames(tissue_tpm) <- gene_trans[GTEX_data$X0, ]$ensembl_gene_id
  colnames(tissue_tpm) <- gsub("\\.+",".", colnames(GTEX_data[-1]))
  
  tissue_tpm$Gene.Type <- genes_infor[rownames(tissue_tpm), "Gene.Type"]
  tissue_tpm$Gene.ID <- rownames(tissue_tpm)
  
  all_tissue_data <- tissue_tpm[which(!is.na(tissue_tpm$Gene.Type)), ]
  colnames(all_tissue_data) <- firstup(colnames(all_tissue_data))
  
  all_data_df <- reshape2::melt(all_tissue_data, id.vars = c("Gene.ID", "Gene.Type"))
  all_data_df$value <- log2(all_data_df$value + 1)
  
  tissue <- as.character(unique(all_data_df$variable))
  
  all_data_non_cancer <- all_tissue_data[all_tissue_data$Gene.Type == "Non-Cancer", tissue]
  all_data_TSG <- all_tissue_data[all_tissue_data$Gene.Type == "TSG", tissue]
  
  all_data_df$len_3UTR_log <- log10(genes_infor[all_data_df$Gene.ID, ]$hg38_3UTR_length)
  
  colnames(all_data_df) <- c("Gene.ID", "Gene.Type", "cancer", "norm_exp", "len_3UTR_log")
  
  return(all_data_df)
}

# Load and process GTEX data
GTEX_data <- read.csv("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_alternative_median_transcript_tpm.csv")
GTEX_data <- GTEX_data[!(duplicated(GTEX_data$X0) | apply(GTEX_data[, -1], 1, max) == "0.0"), ]
tissue_tpm <- as.data.frame(apply(GTEX_data[, -1], 2, as.numeric))
rownames(tissue_tpm) <- GTEX_data[, 1]
colnames(tissue_tpm) <- gsub("\\.+",".", colnames(GTEX_data[-1]))

all_data_df <- reshape2::melt(as.matrix(tissue_tpm))
colnames(all_data_df) <- c("Gene.ID", "cancer", "norm_exp")

alternative_gene_pos$Gene.Type <- genes_infor[alternative_gene_pos$gene_id, "Gene.Type"]
GTEXmedian_tpm_infor_alter <- merge(all_data_df, alternative_gene_pos, by.x = "Gene.ID", by.y = "trans_ID")

# Plot correlation of expression and G-run frequency for alternative isoforms
p_cor_alter <- PlotMutipleCor(GTEXmedian_tpm_infor_alter, "CC_density", "GG_density", "Template", "Non-template", "norm_exp", "Correlation of expression and G-run frequency in different strand in TSS downstream", "spearman") +
  theme(legend.position = "bottom", legend.box = "vertical") +
  ylab("GTEX Cohort") +
  labs(fill = "Strand") +
  xlab("Spearman's rho (expression vs G-run frequency)") +
  ggtitle("Isoforms with alternative TSS")

# Plot scatter plots for canonical and alternative isoforms
sp1_cpd_exp1GTEX <- plot_scatter_2speaman(GTEXmedian_tpm_infor[GTEXmedian_tpm_infor$cancer == "Colon.Sigmoid", ], "GG_density", "norm_exp", "G-run frequency in non-template strand", "Expression in GTEX-Sigmoid\nColon normal samples", 0.01, 15, "black") +
  theme(strip.background = element_blank()) +
  guides(colour = "none", fill = "none") +
  scale_color_manual(values = c("black")) +
  ggtitle("Canonical isoforms") +
  theme_hd() +
  ylim(0, 16)

sp1_cpd_exp1GTEX_alter <- plot_scatter_2speaman(GTEXmedian_tpm_infor_alter[GTEXmedian_tpm_infor_alter$cancer == "Colon.Sigmoid", ], "GG_density", "norm_exp", "G-run frequency in non-template strand", "Expression in GTEX-Sigmoid\nColon normal samples", 0.01, 15, "black") +
  theme(strip.background = element_blank()) +
  guides(colour = "none", fill = "none") +
  scale_color_manual(values = c("black")) +
  ggtitle("Isoforms with alternative TSS") +
  theme_hd() +
  ylim(0, 16)

# Arrange plots for Supplementary Figure 3
pall3s <- ggarrange(p_alternative + theme(axis.text.x = element_text(angle = 45, hjust = 1)), 
                    sp1_cpd_exp1GTEX_alter, 
                    sp1_cpd_exp1GTEX, 
                    nrow = 3, ncol = 1, 
                    labels = c(letters[1:3]), 
                    font.label = list(size = 18, color = "black", face = "bold"))

pall3s2 <- ggarrange(pall3s, p_cor_alter, nrow = 1, ncol = 2, 
                     labels = c("", letters[4]), 
                     font.label = list(size = 18, color = "black", face = "bold"), 
                     widths = c(1, 1.5))

# Save Supplementary Figure 3 as a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure3.tiff", 
     width = 15, height = 15, res = 300, units = "in", compression = "lzw")
print(pall3s2)
dev.off()




tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure4.tiff",width = 16,height = 16,res=300,units="in",compression = "lzw")
types=sort(as.character(unique(genes_infor$G4_predict_type)))
p2_exp_all=ggplot(data = all_data_tsg[all_data_tsg$cancer!="COAD"& !is.na(all_data_tsg$Gene.Type),],aes(y=norm_exp,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Expression ")+theme_hd_minimal()+theme(axis.text.x =element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))+facet_wrap(~cancer)
print(p2_exp_all)
dev.off()



tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure5.tiff",width = 16,height = 16,res=300,units="in",compression = "lzw")
p3_signal_all=ggplot(data = genes_chromatin[genes_chromatin$cancer!="colon"& !is.na(genes_chromatin$Gene.Type),],aes(y=signal,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("DHS signal in TSS downstream  ")+theme_hd_minimal()+theme(axis.text.x = element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))+facet_wrap(~cancer)
print(p3_signal_all)
dev.off()
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure6.tiff",width = 16,height = 16,res=300,units="in",compression = "lzw")
p_meth_all<-ggplot(data = aver_meth_gene[aver_meth_gene$cancer!="COAD",],aes(y=norm_meth,x=G4_predict_type,fill=G4_predict_type))+geom_boxplot()+geom_signif(comparisons = getComp(types),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Methylation in TSS downstream ")+theme_hd_minimal()+theme(axis.text.x = element_blank(),legend.position = "bottom",axis.title.x = element_blank())+scale_fill_manual("Predicted G4s in different strand",values = c("#CDCDC1","#00B2EE", "#FFC125", "#8FBC8F"))+facet_wrap(~cancer)
print(p_meth_all)
dev.off()

pos_all_downMask <- getFreqmaskCG(gene_pos_data, 0, 1000, Hsapiens, "+")
pos_all_downMask300 <- getFreqmaskCG(gene_pos_data, 0, 300, Hsapiens, "+")

# Add masked nucleotide density data to genes_infor
genes_infor$CC_densitydownMask <- pos_all_downMask[genes_infor$ensembl_gene_id, "CC_density"]
genes_infor$GG_densitydownMask <- pos_all_downMask[genes_infor$ensembl_gene_id, "GG_density"]

p_masked<-plotSubSigParied(genes_infor[!is.na(genes_infor$age_type2)& genes_infor$Gene.Type=="Non-Cancer",], "CC_densitydownMask","GG_densitydownMask","age_type2","",cols1=c("#40A4DE", "#F07F16C7"),"Template ","Non-template","Age groups(Million years)","G-run frequency (CpG masked)",cols_x="black")+labs(fill="Strand")+ggtitle("TSS downstream")

GTEXmedian_tpm_infor<-merge(GTEXmedian_tpm,genes_infor,by.x="Gene.ID",by.y="ensembl_gene_id")
exp_masked=PlotMutipleCor(GTEXmedian_tpm_infor,"CC_densitydownMask","GG_densitydownMask","Template","Non-template","norm_exp","Correlation of expression and frequency of G-runs\nin different strand in TSS downstream(CpG masked)","spearman")+theme(legend.position = "bottom",legend.box ="vertical")+ylab("GTEX Cohort")+labs(fill="Strand")+xlab("Spearman's rho")

all_data_tsg$CC_densitydownMask=genes_infor[all_data_tsg$ensembl_gene_id,"CC_densitydownMask"]
aver_meth_gene$CC_densitydownMask=genes_infor[aver_meth_gene$ensembl_gene_id,"CC_densitydownMask"]
aver_meth_gene$GG_densitydownMask=genes_infor[aver_meth_gene$ensembl_gene_id,"GG_densitydownMask"]

meth_masked=PlotMutipleCor(aver_meth_gene,"CC_densitydownMask","GG_densitydownMask","Template strand","Non-template strand","norm_meth","Correlation of promoter methylation and frequency of G-runs\nin different strands in TSS downstream(CpG masked)","spearman")+theme(legend.position = "bottom",legend.box ="vertical")+xlab("Spearman's rho")


pall0=ggarrange(p_masked, meth_masked,nrow=2,ncol=1,labels=c(letters[c(1,3)]),font.label = list(size = 18, color = "black", face = "bold"))

pall2<-ggarrange(pall0,exp_masked,nrow=1,ncol=2,labels=c("",letters[2]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1.5))

tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure15.tiff",width = 17,height = 17,res=300,units="in",compression = "lzw")
print(pall2)
dev.off()



result <- plot_regression_results(GTEXmedian_tpm_infor, "norm_exp")
p_anova<-result$plot0+ggtitle("Contribution of G-runs frequency to\nexpression under control of CpG O/E")


all_data_tsg$AA_density<-genes_infor[all_data_tsg$ensembl_gene_id,"AA_density"]
all_data_tsg$CpG_OE<-genes_infor[all_data_tsg$ensembl_gene_id,"CpG_OE"]
all_data_tsg$GG_density<-genes_infor[all_data_tsg$ensembl_gene_id,"GG_density"]

result_tcga <- plot_regression_results(all_data_tsg, "norm_exp")
p_tcga_anova<-result_tcga$plot0+ggtitle("Contribution of G-runs frequency to\nexpression under control of CpG O/E")
pos_all_down<-getFreq(gene_pos_data,0,1000,Hsapiens,"+")
genes_infor$AA_density<-pos_all_down[,"AA_density"]
genes_infor$CpG_OE<-pos_all_down[,"CpG_OE"]
pos_all_down300<-getFreq(gene_pos_data,0,300,Hsapiens,"+")
genes_infor$down300_GG_density<-pos_all_down300[,"GG_density"]
genes_infor$down300_AA_density<-pos_all_down300[,"AA_density"]
genes_infor$down300_CpG_OE<-pos_all_down300[,"CpG_OE"]
genes_infor$Gene.Type2<-as.factor(genes_infor$Gene.Type)

all_data_tsg$down300_GG_density<-genes_infor[all_data_tsg$ensembl_gene_id,"down300_GG_density"]
all_data_tsg$down300_AA_density<-genes_infor[all_data_tsg$ensembl_gene_id,"down300_AA_density"]
all_data_tsg$down300_CpG_OE<-genes_infor[all_data_tsg$ensembl_gene_id,"down300_CpG_OE"]
all_data_tsg_filtered=filter_low_data(all_data_tsg )
all_data_tsg_filtered$diff_exp=abs(all_data_tsg_filtered$diff_exp)
result_tcga_df <- plot_regression_results(all_data_tsg_filtered, "diff_exp","down300_")
p_df_anova<-result_tcga_df$plot0+ggtitle("Contribution of G-runs frequency to the abs(log2FC) of \nexpression of Down-TSGs under control of CpG O/E")

load("./data/processing_validation_datasets_new.RData")

filter_low_data<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=(quantile(temp_data$norm_exp,0.2,na.rm=T)+quantile(temp_data$tumor_exp,0.2,na.rm=T))/2
    print(quantile_cutoff)
    temp_data=temp_data[which(temp_data$Gene.Type=="TSG" & temp_data$diff_exp<0  & !(temp_data$norm_exp < quantile_cutoff & temp_data$tumor_exp < quantile_cutoff )),]  
    result=rbind(result,temp_data)
  }
  
  return(result)
}

all_data_other_cohort$down300_GG_density<-genes_infor[all_data_other_cohort$ensembl_gene_id,"down300_GG_density"]
all_data_other_cohort$down300_AA_density<-genes_infor[all_data_other_cohort$ensembl_gene_id,"down300_AA_density"]
all_data_other_cohort$down300_CpG_OE<-genes_infor[all_data_other_cohort$ensembl_gene_id,"down300_CpG_OE"]
all_data_other_cohort_filtered=filter_low_data(all_data_other_cohort  )
all_data_other_cohort_filtered$diff_exp<-abs(all_data_other_cohort_filtered$diff_exp)
result_other_df <- plot_regression_results(all_data_other_cohort_filtered, "diff_exp","down300_")
p_df_other_anova<-result_other_df$plot0+ggtitle("Contribution of G-runs frequency to the abs(log2FC) of\nexpression of Down-TSGs under control of CpG O/E")

pallS1=ggarrange(p_tcga_anova+theme(legend.position = "none"),p_df_anova+theme(legend.position = "none"),p_df_other_anova,nrow=3,ncol=1,labels=letters[2:4],font.label = list(size = 18, color = "black", face = "bold"))
pallS2=ggarrange(p_anova+
                   theme(legend.position = "none"),pallS1,nrow=1,ncol=2,labels=c(letters[1],""),font.label = list(size = 18, color = "black", face = "bold"))

tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure18.tiff",width = 20,height = 16,res=300,units="in",compression = "lzw")
print(pallS2)
dev.off()





save.image("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/rdata/figure2.rdata.RData")



