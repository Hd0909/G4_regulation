# Set working directory to the project folder
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new")

# Load gene information data
genes_infor <- readRDS("./data/genes_infor_downloaded_from_ensembl105_gene_infor.rds")

# Load all gene data
all_data_genes <- readRDS("./data/genes_infor_downloaded_from_ensembl105_all_data_genes.rds")

# Load other species data
load(file="./data/other_species_all_data.rdata")

# Source custom functions for evolution analysis
source('functions/functions_evolution_5UTR.r')

# Load promoter window data
load(file="./data/result_promoter_window_ACTG_1000_1000.rdata")

# Prepare gene position data
gene_pos <- genes_infor[, c("genes_main_transcript", "chromosome_name", "transcript_start", "transcript_end", "strand", "ensembl_gene_id")]
colnames(gene_pos) <- c("trans_ID", "seqnames", "start", "end", "strand", "gene_id")

# Convert strand values to '+' and '-'
gene_pos$strand[gene_pos$strand == "1"] <- "+"
gene_pos$strand[gene_pos$strand == "-1"] <- "-"

# Add 'chr' prefix to chromosome names
gene_pos$seqnames <- paste("chr", gene_pos$seqnames, sep="")

# Convert gene position data to GRanges object
gene_pos_data <- makeGRangesFromDataFrame(gene_pos, keep.extra.columns=TRUE)

# Calculate nucleotide density upstream and downstream of TSS
pos_all_up <- getFreq(gene_pos_data, -1000, 0, Hsapiens, "+")
pos_all_down <- getFreq(gene_pos_data, 0, 1000, Hsapiens, "+")

# Add nucleotide density data to genes_infor
genes_infor$CC_density <- pos_all_down[genes_infor$ensembl_gene_id, "CC_density"]
genes_infor$GG_density <- pos_all_down[genes_infor$ensembl_gene_id, "GG_density"]
genes_infor$AA_density <- pos_all_down[genes_infor$ensembl_gene_id, "AA_density"]
genes_infor$CpG_OE <- pos_all_down[genes_infor$ensembl_gene_id, "CpG_OE"]

# Prepare data for plotting
x <- rbind(cbind(type="TSS_upstream", melt(pos_all_up[1:13])), cbind(type="TSS_downstream", melt(pos_all_down[1:13])))
x$variable <- gsub("_density", "", as.character(x$variable))
x$type <- factor(x$type, levels=c("TSS_upstream", "TSS_downstream"))
x$Gene.Type <- genes_infor[x$Group.1, "Gene.Type"]

# Plot G+C content
pGC_content <- ggplot(x[x$variable == "GC_content", ], aes(x=type, y=value)) +
  geom_signif(comparisons=getComp(unique(as.character(x$type))), map_signif_level=TRUE, col="black", test='wilcox.test', step_increase=0.1, tip_length=0.01) +
  geom_boxplot(alpha=I(0.7), outlier.colour=NA, aes(fill=type)) +
  ylab("G+C content") +
  theme_hd_minimal() +
  guides(fill="none") +
  scale_fill_manual(values=c("#7FFF00", "#458B00")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1))

# Plot G+C content (alternative)
p33 <- ggplot(x[x$variable == "GC_content", ], aes(x=type, y=value)) +
  geom_signif(comparisons=getComp(unique(as.character(x$type))), map_signif_level=TRUE, col="black", test='wilcox.test', step_increase=0.1, tip_length=0.01) +
  geom_boxplot(alpha=I(0.7), outlier.colour=NA, aes(fill=type)) +
  ylab("G+C content") +
  theme_hd_minimal() +
  guides(fill="none") +
  scale_fill_manual(values=c("#7FFF00", "#458B00")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1))


# Filter out Drosophila melanogaster
all_data <- all_data[!all_data$species %in% c("Drosophila.melanogaster"), ]
all_data$species <- factor(all_data$species, levels=levels(all_data$species)[-9])
all_data$species <- factor(all_data$species, levels=levels(all_data$species)[10:1])

# Plot species-specific data
p_spe <- plotSubSigParied(all_data, "CC_density", "GG_density", "species", "", cols1=c("#40A4DE", "#F07F16C7"), "Template ", "Non-template", "species", "G-run frequency\nin TSS downstream", cols_x="black") +
  xlab("Species") +
  ylim(0, 0.2) +
  labs(fill="Strand")

# Plot age-specific data
p_age <- plotSubSigParied(genes_infor[!is.na(genes_infor$age_type2), ], "CC_density", "GG_density", "age_type2", "", cols1=c("#40A4DE", "#F07F16C7"), "Template ", "Non-template", "Age groups (Million years)", "G-run frequency\nin TSS downstream", cols_x="black") +
  ylim(0, 0.18) +
  labs(fill="Strand")

# Calculate sliding window data
a_site <- SlidingWindow(median, 1:result_window_1000_1000$length_x, window=50, step=10)

# Plot density of nucleotides around TSS
bm_4 <- plot_density_TSS7(result_window_1000_1000, a_site, "density", result_window_1000_1000$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id, "", genes_infor, result_window_1000_1000$genes_promoter_df)

# Load hallmark gene sets
library(readgmt)
gmt <- read_gmt("./data/h.all.v2022.1.Hs.symbols.gmt", tidy=T)
term2gene <- as.data.frame(gmt)
term2gene$gene_set <- gsub("HALLMARK_", "", term2gene$gene_set)

# Merge with gene information
x <- merge(term2gene, genes_infor, by.x="gene", by.y="hgnc_symbol")

# Calculate mean density for each gene set
temp_x <- aggregate(x[, c("CC_density", "GG_density")], list(gene_set=x$gene_set), mean)
x$gene_set <- factor(x$gene_set, levels=temp_x$gene_set[order(temp_x$GG_density - temp_x$CC_density)])

# Plot gene set-specific data
terms_x <- plotSubSigParied(x[!is.na(x$age_type2), ], "GG_density", "CC_density", "gene_set", "", cols1=c("#F07F16C7", "#40A4DE"), "Non-template", "Template ", "HallMark Geneset", "G-run frequency\nin TSS downstream", cols_x="black") +
  labs(fill="Strand") +
  coord_flip()

# Arrange plots
pall1 <- ggarrange(bm_4 + theme_hd() + scale_color_manual(labels=c("A", "C", "G", "T"), values=c("#104E8B", "#B23AEE", "#CD2626", "#458B00")) + ylab("Nucleotide frequency") + xlab("Distance from TSS (bp)") + theme(axis.text.x=element_text(angle=90)), p33, nrow=1, ncol=2, labels=c(letters[1:2]), font.label=list(size=18, color="black", face="bold"), widths=c(1, 1))
pall2 <- ggarrange(pall1, p_spe + guides(fill="none"), p_age + guides(fill="none"), nrow=3, ncol=1, labels=c("", letters[3:4]), font.label=list(size=18, color="black", face="bold"))
pall3 <- ggarrange(pall2, terms_x, nrow=1, ncol=2, labels=c("", letters[5]), font.label=list(size=18, color="black", face="bold"), common.legend=T, legend="bottom", widths=c(1, 1.5))

# Save figure
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/all_figure1.tiff", width=15, height=18, res=300, units="in", compression="lzw")
print(pall3)
dev.off()

# Load additional promoter window data
result_window_1000_1000More <- readRDS("./data/result_promoter_window_ACTG_1000_1000more.rdata")

# Calculate density ratios
result_window_1000_1000More$CC_densityOE <- result_window_1000_1000More$CC_density / result_window_1000_1000$C_density
result_window_1000_1000More$CCC_densityOE <- result_window_1000_1000More$CCC_density / result_window_1000_1000$C_density
result_window_1000_1000More$GG_densityOE <- result_window_1000_1000More$GG_density / result_window_1000_1000$G_density
result_window_1000_1000More$GGG_densityOE <- result_window_1000_1000More$GGG_density / result_window_1000_1000$G_density

# Calculate sliding window data
a_site <- SlidingWindow(median, 1:result_window_1000_1000More$length_x, window=50, step=10)

# Plot density of G-runs in non-template strand
bm_1 <- plot_density_TSS7(result_window_1000_1000More[c("length_x", "genes_promoter_df", "CC_densityOE", "GG_densityOE")], a_site, "G-run frequency in non-template strand", result_window_1000_1000More$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type != "Oncogene"], "", genes_infor, result_window_1000_1000More$genes_promoter_df) +
  xlab("Distance from TSS (bp)") +
  theme_hd() +
  scale_color_manual(values=c("#40A4DE", "#F07F16C7"), labels=c("Template ", "Non-template")) +
  ylab("Ratio of G-runs count to G Count") +
  ggtitle("G-runs: two or more Gs") +
  scale_fill_manual(values=c("#40A4DE", "#F07F16C7"))

# Plot density of G-runs in non-template strand (three or more Gs)
bm_2 <- plot_density_TSS7(result_window_1000_1000More[c("length_x", "genes_promoter_df", "CCC_densityOE", "GGG_densityOE")], a_site, "G-run frequency in non-template strand", result_window_1000_1000More$genes_promoter_df$gene_id %in% genes_infor$ensembl_gene_id[genes_infor$Gene.Type != "Oncogene"], "", genes_infor, result_window_1000_1000More$genes_promoter_df) +
  xlab("Distance from TSS (bp)") +
  theme_hd() +
  scale_color_manual(values=c("#40A4DE", "#F07F16C7"), labels=c("Template ", "Non-template")) +
  ylab("Ratio of G-runs count to G Count") +
  ggtitle("G-runs: three or more Gs") +
  scale_fill_manual(values=c("#40A4DE", "#F07F16C7"))

# Load random frequency data
pos_all_downrandom <- readRDS("./pos_alldownrandom.rds")

# Calculate observed and random nucleotide frequencies
pos_all_down <- getFreq(gene_pos_data, 0, 1000, Hsapiens, "+")
pos_all_up <- getFreq(gene_pos_data, -1000, 0, Hsapiens, "+")

# Calculate G-runs density ratios
genes_infor$G_ensity <- pos_all_down[genes_infor$ensembl_gene_id, "G_density"]
genes_infor$GG_G_ratio <- genes_infor$GG_density / genes_infor$G_ensity
genes_infor$C_ensity <- pos_all_down[genes_infor$ensembl_gene_id, "C_density"]
genes_infor$CC_C_ratio <- genes_infor$CC_density / genes_infor$C_ensity

# Calculate G-runs length
genes_infor$GG_length <- pos_all_down[genes_infor$ensembl_gene_id, "GG_length"]
genes_infor$CC_length <- pos_all_down[genes_infor$ensembl_gene_id, "CC_length"]
genes_infor$GGG_density<-pos_all_down[genes_infor$ensembl_gene_id,"GGG_density"]
genes_infor$CCC_density<-pos_all_down[genes_infor$ensembl_gene_id,"CCC_density"]

# Plot age-specific G-runs density ratios
p3 <- plotSubSigParied(genes_infor[!is.na(genes_infor$age_type2) & genes_infor$Gene.Type == "Non-Cancer", ], "CCC_density", "GGG_density", "age_type2", "", cols1=c("#40A4DE", "#F07F16C7"), "Template ", "Non-template", "Age groups (Million years)", "G-run frequency in TSS downstream", cols_x="black") +
  ggtitle("G-run: three or more guanine") +
  labs(fill="Strand")

# Plot age-specific G-runs length
p_len <- plotSubSigParied(genes_infor[!is.na(genes_infor$age_type2) & genes_infor$Gene.Type == "Non-Cancer", ], "CC_length", "GG_length", "age_type2", "", cols1=c("#40A4DE", "#F07F16C7"), "Template ", "Non-template", "Age groups (Million years)", "Average size of G-runs", cols_x="black") +
  labs(fill="Strand") +
  ggtitle("TSS downstream")

# Plot age-specific G-runs density ratios
p_ratio <- plotSubSigParied(genes_infor[!is.na(genes_infor$age_type2) & genes_infor$Gene.Type == "Non-Cancer", ], "CC_C_ratio", "GG_G_ratio", "age_type2", "", cols1=c("#40A4DE", "#F07F16C7"), "Template ", "Non-template", "Age groups (Million years)", "Ratio of G-runs count to G count", cols_x="black") +
  labs(fill="Strand") +
  ggtitle("TSS downstream")

# Prepare data for nucleotide frequency plot
x <- rbind(cbind(type="TSS_upstream", melt(pos_all_up[1:13])), cbind(type="TSS_downstream", melt(pos_all_down[1:13])))
x$variable <- gsub("_density", "", as.character(x$variable))
x$type <- factor(x$type, levels=c("TSS_upstream", "TSS_downstream"))
x$Gene.Type <- genes_infor[x$Group.1, "Gene.Type"]

# Plot nucleotide frequency
p44 <- ggplot(x[x$type == "TSS_downstream" & x$variable %in% c("G", "C"), ], aes(x=variable, y=value)) +
  geom_signif(comparisons=getComp(c("G", "C")), map_signif_level=TRUE, col="black", test='wilcox.test', step_increase=0.1, tip_length=0.01) +
  geom_boxplot(alpha=I(0.7), outlier.colour=NA, aes(fill=variable)) +
  ylab("Nucleotide Frequency") +
  theme_hd_minimal() +
  guides(fill="none") +
  scale_fill_manual(values=c("#7FFF00", "#458B00")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1)) +
  ggtitle("TSS downstream")

# Plot age-specific nucleotide density
pt <- ggplot(genes_infor[!is.na(genes_infor$age_type2), ], aes(x=age_type2, y=CC_density, fill=age_type2)) +
  geom_boxplot(show.legend=FALSE) +
  geom_signif(comparisons=getComp(levels(genes_infor$age_type2[!is.na(genes_infor$age_type2)])), test='wilcox.test', map_signif_level=TRUE, col="black", step_increase=0.1) +
  theme_hd() +
  scale_fill_manual(values=colorRampPalette(c("white", "#1C84C9"))(6)[2:6]) +
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  ylab("G-run frequency\nin template strand") +
  xlab("Age groups (Million years)") +
  ggtitle("All genes")
gene_age_infor <- read.delim("./data/gene_age_infor.txt")
rownames(gene_age_infor)<-gene_age_infor$Age.class
genes_infor$age_time<-gene_age_infor[genes_infor$Age_class,4]
genes_infor$age_time2<-factor(paste(genes_infor$age_time,"",sep=""),levels=paste(sort(unique(genes_infor$age_time)),"",sep=""))
genes_infor$age_time2_num <- as.numeric(as.character(genes_infor$age_time2))


kendall_cor <- cor.test(genes_infor$age_time2_num, genes_infor$GG_density, method = "kendall")


gene_counts <- table(genes_infor$age_time2)


gene_counts_df <- data.frame(
  age_time2 = names(gene_counts),
  Gene_Count = as.numeric(gene_counts)
)


p <- ggplot(genes_infor[!is.na(genes_infor$Age_class), ], 
            aes(x = age_time2, y = GG_density)) +
  geom_boxplot() +
  theme_hd() +
  xlab("Gene ages: Million years") +
  ylab("G-run frequency on non-template strand") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = -Inf, y = Inf, label = paste("τ =", round(kendall_cor$estimate, 2), 
                                                    "\tp < 2.2e-16"), 
           hjust = -0.2, vjust = 1.2, size = 4, color = "black") +
  geom_text(data = gene_counts_df, aes(x = age_time2, y = 0.001, label = Gene_Count),
            vjust = 0.5, size = 3.5, color = "black", angle = 90) +
  geom_line(data = gene_counts_df, aes(x = age_time2, y = Gene_Count / 18000, group = 1), 
            color = "lightgray", size = 0.5) +
  geom_point(data = gene_counts_df, aes(x = age_time2, y = Gene_Count / 18000), 
             color = "red", size = 1) +
  scale_y_continuous(
    name = "G-run frequency on non-template strand",
    sec.axis = sec_axis(~ . * 18000, name = "Gene Count")
  )+geom_vline(xintercept = c(7,14,18,23),linetype="dashed",color="blue",size=0.5)

# 打印图形
print(p)
pnt=ggplot(genes_infor[!is.na(genes_infor$age_type2),],aes(x=age_type2,y=GG_density,fill=age_type2))+geom_boxplot(show.legend = FALSE)+geom_signif(comparisons=getComp(levels(genes_infor$age_type2[!is.na(genes_infor$age_type2)])), test='wilcox.test', map_signif_level=TRUE,col="black",step_increase=0.1) +theme_hd()+scale_fill_manual(values=colorRampPalette(c("white" ,"#CF8310"))(6)[2:6])+theme(axis.text.x = element_text(angle = 30,hjust = 1))+ylab("G-run frequency\non non-template strand")+xlab("Age groups(Million years)")+ggtitle("All genes")
# Arrange plots
p_age_all0 <- ggarrange(p44, bm_1, bm_2, p, nrow=1, ncol=4, labels=c(letters[1:4]), font.label=list(size=18, color="black", face="bold"), widths=c(1.2, 2, 2, 3.3))
p_age_all1 <- ggarrange(p_len, pt, pnt, p_ratio, nrow=1, ncol=4, labels=c(letters[5:8]), font.label=list(size=18, color="black", face="bold"))
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure1.tiff", width=18, height=16, res=300, units="in", compression="lzw")
print(ggarrange(p_age_all0, p_age_all1, nrow=2, ncol=1))
dev.off()
