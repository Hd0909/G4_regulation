# Set the working directory to the specified path
setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new")

# Load necessary functions and data
source('functions/functions_evolution_5UTR.r')  # Load custom functions
load(file="./data/all_meth_exp_signal_promoter.rdata")  # Load methylation and expression data
load(file="./data/predict_G4_up_down.rdata")  # Load G4 prediction data

# Load and process data for G-run analysis
all_cg <- readRDS(file="./data/aver_cg_multiple_pos_continues2.rds")  # Load G-run data
cor_dhs <- readRDS(file="./data/aver_dhs_multiple_pos_cor_all.rds")  # Load DHS correlation data
cor_meth <- readRDS(file="./data/aver_meth_multiple_pos_cor_all.rds")  # Load methylation correlation data
cor_exp <- readRDS(file="./data/aver_exp_multiple_pos_cor_all.rds")  # Load expression correlation data

# Extract and process position information for G-run data
all_cg$pos1 <- str_split_fixed(all_cg$pos, " ", 2)[, 1]  # Extract the first part of the position string
temp1 <- str_split_fixed(all_cg$pos1, "_", 2)  # Split the position string by underscore
neg <- ifelse(as.numeric(temp1[, 2]) > as.numeric(temp1[, 1]), 1, -1)  # Determine the direction (positive or negative)
all_cg$position <- (as.numeric(temp1[, 1]) / 2 + as.numeric(temp1[, 2]) / 2 + 50) * neg  # Calculate the midpoint position

# Melt the data for plotting
all_cg_melt <- melt(all_cg, id.vars = c("pos1", "pos", "position", "ensembl_gene_id", "ensembl_trans_id"))
all_cg_melt$value <- as.numeric(all_cg_melt$value)  # Convert the value column to numeric

# Plot the frequency of G-runs
p1_cg_sig <- plotSubSigcg(all_cg_melt[all_cg_melt$position <= 1000 & all_cg_melt$position >= -1000, ], 
                          "position", "value", "variable", "", 
                          cols1 = c("#F07F16C7", "#40A4DE"), "", "G-run Frequency") +
  labs(fill = "Gene.Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  xlab("Distance from TSS (100bp segments)")

# Plot the correlation of expression and G-run frequency
p_cor_exp <- ggplot(cor_exp[cor_exp$position_gg <= 1000 & cor_exp$position_gg >= -1000, ], 
                    aes(x = as.factor(position_gg), y = cor)) +
  geom_boxplot(aes(fill = type_x, alpha = I(0.1)), show.legend = FALSE) +
  geom_jitter(width = 0.25, aes(col = type_x)) +
  theme_hd_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  ylab("Spearman's correlation coefficient") +
  xlab("Distance from TSS (100bp segments)") +
  ggtitle("Expression  VS Gruns in different regions") +
  scale_color_manual(values = c("#F07F16C7", "#40A4DE")) +
  labs(col = "Strands")

# Plot the correlation of promoter methylation and G-run frequency
p_cor_meth <- ggplot(cor_meth[cor_meth$position_gg <= 1000 & cor_meth$position_gg >= -1000, ], 
                     aes(x = as.factor(position_gg), y = cor)) +
  geom_boxplot(aes(fill = type_x, alpha = I(0.1)), show.legend = FALSE) +
  geom_jitter(width = 0.25, aes(col = type_x)) +
  theme_hd_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  ylab("Spearman's correlation coefficient") +
  xlab("Distance from TSS (100bp segments)") +
  ggtitle("Promoter methylation vs. G-run Frequency in each segment") +
  scale_color_manual(values = c("#F07F16C7", "#40A4DE")) +
  labs(col = "Strands")

# Plot the correlation of promoter DHS signal and G-run frequency
p_cor_dhs <- ggplot(cor_dhs[cor_dhs$position_gg <= 1000 & cor_dhs$position_gg >= -1000, ], 
                    aes(x = as.factor(position_gg), y = cor)) +
  geom_boxplot(aes(fill = type_x, alpha = I(0.1)), show.legend = FALSE) +
  geom_jitter(width = 0.25, aes(col = type_x)) +
  theme_hd_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  ylab("Spearman's correlation coefficient") +
  xlab("Distance from TSS (100bp segments)") +
  ggtitle("Promoter DHS signal vs. G-run Frequency in each segment") +
  scale_color_manual(values = c("#F07F16C7", "#40A4DE")) +
  labs(col = "Strands")

# Arrange all plots into a single figure
pall4 <- ggarrange(p1_cg_sig + guides(color = "none"), p_cor_exp, p_cor_meth, p_cor_dhs, 
                   nrow = 2, ncol = 2, 
                   labels = c(letters[1:4]), 
                   font.label = list(size = 18, color = "black", face = "bold"), 
                   common.legend = TRUE, legend = "bottom")

# Save the combined plot to a TIFF file
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/all_figure3.tiff", 
     width = 12, height = 12, res = 300, units = "in", compression = "lzw")
print(pall4)
dev.off()

# Repeat the process for data with a 50 bp window
all_cg50 <- readRDS(file="./data/aver_cg_multiple_pos_continues50.rds")  # Load G-run data with 50 bp window
cor_dhs50 <- readRDS(file="./data/aver_dhs_multiple_pos_cor_all50.rds")  # Load DHS correlation data with 50 bp window
cor_meth50 <- readRDS(file="./data/aver_meth_multiple_pos_cor_all50.rds")  # Load methylation correlation data with 50 bp window
cor_exp50 <- readRDS(file="./data/aver_exp_multiple_pos_cor_all50.rds")  # Load expression correlation data with 50 bp window

# Extract and process position information for G-run data with 50 bp window
all_cg50$pos1 <- str_split_fixed(all_cg50$pos, " ", 2)[, 1]  # Extract the first part of the position string
temp150 <- str_split_fixed(all_cg50$pos1, "_", 2)  # Split the position string by underscore
neg50 <- ifelse(as.numeric(temp150[, 2]) > as.numeric(temp150[, 1]), 1, -1)  # Determine the direction (positive or negative)
all_cg50$position <- (as.numeric(temp150[, 1]) / 2 + as.numeric(temp150[, 2]) / 2 + 25) * neg50  # Calculate the midpoint position

# Melt the data for plotting
all_cg_melt50 <- melt(all_cg50, id.vars = c("pos1", "pos", "position", "ensembl_gene_id", "ensembl_trans_id"))
all_cg_melt50$value <- as.numeric(all_cg_melt50$value)  # Convert the value column to numeric

# Plot the frequency of G-runs for 50 bp window data
p1_cg_sig50 <- plotSubSigcg(all_cg_melt50[all_cg_melt50$position <= 1000 & all_cg_melt50$position >= -1000, ], 
                            "position", "value", "variable", "", 
                            cols1 = c("#F07F16C7", "#40A4DE"), "", "G-run Frequency") +
  labs(fill = "Gene.Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  xlab("Distance from TSS (100bp segments)")

# Plot the correlation of expression and G-run frequency for 50 bp window data
p_cor_exp50 <- ggplot(cor_exp50[cor_exp50$position_gg <= 1000 & cor_exp50$position_gg >= -1000, ], 
                      aes(x = as.factor(position_gg), y = cor)) +
  geom_boxplot(aes(fill = type_x, alpha = I(0.1)), show.legend = FALSE) +
  geom_jitter(width = 0.25, aes(col = type_x)) +
  theme_hd_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  ylab("Spearman's Correlation Coefficient") +
  xlab("Distance from TSS (100bp segments)") +
  ggtitle("Expression vs. G-run Frequency in each segment") +
  scale_color_manual(values = c("#F07F16C7", "#40A4DE")) +
  labs(col = "Strands")

# Plot the correlation of promoter methylation and G-run frequency for 50 bp window data
p_cor_meth50 <- ggplot(cor_meth50[cor_meth50$position_gg <= 1000 & cor_meth50$position_gg >= -1000, ], 
                       aes(x = as.factor(position_gg), y = cor)) +
  geom_boxplot(aes(fill = type_x, alpha = I(0.1)), show.legend = FALSE) +
  geom_jitter(width = 0.25, aes(col = type_x)) +
  theme_hd_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  ylab("Spearman's Correlation Coefficient") +
  xlab("Distance from TSS (100bp segments)") +
  ggtitle("Promoter methylation vs. G-run Frequency in each segment") +
  scale_color_manual(values = c("#F07F16C7", "#40A4DE")) +
  labs(col = "Strands")

# Plot the correlation of promoter DHS signal and G-run frequency for 50 bp window data
p_cor_dhs50 <- ggplot(cor_dhs50[cor_dhs50$position_gg <= 1000 & cor_dhs50$position_gg >= -1000, ], 
                      aes(x = as.factor(position_gg), y = cor)) +
  geom_boxplot(aes(fill = type_x, alpha = I(0.1)), show.legend = FALSE) +
  geom_jitter(width = 0.25, aes(col = type_x)) +
  theme_hd_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  ylab("Spearman's Correlation Coefficient") +
  xlab("Distance from TSS (100bp segments)") +
  ggtitle("Promoter DHS signal vs. G-run Frequency in each segment") +
  scale_color_manual(values = c("#F07F16C7", "#40A4DE")) +
  labs(col = "Strands")

# Arrange all plots into a single figure for 50 bp window data
pall4 <- ggarrange(p1_cg_sig50 + guides(color = "none"), p_cor_exp50, p_cor_meth50, p_cor_dhs50, 
                   nrow = 2, ncol = 2, 
                   labels = c(letters[1:4]), 
                   font.label = list(size = 18, color = "black", face = "bold"), 
                   common.legend = TRUE, legend = "bottom")

# Save the combined plot to a TIFF file for 50 bp window data
tiff("/media/huangdan/hardisk0/HD/my_data/HD_evolution/UTR5_new/all_figures/supplementary_figure7.tiff", 
     width = 17, height = 17, res = 300, units = "in", compression = "lzw")
print(pall4)
dev.off()