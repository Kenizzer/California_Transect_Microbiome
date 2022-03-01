# Analysis of grafted grapevines across time and space.
# Samples were collected in summer 2018 and 2019
# Samples consist of leaf, berry, and root compartments
# Code by: Joel F. Swift

# For the data filtering see: Data_filtering_normalization.R
# For general analysis and figures see: CA_transect_analysis.R
# For machine learning see: CA_transect_machine_learning.R
# For DEseq2 diff. abund. see: CA_transect_differential_abundance.R

#Packages w/ version numbers.
library('DESeq2'); packageVersion('DESeq2')
library('tidyverse'); packageVersion('tidyverse')
library('ggpubr'); packageVersion('ggpubr')
library('phyloseq'); packageVersion('phyloseq')
library('qiime2R'); packageVersion('qiime2R')
library('zinbwave'); packageVersion('zinbwave')
library('BiocParallel'); packageVersion('BiocParallel')

# Theme set and Color Palettes
theme_set(theme_pubr())
rootstock_palette <- c('#1b9e77', '#f0a4af', '#7570b3')
scion_palette <- c('#ed254e', '#0e79b2')
site_palette <- c('#e6ab02', '#281c39', '#12664c')
compartment_palette <- c("#5a1991", "#139d08", "#5c3c0d") #https://lospec.com/palette-list/famicube
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#AA4499", "#332288", "#117733", 
                             "#661100", "#999933", "#44AA99", "#882255", "#6699CC", "#888888")
#Functions

# Given a phyloseq object, return a dataframe of the sample metadata 
# From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}



# Negative bionomial from DESeq2 with pseudocount
# load non-normalized dataset generated from Data_filtering_normalization.R
phy_filt <- readRDS("../Data_files/phyloseq_16s_no_soil_filtered_dataset.rds")
# Remove taxa that are not found with greater than 25 read count in 25 samples 
phy_filt <- filter_taxa(phy_filt, function(x) sum(x > 25) > (25), TRUE) 
# In order for DESEQ2 to run there can be no zeros in the OTU matrix, so I add 1 to each count
otu_table(phy_filt) <- otu_table(phy_filt) + 1

# 2 Breaks Brix 
phy_filt@sam_data$brix_2_breaks <- cut(phy_filt@sam_data$brix, breaks = c(-Inf, 7, Inf), labels = c("Pre-ripening", "Ripening")) # Check if I want to use these terms
# Create Deseq2 object with study design formula
DS2_mod <- phyloseq_to_deseq2(phy_filt, ~ rootstock + scion + plant_body_site + brix_2_breaks + year + site)
#geoMeans <- apply(counts(DS2), 1, gm_mean)
#DS2 <- estimateSizeFactors(DS2, geoMeans = geoMeans)
DS2_mod <- estimateSizeFactors(DS2_mod)
DS2_mod <- estimateDispersions(DS2_mod, fitType = "local")
DS2_nb_mod <- nbinomWaldTest(DS2_mod, maxit = 500)
# Results table contrasts
resultsNames(DS2_nb_mod)
# Main effects
## rootstock
x <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name =  "rootstock_Freedom_vs_1103.Paulsen")),padj<0.05)
y <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name =  "rootstock_Teleki.5C_vs_1103.Paulsen")),padj<0.05)
z <- as.data.frame(subset(results(DS2_nb_mod, contrast=c('rootstock', 'Freedom', 'Teleki.5C'), lfcThreshold=0, alpha=0.05), padj<0.05))
# Add ASV column
x$ASV <- rownames(x)
y$ASV <- rownames(y)
z$ASV <- rownames(z)
R_lf2c <- rbind(x,y,z)
R_lf2c <- R_lf2c[!duplicated(R_lf2c$ASV),]
## compartment
x <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name =  "plant_body_site_leaf_vs_berry")),padj<0.05)
y <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name =  "plant_body_site_root_vs_berry")),padj<0.05)
z <- as.data.frame(subset(results(DS2_nb_mod, contrast=c('plant_body_site', 'root', 'leaf'), lfcThreshold=0, alpha=0.05), padj<0.05))
# Add ASV column
x$ASV <- rownames(x)
y$ASV <- rownames(y)
z$ASV <- rownames(z)
C_lf2c <- rbind(x,y,z)
C_lf2c <- C_lf2c[!duplicated(C_lf2c$ASV),]
## brix
B_lf2c <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name ="brix_2_breaks_Ripening_vs_Pre.ripening")),padj<0.05)
## scion
Sc_lf2c <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name ="scion_chardonnay_vs_cabernet.sauvignon")),padj<0.05)
## year
Y_lf2c <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name ="year_2019_vs_2018")),padj<0.05)
## site
x <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name =  "site_livingston_vs_liberty")),padj<0.05)
y <- subset(as.data.frame(results(DS2_nb_mod, cooksCutoff = FALSE, name =  "site_ripperdan_vs_liberty")),padj<0.05)
z <- as.data.frame(subset(results(DS2_nb_mod, contrast=c('site', 'livingston', 'ripperdan'), lfcThreshold=0, alpha=0.05), padj<0.05))
# Add ASV column
x$ASV <- rownames(x)
y$ASV <- rownames(y)
z$ASV <- rownames(z)
Si_lf2c <- rbind(x,y,z)
Si_lf2c <- Si_lf2c[!duplicated(Si_lf2c$ASV),]

# Extract Log2FoldChange column, Main Effects
R_lf2c <- as.data.frame(abs(R_lf2c$log2FoldChange))
Sc_lf2c <- as.data.frame(abs(Sc_lf2c$log2FoldChange))
C_lf2c <- as.data.frame(abs(C_lf2c$log2FoldChange))
B_lf2c <- as.data.frame(abs(B_lf2c$log2FoldChange))
Y_lf2c <- as.data.frame(abs(Y_lf2c$log2FoldChange))
Si_lf2c <- as.data.frame(abs(Si_lf2c$log2FoldChange))
# Rename 1st column to Log2FoldChange
colnames(R_lf2c)[1] <- "Log2FoldChange"
colnames(Sc_lf2c)[1] <- "Log2FoldChange"
colnames(C_lf2c)[1] <- "Log2FoldChange"
colnames(B_lf2c)[1] <- "Log2FoldChange"
colnames(Y_lf2c)[1] <- "Log2FoldChange"
colnames(Si_lf2c)[1] <- "Log2FoldChange"
# Add factor column
R_lf2c["Factor"] = "Rootstock"
Sc_lf2c["Factor"] = "Scion"
C_lf2c["Factor"] = "Compartment"
B_lf2c["Factor"] = "Sugar Content"
Y_lf2c["Factor"] = "Year"
Si_lf2c["Factor"] = "Site"

# Number of differential abundant ASVS
length(R_lf2c$Log2FoldChange) # 676
length(Sc_lf2c$Log2FoldChange) # 464
length(C_lf2c$Log2FoldChange) # 903
length(B_lf2c$Log2FoldChange) # 322
length(Y_lf2c$Log2FoldChange) # 411
length(Si_lf2c$Log2FoldChange) # 829

# Mean and SE of absolute log2fold change of significant ASVs
mean(R_lf2c$Log2FoldChange) # 0.886
mean(Sc_lf2c$Log2FoldChange) # 0.723
mean(C_lf2c$Log2FoldChange) # 2.515
mean(B_lf2c$Log2FoldChange) # 0.777
mean(Y_lf2c$Log2FoldChange) # 0.864
mean(Si_lf2c$Log2FoldChange) # 1.278
se <- function(x) sqrt(var(x)/length(x))
se(R_lf2c$Log2FoldChange)
se(Sc_lf2c$Log2FoldChange)
se(C_lf2c$Log2FoldChange)
se(B_lf2c$Log2FoldChange)
se(Y_lf2c$Log2FoldChange)
se(Si_lf2c$Log2FoldChange)
# Bind by row all dfs
Log2fold_by_factor_DS2_NB.DE <- rbind(R_lf2c, Sc_lf2c, C_lf2c, B_lf2c, Y_lf2c, Si_lf2c)
# Tukey test for significant difference in means 
TukeyHSD(aov(Log2FoldChange ~ Factor, data = Log2fold_by_factor_DS2_NB.DE), conf.level = 0.95)
# Violin plot
# Make dataframe to hold significance letters and plot
labels_df <- tibble(Factor = levels(as.factor(Log2fold_by_factor_DS2_NB.DE$Factor)), Mlog2=max(Log2fold_by_factor_DS2_NB.DE$Log2FoldChange) * 1.4)
Log2fold_violin_plot_2_breaks <- ggplot(Log2fold_by_factor_DS2_NB.DE, aes(x=Factor, y=Log2FoldChange)) + 
  geom_violin(trim = FALSE, scale = "width", fill = "grey") +
  ylim(NA,10) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", position = position_dodge(width = 0.9)) +
  xlab("Factor") + labs(y=expression(paste(Log[2]," fold change")), x=("Source of variation"))
# Add significance letters after in order to preserve stat_summary
Log2fold_violin_plot <- Log2fold_violin_plot_2_breaks + geom_text(data=labels_df, aes(x = Factor, y = Mlog2, label=c("a","b","c","d", "bc", "bc")), size = 6)
ggsave("Log2fold_violin_plot.png", Log2fold_violin_plot, height = 6, width = 10, bg = 'white', dpi = 1200)
ggsave("Log2fold_violin_plot.pdf", Log2fold_violin_plot, height = 6, width = 10)
ggsave("Log2fold_violin_plot.svg", Log2fold_violin_plot, height = 6, width = 10)
