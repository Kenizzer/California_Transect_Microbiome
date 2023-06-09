# Processing of soil samples
# Code by: Joel F. Swift

# For the data filtering see: Data_filtering_normalization.R
# For general analysis and figures see: CA_transect_analysis.R
# For machine learning see: CA_transect_machine_learning.R
# For DEseq2 diff. abund. see: CA_transect_differential_abundance.R

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('ggpubr'); packageVersion('ggpubr')
library('car'); packageVersion('car')
library('vegan'); packageVersion('vegan')
library('phyloseq'); packageVersion('phyloseq')
library('qiime2R'); packageVersion('qiime2R')
library('emmeans'); packageVersion('emmeans')
library('lubridate'); packageVersion('lubridate')
library('rebus'); packageVersion('rebus')
library('viridis'); packageVersion('viridis')
library('corrplot'); packageVersion('corrplot')
library('factoextra'); packageVersion('factoextra')
library('data.table'); packageVersion('data.table')
library('patchwork'); packageVersion("patchwork")

# Theme set and Color Palettes
theme_set(theme_pubr())
rootstock_palette <- c('#1b9e77', '#f0a4af', '#7570b3')
scion_palette <- c('#ed254e', '#0e79b2')
site_palette <- c('#e6ab02', '#281c39', '#12664c')
compartment_palette <- c("#5a1991", "#139d08", "#5c3c0d") #https://lospec.com/palette-list/famicube
safe_colorblind_palette <- c("Acidobacteriota" = "#88CCEE", "Actinobacteriota" = "#CC6677",
                             "Bacteroidota" = "#DDCC77", "Chloroflexi" = "#AA4499",
                             "Deinococcota" = "#332288","Firmicutes" = "#117733",
                             "Myxococcota" = "#661100", "Planctomycetota" = "#999933",
                             "Proteobacteria" = "#44AA99", "Verrucomicrobiota" = "#882255",
                             "Desulfobacterota" = "#888888","Crenarchaeota" = "#D55E00",
                             "Other" = "#6699CC") 
# Functions
# Given a phyloseq object, return a dataframe of the sample metadata 
# From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

PLOT_PCoA <- function(plot_data, distance_matrix, axis1, axis2, split_by){
  if (split_by == 'site'){
    temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
    temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
    Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
      geom_point(aes(fill=lower_political, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
      scale_shape_manual(values=c(21, 24, 23, 22)) +  # Modified to work with soil (4 shapes not 3)
      scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
      labs(shape= "Compartment", color= "Site") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
      guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1))) +
      theme(legend.position="right")
    return(Plot)
  } else if(split_by == 'rootstock'){
    temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
    temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
    Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
      geom_point(aes(fill=rootstock, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
      scale_shape_manual(values=c(22, 24, 21, 23)) +  # Modified to work with soil (4 shapes not 3)
      scale_fill_manual(name = "Rootstock", values=rootstock_palette) +
      labs(shape= "Compartment", color= "Site") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
      guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))
    return(Plot)
  } else if(split_by == 'year'){
    temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
    temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
    Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
      geom_point(aes(fill=year, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
      scale_shape_manual(values=c(22, 24, 21, 23)) +  # Modified to work with soil (4 shapes not 3)
      labs(shape= "Compartment", color= "Site") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
      guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))
    return(Plot)
  }    else if(split_by == 'block'){
    temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
    temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
    Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
      geom_point(aes(fill=block, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
      scale_shape_manual(values=c(22, 24, 21, 23)) +  # Modified to work with soil (4 shapes not 3)
      labs(shape= "Compartment", color= "Site") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
      guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "black")))
    return(Plot)
  }
}

is_gt <- function(object, dist, threshold){
  samples <- rownames(object)[dist > threshold | dist < -threshold]
  return(samples)
}

# Function from the R cookbook
# From: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
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
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

# load dataset generated from Data_filtering_normalization.R
phy_soil_vst <- readRDS("../../Data_files/phyloseq_16s_filtered_vst_dataset.rds")
phy_soil_only_vst <- subset_samples(phy_soil_vst, plant_body_site=="soil")



# Figure 2 panel A
Full_df <- readRDS("dataframe_for_analysis.RDS")
Full_df$Year <-c("2018", "2018", "2018", "2018", "2018", "2018", "2018", "2018", "2018", "2018", "2018", "2019", "2019", "2019", "2019", "2019", "2019", "2019", "2019", "2019", "2019", "2019", "2019", "2019", "2019")
#Rename the rownames to first column called Sample
setDT(Full_df, keep.rownames = "Sample")
# Make the data long form for GGplot
Full_df_long <- reshape2::melt(Full_df)
Full_df_long$Site <- as.factor(Full_df_long$Site)
Full_df_long$Site <- recode_factor(Full_df_long$Site, Ripperdan = "Madera", Livingston ="Merced", Liberty = "San Joaquin")
# Generate plot for means of sand/silt/clay by site irrespective of year.
# Generate means
x1 <- aggregate(Clay ~ Site, Full_df, mean)
x2 <- aggregate(Silt ~ Site, Full_df, mean)
x3 <- aggregate(Sand ~ Site, Full_df, mean)
# Create the dataframe
Site_means <- as.data.table(cbind(x1, Silt=x2$Silt, Sand=x3$Sand))
Site_means$Site <- c("San Joaquin", "Merced", "Madera")
# Make the data long form for GGplot2
Site_means_long <- melt(Site_means)
# PLOT
Panel_A_stacked <- ggplot(Site_means_long, aes(fill=variable, y=value, x=Site)) +
  geom_bar(position="fill", stat="identity", color = "black") +
  theme(legend.title = element_blank(), legend.position = "right") +
  ylab("Proportion") +
  xlab(element_blank())

Panel_A <- ggplot(Full_df_long, aes(fill=variable, y=value, x=Sample)) +
  geom_bar(position="fill", stat="identity", color = "black") +
  theme(legend.title = element_blank(), legend.position = "right") +
  ylab("Proportion") +
  xlab(element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(~Site, scales = "free")  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# Figure 2 panel B
Ark <- as.data.frame(read_tsv("Soil_elements.tsv", col_names = TRUE))
# Remove Mount Vernon samples
Ark_CA <- Ark[10:37,]
# PCAs With Mt Vernon DDPSC ionomics
d <- Ark_CA %>% dplyr::select(pH:B)
rownames(d) <- Ark_CA$sample
d <- as.data.frame(scale(d, scale=T, center=T))
x <- apply(d, 2, is_gt, object=Ark_CA, threshold=5)
outliers <- unique(unlist(x))
data_clean <- Ark_CA[!(Ark_CA$sample %in% outliers),]
d_clean <- data_clean %>% dplyr::select(pH:B)
d_cleanScale <- scale(d_clean, scale=T, center=T)
pca <- prcomp(d_cleanScale)
x <- summary(pca)

# PLOT
Panel_B <- fviz_pca_biplot(pca, repel = TRUE, geom = "point", fill.ind = Ark_CA$Site, pointshape = 21, pointsize = 5, alpha.ind = 0.8, col.var = "black",
                label = "var", palette = site_palette,invisible = "quali", title = NULL) +
                labs(color = "Site", fill = "Site") + xlab("PC1 (52.6%)") + ylab("PC2 (13.8%)")


# Figure 2 panel C
# Get data from phyloseq object.
phylum_lvl <- tax_glom(phy_soil_only_vst, taxrank = "Phylum") # 42 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
phylum_lvl <- transform_sample_counts(phylum_lvl, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
phylum_lvl <- fantaxtic::get_top_taxa(physeq_obj = phylum_lvl, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Melt data frame with phyloseq function for plotting.
phylum_lvl <- psmelt(phylum_lvl)

# Reorder levels to put other to the end, other in alphabetical
phylum_lvl$Phylum <- forcats::fct_relevel(as.factor(phylum_lvl$Phylum), "Other", after = Inf)
summary(phylum_lvl$Phylum)
# PLOT
levels(phylum_lvl$lower_political) <- c("Madera", "Merced", "San Joaquin")
Panel_C <- ggplot(phylum_lvl, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(position="fill", stat = "identity") +
  ylab("Relative abundance") +
  facet_grid(~lower_political, scale = "free") +
  scale_fill_manual(values=safe_colorblind_palette)  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')
# Stacked barplot
summary(phy_soil_only_vst@sam_data)
# madera 8, merced 8, san joaquin 12
site_phylum <- phylum_lvl
site_phylum[site_phylum$lower_political == "Madera",]$Abundance <- site_phylum[site_phylum$lower_political == "Madera",]$Abundance / 8
site_phylum[site_phylum$lower_political == "Merced",]$Abundance <- site_phylum[site_phylum$lower_political == "Merced",]$Abundance / 8
site_phylum[site_phylum$lower_political == "San Joaquin",]$Abundance <- site_phylum[site_phylum$lower_political == "San Joaquin",]$Abundance / 12
#PLOT
Panel_C_stacked <- ggplot(site_phylum, aes(x=lower_political, y=Abundance, fill=Phylum)) +
                          geom_bar(stat="identity") +
                          scale_fill_manual(values=safe_colorblind_palette) + 
                          theme(legend.text.align = 0) +
                          theme(legend.position = 'right')

# Figure 2 panel D
# PCoA W/ only soil samples 
# Calculate bray-curtis
out.bray <- ordinate(phy_soil_only_vst, method = "MDS", distance = "bray")
# Calculate distance matrixes for use in Vegan
out.dist.bray <- phyloseq::distance(phy_soil_only_vst, method = "bray")
# get variance explained by axes 1-3
sum(out.bray$values$Eigenvalues[1:3])/sum(out.bray$values$Eigenvalues) #43%
# Plot PCoAs 1x2-3 for bray by site
Panel_D <- PLOT_PCoA(phy_soil_only_vst, out.bray, 1, 2, split_by = 'site')
#pcoa_bray_1_3 <- PLOT_PCoA(phy_soil_only_vst, out.bray, 1, 3, split_by = 'site')
#Bray123_pcoa <- ggarrange(pcoa_bray_1_2, pcoa_bray_1_3, common.legend = TRUE, legend = "right", labels = c("Bray-Curtis Dissimilarity"," "))


# Assembled plot
#ggarrange(Panel_A, Panel_B, Panel_C, Panel_D, labels = "AUTO") # This doesnt work to align the plots
FigureX_multipanel <- (Panel_A | Panel_B) / (Panel_C | Panel_D) + plot_annotation(tag_levels = 'A')

ggsave("Figure_2_soil.pdf", FigureX_multipanel, dpi = 600, width = 16, height = 10)
ggsave("Figure_2_soil.svg", FigureX_multipanel, dpi = 600, width = 16, height = 10)


# Linear modeling
# Soil texture
# lms for % sand silt clay
silt_model <- lm(Silt ~ Site + Year, Full_df)
sand_model <- lm(Sand ~ Site + Year, Full_df)
clay_model <- lm(Clay ~ Site + Year, Full_df)
# anova
Anova(silt_model, type = "II")
Anova(sand_model, type = "II")
Anova(clay_model, type = "II")
# posthocs
as.data.frame(pairs(emmeans(clay_model, ~Site)))

summarySE(Full_df, measurevar = 'Clay', groupvars = 'Site' )
# Soil elemental composition
Soil_ele_lm.df <- cbind(as.data.frame(pca$x), Ark_CA)
# lms
Soil_ele_pc1_model <- lm(PC1 ~ Site + Year, Soil_ele_lm.df)
Soil_ele_pc2_model <- lm(PC2 ~ Site + Year, Soil_ele_lm.df)
# anova
Anova(Soil_ele_pc1_model, type = "II")
Anova(Soil_ele_pc2_model, type = "II")
# posthocs
as.data.frame(pairs(emmeans(Soil_ele_pc1_model, ~Site)))

# Alpha diversity
Alpha_div_metrics <- function(phyloseq_obj){
  SH <- as.data.frame(vegan::diversity(t(otu_table(phyloseq_obj)), index="shannon"))
  IS <- as.data.frame(vegan::diversity(t(otu_table(phyloseq_obj)), index="invsimpson"))
  CH <- as.data.frame(t(estimateR(round(t(otu_table(phyloseq_obj))))))
  CH <- as.data.frame(CH$S.chao1)
  alpha_div <- cbind(SH, IS, CH)
  colnames(alpha_div) <-c("Shannon", "Invsimpson", "Chao1")
  R <- picante::pd(samp = t(otu_table(phyloseq_obj)), tree = phy_tree(phyloseq_obj), include.root = TRUE)
  alpha_div <- cbind(alpha_div, "Faithpd" = R$PD, as.data.frame(phyloseq_obj@sam_data))
  return(alpha_div)
}

Alpha_soil.df <- Alpha_div_metrics(phy_soil_only_vst)
# lms
Alpha_soil_faith_model <- lm(Faithpd ~ lower_political + year, Alpha_soil.df)
Alpha_soil_chao1_model <- lm(Chao1 ~ lower_political + year, Alpha_soil.df)
# anova
Anova(Alpha_soil_faith_model, type = "II")
Anova(Alpha_soil_chao1_model, type = "II")
# posthocs
as.data.frame(pairs(emmeans(Alpha_soil_faith_model, ~ year)))
as.data.frame(pairs(emmeans(Alpha_soil_chao1_model, ~ year)))


# beta diversity
### Adonis2 ###
# Create sample data dataframe from phyloseq object to use in vegan/adonis
phy_vst_soil_only_metadata <- pssd2veg(phy_soil_only_vst)
adonis2(out.dist.bray ~ site + year, data = phy_vst_soil_only_metadata, by = "margin", permutations = 1000)

# Testing PC1 (soil elemental comp) with beta diversity
metadata_w_PCs <- cbind(phy_vst_soil_only_metadata, pca$x)
adonis2(out.dist.bray ~ PC1 + PC2 + year, data = metadata_w_PCs, by = "margin", permutations = 10000) 
# PC1 non-significant p = 0.083
# PC2 non-significant p = 0.095

# Testing individual elements and pH with beta diversity
metadata_w_elements <- cbind(phy_vst_soil_only_metadata, Ark_CA)
adonis2(out.dist.bray ~ pH + P + K + Ca + Mg + S + Na + Fe + Mn + Zn + Cu + B + year, data = metadata_w_elements, by = "margin", permutations = 1000) 
# only Mg is significant, Mg loads strongly negative on PC1 and very little on to PC2