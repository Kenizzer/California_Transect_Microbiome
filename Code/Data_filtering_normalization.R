# Data pre-processing, filtering and normalization
# Code by: Joel F. Swift

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('phyloseq'); packageVersion('phyloseq')
library('qiime2R'); packageVersion('qiime2R')
library('DESeq2'); packageVersion('DESeq2')
library('decontam'); packageVersion('decontam')
library('lubridate'); packageVersion('lubridate')

# Load QIIME2 objects into a phyloseq class object.
phyloseq_16s <- qza_to_phyloseq(features = '../Qiime_files/table-no-mitochondria-no-chloroplast.qza', tree = '../Qiime_files/rooted-tree.qza', taxonomy = '../Qiime_files/taxonomy.qza', metadata = '../Metadata/Metadata_all_plates.tsv')

# Sanity check make sure factors and levels are correct.
summary(phyloseq_16s@sam_data)

#### Use decontam to remove contaminant ASVs and then remove control samples ####
sample_data(phyloseq_16s)$is.neg <- sample_data(phyloseq_16s)$Sample_or_Control == "Control" # Column in metadata to denote controls and real samples
contamdf.prev <- isContaminant(phyloseq_16s, method="prevalence", neg="is.neg", threshold=0.5) # default threshold
table(contamdf.prev$contaminant) # 183 flagged as contaminants
To_remove_ASVs <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])
Taxa <- taxa_names(phyloseq_16s)
Taxa <- Taxa[!(Taxa %in% To_remove_ASVs)]
phyloseq_16s <- prune_taxa(Taxa, phyloseq_16s)
phyloseq_16s #Sanity check started with 45515 ASVs - 183 contaminates = 45332
rm(contamdf.prev, Taxa)
# Remove control samples (and soil)
controls_to_remove <- c('NEG1', 'NEG2', 'NEG3', 'NEG4', 'NEG5', 'NEG6', 'Neg7', 'Neg8', 'Neg9', 'Neg10', 'Neg11', 'Neg12', 'Neg13', 'Neg14', 'POS1', 'POS2', 'POS3', 'Pos4', 'Pos5', 'Pos6', 'Pos7')
soil_to_remove <- c('CA.soil.1', 'CA.soil.2', 'CA.soil.3', 'CA.soil.4', 'CA.soil.5', 'CA.soil.6', 'CA.soil.7', 'CA.soil.8', 'CA.soil.9', 'CA.soil.10', 'CA.soil.11', 'CA.soil.12', 'CA.soil.13', 'CA.soil.14', 'CA.soil.15', 'CA.soil.16', 'CA.soil.17', 'CA.soil.18', 'CA.soil.19', 'CA.soil.20', 'CA.soil.21', 'CA.soil.22', 'CA.soil.23', 'CA.soil.24', 'CA.soil.25', 'CA.soil.26', 'CA.soil.27', 'CA.soil.28')
Sample_names <- sample_names(phyloseq_16s) # vector of samples
Sample_names <- Sample_names[!(Sample_names %in% controls_to_remove)] # remove controls from vector of samples
phyloseq_16s <- prune_samples(Sample_names, phyloseq_16s) # Sanity check 672 - 21 controls = 651
Sample_names <- sample_names(phyloseq_16s) # vector of samples, this time minus controls
Sample_names <- Sample_names[!(Sample_names %in% soil_to_remove)] # 651 - 28 soil samples = 623
phyloseq_16s_no_soil <- prune_samples(Sample_names, phyloseq_16s) 
# Add factor for year 
phyloseq_16s@sam_data$date_collect <- mdy(phyloseq_16s@sam_data$date_collect)
phyloseq_16s@sam_data$year <- as.factor(year(phyloseq_16s@sam_data$date_collect))
phyloseq_16s_no_soil@sam_data$date_collect <- mdy(phyloseq_16s_no_soil@sam_data$date_collect)
phyloseq_16s_no_soil@sam_data$year <- as.factor(year(phyloseq_16s_no_soil@sam_data$date_collect)) 
# Add sequence depth
phyloseq_16s@sam_data$sequence_depth <- sample_sums(phyloseq_16s)
phyloseq_16s_no_soil@sam_data$sequence_depth <- sample_sums(phyloseq_16s_no_soil)


#### Filter ASVs based on occurrence to remove singletons ####
sum(phyloseq_16s@otu_table) # Total read depth 28424778 
mean(taxa_sums(phyloseq_16s)) # ASV mean read depth 627
ASVs_to_keep <- apply(X = otu_table(phyloseq_16s), 
                      MARGIN = ifelse(taxa_are_rows(phyloseq_16s), yes = 1, no = 2), 
                      FUN = function(x){sum(x > 0)}) >= 5 # remove ASV not found in greater than or equal to 5 samples
ASVs_to_keep <- ASVs_to_keep[ASVs_to_keep == TRUE]
ASVs_to_keep <- names(ASVs_to_keep)
phyloseq_16s_filtered <- prune_taxa(ASVs_to_keep, phyloseq_16s) #filtered dataset
sum(phyloseq_16s_filtered@otu_table) # New total read depth 25380394 & 8838 ASVs
# Dataset without soil
sum(phyloseq_16s_no_soil@otu_table) # Total read depth 27528041 
mean(taxa_sums(phyloseq_16s_no_soil)) # ASV mean read depth 607
ASVs_to_keep <- apply(X = otu_table(phyloseq_16s_no_soil), 
                      MARGIN = ifelse(taxa_are_rows(phyloseq_16s_no_soil), yes = 1, no = 2), 
                      FUN = function(x){sum(x > 0)}) >= 5 # remove ASV not found in greater than or equal to 5 samples
ASVs_to_keep <- ASVs_to_keep[ASVs_to_keep == TRUE]
ASVs_to_keep <- names(ASVs_to_keep)
phyloseq_16s_no_soil_filtered <- prune_taxa(ASVs_to_keep, phyloseq_16s_no_soil) #filtered dataset
sum(phyloseq_16s_no_soil_filtered@otu_table) # New total read depth 24520890 & 7981 ASVs
mean(sample_sums(phyloseq_16s_no_soil_filtered@otu_table)) # Mean sample read count depth 39359

#### Remove samples with less than 1000 sequences ####
Samples_to_remove <- sample_sums(phyloseq_16s_filtered) > 1000
length(Samples_to_remove[Samples_to_remove == FALSE]) # 29 samples of 651 = 4.4%
Samples_to_remove <- labels(Samples_to_remove[Samples_to_remove == FALSE])
Sample_names <- sample_names(phyloseq_16s_filtered)
Sample_names <- Sample_names[!(Sample_names %in% Samples_to_remove)]
phyloseq_16s_filtered <- prune_samples(Sample_names, phyloseq_16s_filtered) #622 samples
mean(sample_sums(phyloseq_16s_filtered@otu_table)) # Mean sample read count depth 40774
sum(phyloseq_16s_filtered@otu_table) # 25361797
#Save for use in DESeq later (CA_transect_differential_abundance.R)
saveRDS(phyloseq_16s_filtered, "phyloseq_16s_filtered_dataset.rds")

# Dataset without soil
Samples_to_remove <- sample_sums(phyloseq_16s_no_soil_filtered) > 1000
length(Samples_to_remove[Samples_to_remove == FALSE]) # 29 samples of 623 = 4.7%
Samples_to_remove <- labels(Samples_to_remove[Samples_to_remove == FALSE])
Sample_names <- sample_names(phyloseq_16s_no_soil_filtered)
Sample_names <- Sample_names[!(Sample_names %in% Samples_to_remove)]
phyloseq_16s_no_soil_filtered <- prune_samples(Sample_names, phyloseq_16s_no_soil_filtered) #594 samples
mean(sample_sums(phyloseq_16s_no_soil_filtered@otu_table)) # Mean sample read count depth 41250
sum(phyloseq_16s_no_soil_filtered@otu_table) # 24502564
#Save for use in DESeq later (CA_transect_differential_abundance.R)
saveRDS(phyloseq_16s_no_soil_filtered, "phyloseq_16s_no_soil_filtered_dataset.rds")


#### Normalization ####
# Transform/normalize data using DESeq2 vst 
summary(phyloseq_16s_filtered@sam_data)
# Function from P.J. McMurdie and S. Holmes (bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html)
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

dds <- phyloseq_to_deseq2(phyloseq_16s_filtered, ~ plant_body_site + year + site) # have to drop rootstock and scion as the soil samples contain NAs for those columns. 
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- estimateDispersions(dds)
vst <- varianceStabilizingTransformation(dds, blind=FALSE,fitType = 'local')
dim(vst) #Same dim as otu table
otu_vst <- as(assay(vst),'matrix') %>% otu_table(.,taxa_are_rows = TRUE) # convert to an otu_table class
phyloseq_16s_filtered_vst <- phyloseq_16s_filtered
otu_table(phyloseq_16s_filtered_vst) <- otu_table(otu_vst, taxa_are_rows = TRUE)
phyloseq_16s_filtered_vst <- transformSampleCounts(phyloseq_16s_filtered_vst,function(x) ifelse(x<0,0,x)) # If negative number after VST set = 0 else leave as x (see https://github.com/joey711/phyloseq/issues/445)
# Save a copy to load later
saveRDS(phyloseq_16s_filtered_vst, "phyloseq_16s_filtered_vst_dataset.rds")
# Dataset without soil
dds <- phyloseq_to_deseq2(phyloseq_16s_no_soil_filtered, ~ rootstock + scion + plant_body_site + year + site + brix)
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- estimateDispersions(dds)
vst <- varianceStabilizingTransformation(dds, blind=FALSE,fitType = 'local')
dim(vst) #Same dim as otu table
otu_vst <- as(assay(vst),'matrix') %>% otu_table(.,taxa_are_rows = TRUE) # convert to an otu_table class
phyloseq_16s_no_soil_filtered_vst <- phyloseq_16s_no_soil_filtered
otu_table(phyloseq_16s_no_soil_filtered_vst) <- otu_table(otu_vst, taxa_are_rows = TRUE)
phyloseq_16s_no_soil_filtered_vst <- transformSampleCounts(phyloseq_16s_no_soil_filtered_vst,function(x) ifelse(x<0,0,x)) # If below negative number after VST 0 else leave as x (see https://github.com/joey711/phyloseq/issues/445)
# Save a copy to load later
saveRDS(phyloseq_16s_no_soil_filtered_vst, "phyloseq_16s_no_soil_filtered_vst_dataset.rds")