# California transect microbiomes

## Abridged abstract

Coming Soon!


Table of markers used within the study.

| Marker   | F Primer| F Primer Sequence       | R Primer| R Primer Sequence    | Size (bp) | Citation                                       |
|:--------:|:-------:|:-----------------------:|:-------:|:--------------------:|:---------:|:----------------------------------------------:|
| 16s      | 515F    | GTGYCAGCMGCCGCGGTAA     | 806R    | GGACTACNVGGGTWTCTAAT | 390       | Parada *et al.* 2016 and Apprill *et al.* 2015 |

## Experimental and Sampling Design

Multiple scion/rootstock combinations were sampled across vineyards in the 2018 and 2019 growing seasons in the Central Valley of California. **A)** Each of the colored counties represents a vineyard that was sampled. Within each county box (San Joaquin, Merced, and Madera) the scion/rootstock combinations are depicted, corresponding to the legend in the bottom left. The sugar content (measured in degrees Brix) was collected for vines across the **B)** 2018 and **C)** 2019 growing seasons. Grey shading represents collection windows for microbiome sampling.

![Image of experimental and sampling design](https://github.com/Kenizzer/California_Transect_Microbiome/blob/master/Experimenetal_design_image/Figure1.png)

## Processing 

### QIIME2 
***See Qiime_analysis.md***

**QIIME2** (Bolyen et al. 2019) was used to demultiplex samples according to barcode sequence. The DADA2 plugin (Callahan et al. 2016) was used to denoise, dereplicate, and filter chimeric sequences on each sequencing run individually for accurate error model generation. Afterward, the resulting amplicon sequence variant (ASV) tables and catalogs of representative sequences for each sequence plate were merged. A Naive Bayes classifier pre-trained on the SILVA v.138 16S database (Yilmaz et al. 2014; Bokulich et al. 2018, 2021) was used for taxonomic classification. ASVs not assigned to a phylum were removed along with ASVs assigned to chloroplasts and mitochondria. 

### Filtering
***See Data_filtering_normalization.r***

The package **Decontam** (Davis et al. 2017) was used to remove contaminants using the prevalence-based detection method with a threshold of 0.5, removing contaminant ASVs more prevalent in the negative controls than real samples. The data set was filtered to remove singletons by retaining only ASVs present in five or more samples and removing samples with a read count less than 1,000. ASV counts were then normalized by applying a variance stabilizing transformation from the package **DESeq2** (Love et al. 2014) with a model containing all of the main effects.


## Bacterial Analysis

### Alpha and beta diversity

***See CA_transect_analysis.r***
 Alpha diversity statistics were calculated using **vegan** (Anderson, M. J. 2001; Oksanen *et al.* 2019) and **picante** (Kembel et al. 2010) in the case of Faith’s Phylogenetic distance (Faith 1992). Linear mixed models were fit via the lmerTest package v3.1.3 (Kuznetsova et al. 2017) to assess the effects of the experimental design on alpha diversity indices: `response ~ Rootstock (R) + Scion (S) + Compartment (C) + Sugar Content (Su) + Year + Site + R×S + C×Su`. PERMANOVA analyses were conducted using the function Adonis from the package **vegan**. For each plant compartment (berries, leaves, roots), a model with Bray-Curtis dissimilarity as the response and fit with all factors as marginal fixed effects with 1,000 permutations per model.

### Differential abundance
***See CA_transect_differential_abundance***

Differential abundance analysis was conducted using **DESeq2** (Love et al. 2014). Using raw reads, we removed ASVs not represented by a depth of at least 25 reads in more than 25 samples. DESeq2 was fit using a model composed of each factor as a fixed main effect. For each factor, we extracted all of the ASVs that exhibiting significant differential abundance after P-Value correction (Benjamini and Hochberg 1995), as well as their fold changes (Log2) to generate summary plots.

### Machine Learning
***See Code/Machine_learning***

 We used **ranger’s** (Wright and Ziegler 2017) implementation of the random forest algorithm in the package **caret** (Kuhn 2008). For training the random forest classifier, the dataset was randomly split into a train set (80%) and a test set (20%). Optimal hyperparameters for each classifier were determined using a grid search over the number of trees (1-501), minimum node size (1, 5, 10), and number of features available at each node (10-100% of the ASVs). For each combination in the grid, the performance of the classifier on out-of-bag samples was assessed with 10-fold cross validation. Classifiers were then trained to predict each of the categorical factors (rootstock, scion, compartment, year, and site). Tile plots were used to visualize the output confusion matrix results. We determined the relative importance of phyla in classification accuracy per factor, as well as ASVs that contributing considerably to classifier accuracy (i.e., high gini importance), top ASVs were further examined using boxplots.

### Soil Analysis
***See Code/Soil_Analysis***

For soil texture, hydrometer readings were processed using the package **envalysis** (Steinmetz 2021) to obtain percentages of sand, silt, and clay. Independent linear models for sand, silt, and clay were fit with collection site and year as main effects. For elemental composition analysis, concentrations values greater than five standard deviations from the mean were removed. A biplot was generated using the package **factoextra** (Kassambara and Mundt 2020) to visualize clustering of soil samples by collection site along with the loadings of the principal component analysis (PCA). A linear model was fit to the first two principal components (PC) with collection site and year as main effects. Each of the linear models was assessed via a type-2 ANOVA framework. Soil microbiota samples were processed similar to the plant compartment samples above.

### Brix Analysis
***See CA_transect_analysis.R (lines 224-269)***

For each vine, we measured the total soluble solids (measure in °Bx) using a hand-held refractometer (ATAGO) by selecting berries from a damage-free representative cluster. We modeled the sugar content of the berries of the vines across the growing season. Using a linear model, we assessed the effects from the experimental design on the sugar content (°Brix) using the following model: `Sugar Content ~ Rootstock + Site + Scion + Collection Week + Scion×Collection Week`. The package **car** (Fox and Weisberg 2019) was used to assess the model under a type-3 ANOVA framework.


**Citations**
---
1.  Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 2019; 37: 852–857. 
2.  Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 2016; 13: 581–583. 
3.  Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, et al. The SILVA and “All-species Living Tree Project (LTP)” taxonomic frameworks. Nucleic Acids Res 2014; 42: D643–D648. 
4.  Bokulich NA, Kaehler BD, Rideout JR, Dillon M, Bolyen E, Knight R, et al. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. Microbiome 2018; 6: 1–17. 
5.  Bokulich NA, Robeson M, Dillon M, Ziemski M, Kaehler B, O’Rourke D. RESCRIPt: 2021.8.0.dev0. Zenodo. . 
6.  Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Simple Stat Identif Remov Contam Seq marker-gene metagenomics data 2017; 221499. 
7.  Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 2014; 15: 550. 
8.  Anderson MJ. A new method for non-parametric multivariate analysis of variance. Austral Ecol 2001; 26: 32–46. 
9.  Oksanen J, Guillaume Blanchet F, Friendly M, Kindt R, Legendre P, McGlinn D, et al. vegan: Community Ecology Package. 2019. 
10.   Kembel SW, Cowan PD, Helmus MR, Cornwell WK, Morlon H, Ackerly DD, et al. Picante: R tools for integrating phylogenies and ecology. Bioinformatics 2010; 26: 1463–1464. 
11.   Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv 1992; 61: 1–10. 
12.   Kuznetsova A, Brockhoff PB, Christensen RHB. lmerTest Package: Tests in Linear Mixed Effects Models. J Stat Softw 2017; 82. 
13.   Benjamini Y, Hochberg Y. Controlling the False Discovery Rate : A Practical and Powerful Approach to Multiple Testing. J R Stat Soc 1995; 57: 289–300. 
14.   Wright MN, Ziegler A. ranger : A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 2017; 77. 
15.   Kuhn M. Building Predictive Models in R Using the caret Package. J Stat Softw 2008; 28: 159–160. 
16.   Steinmetz Z. envalysis: Miscellaneous Functions for Environmental Analyses. 2021. CRAN. 
17.   Kassambara A, Mundt F. factoextra: Extract and Visualize the Results of Multivariate Data Analyses. 2020. CRAN. 
18.   Fox J, Weisberg S. An R Companion to Applied Regression, 3rd ed. 2019. Sage, Thousands Oaks, CA.