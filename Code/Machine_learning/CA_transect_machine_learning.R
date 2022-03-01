# Analysis of grafted grapevines across time and space.
# Samples were collected in summer 2018 and 2019
# Samples consist of leaf, berry, and root compartments
# Code by: Joel F. Swift

# For the data filtering see: Data_filtering_normalization.R
# For general analysis and figures see: CA_transect_analysis.R
# For machine learning see: CA_transect_machine_learning.R
# For DEseq2 diff. abund. see: CA_transect_differential_abundance.R

library('tidyverse'); packageVersion('tidyverse')
library("phyloseq"); packageVersion('phyloseq')
library("ggpubr"); packageVersion('ggpubr')
library("vegan"); packageVersion('vegan')
library("MASS"); packageVersion('MASS')
library("scales"); packageVersion('scales')
library("picante"); packageVersion('picante')
library("caret"); packageVersion('caret')
library("AppliedPredictiveModeling"); packageVersion('AppliedPredictiveModeling')
library("ranger"); packageVersion('ranger')
library("e1071"); packageVersion('e1071')
library("randomForest"); packageVersion('randomForest')
library("alluvial"); packageVersion('alluvial')
library("matrixStats"); packageVersion("matrixStats")

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
# Set seed for analysis
set.seed(1154829343)

# Functions

# Function to return metadata df from phyloseq object
pssd2veg <- function(physeq) {
  # From a phyloseq object return a dataframe of the sample metadata for use in vegan
  # From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# Function to plot or run a linear model on a ASV from a phyloseq object
plot_deseq2_DiffAbunMicob <- function(physeq_obj, ASV_number, FACTOR1 = "rootstock", FACTOR2 = FALSE){
  temp_sample_tab <- pssd2veg(physeq_obj)
  otu_matrix <- as(otu_table(physeq_obj), "matrix")
  TEMP <- data.frame(ASV_count = otu_matrix[ASV_number,])
  TEMP2<- cbind(temp_sample_tab, TEMP)
  # Main effects
  if (isFALSE(FACTOR2) && (FACTOR1) == "rootstock"){
    x <- ggplot(TEMP2, aes(plant_body_site, ASV_count, fill= rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=rootstock_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Compartment") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if (isFALSE(FACTOR2) && (FACTOR1) == "scion"){
    x <- ggplot(TEMP2, aes(plant_body_site, ASV_count, fill= scion)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Scion", values=scion_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Compartment") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if (isFALSE(FACTOR2) && (FACTOR1) == "site"){
    x <- ggplot(TEMP2, aes(plant_body_site, ASV_count, fill= site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Site", values=site_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Compartment") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if (isFALSE(FACTOR2) && (FACTOR1) == "year"){
    x <- ggplot(TEMP2, aes(plant_body_site, ASV_count, fill= year)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Year", values=scion_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Compartment") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if (isFALSE(FACTOR2) && (FACTOR1) == "compartment"){
    x <- ggplot(TEMP2, aes(plant_body_site, ASV_count, fill= plant_body_site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Compartment", values=compartment_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Compartment") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if (isFALSE(FACTOR2) && (FACTOR1) == "sugar content"){
    x <- ggplot(TEMP2, aes(brix_2_breaks, ASV_count, fill= brix_2_breaks)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Sugar Content", values=scion_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Sugar Content") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  }  
  # Joint predictions
  if ((FACTOR1) == "rootstock" && (FACTOR2) == "site" || (FACTOR1) == "site" && (FACTOR2) == "rootstock"){
    x <- ggplot(TEMP2, aes(rootstock, ASV_count, fill= site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Site", values=site_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Rootstock") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "rootstock" && (FACTOR2) == "compartment" || (FACTOR1) == "compartment" && (FACTOR2) == "rootstock"){
    x <- ggplot(TEMP2, aes(rootstock, ASV_count, fill= plant_body_site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Compartment", values=compartment_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Rootstock") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "rootstock" && (FACTOR2) == "scion" || (FACTOR1) == "scion" && (FACTOR2) == "rootstock" ){
    x <- ggplot(TEMP2, aes(scion, ASV_count, fill= rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=rootstock_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Scion") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "rootstock" && (FACTOR2) == "year" || (FACTOR1) == "year" && (FACTOR2) == "rootstock" ){
    x <- ggplot(TEMP2, aes(year, ASV_count, fill= rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=rootstock_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Year") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "scion" && (FACTOR2) == "site" || (FACTOR1) == "site" && (FACTOR2) == "scion"){
    x <- ggplot(TEMP2, aes(scion, ASV_count, fill= site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Site", values=site_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Scion") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "scion" && (FACTOR2) == "compartment" || (FACTOR1) == "compartment" && (FACTOR2) == "scion"){
    x <- ggplot(TEMP2, aes(scion, ASV_count, fill= plant_body_site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Compartment", values=compartment_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Scion") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "scion" && (FACTOR2) == "year" || (FACTOR1) == "year" && (FACTOR2) == "scion"){
    x <- ggplot(TEMP2, aes(year, ASV_count, fill= scion)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Scion", values=scion_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Year") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "year" && (FACTOR2) == "compartment" || (FACTOR1) == "compartment" && (FACTOR2) == "scion"){
    x <- ggplot(TEMP2, aes(year, ASV_count, fill= plant_body_site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Compartment", values=compartment_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Year") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  } else if((FACTOR1) == "year" && (FACTOR2) == "site" || (FACTOR1) == "site" && (FACTOR2) == "scion"){
    x <- ggplot(TEMP2, aes(year, ASV_count, fill= site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Site", values=site_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Year") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  }else if((FACTOR1) == "compartment" && (FACTOR2) == "site" || (FACTOR1) == "site" && (FACTOR2) == "compartment"){
    x <- ggplot(TEMP2, aes(site, ASV_count, fill= plant_body_site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Compartment", values=compartment_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Site") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  }else if((FACTOR1) == "sugar content" && (FACTOR2) == "year" || (FACTOR1) == "year" && (FACTOR2) == "sugar content"){
    x <- ggplot(TEMP2, aes(brix_2_breaks, ASV_count, fill= year)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Year", values=scion_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Sugar content") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  }else if((FACTOR1) == "sugar content" && (FACTOR2) == "site" || (FACTOR1) == "site" && (FACTOR2) == "sugar content"){
    x <- ggplot(TEMP2, aes(brix_2_breaks, ASV_count, fill= site)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Site", values=site_palette) + scale_y_continuous(name="VS Transformed Abundance") + xlab("Sugar content") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(paste("ASV", aes_string(ASV_number)))
  }
  return(x)
}

# Function to plot confusion matrix using ggtile plot from a confussion matrix object
# By user: Enrique Perez Herrero 
# on https://stackoverflow.com/questions/46063234/how-to-produce-a-confusion-matrix-and-find-the-misclassification-rate-of-the-na%C3%AF
ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  
  d <- as.data.frame.matrix(m$table)
  drn <- colnames(d)
  drr <- rownames(d)
  d <- as.data.frame(t(t(d)/colSums(d))) # Standardize the shading in the ggplot by # samples for a given factor level in the reference. 
  d <- d %>% gather(x, value)
  Y <- cbind(as.data.frame(m$table), Proportion = d$value)
  Y$Reference <- fct_rev(Y$Reference) # Added this line to get a downward diagonal 
  p <-
    ggplot(data = Y, aes(x = Reference, y = Prediction, fill= Proportion)) +
    geom_tile( colour = "white") +
    scale_fill_gradient(low = "white", high = "#14A02E", na.value = "white", limits=c(0,1)) +
    ggtitle(mytitle) +
    theme(legend.position = "right", axis.text.x = element_text(angle = 60, hjust = 1)) +
    guides(fill = guide_colorbar(frame.colour = "black", ticks = FALSE))
  return(p)
}

# Define a custom Ranger RF model that saves the per fold info to a seperate folder.
ranger_funcs <- getModelInfo("ranger", regex = FALSE)[[1]]
ranger_funcs$fit <- function (x, y, wts, param, lev, last, classProbs, ...) 
{
  if ((!is.data.frame(x)) || dplyr::is.tbl(x)) 
    x <- as.data.frame(x, stringsAsFactors = TRUE)
  x$.outcome <- y
  if (!is.null(wts)) {
    out <- ranger::ranger(dependent.variable.name = ".outcome", 
                          data = x, mtry = min(param$mtry, ncol(x)), min.node.size = param$min.node.size, 
                          splitrule = as.character(param$splitrule), write.forest = TRUE, 
                          probability = classProbs, case.weights = wts, ...)
  }
  else {
    out <- ranger::ranger(dependent.variable.name = ".outcome", 
                          data = x, mtry = min(param$mtry, ncol(x)), min.node.size = param$min.node.size, 
                          splitrule = as.character(param$splitrule), write.forest = TRUE, 
                          probability = classProbs, ...)
  }
  if (!last) 
    out$y <- y
  save(out, file = paste("./ranger/ranger", param$mtry, param$splitrule, param$min.node.size, format(Sys.time(), "%H_%M_%S.RData"), sep = "_"))
  out
}

# Machine learning main function
MachineLearning_RF_ranger <- function(PHYSEQ_OBJ_1, GROUPING, TREES) {
  # Remove ASV Table and meta data from phyloseq objects
  ASV.df <- as.data.frame(otu_table(PHYSEQ_OBJ_1))
  ASV_metadata.df <- as.data.frame(sample_data(PHYSEQ_OBJ_1))
  # Format ASV table to be used for machine learning applications and make metadata df
  ASV.df <- t(ASV.df)
  ASV_meta.df <- data.frame(Sample = rownames(ASV_metadata.df), Year = ASV_metadata.df$year, Scion = ASV_metadata.df$scion, 
                            Rootstock = ASV_metadata.df$rootstock, Compartment = ASV_metadata.df$plant_body_site, Site = ASV_metadata.df$site,
                            Sugar_content = ASV_metadata.df$brix_2_breaks,
                            Compartment_Rootstock = paste(ASV_metadata.df$plant_body_site, ASV_metadata.df$rootstock, sep = "_"))
  ASV_prefiltered.df <- cbind(ASV.df, ASV_meta.df)
  train_index <- as.data.frame(ASV_prefiltered.df %>% sample_n(475))
  rownames(train_index) <- train_index$Sample
  train_index <- match(rownames(train_index), rownames(ASV_prefiltered.df))
  train_x <- as.data.frame(ASV.df[train_index, ])
  test_y <- as.data.frame(ASV.df[-train_index, ])
  # Train set, 474
  train_x$Sample <- rownames(train_x)
  Training_meta.df <- merge(train_x, ASV_meta.df, by = 'Sample')
  train_x <- subset(Training_meta.df, select = -c(Compartment, Site, Scion, Year, Rootstock, Compartment_Rootstock, Sugar_content))
  rownames(train_x) <- train_x$Sample
  train_x <- subset(train_x, select = -c(Sample))
  Training_meta.df <- subset(Training_meta.df, select = c(Compartment, Site, Scion, Year, Rootstock, Compartment_Rootstock, Sugar_content))
  rownames(Training_meta.df) <- Training_meta.df$Sample 
  # Test set, 120 samples
  test_y$Sample <- rownames(test_y)
  Testing_meta.df <- merge(test_y, ASV_meta.df, by = "Sample")
  test_y <- subset(Testing_meta.df, select = -c(Compartment, Site, Scion, Year, Rootstock, Compartment_Rootstock, Sugar_content))
  rownames(test_y) <- test_y$Sample
  test_y <- subset(test_y, select = -c(Sample))
  Testing_meta.df <- subset(Testing_meta.df, select = c(Compartment, Site, Scion, Year, Rootstock, Compartment_Rootstock, Sugar_content))
  rownames(Testing_meta.df) <- Testing_meta.df$Sample 
  # Training model
  Training_grid <- expand.grid(.mtry = seq(10, length(train_x), round(length(train_x)*0.1)), .splitrule= "gini",
                               .min.node.size = c(1, 5, 10))
  train_control <- trainControl(method="cv", number=10)
  RF_CM <- list()
  RF_CM[["RF_model"]] <- train(x = train_x, y = Training_meta.df[[GROUPING]], method = ranger_funcs, importance = "impurity",
                               tuneGrid = Training_grid, trControl = train_control, num.trees = TREES)
  RF_prediction_3 <- predict(RF_CM[["RF_model"]], test_y)
  RF_CM[["CMatrix"]] <- confusionMatrix(RF_prediction_3, as.factor(Testing_meta.df[[GROUPING]]), mode = "everything")
  RF_CM[["CMatrixPLOT"]] <- ggplotConfusionMatrix(RF_CM[["CMatrix"]])
  RF_CM[["VarImporance"]] <- varImp(RF_CM[["RF_model"]])
  return(RF_CM)
}

# Function extract taxonomy for ASVs with importance data from either a optimal model or per fold model
Get_taxa_importance <- function(RF_MODEL = NULL, FOLD_MODEL = NULL, TOP_N = 100){
  if (!is.null(RF_MODEL)){
    # Make vector of ASV numbers with high importance (specified or limited by TOP_N)
    X <- c(order(RF_MODEL[["VarImporance"]]$importance, decreasing=TRUE)[1:TOP_N])
    my_vect <- c()
    # For loop to go through ASV numbers and get ASV_hashs
    for (i in X) {
      Z <- rownames(RF_MODEL[["VarImporance"]]$importance)[i]
      my_vect <- append(my_vect, c(Z))  
    }
    #return(my_vect) # for testing XXX
    # For loop to get taxonomic assignments of the ASV_Hashs
    taxa_df <- data.frame(Kingdom=character(), Phylum=character(), Class=character(), Order=character(), Family=character(), Genus=character(), Species=character(), stringsAsFactors=FALSE)
    for (i in my_vect){
      taxa_df[i, ] <- c(phy_vst@tax_table[i,])
    }
    #return(taxa_df) # for testing XXX
    # Connect these back to their importance values
    final_df <- merge(taxa_df,RF_MODEL[["VarImporance"]]$importance, by="row.names")
    colnames(final_df) <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Importance")
    # Fix first column name, coerce all columns but importance to factor()
    cols <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    final_df[,cols] <- lapply(final_df[cols], factor) 
    return(final_df)
  } else if (!is.null(FOLD_MODEL)){
    X <- as.data.frame(rescale(FOLD_MODEL$variable.importance, to = c(0,100))) # Rescale to 0 to 100 like caret does for varIMP
    X <- c(order(X, decreasing=TRUE)[1:TOP_N])
    my_vect <- c()
    # For loop to go through ASV numbers and get ASV_hashs
    for (i in X) {
      Z <- rownames(as.data.frame(FOLD_MODEL$variable.importance))[i]
      my_vect <- append(my_vect, c(Z))  
      #return(my_vect) # for testing XXX
    }
    taxa_df <- data.frame(Kingdom=character(), Phylum=character(), Class=character(), Order=character(), Family=character(), Genus=character(), Species=character(), stringsAsFactors=FALSE)
    for (i in my_vect){
      taxa_df[i, ] <- c(phy_vst@tax_table[i,])
    }
    #return(taxa_df) # for testing XXX
    # Connect these back to their importance values
    final_df <- merge(taxa_df, as.data.frame(rescale(FOLD_MODEL$variable.importance, to = c(0,100))), by="row.names")
    colnames(final_df) <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Importance")
    # Fix first column name, coerce all columns but importance to factor()
    cols <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    final_df[,cols] <- lapply(final_df[cols], factor) 
    return(final_df)
  }
}

# Function to obtain importance values with taxonomy for ASVs from a list of Ranger random forest fold outputs from caret.
# This should take ~3 mins per list (its getting almost 8k values from 10 folds).
Extract_importance <- function(list_of_imp = NULL){
  VarI_list<- c()
  for (i in seq(1,10,1)){
    X <- list_of_imp[[i]]$out
    X <- Get_taxa_importance(FOLD_MODEL = X, TOP_N = 7981)
    VarI_list[[i]]<- X
  }
  return(VarI_list)
}

# Merge importance columns from all entries
merge_imp <- function(obj){
  X <- cbind(obj[[1]],
             "Importance_02" = obj[[2]][,c(9)],
             "Importance_03" = obj[[3]][,c(9)],
             "Importance_04" = obj[[4]][,c(9)],
             "Importance_05" = obj[[5]][,c(9)],
             "Importance_06" = obj[[6]][,c(9)],
             "Importance_07" = obj[[7]][,c(9)],
             "Importance_08" = obj[[8]][,c(9)],
             "Importance_09" = obj[[9]][,c(9)],
             "Importance_10" = obj[[10]][,c(9)])
  return(X)
}

#------------------------------------------------------------------------------------#
# load dataset generated from Data_filtering_normalization.R
phy_vst <- readRDS("phyloseq_16s_no_soil_filtered_vst_dataset.rds")
# 2 Breaks Brix 
phy_vst@sam_data$brix_2_breaks <- cut(phy_vst@sam_data$brix, breaks = c(-Inf, 7, Inf), labels = c("Pre-ripening", "Ripening")) # Check if I want to use these terms
summary(phy_vst@sam_data)

# Running with 10 CV on 80:20 dataspilt with optimal # trees from CA_ML_testing.R
###### THIS CODE BLOCK WAS RUN ON THE ARIES CLUSTER IN JUPYTER #######
# Machine_learning_chap2
#RF_CM_SITE <- MachineLearning_RF_ranger(phy_vst, "Site", 237)
#saveRDS(RF_CM_SITE, "RF_CM_CV10_237tree_Site.rand_split.rds")
#RF_CM_ROOTSTOCK <- MachineLearning_RF_ranger(phy_vst, "Rootstock", 255)
#saveRDS(RF_CM_ROOTSTOCK, "RF_CM_CV10_255tree_Rootstock.rand_split.rds")
#RF_CM_COMPARTMENT <- MachineLearning_RF_ranger(phy_vst, "Compartment", 63)
#saveRDS(RF_CM_COMPARTMENT, "RF_CM_CV10_63tree_Compartment.rand_split.rds")
#RF_CM_YEAR <- MachineLearning_RF_ranger(phy_vst, "Year", 477)
#saveRDS(RF_CM_YEAR, "RF_CM_CV10_477tree_Year.rand_split.rds")
#RF_CM_SCION <- MachineLearning_RF_ranger(phy_vst, "Scion", 433)
#saveRDS(RF_CM_SCION, "RF_CM_CV10_433tree_Scion.rand_split.rds")
#RF_CM_Sugar <- MachineLearning_RF_ranger(phy_vst, "Sugar_content", 267)
#saveRDS(RF_CM_Sugar, "RF_CM_CV10_267tree_Sugar.rand_split.rds")

#### Main effects  ####
# Loading completed ML runs from the cluster
RF_root <- readRDS("ML_testing/RF_CM_CV10_255tree_Rootstock.rand_split.rds")
RF_scio <- readRDS("ML_testing/RF_CM_CV10_433tree_Scion.rand_split.rds")
RF_site <- readRDS("ML_testing/RF_CM_CV10_237tree_Site.rand_split.rds")
RF_year <- readRDS("ML_testing/RF_CM_CV10_477tree_Year.rand_split.rds")
RF_comp <- readRDS("ML_testing/RF_CM_CV10_63tree_Compartment.rand_split.rds")
RF_Sugr <- readRDS("ML_testing/RF_CM_CV10_267tree_Sugar.rand_split.rds")
# Confusion Matrix plot for supplement
a <- RF_comp$CMatrixPLOT
b <- RF_root$CMatrixPLOT
c <- RF_site$CMatrixPLOT
d <- RF_year$CMatrixPLOT
e <- RF_scio$CMatrixPLOT
f <- RF_Sugr$CMatrixPLOT
CM_PLOT <- ggarrange(a,b,c,d,e,f  , align = "hv", common.legend = TRUE, legend = 'right', labels = c("A","B","C","D","E"))
ggsave("Confusion_matrix_plots.svg", CM_PLOT, width = 12, height = 8)
ggsave("Confusion_matrix_plots.pdf", CM_PLOT, width = 12, height = 8)

# Plotting ASVs that aid in predictions of the main effects
# Load a list of files to get mean and standard error for the OOB samples
# Rootstock
files <- list.files(path = "./ML_testing/ranger_rootstock/", pattern = ".*1606_gini_10.*")
files <- paste("./ML_testing/ranger_rootstock/", files, sep = "")
RF_list_root <- lapply(files, function(x) mget(load(x)))
RF_list_root <- Extract_importance(RF_list_root)
RF_list_root <- merge_imp(RF_list_root)
RF_list_root[19] <- as.data.frame(rowSds(as.matrix(RF_list_root[,9:18])))
RF_list_root[20] <- as.data.frame(rowMeans((RF_list_root[,9:18])))
colnames(RF_list_root)[19] <- "SD"
colnames(RF_list_root)[20] <- "mean"
RF_list_root <- plyr::join(data.frame('ASV' = rownames(phy_vst@tax_table)), RF_list_root) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_root$ASV_numb <- paste("ASV", rownames(RF_list_root))
# Scion
files <- list.files(path = "./ML_testing/ranger_scion/", pattern = ".*4000_gini_5.*")
files <- paste("./ML_testing/ranger_scion/", files, sep = "")
RF_list_scio <- lapply(files, function(x) mget(load(x)))
RF_list_scio <- Extract_importance(RF_list_scio)
RF_list_scio <- merge_imp(RF_list_scio)
RF_list_scio[19] <- as.data.frame(rowSds(as.matrix(RF_list_scio[,9:18])))
RF_list_scio[20] <- as.data.frame(rowMeans((RF_list_scio[,9:18])))
colnames(RF_list_scio)[19] <- "SD"
colnames(RF_list_scio)[20] <- "mean"
RF_list_scio <- plyr::join(data.frame('ASV' = rownames(phy_vst@tax_table)), RF_list_scio) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_scio$ASV_numb <- paste("ASV", rownames(RF_list_scio))
# Site
files <- list.files(path = "./ML_testing/ranger_site/", pattern = ".*808_gini_5.*")
files <- paste("./ML_testing/ranger_site/", files, sep = "")
RF_list_site <- lapply(files, function(x) mget(load(x)))
RF_list_site <- Extract_importance(RF_list_site)
RF_list_site <- merge_imp(RF_list_site)
RF_list_site[19] <- as.data.frame(rowSds(as.matrix(RF_list_site[,9:18])))
RF_list_site[20] <- as.data.frame(rowMeans((RF_list_site[,9:18])))
colnames(RF_list_site)[19] <- "SD"
colnames(RF_list_site)[20] <- "mean"
RF_list_site <- plyr::join(data.frame('ASV' = rownames(phy_vst@tax_table)), RF_list_site) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_site$ASV_numb <- paste("ASV", rownames(RF_list_site))
# Year
files <- list.files(path = "./ML_testing/ranger_year/", pattern = ".*1606_gini_10.*")
files <- paste("./ML_testing/ranger_year/", files, sep = "")
RF_list_year <- lapply(files, function(x) mget(load(x)))
RF_list_year <- Extract_importance(RF_list_year)
RF_list_year <- merge_imp(RF_list_year)
RF_list_year[19] <- as.data.frame(rowSds(as.matrix(RF_list_year[,9:18])))
RF_list_year[20] <- as.data.frame(rowMeans((RF_list_year[,9:18])))
colnames(RF_list_year)[19] <- "SD"
colnames(RF_list_year)[20] <- "mean"
RF_list_year <- plyr::join(data.frame('ASV' = rownames(phy_vst@tax_table)), RF_list_year) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_year$ASV_numb <- paste("ASV", rownames(RF_list_year))
# Compartment
files <- list.files(path = "./ML_testing/ranger_compartment/", pattern = ".*1606_gini_10.*")
files <- paste("./ML_testing/ranger_compartment/", files, sep = "")
RF_list_comp <- lapply(files, function(x) mget(load(x)))
RF_list_comp <- Extract_importance(RF_list_comp)
RF_list_comp <- merge_imp(RF_list_comp)
RF_list_comp[19] <- as.data.frame(rowSds(as.matrix(RF_list_comp[,9:18])))
RF_list_comp[20] <- as.data.frame(rowMeans((RF_list_comp[,9:18])))
colnames(RF_list_comp)[19] <- "SD"
colnames(RF_list_comp)[20] <- "mean"
RF_list_comp <- plyr::join(data.frame('ASV' = rownames(phy_vst@tax_table)), RF_list_comp) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_comp$ASV_numb <- paste("ASV", rownames(RF_list_comp))
# Sugar Content
files <- list.files(path = "./ML_testing/ranger_sugar/", pattern = ".*5596_gini_10.*")
files <- paste("./ML_testing/ranger_sugar/", files, sep = "")
RF_list_sugar <- lapply(files, function(x) mget(load(x)))
RF_list_sugar <- Extract_importance(RF_list_sugar)
RF_list_sugar <- merge_imp(RF_list_sugar)
RF_list_sugar[19] <- as.data.frame(rowSds(as.matrix(RF_list_sugar[,9:18])))
RF_list_sugar[20] <- as.data.frame(rowMeans((RF_list_sugar[,9:18])))
colnames(RF_list_sugar)[19] <- "SD"
colnames(RF_list_sugar)[20] <- "mean"
RF_list_sugar <- plyr::join(data.frame('ASV' = rownames(phy_vst@tax_table)), RF_list_sugar) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_sugar$ASV_numb <- paste("ASV", rownames(RF_list_sugar))

# Subplots for each main effect
a <- ggplot(RF_list_comp[RF_list_comp$mean > 25, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  theme(axis.title.x = element_blank()) + ylab("ASV") + xlim(c(0,100))

b <- ggplot(RF_list_root[RF_list_root$mean > 25, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Relative Decrease in Gini ") + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + xlim(c(0,100))

c <- ggplot(RF_list_site[RF_list_site$mean > 25, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Relative Decrease in Gini ") + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + xlim(c(0,100))

d <- ggplot(RF_list_year[RF_list_year$mean > 25, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Relative Decrease in Gini ") + ylab("ASV") + xlim(c(0,100))

e <- ggplot(RF_list_scio[RF_list_scio$mean > 25, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Relative Decrease in Gini ") + theme(axis.title.y = element_blank()) + xlim(c(0,100))

f <- ggplot(RF_list_sugar[RF_list_sugar$mean > 25, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Relative Decrease in Gini ") + theme(axis.title.y = element_blank()) + xlim(c(0,100))

# Combination plot
ggarrange(a,b,c,d,e,f, common.legend = TRUE, labels = c("Compartment","Rootstock", "Site", "Year", "Scion", "Sugar"), align = "hv", label.y = 1, legend = "right")
# Supplement or main? XXX
ASVs_deciding_ML <- ggarrange(a,b,c,d,e,f, common.legend = TRUE, labels = c("AUTO"), align = "hv", label.y = 1, legend = "right")
ggsave("ASVs_ML_main_effects.pdf", ASVs_deciding_ML, height = 8, width = 14)
ggsave("ASVs_ML_main_effects.svg", ASVs_deciding_ML, height = 8, width = 14)
rm(RF_list_comp, RF_list_root, RF_list_scio, RF_list_site, RF_list_year) # clean env

# Getting taxonomy for main effects ASVs from model prediction of the test set
RF_comp.df <- Get_taxa_importance(RF_comp, TOP_N = 7981)
RF_site.df <- Get_taxa_importance(RF_site, TOP_N = 7981)
RF_root.df <- Get_taxa_importance(RF_root, TOP_N = 7981)
RF_scio.df <- Get_taxa_importance(RF_scio, TOP_N = 7981)
RF_year.df <- Get_taxa_importance(RF_year, TOP_N = 7981)
RF_sugr.df <- Get_taxa_importance(RF_Sugr, TOP_N = 7981)

# Add a factor to each prediction and convert to relative importance
a <- data.frame(RF_comp.df, Factor = "Compartment")
a$Importance <- a$Importance / sum(a$Importance)

b <- data.frame(RF_root.df, Factor = "Rootstock")
b$Importance <- b$Importance / sum(b$Importance)

c <- data.frame(RF_site.df, Factor = "Site")
c$Importance <- c$Importance / sum(c$Importance)

d <- data.frame(RF_year.df, Factor = "Year")
d$Importance <- d$Importance / sum(d$Importance)

e <- data.frame(RF_scio.df, Factor = "Scion")
e$Importance <- e$Importance / sum(e$Importance)

f <- data.frame(RF_sugr.df, Factor = "Sugar Content")
f$Importance <- f$Importance / sum(f$Importance)

g <- rbind(a, b, c, d, e, f)



# Check which class of proteobacteria make up the barchart
# Proteobacteria was the dominant phyla that assisted predictions so I broke up the classes here.
sum(as.numeric(a[a$Phylum == "Proteobacteria",]$Importance), na.rm=TRUE) # 54.6%
sum(as.numeric(b[b$Phylum == "Proteobacteria",]$Importance), na.rm=TRUE) # 44.3%
sum(as.numeric(c[c$Phylum == "Proteobacteria",]$Importance), na.rm=TRUE) # 41.0%
sum(as.numeric(d[d$Phylum == "Proteobacteria",]$Importance), na.rm=TRUE) # 39.7%
sum(as.numeric(e[e$Phylum == "Proteobacteria",]$Importance), na.rm=TRUE) # 46.9%
sum(as.numeric(f[f$Phylum == "Proteobacteria",]$Importance), na.rm=TRUE) # 42.3%

sum(as.numeric(a[a$Phylum == "Actinobacteriota",]$Importance), na.rm=TRUE) # 14.0%
sum(as.numeric(b[b$Phylum == "Actinobacteriota",]$Importance), na.rm=TRUE) # 20.1%
sum(as.numeric(c[c$Phylum == "Actinobacteriota",]$Importance), na.rm=TRUE) # 17.6%
sum(as.numeric(d[d$Phylum == "Actinobacteriota",]$Importance), na.rm=TRUE) # 24.9%
sum(as.numeric(e[e$Phylum == "Actinobacteriota",]$Importance), na.rm=TRUE) # 19.7%
sum(as.numeric(f[f$Phylum == "Actinobacteriota",]$Importance), na.rm=TRUE) # 23.4%
# compartment
sum(as.numeric(a[a$Class == "Alphaproteobacteria",]$Importance), na.rm=TRUE) # 32.6%
sum(as.numeric(a[a$Class == "Gammaproteobacteria",]$Importance), na.rm=TRUE) # 21.9%
# rootstock
sum(as.numeric(b[b$Class == "Alphaproteobacteria",]$Importance), na.rm=TRUE) # 20.1%
sum(as.numeric(b[b$Class == "Gammaproteobacteria",]$Importance), na.rm=TRUE) # 24.1%
# site
sum(as.numeric(c[c$Class == "Alphaproteobacteria",]$Importance), na.rm=TRUE) # 17.8%
sum(as.numeric(c[c$Class == "Gammaproteobacteria",]$Importance), na.rm=TRUE) # 23.1%
# year
sum(as.numeric(d[d$Class == "Alphaproteobacteria",]$Importance), na.rm=TRUE) # 15.9%
sum(as.numeric(d[d$Class == "Gammaproteobacteria",]$Importance), na.rm=TRUE) # 23.7%
# scion
sum(as.numeric(e[e$Class == "Alphaproteobacteria",]$Importance), na.rm=TRUE) # 21.4%
sum(as.numeric(e[e$Class == "Gammaproteobacteria",]$Importance), na.rm=TRUE) # 25.5%
# Sugar content
sum(as.numeric(f[f$Class == "Alphaproteobacteria",]$Importance), na.rm=TRUE) # 15.8%
sum(as.numeric(f[f$Class == "Gammaproteobacteria",]$Importance), na.rm=TRUE) # 26.4%

# Pick phyla to retain and collapse rest to other (based on ASV phyla with ASV of over 25 relative gini importance; see figure above)
phyla_to_keep <- c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", 
                   "Chloroflexi", "Crenarchaeota", "Deinococcota", "Desulfobacterota", 
                   "Firmicutes", "Myxococcota", "Planctomycetota", 
                   "Proteobacteria", "Verrucomicrobiota")
g$Phylum <- forcats::fct_other(g$Phylum, keep = phyla_to_keep, other_level = "Other")
g <- as.data.frame(g %>% group_by(Factor, Phylum) %>% summarise(Importance = sum(Importance))) # Collapse to importance per phyla to aid in use in inkscape

g$Factor <- factor(g$Factor, levels = c("Compartment", "Rootstock", "Site", "Year", "Scion", "Sugar Content"))

h <- ggplot(g, aes(x = Factor, y = Importance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=safe_colorblind_palette) +
  ylab("Relative Importance to Classifier") +
  theme(legend.position = "right", axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Plot showing the relative importance of bacterial phyla in ML decisions.
# Save Supplement or main?!? XXX
ggsave("Relative importance of phyla in classification.pdf", h, height = 8, width = 8)
ggsave("Relative importance of phyla in classification.svg", h, height = 8, width = 8)


rm(a,b,c,d,e,f,g,RF_root.df,RF_comp.df,RF_scio.df,RF_site.df,RF_year.df,plots,Legends,Legend1,Legend2,to_save) # clean env


# Visualizing the boxplots of all the ASVs that contribute to predictions above 25%
# For looking at the important ASVs ROOTSTOCK
most_imp_ASV <- rownames(RF_root[["VarImporance"]]$importance)[RF_root[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 664, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 1147, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 1312, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 1340, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2465, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2470, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2592, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 2983, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 3017, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 3050, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 3248, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 3938, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 6460, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 7495, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 7791, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "rootstock")
plot_deseq2_DiffAbunMicob(phy_vst, 7892, FACTOR1 = "rootstock")

# For looking at the important ASVs SITE
most_imp_ASV <- rownames(RF_site[["VarImporance"]]$importance)[RF_site[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 24, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 28, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 48, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 50, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 233, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 407, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1312, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2461, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2470, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2580, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3050, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3194, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 4018, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 6527, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "site")

# For looking at the important ASVs SCION
most_imp_ASV <- rownames(RF_scio[["VarImporance"]]$importance)[RF_scio[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 664, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1170, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1312, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1340, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1342, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2461, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2465, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2470, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2592, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 3261, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7876, FACTOR1 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7889, FACTOR1 = "scion")


# For looking at the important ASVs COMPARTMENT
most_imp_ASV <- rownames(RF_comp[["VarImporance"]]$importance)[RF_comp[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 35, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 238, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 670, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 4940, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 5270, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7436, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7716, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7980, FACTOR1 = "compartment")

# For looking at the important ASVs YEAR
most_imp_ASV <- rownames(RF_year[["VarImporance"]]$importance)[RF_year[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1920, FACTOR1 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2115, FACTOR1 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2116, FACTOR1 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "year")

# For looking at the important ASVs Sugar Content
most_imp_ASV <- rownames(RF_Sugr[["VarImporance"]]$importance)[RF_Sugr[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "sugar content")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "sugar content")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "sugar content")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "sugar content") # + geom_point(aes(colour = year)) #something interesting with this one and year
plot_deseq2_DiffAbunMicob(phy_vst, 2599, FACTOR1 = "sugar content")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "sugar content")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "sugar content")
plot_deseq2_DiffAbunMicob(phy_vst, 3216, FACTOR1 = "sugar content")
plot_deseq2_DiffAbunMicob(phy_vst, 4114, FACTOR1 = "sugar content")

# Example boxplots for supplemental 

# Compartment/Site ASVs
crenarch_plot_1 <- plot_deseq2_DiffAbunMicob(phy_vst, 35, FACTOR1 = "site")
crenarch_plot_2 <- plot_deseq2_DiffAbunMicob(phy_vst, 48, FACTOR1 = "site")
crenarch_plot_3 <- plot_deseq2_DiffAbunMicob(phy_vst, 28, FACTOR1 = "site")

# Year ASV
Actino_1 <- plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "year")

# Scion/Rootstock
Proteo_1 <- plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "scion", FACTOR2 = "rootstock")

# Sugar content
Proteo_2 <- plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "sugar content", FACTOR2 = "year")
Proteo_3 <- plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "sugar content", FACTOR2 = "site") 

ASVs_of_interest_plot <- ggarrange(crenarch_plot_1,crenarch_plot_2, Actino_1, Proteo_1, Proteo_2, Proteo_3, nrow = 2,  ncol = 3, labels = "AUTO")

ggsave("ASV_boxplots_of_ML_importance.pdf", ASVs_of_interest_plot, width = 18, height = 10)
ggsave("ASV_boxplots_of_ML_importance.svg", ASVs_of_interest_plot, width = 18, height = 10)


#### Joint Predictions ####

# Load joint models for confusion matrix plot for supplement
RF_CxSi <- readRDS("ML_testing/RF_CM_CV10_400tree_Compartment_by_Site.rand_split.rds")
RF_RxC  <- readRDS("ML_testing/RF_CM_CV10_400tree_Rootstock_by_Compartment.rand_split.rds")
RF_RxS  <- readRDS("ML_testing/RF_CM_CV10_400tree_Rootstock_by_Scion.rand_split.rds")
RF_RxSi <- readRDS("ML_testing/RF_CM_CV10_400tree_Rootstock_by_Site.rand_split.rds")
RF_RxY  <- readRDS("ML_testing/RF_CM_CV10_400tree_Rootstock_by_Year.rand_split.rds")
RF_SxC  <- readRDS("ML_testing/RF_CM_CV10_400tree_Scion_by_Compartment.rand_split.rds")
RF_SxSi <- readRDS("ML_testing/RF_CM_CV10_400tree_Scion_by_Site.rand_split.rds")
RF_SxY  <- readRDS("ML_testing/RF_CM_CV10_400tree_Scion_by_Year.rand_split.rds")
RF_YxC  <- readRDS("ML_testing/RF_CM_CV10_400tree_Year_by_Compartment.rand_split.rds")
RF_YxSi <- readRDS("ML_testing/RF_CM_CV10_400tree_Year_by_Site.rand_split.rds")

a <- RF_CxSi$CMatrixPLOT
b <- RF_RxC$CMatrixPLOT
c <- RF_RxS$CMatrixPLOT
d <- RF_RxSi$CMatrixPLOT
e <- RF_RxY$CMatrixPLOT
f <- RF_SxC$CMatrixPLOT
g <- RF_SxSi$CMatrixPLOT
h <- RF_SxY$CMatrixPLOT
i <- RF_YxC$CMatrixPLOT
j <- RF_YxSi$CMatrixPLOT

pdf("Confusion_matrix_plots_JOINT_PRED.pdf")
a
b
c
d
e
f
g
h
i
j
dev.off()

# Clean up the environment
rm(a,b,c,d,e,f,g,h,i,j)

# Visualizing the boxplots of all the ASVs that contribute to predictions above 10%
#Compartment_by_Site
most_imp_ASV <- rownames(RF_CxSi[["VarImporance"]]$importance)[RF_CxSi[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 48, FACTOR1 = "compartment", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 407, FACTOR1 = "compartment", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "compartment", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "compartment", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "compartment", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "compartment", FACTOR2 = "site")

#Rootstock_by_Compartment
most_imp_ASV <- rownames(RF_RxC[["VarImporance"]]$importance)[RF_RxC[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 35, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 603, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7716, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7892, FACTOR1 = "rootstock", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7980, FACTOR1 = "rootstock", FACTOR2 = "compartment")


#Rootstock_by_Scion
most_imp_ASV <- rownames(RF_RxS[["VarImporance"]]$importance)[RF_RxS[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2470, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 3261, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "rootstock", FACTOR2 = "scion")
plot_deseq2_DiffAbunMicob(phy_vst, 7889, FACTOR1 = "rootstock", FACTOR2 = "scion")

#Rootstock_by_Site
most_imp_ASV <- rownames(RF_RxSi[["VarImporance"]]$importance)[RF_RxSi[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 16, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 48, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1312, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2461, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2470, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3126, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3261, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3961, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "rootstock", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "rootstock", FACTOR2 = "site")


#Rootstock_by_Year
most_imp_ASV <- rownames(RF_RxY[["VarImporance"]]$importance)[RF_RxY[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2983, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "rootstock", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "rootstock", FACTOR2 = "year")

#Scion_by_Compartment
most_imp_ASV <- rownames(RF_SxC[["VarImporance"]]$importance)[RF_SxC[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 35, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 238, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1312, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 4940, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7436, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7716, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "scion", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7980, FACTOR1 = "scion", FACTOR2 = "compartment")

#Scion_by_Site
most_imp_ASV <- rownames(RF_SxSi[["VarImporance"]]$importance)[RF_SxSi[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 16, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 48, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1312, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2461, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2470, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3126, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3261, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3961, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "scion", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "scion", FACTOR2 = "site")

#Scion_by_Year
most_imp_ASV <- rownames(RF_SxY[["VarImporance"]]$importance)[RF_SxY[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 1920, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 3582, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "scion", FACTOR2 = "year")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "scion", FACTOR2 = "year")

#Year_by_Compartment
most_imp_ASV <- rownames(RF_YxC[["VarImporance"]]$importance)[RF_YxC[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 35, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 1920, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2115, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 4940, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7716, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "year", FACTOR2 = "compartment")
plot_deseq2_DiffAbunMicob(phy_vst, 7980, FACTOR1 = "year", FACTOR2 = "compartment")

#Year_by_Site
most_imp_ASV <- rownames(RF_YxSi[["VarImporance"]]$importance)[RF_YxSi[["VarImporance"]]$importance > 25.0]
most_imp_number <- match(most_imp_ASV, rownames(phy_vst@tax_table))
plot_deseq2_DiffAbunMicob(phy_vst, 48, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 407, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 627, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 661, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 816, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 843, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 844, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 953, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1054, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1091, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1171, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1172, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 1920, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2073, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2389, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2461, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2470, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2473, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2802, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2949, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2983, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 2986, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3050, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3244, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 3262, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 4114, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7569, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7636, FACTOR1 = "year", FACTOR2 = "site")
plot_deseq2_DiffAbunMicob(phy_vst, 7878, FACTOR1 = "year", FACTOR2 = "site")