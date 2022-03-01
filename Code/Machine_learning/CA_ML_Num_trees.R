# Analysis of grafted grapevines across time and space.
# Samples were collected in summer 2018 and 2019
# Samples consist of leaf, berry, and root compartments
# Code by: Joel F. Swift

# For the data filtering see: Data_filtering_normalization.R
# For general analysis and figures see: CA_transect_analysis.R
# For machine learning see: CA_transect_machine_learning.R
# For testing ML parameters see: CA_ML_testing.R
# For DEseq2 diff. abund. see: CA_transect_differential_abundance.R

library("reshape2"); packageVersion("reshape2")
library("tidyverse"); packageVersion("tidyverse")
library("ggpubr"); packageVersion("ggpubr")
theme_set(theme_pubr())


delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
  x.list[unlist(lapply(x.list, length) != 0)]
} 


# Load dataset from Jupyter notebook
# Constructed by running a random forest model given a specific number of trees
# 1 to 501 by 2.
load("NumTree_OBBerror_6models.Rda")

datalist <- delete.NULLs(datalist)
df_all <- data.frame(matrix(unlist(datalist), nrow=length(datalist), byrow=T))
# Add column names, add a column to contain the number of trees used
colnames(df_all) <- c("Site_Model", "Rootstock_Model", "Compartment_Model", "Year_Model", "Scion_Model", "Sugar_content_Model")
df_all$trees <- seq(1,501,2)

# Look at quantiles and identify the number of trees that results
# in the lowest OOB Error Estimate (i.e. the highest accuracy).
summary(df_all)

seq(1,501,2)[119]


which.min(df_all$Site_Model) #119
seq(1,501,2)[119] # Trees = 237

which.min(df_all$Rootstock_Model) # 128
seq(1,501,2)[128] # Trees = 255

which.min(df_all$Compartment_Model) #32
seq(1,501,2)[32] # Trees = 63

which.min(df_all$Year_Model) #239
seq(1,501,2)[239] # Trees = 477

which.min(df_all$Scion_Model) #217
seq(1,501,2)[217] # Trees = 433

which.min(df_all$Sugar_content_Model) #134
seq(1,501,2)[134] # Trees = 267

# Melt df into format for use in ggplot2
# Convert trees column from character to numeric
new_df <- melt(df_all, id=c("trees"))
new_df$trees <- as.numeric(new_df$trees)

# Plot
p <- ggplot(new_df) + 
  geom_line(aes(x=trees, y=value, color = variable, group = variable), size = 1) +
  scale_color_manual("Model", values = c("red", "blue", "purple", "black", "green", "orange")) +
  scale_x_continuous(breaks = seq(1,501,100), expand = c(0,0)) +
  xlab("Number of Trees") +
  ylab("OOB Error Estimate") +
  geom_vline(xintercept=237, linetype='dashed', col = 'red') +
  geom_vline(xintercept=255, linetype='dashed', col = 'blue') +
  geom_vline(xintercept=63, linetype='dashed', col = 'purple') +
  geom_vline(xintercept=477, linetype='dashed', col = 'black') +
  geom_vline(xintercept=433, linetype='dashed', col = 'green') +
  geom_vline(xintercept=267, linetype='dashed', col = 'orange') +
  theme(legend.position = "right")

# Save plot
ggsave(filename = "Number_of_trees_plot.pdf", plot = p, units = "in", height = 6, width = 12)
ggsave(filename = "Number_of_trees_plot.png", plot = p, units = "in", height = 6, width = 12)