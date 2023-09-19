# Analysis of grafted grapevines across time and space.
# Samples were collected in summer 2018 and 2019
# Samples consist of leaf, berry, and root compartments
# Code by: Joel F. Swift

# For the data filtering see: Data_filtering_normalization.R
# For general analysis and figures see: CA_transect_analysis.R
# For machine learning see: CA_transect_machine_learning.R
# For DEseq2 diff. abund. see: CA_transect_differential_abundance.R

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('ggpubr'); packageVersion('ggpubr')
library('car'); packageVersion('car')
library('lme4'); packageVersion('lme4')
library('lmerTest'); packageVersion('lmerTest')
library('vegan'); packageVersion('vegan')
library('phyloseq'); packageVersion('phyloseq')
library('qiime2R'); packageVersion('qiime2R')
library('DESeq2'); packageVersion('DESeq2')
library('emmeans'); packageVersion('emmeans')
library('lubridate'); packageVersion('lubridate')
library('rebus'); packageVersion('rebus')
library('viridis'); packageVersion('viridis')
library('corrplot'); packageVersion('corrplot')
library('factoextra'); packageVersion('factoextra')

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
##### Functions #####
# Calculate shannon, inverse simpson, and faith's phylogeny diversity metrics from a phyloseq object
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

# Function by ZNH to id outliers past a defined distance.
is_gt <- function(object, dist, threshold){
    samples <- rownames(object)[dist > threshold | dist < -threshold]
    return(samples)
}

# Plot PCoA using ggplot from phyloseq ordination
PLOT_PCoA <- function(plot_data, distance_matrix, axis1, axis2, split_by){
    if (split_by == 'site'){
        temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
        temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
        Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
            geom_point(aes(fill=lower_political, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
            scale_shape_manual(values=c(22, 24, 21), labels = c("Berry", "Leaf", "Root")) +  
            scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
            labs(shape= "Compartment", color= "Site") +
            xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
            ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
            guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))
        return(Plot)
    } else if(split_by == 'rootstock'){
        temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
        temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
        Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
            geom_point(aes(fill=rootstock, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
            scale_shape_manual(values=c(22, 24, 21)) +  
            scale_fill_manual(name = "Rootstock", values=rootstock_palette) +
            labs(shape= "Compartment", color= "Site") +
            xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
            ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
            guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))
        return(Plot)
    } else if(split_by == 'scion'){
        temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
        temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
        Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
            geom_point(aes(fill=scion, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
            scale_shape_manual(values=c(22, 24, 21)) +  
            scale_fill_manual(name = "Scion", values=scion_palette) +
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
            scale_shape_manual(values=c(22, 24, 21)) +  
            #scale_fill_manual(name = "Scion", values=scion_palette) +
            labs(shape= "Compartment", color= "Site") +
            xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
            ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
            guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))
        return(Plot)
    }  else if(split_by == 'col_week'){
        temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
        temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
        Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
            geom_point(aes(fill=col_week, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
            scale_shape_manual(values=c(22, 24, 21)) +  
            #scale_fill_manual(name = "Scion", values=scion_palette) +
            labs(shape= "Compartment", color= "Site") +
            xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
            ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
            guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))
        return(Plot)
    }   else if(split_by == 'block'){
        temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
        temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
        Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
            geom_point(aes(fill=block, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
            scale_shape_manual(values=c(22, 24, 21)) +  
            #scale_fill_manual(name = "Scion", values=scion_palette) +
            labs(shape= "Compartment", color= "Site") +
            xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
            ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
            guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "black")))
        return(Plot)
    }
}

# Given a phyloseq object, return a dataframe of the sample metadata 
# From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
pssd2veg <- function(physeq) {
    sd <- sample_data(physeq)
    return(as(sd,"data.frame"))
}

# load dataset generated from Data_filtering_normalization.R
phy_vst <- readRDS("../Data_files/phyloseq_16s_no_soil_filtered_vst_dataset.rds")

##### 1.0) ALPHA DIVERSITY ANALYSIS #####
alpha_diversity.df <- Alpha_div_metrics(phy_vst)
# Lmer models
fphd_16s_fit.mod <- lmerTest::lmer(Faithpd ~ rootstock*scion + plant_body_site*brix + year + site + (1|sequence_depth), data = alpha_diversity.df)
simpI_16s_fit.mod <- lmerTest::lmer(Invsimpson ~ rootstock*scion + plant_body_site*brix + year + site + (1|sequence_depth), data = alpha_diversity.df)
chao1_16s_fit.mod <- lmerTest::lmer(Chao1 ~ rootstock*scion + plant_body_site*brix + year + site + (1|sequence_depth), data = alpha_diversity.df)
# anovas, paste these into excel
anova(fphd_16s_fit.mod)
anova(simpI_16s_fit.mod)
anova(chao1_16s_fit.mod)
# check the random effects
rand(fphd_16s_fit.mod)
rand(simpI_16s_fit.mod)
rand(chao1_16s_fit.mod)
# posthocs
pairs(emmeans(fphd_16s_fit.mod, ~ rootstock))
pairs(emmeans(fphd_16s_fit.mod, ~ plant_body_site))
pairs(emmeans(fphd_16s_fit.mod, ~ year))
pairs(emmeans(fphd_16s_fit.mod, ~ site))
pairs(emmeans(fphd_16s_fit.mod, ~ rootstock|scion))
pairs(emmeans(fphd_16s_fit.mod, ~ rootstock|scion))


pairs(emmeans(simpI_16s_fit.mod, ~ plant_body_site))
pairs(emmeans(simpI_16s_fit.mod, ~ year))
pairs(emmeans(simpI_16s_fit.mod, ~ site))
pairs(emmeans(simpI_16s_fit.mod, ~ rootstock|scion))

pairs(emmeans(chao1_16s_fit.mod, ~ rootstock))
pairs(emmeans(chao1_16s_fit.mod, ~ plant_body_site))
pairs(emmeans(chao1_16s_fit.mod, ~ year))
pairs(emmeans(chao1_16s_fit.mod, ~ site))
pairs(emmeans(chao1_16s_fit.mod, ~ rootstock|scion))

##### Supplemental Figure 2 #####
# Discretizing sugar content values
Brix_values <- as.data.frame(phy_vst@sam_data)

Brix_spilt_plot <- ggplot(Brix_values[Brix_values$plant_body_site == "berry",], aes(x = 1:184, y = brix)) + 
  geom_point()  + geom_hline(yintercept = 7, color = "red", linetype = 'dashed') +
  xlab("Index") + ylab("Sugar content (\u00B0Bx)")

ggsave("Brix_discret_split_plot.svg", Brix_spilt_plot, bg = 'white', dpi = 600, height = 8, width = 6)
ggsave("Brix_discret_split_plot.pdf", Brix_spilt_plot, bg = 'white', dpi = 600, height = 8, width = 6)

###### Figure 1 B-D ######
# Load Brix data set from csv file
Brixdf <- read.csv("../Data_files/Brix_data_CA_2018-2019.csv")
# Convert year to a factor
Brixdf$Year <- as.factor(Brixdf$Year)
# Convert date strings to an R accepted format year-month-day
Brixdf$Sample.Date <- dmy(Brixdf$Sample.Date)

# 2018/2019 scion panels
a <- ggplot(Brixdf[Brixdf$Year == "2018",], aes(x=Sample.Date, y=brix, color = Scion)) + 
    annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-06-19"), xmax=as.Date("2018-06-21"), ymin=-Inf, ymax=Inf) + 
    annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-10"), xmax=as.Date("2018-07-13"), ymin=-Inf, ymax=Inf) + 
    annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-31"), xmax=as.Date("2018-08-02"), ymin=-Inf, ymax=Inf) +
    geom_point(aes(fill = Scion), shape = 21, size = 3, color = 'black', alpha = 0.8) +
    geom_smooth(method = "loess", level=0.95, size = 1, se=FALSE) +
    scale_color_manual(values = scion_palette) +
    scale_fill_manual(values = scion_palette) +
    scale_x_date(limits = as.Date(c("2018-06-15","2018-08-10"))) +
    #\u00B0 makes the degree symbol#
    ylab("Sugar content (\u00B0Bx)") +
    ylim(3, 23) +
    theme(legend.position = "right", axis.title.x=element_blank()) +
    ggtitle("2018")
b <- ggplot(Brixdf[Brixdf$Year == "2019",], aes(x=Sample.Date, y=brix, color = Scion)) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-06-25"), xmax=as.Date("2019-06-27"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-07-23"), xmax=as.Date("2019-07-25"), ymin=-Inf, ymax=Inf) + 
  geom_point(aes(fill = Scion), shape = 21, size = 3, color = 'black', alpha = 0.8) +
  geom_smooth(method = "loess", level=0.95, size = 1, se=FALSE) +
  scale_color_manual(values = scion_palette) +
  scale_fill_manual(values = scion_palette) +
  scale_x_date(limits = as.Date(c("2019-06-15","2019-08-10"))) +
  #\u00B0 makes the degree symbol#
  ylab("Sugar content (\u00B0Bx)") +
  ylim(3, 23) +
  theme(legend.position = "right",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank()) +
  ggtitle("2019")

part1 <- ggarrange(a,b, ncol = 2, align = 'hv', common.legend = TRUE)

# 2018/2019 site panels
c <- ggplot(Brixdf[Brixdf$Year == "2018",], aes(x=Sample.Date, y=brix, color = Sample.Location)) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-06-19"), xmax=as.Date("2018-06-21"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-10"), xmax=as.Date("2018-07-13"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-31"), xmax=as.Date("2018-08-02"), ymin=-Inf, ymax=Inf) +
  geom_point(aes(fill = Sample.Location), shape = 21, size = 3, color = 'black', alpha = 0.8) +
  geom_smooth(method = "loess", level=0.95, size = 1, se=FALSE) +
  scale_color_manual(name = "Site", values = site_palette) +
  scale_fill_manual(name = "Site", values = site_palette) +
  scale_x_date(limits = as.Date(c("2018-06-15","2018-08-10"))) +
  #\u00B0 makes the degree symbol#
  ylab("Sugar content (\u00B0Bx)") +
  ylim(3, 23) +
  theme(legend.position = "right", axis.title.x=element_blank()) +
  ggtitle("2018")
d <- ggplot(Brixdf[Brixdf$Year == "2019",], aes(x=Sample.Date, y=brix, color = Sample.Location)) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-06-25"), xmax=as.Date("2019-06-27"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-07-23"), xmax=as.Date("2019-07-25"), ymin=-Inf, ymax=Inf) + 
  geom_point(aes(fill = Sample.Location), shape = 21, size = 3, color = 'black', alpha = 0.8) +
  geom_smooth(method = "loess", level=0.95, size = 1, se=FALSE) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  scale_x_date(limits = as.Date(c("2019-06-15","2019-08-10"))) +
  #\u00B0 makes the degree symbol#
  ylab("Sugar content (\u00B0Bx)") +
  ylim(3, 23) +
  theme(legend.position = "right",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank()) +
  ggtitle("2019")

part2 <- ggarrange(c,d, ncol = 2, align = 'hv', common.legend = TRUE)

# Combine and save
fig1 <- ggarrange(part1, part2, nrow = 2)
ggsave("figure1_brix_line_plots.svg", fig1, height = 9, width = 9)
ggsave("figure1_brix_line_plots.pdf", fig1, height = 9, width = 9)

# Brix statistics | Berry only
alpha_berries <- alpha_diversity.df[alpha_diversity.df$plant_body_site == 'berry',]
brix_berry.mod <- lm(brix ~ rootstock + lower_political + scion*col_week, data = alpha_berries) # Model brix data with rootstock, scion, and scion*col_week
car::Anova(brix_berry.mod , type = "III") # anova type III because of interaction
pairs(emmeans(brix_berry.mod , ~ lower_political)) # San joaquin lower than others.
pairs(emmeans(brix_berry.mod , ~ scion|col_week)) # chard and cab diverge later in season
shan_berry_brix.mod <- lmerTest::lmer(Shannon ~ rootstock + scion + site + brix + brix:year + (1|sequence_depth) + (1|seq_plate), data = alpha_berries)
performance::check_collinearity(shan_berry_brix.mod)
anova(shan_berry_brix.mod)
pairs(emtrends(shan_berry_brix.mod, 'year', var = "brix"))
# Roots
alpha_root <- alpha_diversity.df[alpha_diversity.df$plant_body_site == 'root',]
ggplot(alpha_root, aes(x=brix, y=Shannon, fill=year)) + geom_point(shape=21, size =3) + geom_smooth(method = 'lm') #+ facet_wrap(~col_week)
# Leaves
alpha_leaf <- alpha_diversity.df[alpha_diversity.df$plant_body_site == 'leaf',]
ggplot(alpha_leaf, aes(x=brix, y=Shannon, fill=year)) + geom_point(shape=21, size =3) + geom_smooth(method = 'lm') #+ facet_wrap(~col_week)


####### Figure 3 ######
# alpha diversity (Faith across compartments)
# Make a tibble to store where to place significance letters from the Tukey test
labels_df <- tibble(plant_body_site=levels(alpha_diversity.df$plant_body_site), rootstock=levels(alpha_diversity.df$rootstock), Mfaith=max(alpha_diversity.df$Faithpd) * 1.2)
TukeyHSD(aov(Faithpd ~ plant_body_site, data = alpha_diversity.df), conf.level = 0.95)
fig3A <- ggplot(alpha_diversity.df, aes(plant_body_site, Faithpd, fill = plant_body_site)) +
    geom_jitter(width = 0.25, color = "black") +
    xlab ("Compartment") +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
    theme(legend.position = "top") +
    ylab("Faith's phylogenetic diversity") +
    scale_fill_manual(name="Compartment", breaks=c("berry", "leaf", "root"),
                      labels=c("Berry","Leaf", "Root"), values = compartment_palette) +
    scale_x_discrete(breaks=c("berry", "leaf", "root"), labels=c("Berry", "Leaf", "Root")) +
    geom_text(data=labels_df, aes(plant_body_site, Mfaith, label=c("a","a","b")), size = 6) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# Stats mean faith for compartments for main text
mean(alpha_diversity.df[alpha_diversity.df$plant_body_site == 'berry',]$Faithpd) #4.633959
mean(alpha_diversity.df[alpha_diversity.df$plant_body_site == 'leaf',]$Faithpd) #5.709501
mean(alpha_diversity.df[alpha_diversity.df$plant_body_site == 'root',]$Faithpd) #61.06594
# Boxplot of scion X rootstock X compartment faith
facet_label <- c("Cabernet Sauvignon", "Chardonnay")
names(facet_label) <- c("cabernet sauvignon", "chardonnay")

TukeyHSD(aov(Faithpd ~ rootstock*scion, data = alpha_diversity.df[alpha_diversity.df$plant_body_site == "root",]), conf.level = 0.95)

fig3B <- ggplot(alpha_diversity.df[alpha_diversity.df$plant_body_site == "root",], aes(x=rootstock, y=Faithpd,  fill=rootstock)) +
    geom_jitter(width = 0.2, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") + xlab ("Rootstock") + ylab("Faith's phylogenetic diversity") +
    scale_fill_manual(name = "Rootstock", values=rootstock_palette) +
    geom_text(data=labels_df, aes(rootstock, Mfaith, label=c("ab","a","a","a","a","b")), size = 6) + ylim(0, NA) +
    facet_wrap(~scion, labeller = labeller(scion = facet_label)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


fig3 <- ggarrange(fig3A,fig3B, labels = c('A', 'B'))
# save
ggsave("figure3_faith_by_compartment_and_rootsock-scion.svg", fig3, height = 6, width = 10, bg = 'white', dpi = 1200)
ggsave("figure3_faith_by_compartment_and_rootsock-scion.pdf", fig3, height = 6, width = 10)
# Stats 
# We found that Teleki 5C showed lower alpha diversity (faith) than the other rootstocks,
# this was particularly apparent in the root samples. Below are some summary stats for the 
# rootstocks. It seems that at Merced in particular the values are much lower than the mean 
# for other rootstocks and even the scion/rootstock combo at another site.
summarySE(alpha_diversity.df, measurevar = "Faithpd", groupvars = c("plant_body_site", "rootstock", "scion"))
summarySE(alpha_diversity.df, measurevar = "Faithpd", groupvars = c("plant_body_site", "rootstock", "scion", "site"))

##### 2.0) BETA DIVERSITY ANALYSIS #####
# Calculate bray-curtis
out.bray <- ordinate(phy_vst, method = "MDS", distance = "bray")
# Calculate distance matrixes for use in Vegan
out.dist.bray <- phyloseq::distance(phy_vst, method = "bray")
# get variance explained by axes 1-3
sum(out.bray$values$Eigenvalues[1:3])/sum(out.bray$values$Eigenvalues) #27%
# Plot PCoAs 1x2-3 for bray by site
pcoa_bray_1_2 <- PLOT_PCoA(phy_vst, out.bray, 1, 2, split_by = 'site')
pcoa_bray_1_3 <- PLOT_PCoA(phy_vst, out.bray, 1, 3, split_by = 'site')
Bray123_pcoa <- ggarrange(pcoa_bray_1_2, pcoa_bray_1_3, common.legend = TRUE, legend = "right", labels = c("Bray-Curtis Dissimilarity"," "))
ggsave("Bray-Curtis PCoA by site.pdf", Bray123_pcoa, height = 8, width = 16)

### Adonis2 ###
# Create sample data dataframe from phyloseq object to use in vegan/adonis
phy_vst_metadata <- pssd2veg(phy_vst)
# Adonis Full model
# This shows significance for all terms but is kind of overpowered by the root samples
# Running adonis2 with compartment separate as this is the largest factor in the dataset.
# adonis2(out.dist.bray ~ rootstock*scion + plant_body_site*brix + year  + site, data = phy_vst_metadata, by = "margin", permutations = 1000)
# Adonis split by compartment
#Berry Adonis
phy_vst_berry <- subset_samples(phy_vst, plant_body_site == 'berry')
out.dist.bray.berry <- phyloseq::distance(phy_vst_berry, method = "bray")
phy_vst_metadata_berry <- pssd2veg(phy_vst_berry)
adonis2(out.dist.bray.berry ~ rootstock + scion + brix + year + site, data = phy_vst_metadata_berry, by = "margin", permutations = 1000)
# Leaf Adonis
phy_vst_leaf <- subset_samples(phy_vst, plant_body_site == 'leaf')
out.dist.bray.leaf <- phyloseq::distance(phy_vst_leaf, method = "bray")
phy_vst_metadata_leaf <- pssd2veg(phy_vst_leaf)
adonis2(out.dist.bray.leaf ~ rootstock + scion + brix + year + site, data = phy_vst_metadata_leaf, by = "margin", permutations = 1000)
# Root Adonis
phy_vst_root <- subset_samples(phy_vst, plant_body_site == 'root')
out.dist.bray.root <- phyloseq::distance(phy_vst_root, method = "bray")
phy_vst_metadata_root <- pssd2veg(phy_vst_root)
adonis2(out.dist.bray.root ~ rootstock + scion + brix + year + site, data = phy_vst_metadata_root,by = "margin", permutations = 1000)

# Average Bray-Curtis distance
se <- function(x) sqrt(var(x)/length(x))

mean(out.dist.bray)
mean(out.dist.bray.root)
mean(out.dist.bray.berry)
mean(out.dist.bray.leaf)

se(out.dist.bray)
se(out.dist.bray.root)
se(out.dist.bray.berry)
se(out.dist.bray.leaf)

##### 3.0) Taxonomic Barplots #####
# Get data from phyloseq object.
phylum_lvl <- tax_glom(phy_vst, taxrank = "Phylum") # 42 taxa
class_lvl <- tax_glom(phy_vst, taxrank = "Class") # 107 taxa
order_lvl <- tax_glom(phy_vst, taxrank = "Order") # 265 taxa
family_lvl <- tax_glom(phy_vst, taxrank = "Family") # 430 taxa
genus_lvl <- tax_glom(phy_vst, taxrank = "Genus") # 800 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
phylum_lvl <- transform_sample_counts(phylum_lvl, function(x) x/sum(x))
class_lvl <- transform_sample_counts(class_lvl, function(x) x/sum(x))
order_lvl <- transform_sample_counts(order_lvl, function(x) x/sum(x))
family_lvl <- transform_sample_counts(family_lvl, function(x) x/sum(x))
genus_lvl <- transform_sample_counts(genus_lvl, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
phylum_lvl<- fantaxtic::get_top_taxa(physeq_obj = phylum_lvl, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
class_lvl<- fantaxtic::get_top_taxa(physeq_obj = class_lvl, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
order_lvl<- fantaxtic::get_top_taxa(physeq_obj = order_lvl, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
family_lvl<- fantaxtic::get_top_taxa(physeq_obj = family_lvl, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
genus_lvl<- fantaxtic::get_top_taxa(physeq_obj = genus_lvl, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Melt data frame with phyloseq function for plotting.
phylum_lvl <- psmelt(phylum_lvl)
class_lvl <- psmelt(class_lvl)
order_lvl <- psmelt(order_lvl)
family_lvl <- psmelt(family_lvl)
genus_lvl <- psmelt(genus_lvl)
# Reorder levels to put other to the end, other in alphabetical
phylum_lvl$Phylum <- forcats::fct_relevel(as.factor(phylum_lvl$Phylum), "Other", after = Inf)
class_lvl$Class <- forcats::fct_relevel(as.factor(class_lvl$Class), "Other", after = Inf)
order_lvl$Order <- forcats::fct_relevel(as.factor(order_lvl$Order), "Other", after = Inf)
family_lvl$Family <- forcats::fct_relevel(as.factor(family_lvl$Family), "Other", after = Inf)
genus_lvl$Genus <- forcats::fct_relevel(as.factor(genus_lvl$Genus), "Other", after = Inf)
# phylum
phylum_compartment_by_year <- ggplot(phylum_lvl, aes(x=Sample, y=Abundance, fill = Phylum)) +
    geom_bar(stat = "identity") +
    ylab("Relative abundance") +
    xlab("sample") +
    facet_wrap(~plant_body_site:year, scale = "free", ncol = 2) +
    scale_fill_manual(values=safe_colorblind_palette) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = 'right')

ggsave("Top 10 phylum by compartment and year.pdf", phylum_compartment_by_year, height = 24, width = 16)

phylum_compartment_by_site<- ggplot(phylum_lvl, aes(x= Sample, y=Abundance, fill = Phylum)) +
    geom_bar(stat = "identity") +
    ylab("Relative abundance") +
    xlab("sample") +
    facet_wrap(~plant_body_site:lower_political, scale = "free", ncol = 3) +
    scale_fill_manual(values=safe_colorblind_palette) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = 'right')

ggsave("Top 10 phylum by compartment and site.pdf", phylum_compartment_by_site, height = 24, width = 24)

ggplot(phylum_lvl[phylum_lvl$plant_body_site == 'root',], aes(x= Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~lower_political:year, scale = "free", ncol = 3) +
  scale_fill_manual(values=safe_colorblind_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')


###### Figure S4 PCoAs for compartments together ######
# load dataset with soil samples included
phy_with_soil_vst <- readRDS("../Data_files/phyloseq_16s_filtered_vst_dataset.rds")
# ordinations, axes 1/2/3
out.bray_all <- ordinate(phy_with_soil_vst, method = "MDS", distance = "bray")
A1.2 <- plot_ordination(phy_with_soil_vst, out.bray_all, axes = c(1,2)) 
A1.3 <- plot_ordination(phy_with_soil_vst, out.bray_all, axes = c(1,3)) 

a <- A1.2 + geom_point(aes(fill=plant_body_site, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
  scale_shape_manual(name = "Compartment", values = c(22,24,21,23), labels = c("Berry", "Leaf", "Root", "Soil")) +
  scale_fill_manual(name = "Compartment", values= c("#5a1991", "#139d08", "#5c3c0d", "grey"), labels = c("Berry", "Leaf", "Root", "Soil")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y)))

b <- A1.3 + geom_point(aes(fill=plant_body_site, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
  scale_shape_manual(name = "Compartment", values = c(22,24,21,23), labels = c("Berry", "Leaf", "Root", "Soil")) +
  scale_fill_manual(name = "Compartment", values= c("#5a1991", "#139d08", "#5c3c0d", "grey"), labels = c("Berry", "Leaf", "Root", "Soil")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 3", sub(".*\\ ", "", temp$labels$y)))

part1 <- ggarrange(a,b, legend = 'right', common.legend = TRUE, align = 'hv', labels = "AUTO")

c <- A1.2 + geom_point(aes(fill=lower_political, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
  scale_shape_manual(name = "Compartment", values = c(22,24,21,23), labels = c("Berry", "Leaf", "Root", "Soil")) +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)),
         shape = "none")

d <- A1.3 + geom_point(aes(fill=lower_political, shape=plant_body_site), size = 5, alpha = 0.80, color = "black") +
  scale_shape_manual(name = "Compartment", values = c(22,24,21,23), labels = c("Berry", "Leaf", "Root", "Soil")) +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 3", sub(".*\\ ", "", temp$labels$y))) +
  guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)),
         shape = "none")

part2 <- ggarrange(c,d, legend = 'right', common.legend = TRUE, align = 'hv', labels = c('C', 'D'))

# combine plots and save
PCoAs_compartments_together <- ggarrange(part1,part2, nrow = 2)
ggsave("FigureS4_experimental_factors_pcoas.svg", PCoAs_compartments_together, width = 12, height = 12)


###### Figure S5 PCoAs for each experimental factor by compartment  ###### 
# PCoA plots for each microbiome separately
# Phyloseq data objects subset from above
phy_vst_leaf
phy_vst_berry
phy_vst_root
# Ordinations based off the subset objects
out.bray_leaf <- ordinate(phy_vst_leaf, method = "MDS", distance = "bray")
out.bray_berry <- ordinate(phy_vst_berry, method = "MDS", distance = "bray")
out.bray_root <- ordinate(phy_vst_root, method = "MDS", distance = "bray")

# Berry PCoAs
temp <- plot_ordination(phy_vst_berry, out.bray_berry, axes = c(1,2)) 
# in order rootstock, scion, year, site, brix
a <-  temp +
  geom_point(aes(fill=rootstock, shape=plant_body_site), shape = 22, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Rootstock", values=rootstock_palette) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y)))
b <-  temp +
  geom_point(aes(fill=scion, shape=plant_body_site), shape = 22, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Scion", values=scion_palette, labels = c("Cabernet Sauvignon", "Chardonnay")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y)))
c <-  temp +
  geom_point(aes(fill=year, shape=plant_body_site), shape = 22, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Year", values = c("grey50", "white")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y)))
d <-  temp +
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 22, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y)))
e <-  temp +
  geom_point(aes(fill=brix, shape=plant_body_site), shape = 22, size = 5, alpha = 0.95, color = "black") +
  scale_fill_gradient(low = "#c8ed8c", high = "#722F37", name = "Sugar content (\u00B0Bx)") +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y)))

berry_pcoas <- ggarrange(a,b,c,d,e, nrow = 1, align = 'hv', labels = c("A"))

# Leaf PCoAs
temp <- plot_ordination(phy_vst_leaf, out.bray_leaf, axes = c(1,2)) 
# in order rootstock, scion, year, site, brix
a <-  temp +
  geom_point(aes(fill=rootstock, shape=plant_body_site), shape = 24, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Rootstock", values=rootstock_palette) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
b <-  temp +
  geom_point(aes(fill=scion, shape=plant_body_site), shape = 24, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Scion", values=scion_palette, labels = c("Cabernet Sauvignon", "Chardonnay")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
c <-  temp +
  geom_point(aes(fill=year, shape=plant_body_site), shape = 24, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Year", values = c("grey50", "white")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
d <-  temp +
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 24, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
e <-  temp +
  geom_point(aes(fill=brix, shape=plant_body_site), shape = 24, size = 5, alpha = 0.95, color = "black") +
  scale_fill_gradient(low = "#c8ed8c", high = "#722F37", name = "Sugar content (\u00B0Bx)") +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")

leaf_pcoas <- ggarrange(a,b,c,d,e, nrow = 1, align = 'hv', labels = c("B"))

# Root PCoAs
temp <- plot_ordination(phy_vst_root, out.bray_root, axes = c(1,2)) 
# in order rootstock, scion, year, site, brix
a <-  temp +
  geom_point(aes(fill=rootstock, shape=plant_body_site), shape = 21, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Rootstock", values=rootstock_palette) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
b <-  temp +
  geom_point(aes(fill=scion, shape=plant_body_site), shape = 21, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Scion", values=scion_palette) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
c <-  temp +
  geom_point(aes(fill=year, shape=plant_body_site), shape = 21, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Year", values = c("grey50", "white")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
d <-  temp +
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 21, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
e <-  temp +
  geom_point(aes(fill=brix, shape=plant_body_site), shape = 21, size = 5, alpha = 0.95, color = "black") +
  scale_fill_gradient(low = "#c8ed8c", high = "#722F37", name = "Sugar content (\u00B0Bx)") +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")

root_pcoas <- ggarrange(a,b,c,d,e, nrow = 1, align = 'hv', labels = c("C"))

# combine together compartment plots and save
PCoAs_all_factors <- ggarrange(berry_pcoas, leaf_pcoas, root_pcoas, nrow = 3)
ggsave("FigureS5_experimental_factors_pcoas.svg", PCoAs_all_factors, width = 22, height = 12)








###### Figure 4 ###### 
### Panel A and B
### Bray curtis PCoA with all compartments
### Shapes are compartment and fill is sites
pcoa_bray_1_2 <- PLOT_PCoA(phy_vst, out.bray, 1, 2, split_by = 'site')
pcoa_bray_1_3 <- PLOT_PCoA(phy_vst, out.bray, 1, 3, split_by = 'site')
Bray123_pcoa <- ggarrange(pcoa_bray_1_2, pcoa_bray_1_3, ncol = 2, common.legend = TRUE, legend = "right", labels = c("A", "B"))
### Panel C
### Taxonomic barplots for root samples split by site
### Bars are relative abundance of the top ten phyla
phylum_lvl <- tax_glom(phy_vst, taxrank = "Phylum") # 42 taxa
phylum_lvl <- subset_samples(phylum_lvl, extraction_num != "492" & extraction_num != "615") # Remove two outlier samples
phylum_lvl <- transform_sample_counts(phylum_lvl, function(x) x/sum(x))
phylum_lvl <- fantaxtic::get_top_taxa(physeq_obj = phylum_lvl, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
phylum_lvl <- psmelt(phylum_lvl)
phylum_lvl$Phylum <- forcats::fct_relevel(as.factor(phylum_lvl$Phylum), "Other", after = Inf)
ylabels <- c('madera'='Madera', 'merced' = 'Merced', 'san joaquin' = "San Joaquin")
phylum_root_site <- ggplot(phylum_lvl[phylum_lvl$plant_body_site == "root", ], aes(x= Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~lower_political, scale = "free", ncol = 3, labeller = as_labeller(ylabels)) +
  scale_fill_manual(values=safe_colorblind_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')
###  panel D
### This is a heat map with each of the top 10 phyla + the other category
### The heat map will show the % variance explained after correcting for multiple
### testing across all of the models together (BH, Benjamini & Hochberg).
# Get data from phyloseq object.
phylum_lvl_panelD <- tax_glom(phy_vst, taxrank = "Phylum") # 42 taxa
phylum_lvl_panelD <- subset_samples(phylum_lvl_panelD, extraction_num != "492" & extraction_num != "615")
# Make relative abundance AKA transform sample counts by diving by the ASV total.
phylum_lvl_panelD <- transform_sample_counts(phylum_lvl_panelD, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
phylum_lvl_panelD <- fantaxtic::get_top_taxa(physeq_obj = phylum_lvl_panelD, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Instead of using psmelt, instead keep data in wide format (taxa as rows)
Phylum_relative_abund <- as.data.frame(t(otu_table(phylum_lvl_panelD)))
colnames(Phylum_relative_abund) <- as.data.frame(phylum_lvl_panelD@tax_table)$Phylum # Name the columns their by phyla
Phylum_w_meta <- cbind(Phylum_relative_abund, pssd2veg(phylum_lvl_panelD))# Attach metadata
# Firmicutes and Deinococcota average relative abundance.
sum(Phylum_w_meta$Deinococcota * 100) / 592 # 1.08%
sum(Phylum_w_meta$Firmicutes * 100) / 592 # 5.45%
# Linear models but only for the root samples
model_ot <- lm(Other ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_my <- lm(Myxococcota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_pr <- lm(Proteobacteria ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_ch <- lm(Chloroflexi ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_de <- lm(Deinococcota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_fi <- lm(Firmicutes ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_act <- lm(Actinobacteriota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_ba <- lm(Bacteroidota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_pl <- lm(Planctomycetota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_ve <- lm(Verrucomicrobiota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
model_aci <- lm(Acidobacteriota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "root",])
# San Joaquin different than others
pairs(emmeans(model_pl,~site))
pairs(emmeans(model_ba,~site))
pairs(emmeans(model_aci,~site))
# Merced different than others
pairs(emmeans(model_my,~site))
pairs(emmeans(model_pr,~site))
# Teleki 5c different than others
pairs(emmeans(model_ch,~rootstock))
pairs(emmeans(model_pl,~rootstock))
pairs(emmeans(model_my,~rootstock))
pairs(emmeans(model_ve,~rootstock))
# Higher actinobacteria in 2019 v 2018
pairs(emmeans(model_act, ~year))
### Heat map of LM results ###
Return_SS_proportion <- function(model){
  SS <- c(anova(model)$Sum)
  output <- (SS / sum(SS)) * 100
  output <- as.data.frame(output) 
  rownames(output) <- c(rownames(anova(model)))
  output <- t(output)
  return(as.data.frame(output))}
# Make list of ion names in order and run the function to generate the % variance
taxa_list <- as.data.frame((phylum_lvl_panelD@tax_table))$Phylum
Percent_var <- bind_rows(Return_SS_proportion(model_ot),
                         Return_SS_proportion(model_my),
                         Return_SS_proportion(model_pr),
                         Return_SS_proportion(model_ch),
                         Return_SS_proportion(model_de),
                         Return_SS_proportion(model_fi),
                         Return_SS_proportion(model_act),
                         Return_SS_proportion(model_ba),
                         Return_SS_proportion(model_pl),
                         Return_SS_proportion(model_ve),
                         Return_SS_proportion(model_aci))
rownames(Percent_var) <- taxa_list
# remove residual 
Percent_var_noint_nores <- Percent_var[, !names(Percent_var) %in% "Residuals"]
# Make dataframe of P values to connect to the % var df above
Get_Pval<- function(model){
  Pval<- as.data.frame(anova(model)$`Pr(>F)`)
  rownames(Pval) <- paste(rownames(anova(model)), "p", sep = "_")
  Pval<- t(Pval)
  return(as.data.frame(Pval))}
# Run function above on all anova and bind output by row
Pval <- bind_rows(Get_Pval(model_ot),
                  Get_Pval(model_my),
                  Get_Pval(model_pr),
                  Get_Pval(model_ch),
                  Get_Pval(model_de),
                  Get_Pval(model_fi),
                  Get_Pval(model_act),
                  Get_Pval(model_ba),
                  Get_Pval(model_pl),
                  Get_Pval(model_ve),
                  Get_Pval(model_aci))
rownames(Pval) <- taxa_list
# Combine
SS_pval_combo_df <- cbind(Percent_var_noint_nores,Pval)
SS_pval_combo_df$taxa <- rownames(SS_pval_combo_df)
# Get variance columns
total_var <- SS_pval_combo_df %>% dplyr::select(taxa, "rootstock":"site")
# Reorganize and rename to format for plotting
total_var <- total_var %>% gather(key=factor, value=var, -taxa)
total_var <- total_var %>% mutate(factor=str_replace(factor, "_var$", ""))
# Get p-value columns
total_p_name <- SS_pval_combo_df %>% dplyr::select("rootstock_p":"site_p")
total_p <- data.frame(t(apply(total_p_name, 1, FUN=p.adjust, method='BH')))
colnames(total_p)<- colnames(total_p_name) #This is hacky but I need to preserve the colnames to correctly join them later on.
total_p$taxa <- SS_pval_combo_df$taxa
# Reorganize and rename to format for plotting
total_p <- total_p %>% gather(key=factor, value=p_value, -taxa)
total_p <- total_p %>% mutate(factor=str_replace(factor,"_p$", ""))
# Join variance and p-value tables back together 
total_var_p <- full_join(total_var, total_p,by=c("taxa", "factor"))
# Have to remove Deinococcota and Firmicutes prior to adding number for Y-axis plotting
# This is because they have no significant effects of the factors in the experiment
total_var_p <- total_var_p[total_var_p$taxa != "Deinococcota", ]
total_var_p <- total_var_p[total_var_p$taxa != "Firmicutes", ]
# Reorder the taxa to be in alphabetical order
total_var_p$element <- forcats::fct_relevel(as.factor(total_var_p$taxa), "Other", after = Inf)
# Taxa need a number for y-axis when plotting, there are 11 taxa - 2 with no sifinicant comparisons
total_var_p <- total_var_p %>% arrange(element) %>% mutate(element_number=rep(1:9, each=5))
# Only plot significant p_values
total_var_p_sig <- total_var_p %>% filter(p_value < 0.05) 
# Set factor levels and make plot
total_var_p_sig$factor <- factor(total_var_p_sig$factor, levels=c('rootstock', 'scion', 'year', 'brix', 'site'))
Lm_heat_map_var_exp <- ggplot(data=total_var_p_sig, aes(x=factor, y=element_number, fill=var)) + 
  geom_tile(aes(fill=var), color="white", size=1)+
  geom_text(aes(label=paste(plyr::desc(-round(var, digits=2)),"%")), color="black", fontface="bold", size = 3.5) +
  labs(x = "Factor", y="Phylum") + 
  theme(axis.text=element_text(size=12, colour="black"),axis.title=element_text(size=14,face="bold", colour="black"),legend.position = "bottom")+
  scale_y_continuous(position = "right", trans="reverse", breaks = seq(1, 9, 1), minor_breaks=NULL, labels=unique(total_var_p_sig[,"taxa"]), expand = c(0, 0)) +
  scale_x_discrete(position = "top", expand = c(0, 0), labels= c('Rootstock', 'Scion', 'Year', 'Brix', 'Site')) +
  scale_fill_viridis(option="plasma", name="% variance explained", direction = -1)
temp <- ggarrange(phylum_root_site, Lm_heat_map_var_exp, ncol = 2, widths = c(0.75,0.25), labels = c("C","D"))
# Save it
Figure4 <- ggarrange(Fig4, temp, nrow = 2, common.legend = FALSE, legend = "right")
ggsave("Figure4_Bray_PCoA_1x2_1x3_Taxonomic_barplot_root_by_site_LM_heatmap.png", Figure4, height = 12, width = 18.7, units = "in")
ggsave("Figure4_Bray_PCoA_1x2_1x3_Taxonomic_barplot_root_by_site_LM_heatmap.pdf", Figure4, height = 12, width = 18.7, units = "in")

####### NEW Panel A/B figure 4 #######
# panel was added with inkscape to preserve space and sizing of other plot elements
# PCoA plots for each microbiome seperately
# Phyloseq data objects subset from above
phy_vst_berry
phy_vst_leaf
phy_vst_root
# Ordinations based off the subset objects
out.bray_berry <- ordinate(phy_vst_berry, method = "MDS", distance = "bray")
out.bray_leaf <- ordinate(phy_vst_leaf, method = "MDS", distance = "bray")
out.bray_root <- ordinate(phy_vst_root, method = "MDS", distance = "bray")
# berry
temp <- plot_ordination(phy_vst_berry, out.bray_berry, axes = c(1,2)) 
a <- temp +
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 22, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
temp <- plot_ordination(phy_vst_berry, out.bray_berry, axes = c(1,3)) 
b <- temp + 
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 22, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 3", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
# leaf
temp <- plot_ordination(phy_vst_leaf, out.bray_leaf, axes = c(1,2)) 
c <- temp +
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 24, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
temp <- plot_ordination(phy_vst_leaf, out.bray_leaf, axes = c(1,3)) 
d <- temp + 
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 24, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 3", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
# root
temp <- plot_ordination(phy_vst_root, out.bray_root, axes = c(1,2)) 
e <- temp +
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 21, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 2", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")
temp <- plot_ordination(phy_vst_root, out.bray_root, axes = c(1,3)) 
f <- temp + 
  geom_point(aes(fill=lower_political, shape=plant_body_site), shape = 21, size = 5, alpha = 0.95, color = "black") +
  scale_fill_manual(name = "Site", values=site_palette, labels = c("Madera", "Merced", "San Joaquin")) +
  xlab(paste("PCoA 1", sub(".*\\ ", "", temp$labels$x))) +
  ylab(paste("PCoA 3", sub(".*\\ ", "", temp$labels$y))) +
  theme(legend.position = "none")

part1 <- ggarrange(a,c,e, nrow = 3, ncol = 1, align = 'hv')
part2 <- ggarrange(b,d,f, nrow = 3, ncol = 1, align = 'hv')
Fig4A.B <- ggarrange(part1, part2, ncol = 2, align = 'hv', labels = c("A", "B"))

ggsave("Figure4_new_panel_A-b.svg", Fig4A.B, width = 10.239, height = 14.091)

####### NEW Panel C figure 4 #######
library(fantaxtic)
# Load data 
phy_vst <- readRDS("../Data_files/phyloseq_16s_no_soil_filtered_vst_dataset.rds")
# Remove two outlier root samples, IDed above 
phy_vst <- subset_samples(phy_vst, extraction_num != "492" & extraction_num != "615")
# make subset root top ten
phy_root_only_vst <- subset_samples(phy_vst, plant_body_site=="root")
# Make Nested dataframe for plotting
top_nested <- nested_top_taxa(phy_root_only_vst,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Class",
                              n_top_taxa = 10, 
                              n_nested_taxa = 3,
                              by_proportion = TRUE) 
# Facet labels capitalized for consistency
ylabels <- c('madera'='Madera', 'merced' = 'Merced', 'san joaquin' = "San Joaquin")
# Nest plot with facet between sites
fig4C <- plot_nested_bar(top_nested$ps_obj,
                                  top_level = "Phylum",
                                  nested_level = "Class",
                                  legend_title = "Phylum and Class",
                                  palette = safe_colorblind_palette) +
  labs(y = "Relative Abuance") +
  facet_wrap(~lower_political, scales = "free_x", labeller = as_labeller(ylabels)) +
  # mostly changes made to fit the pubr theme without breaking legend text formatting from fantaxtic
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), legend.key.size = unit(10, "points"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right',
        strip.background = element_rect(fill = "#F2F2F2", colour = "black", size = 0.7),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(color = "black"),
        strip.text = element_text(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line("black"), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  guides(fill=guide_legend(ncol=2))
# Save to exact dimensions 
ggsave("Figure4_new_panel_C.svg", fig4C, width = 17.664, height = 7.083)

###### Figure S4 ######
#Figure S4A/B
## Berry plots
ylabels <- c('madera'='Madera', 'merced' = 'Merced', 'san joaquin' = "San Joaquin")
Figure4A <- ggplot(phylum_lvl[phylum_lvl$plant_body_site == "berry",], aes(x= Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~lower_political, scale = "free", ncol = 3, labeller = as_labeller(ylabels)) +
  scale_fill_manual(values=safe_colorblind_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')
## Leaf plots
Figure4B <- ggplot(phylum_lvl[phylum_lvl$plant_body_site == "leaf",], aes(x= Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~lower_political, scale = "free", ncol = 3, labeller = as_labeller(ylabels)) +
  scale_fill_manual(values=safe_colorblind_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# Save it
Figure4 <- ggarrange(Figure4A, Figure4B, nrow = 2, common.legend = TRUE,legend = 'right', labels = c('A','B'))
ggsave("FigureS4_Taxonomic_barplot_berry&leaf_by_site_LM_heatmap.png", Figure4, height = 12, width = 18.7, units = "in")
ggsave("FigureS4_Taxonomic_barplot_berry&leaf_by_site_LM_heatmap.pdf", Figure4, height = 12, width = 18.7, units = "in")

### NEW FIGURE 5 Berry and leaf taxabarplots ###
# subset to berries and leaves
phy_berry_leaf_vst <- subset_samples(phy_vst, plant_body_site %in% c("berry", "leaf"))
top_nested <- nested_top_taxa(phy_berry_leaf_vst,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Class",
                              n_top_taxa = 10, 
                              n_nested_taxa = 3,
                              by_proportion = TRUE) 
ylabels <- c('madera'='Madera', 'merced' = 'Merced', 'san joaquin' = "San Joaquin", "berry" = "Berry", "leaf" = "Leaf")
Fig_5 <- plot_nested_bar(top_nested$ps_obj,
                top_level = "Phylum",
                nested_level = "Class",
                legend_title = "Phylum and Class",
                palette = safe_colorblind_palette) +
  labs(y = "Relative Abuance") +
  facet_wrap(~plant_body_site+lower_political, scales = "free_x", labeller = as_labeller(ylabels)) +
  # mostly changes made to fit the pubr theme without breaking legend text formatting from fantaxtic
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), legend.key.size = unit(10, "points"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right',
        strip.background = element_rect(fill = "#F2F2F2", colour = "black", size = 0.7),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(color = "black"),
        strip.text = element_text(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line("black"), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  guides(fill=guide_legend(ncol=2))

ggsave("Figure5_new_figure.svg", Fig_5, width = 20, height = 10)


### Berry top ten phyla lms
### This is a heat map with each of the top 10 phyla + the other category
### The heat map will show the % variance explained after correcting for multiple
### testing across all of the models together (BH, Benjamini & Hochberg).
# Get data from phyloseq object.
phylum_lvl_panelD <- tax_glom(phy_vst, taxrank = "Phylum") # 42 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
phylum_lvl_panelD <- transform_sample_counts(phylum_lvl_panelD, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
phylum_lvl_panelD <- fantaxtic::get_top_taxa(physeq_obj = phylum_lvl_panelD, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Instead of using psmelt, instead keep data in wide format (taxa as rows)
Phylum_relative_abund <- as.data.frame(t(otu_table(phylum_lvl_panelD)))
colnames(Phylum_relative_abund) <- as.data.frame(phylum_lvl_panelD@tax_table)$Phylum # Name the columns their by phyla
Phylum_w_meta <- cbind(Phylum_relative_abund, pssd2veg(phylum_lvl_panelD))# Attach metadata
# Linear models but only for the berry samples
model_ot <- lm(Other ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_my <- lm(Myxococcota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_pr <- lm(Proteobacteria ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_ch <- lm(Chloroflexi ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_de <- lm(Deinococcota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_fi <- lm(Firmicutes ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_act <- lm(Actinobacteriota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_ba <- lm(Bacteroidota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_pl <- lm(Planctomycetota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_ve <- lm(Verrucomicrobiota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
model_aci <- lm(Acidobacteriota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "berry",])
# Posthocs
pairs(emmeans(model_ba, ~ year))
# Make list of ion names in order and run the function to generate the % variance
taxa_list <- as.data.frame((phylum_lvl_panelD@tax_table))$Phylum
Percent_var <- bind_rows(Return_SS_proportion(model_ot),
                         Return_SS_proportion(model_my),
                         Return_SS_proportion(model_pr),
                         Return_SS_proportion(model_ch),
                         Return_SS_proportion(model_de),
                         Return_SS_proportion(model_fi),
                         Return_SS_proportion(model_act),
                         Return_SS_proportion(model_ba),
                         Return_SS_proportion(model_pl),
                         Return_SS_proportion(model_ve),
                         Return_SS_proportion(model_aci))
rownames(Percent_var) <- taxa_list
# remove residual 
Percent_var_noint_nores <- Percent_var[, !names(Percent_var) %in% "Residuals"]
# Make dataframe of P values to connect to the % var df above
Pval <- bind_rows(Get_Pval(model_ot),
                  Get_Pval(model_my),
                  Get_Pval(model_pr),
                  Get_Pval(model_ch),
                  Get_Pval(model_de),
                  Get_Pval(model_fi),
                  Get_Pval(model_act),
                  Get_Pval(model_ba),
                  Get_Pval(model_pl),
                  Get_Pval(model_ve),
                  Get_Pval(model_aci))
rownames(Pval) <- taxa_list
# Combine
SS_pval_combo_df <- cbind(Percent_var_noint_nores,Pval)
SS_pval_combo_df$taxa <- rownames(SS_pval_combo_df)
# Get variance columns
total_var <- SS_pval_combo_df %>% dplyr::select(taxa, "rootstock":"site")
# Reorganize and rename to format for plotting
total_var <- total_var %>% gather(key=factor, value=var, -taxa)
total_var <- total_var %>% mutate(factor=str_replace(factor, "_var" %R% END, ""))
# Get p-value columns
total_p_name <- SS_pval_combo_df %>% dplyr::select("rootstock_p":"site_p")
total_p <- data.frame(t(apply(total_p_name, 1, FUN=p.adjust, method='BH')))
colnames(total_p)<- colnames(total_p_name) #This is hacky but I need to preserve the colnames to correctly join them later on.
total_p$taxa <- SS_pval_combo_df$taxa
# Reorganize and rename to format for plotting
total_p <- total_p %>% gather(key=factor, value=p_value, -taxa)
total_p <- total_p %>% mutate(factor=str_replace(factor,"_p" %R% END, ""))
# Join variance and p-value tables back together 
total_var_p <- full_join(total_var, total_p,by=c("taxa", "factor"))
total_var_p$p_value < 0.05


### Leaf top ten phyla lms
### This is a heat map with each of the top 10 phyla + the other category
### The heat map will show the % variance explained after correcting for multiple
### testing across all of the models together (BH, Benjamini & Hochberg).
# Get data from phyloseq object.
phylum_lvl_panelD <- tax_glom(phy_vst, taxrank = "Phylum") # 42 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
phylum_lvl_panelD <- transform_sample_counts(phylum_lvl_panelD, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
phylum_lvl_panelD <- fantaxtic::get_top_taxa(physeq_obj = phylum_lvl_panelD, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Instead of using psmelt, instead keep data in wide format (taxa as rows)
Phylum_relative_abund <- as.data.frame(t(otu_table(phylum_lvl_panelD)))
colnames(Phylum_relative_abund) <- as.data.frame(phylum_lvl_panelD@tax_table)$Phylum # Name the columns their by phyla
Phylum_w_meta <- cbind(Phylum_relative_abund, pssd2veg(phylum_lvl_panelD))# Attach metadata
# Linear models but only for the berry samples
model_ot <- lm(Other ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_my <- lm(Myxococcota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_pr <- lm(Proteobacteria ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_ch <- lm(Chloroflexi ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_de <- lm(Deinococcota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_fi <- lm(Firmicutes ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_act <- lm(Actinobacteriota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_ba <- lm(Bacteroidota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_pl <- lm(Planctomycetota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_ve <- lm(Verrucomicrobiota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
model_aci <- lm(Acidobacteriota ~ rootstock + scion + year + brix + site, data=Phylum_w_meta[Phylum_w_meta$plant_body_site == "leaf",])
#Post-Hocs
anova(model_act)
anova(model_fi)
pairs(emmeans(model_act, ~ year))
pairs(emmeans(model_fi, ~ site))
# Make list of ion names in order and run the function to generate the % variance
taxa_list <- as.data.frame((phylum_lvl_panelD@tax_table))$Phylum
Percent_var <- bind_rows(Return_SS_proportion(model_ot),
                         Return_SS_proportion(model_my),
                         Return_SS_proportion(model_pr),
                         Return_SS_proportion(model_ch),
                         Return_SS_proportion(model_de),
                         Return_SS_proportion(model_fi),
                         Return_SS_proportion(model_act),
                         Return_SS_proportion(model_ba),
                         Return_SS_proportion(model_pl),
                         Return_SS_proportion(model_ve),
                         Return_SS_proportion(model_aci))
rownames(Percent_var) <- taxa_list
# remove residual 
Percent_var_noint_nores <- Percent_var[, !names(Percent_var) %in% "Residuals"]
# Make dataframe of P values to connect to the % var df above
Pval <- bind_rows(Get_Pval(model_ot),
                  Get_Pval(model_my),
                  Get_Pval(model_pr),
                  Get_Pval(model_ch),
                  Get_Pval(model_de),
                  Get_Pval(model_fi),
                  Get_Pval(model_act),
                  Get_Pval(model_ba),
                  Get_Pval(model_pl),
                  Get_Pval(model_ve),
                  Get_Pval(model_aci))
rownames(Pval) <- taxa_list
# Combine
SS_pval_combo_df <- cbind(Percent_var_noint_nores,Pval)
SS_pval_combo_df$taxa <- rownames(SS_pval_combo_df)
# Get variance columns
total_var <- SS_pval_combo_df %>% dplyr::select(taxa, "rootstock":"site")
# Reorganize and rename to format for plotting
total_var <- total_var %>% gather(key=factor, value=var, -taxa)
total_var <- total_var %>% mutate(factor=str_replace(factor, "_var" %R% END, ""))
# Get p-value columns
total_p_name <- SS_pval_combo_df %>% dplyr::select("rootstock_p":"site_p")
total_p <- data.frame(t(apply(total_p_name, 1, FUN=p.adjust, method='BH')))
colnames(total_p)<- colnames(total_p_name) #This is hacky but I need to preserve the colnames to correctly join them later on.
total_p$taxa <- SS_pval_combo_df$taxa
# Reorganize and rename to format for plotting
total_p <- total_p %>% gather(key=factor, value=p_value, -taxa)
total_p <- total_p %>% mutate(factor=str_replace(factor,"_p" %R% END, ""))
# Join variance and p-value tables back together 
total_var_p <- full_join(total_var, total_p,by=c("taxa", "factor"))
total_var_p$p_value < 0.05


####### Looking at the top ten phyla for each compartment individually. #######
Top_10_by_compartment_phyla <- function(phyloseq_obj, split_by){
  if(split_by == "leaf"){
  X <- subset_samples(physeq = phy_vst, plant_body_site=="leaf")
  } else if(split_by == "berry"){
  X <- subset_samples(physeq = phy_vst, plant_body_site=="berry")
  } else if(split_by == "root"){
  X <- subset_samples(physeq = phy_vst, plant_body_site=="root")
  }
  # Get data from phyloseq object.
  X <- tax_glom(phy_vst, taxrank = "Phylum") # 42 taxa
  # Make relative abundance AKA transform sample counts by diving by the ASV total.
  Y <- transform_sample_counts(X, function(x) x/sum(x))
  # Take only the n number of taxa per taxonomic rank based on relative abundance.
  # Taxa >n will be added to a other label.
  Y <- fantaxtic::get_top_taxa(physeq_obj = Y, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
  # Instead of using psmelt, instead keep data in wide format (taxa as rows)
  Phylum_relative_abund <- as.data.frame(t(otu_table(Y)))
  colnames(Phylum_relative_abund) <- as.data.frame(Y@tax_table)$Phylum # Name the columns their by phyla
  Phylum_w_meta <- cbind(Phylum_relative_abund, pssd2veg(Y))# Attach metadata
  return(Phylum_w_meta)
}
summary(Top_10_by_compartment_phyla(phy_vst, "leaf"))[,1:11]
summary(Top_10_by_compartment_phyla(phy_vst, "berry"))[,1:11]
summary(Top_10_by_compartment_phyla(phy_vst, "root"))[,1:11]
# Conclusion, at the phyla level we do not see the 10 ten change if we select only
# within a single compartment at a time.