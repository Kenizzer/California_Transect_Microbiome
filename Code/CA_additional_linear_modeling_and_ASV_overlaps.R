# Analysis of grafted grapevines across time and space.
# Samples were collected in summer 2018 and 2019
# Samples consist of leaf, berry, and root compartments
# Code by: Joel F. Swift

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('ggpubr'); packageVersion('ggpubr')
library('car'); packageVersion('car')
library('lme4'); packageVersion('lme4')
library('lmerTest'); packageVersion('lmerTest')
library('vegan'); packageVersion('vegan')
library('phyloseq'); packageVersion('phyloseq')
library('emmeans'); packageVersion('emmeans')
library('lubridate'); packageVersion('lubridate')
library('rebus'); packageVersion('rebus')
library('viridis'); packageVersion('viridis')
library('corrplot'); packageVersion('corrplot')
library('factoextra'); packageVersion('factoextra')

## Functions
# Adjust P values across taxonomic level but separately per term
# Get corrected P values from here
Padj.term.func <- function(factor, modelDF){
  # Check if any are significant, save those their indexes to temp
  temp <- which(p.adjust(modelDF[modelDF$term == factor, ]$p.value, method = "BH") < 0.05)
  temp2 <- filter(modelDF, term == factor)
  # adjust p values by term and add this to a term specific tibble save to temp2
  temp2$p.adj <- p.adjust(modelDF[modelDF$term == factor, ]$p.value, method = "BH")
  temp3 <- temp2[c(temp), ] # index temp2 with temp
  temp4 <- temp3 %>% mutate_if(is.numeric, signif, digits = 4)
  return(as.data.frame(temp4))
  # Returns an empty tibble if none are significant after p.adjust
}


#### Linear modeling ####
# load dataset generated from Data_filtering_normalization.R
phy_vst <- readRDS("../Data_files/phyloseq_16s_no_soil_filtered_vst_dataset.rds")

# split by compartment
phy_vst_root <- subset_samples(phy_vst, plant_body_site == 'root')
phy_vst_leaf <- subset_samples(phy_vst, plant_body_site == 'leaf')
phy_vst_berr <- subset_samples(phy_vst, plant_body_site == 'berry')

##### Root modeling #####
phy_vst_phy <- psmelt(tax_glom(phy_vst_root, "Phylum"))
phy_vst_cla <- psmelt(tax_glom(phy_vst_root, "Class"))
phy_vst_ord <- psmelt(tax_glom(phy_vst_root, "Order"))
phy_vst_fam <- psmelt(tax_glom(phy_vst_root, "Family"))
phy_vst_gen <- psmelt(tax_glom(phy_vst_root, "Genus"))

# Phylum
phy_vst_aovs_phy <- phy_vst_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Class
phy_vst_aovs_cla <- phy_vst_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Order
phy_vst_aovs_ord <- phy_vst_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Family
phy_vst_aovs_fam <- phy_vst_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Genus
phy_vst_aovs_gen <- phy_vst_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))

# Rootstock
a <- Padj.term.func("rootstock", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("rootstock", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("rootstock", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("rootstock", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("rootstock", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Rootstock_all_levels <- rbind(a,b,c,d,e)

# Scion
a <- Padj.term.func("scion", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("scion", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("scion", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("scion", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("scion", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Scion_all_levels <- rbind(a,b,c,d,e)



# Year
a <- Padj.term.func("year", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("year", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("year", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("year", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("year", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Year_all_levels <- rbind(a,b,c,d,e)

# Brix
a <- Padj.term.func("brix", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("brix", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("brix", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("brix", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("brix", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Brix_all_levels <- rbind(a,b,c,d,e)
Brix_all_levels$term <- "Sugar Content"
# Site
a <- Padj.term.func("site", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("site", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("site", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("site", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("site", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Site_all_levels <- rbind(a,b,c,d,e)


All_terms_all_levels <- rbind(Rootstock_all_levels,
                             Scion_all_levels,
                             Year_all_levels,
                             Brix_all_levels,
                             Site_all_levels)

write.csv(All_terms_all_levels, "Supplemental_table_Root_Linear_models_all_terms_and_levels.csv", row.names=FALSE)


rm(a,b,c,d,e, Rootstock_all_levels, Scion_all_levels, Year_all_levels, Brix_all_levels, Site_all_levels, All_terms_all_levels,
   phy_vst_aovs_phy, phy_vst_aovs_ord, phy_vst_aovs_gen, phy_vst_aovs_fam, phy_vst_aovs_cla,
   phy_vst_phy, phy_vst_ord, phy_vst_gen, phy_vst_fam, phy_vst_cla)


##### leaf modeling #####
phy_vst_phy <- psmelt(tax_glom(phy_vst_leaf, "Phylum"))
phy_vst_cla <- psmelt(tax_glom(phy_vst_leaf, "Class"))
phy_vst_ord <- psmelt(tax_glom(phy_vst_leaf, "Order"))
phy_vst_fam <- psmelt(tax_glom(phy_vst_leaf, "Family"))
phy_vst_gen <- psmelt(tax_glom(phy_vst_leaf, "Genus"))

# Phylum
phy_vst_aovs_phy <- phy_vst_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Class
phy_vst_aovs_cla <- phy_vst_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Order
phy_vst_aovs_ord <- phy_vst_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Family
phy_vst_aovs_fam <- phy_vst_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Genus
phy_vst_aovs_gen <- phy_vst_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))

# Rootstock
a <- Padj.term.func("rootstock", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("rootstock", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("rootstock", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("rootstock", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("rootstock", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

#Rootstock_all_levels <- rbind(a,b,c,d,e)
## None significant

# Scion
a <- Padj.term.func("scion", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("scion", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("scion", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("scion", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("scion", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Scion_all_levels <- rbind(a,b,c,d,e)

# Year
a <- Padj.term.func("year", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("year", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("year", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("year", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("year", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Year_all_levels <- rbind(a,b,c,d,e)

# Brix
a <- Padj.term.func("brix", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("brix", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("brix", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("brix", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("brix", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

# Brix_all_levels <- rbind(a,b,c,d,e)
# Brix_all_levels$term <- "Sugar Content"
## None significant

# Site
a <- Padj.term.func("site", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("site", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("site", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("site", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("site", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Site_all_levels <- rbind(a,b,c,d,e)


All_terms_all_levels <- rbind(Scion_all_levels,
                              Year_all_levels,
                              Site_all_levels)

write.csv(All_terms_all_levels, "Supplemental_table_Leaf_Linear_models_all_terms_and_levels.csv", row.names=FALSE)


rm(a,b,c,d,e, Rootstock_all_levels, Scion_all_levels, Year_all_levels, Brix_all_levels, Site_all_levels, All_terms_all_levels,
   phy_vst_aovs_phy, phy_vst_aovs_ord, phy_vst_aovs_gen, phy_vst_aovs_fam, phy_vst_aovs_cla,
   phy_vst_phy, phy_vst_ord, phy_vst_gen, phy_vst_fam, phy_vst_cla)

##### Berry modeling #####
phy_vst_phy <- psmelt(tax_glom(phy_vst_berr, "Phylum"))
phy_vst_cla <- psmelt(tax_glom(phy_vst_berr, "Class"))
phy_vst_ord <- psmelt(tax_glom(phy_vst_berr, "Order"))
phy_vst_fam <- psmelt(tax_glom(phy_vst_berr, "Family"))
phy_vst_gen <- psmelt(tax_glom(phy_vst_berr, "Genus"))

# Phylum
phy_vst_aovs_phy <- phy_vst_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Class
phy_vst_aovs_cla <- phy_vst_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Order
phy_vst_aovs_ord <- phy_vst_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Family
phy_vst_aovs_fam <- phy_vst_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))
# Genus
phy_vst_aovs_gen <- phy_vst_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(Abundance ~ rootstock + scion + year + brix + site, data=data)))) %>%
  summarize(broom::tidy(mod))

# Rootstock
a <- Padj.term.func("rootstock", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("rootstock", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("rootstock", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("rootstock", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("rootstock", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Rootstock_all_levels <- rbind(a,b,c,d,e)

# Scion
a <- Padj.term.func("scion", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("scion", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("scion", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("scion", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("scion", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

#Scion_all_levels <- rbind(a,b,c,d,e)
## None significant

# Year
a <- Padj.term.func("year", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("year", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("year", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("year", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("year", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Year_all_levels <- rbind(a,b,c,d,e)

# Brix
a <- Padj.term.func("brix", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("brix", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("brix", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("brix", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("brix", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Brix_all_levels <- rbind(a,b,c,d,e)
Brix_all_levels$term <- "Sugar Content"

# Site
a <- Padj.term.func("site", phy_vst_aovs_phy)
names(a)[names(a) == 'Phylum'] <- 'Taxon'
a <- cbind(Level = "Phylum", a)

b <- Padj.term.func("site", phy_vst_aovs_cla)
names(b)[names(b) == 'Class'] <- 'Taxon'
b <- cbind(Level = "Class", b)

c <- Padj.term.func("site", phy_vst_aovs_ord)
names(c)[names(c) == 'Order'] <- 'Taxon'
c <- cbind(Level = "Order", c)

d <- Padj.term.func("site", phy_vst_aovs_fam)
names(d)[names(d) == 'Family'] <- 'Taxon'
d <- cbind(Level = "Family", d)

e <- Padj.term.func("site", phy_vst_aovs_gen)
names(e)[names(e) == 'Genus'] <- 'Taxon'
e <- cbind(Level = "Genus", e)

Site_all_levels <- rbind(a,b,c,d,e)


All_terms_all_levels <- rbind(Rootstock_all_levels,
                              Year_all_levels,
                              Brix_all_levels,
                              Site_all_levels)

write.csv(All_terms_all_levels, "Supplemental_table_Berry_Linear_models_all_terms_and_levels.csv", row.names=FALSE)


rm(a,b,c,d,e, Rootstock_all_levels, Scion_all_levels, Year_all_levels, Brix_all_levels, Site_all_levels, All_terms_all_levels,
   phy_vst_aovs_phy, phy_vst_aovs_ord, phy_vst_aovs_gen, phy_vst_aovs_fam, phy_vst_aovs_cla,
   phy_vst_phy, phy_vst_ord, phy_vst_gen, phy_vst_fam, phy_vst_cla)


#### ASV overlap between compartments ####
# load venn diagram package
library(ggvenn)
# load dataset with soil samples 
phy_with_soil <- readRDS("../Data_files/phyloseq_16s_filtered_dataset.rds")
# Soil ASVs
soil_samples <- subset_samples(phy_with_soil, plant_body_site == "soil")
soil_samples <- prune_taxa(taxa_sums(soil_samples) > 0, soil_samples)
soil_asvs <- row.names(soil_samples@tax_table)
# Root ASVs
root_samples <- subset_samples(phy_with_soil, plant_body_site == "root")
root_samples <- prune_taxa(taxa_sums(root_samples) > 0, root_samples)
root_asvs <- row.names(root_samples@tax_table)
# Leaf ASVs
leaf_samples <- subset_samples(phy_with_soil, plant_body_site == "leaf")
leaf_samples <- prune_taxa(taxa_sums(leaf_samples) > 0, leaf_samples)
leaf_asvs <- row.names(leaf_samples@tax_table)
# Berry ASVs
berry_samples <- subset_samples(phy_with_soil, plant_body_site == "berry")
berry_samples <- prune_taxa(taxa_sums(berry_samples) > 0, berry_samples)
berry_asvs <- row.names(berry_samples@tax_table)
# add lists together and give capitalized names for plot labels
with_soil_list <- c(list(Berry = berry_asvs), list(Leaf = leaf_asvs), list(Root = root_asvs), list(Soil = soil_asvs))
Plot_with_soil <- ggvenn(with_soil_list, fill_alpha = 0.5, fill_color = c("#5a1991", "#139d08", "#5c3c0d", "black"),
       text_color = 'white', show_percentage = TRUE, digits = 0, stroke_color = "black", text_size = 5)

# Change percentages less than zero to be <1% (otherwise they are displayed as 0%)
Plot_with_soil[["layers"]][[4]][["data"]][["text"]]
# [1] "1\n(0%)"     "14\n(0%)"    "3105\n(35%)" "18\n(0%)"    "140\n(2%)"   "616\n(7%)"   "2506\n(28%)" "353\n(4%)"   "9\n(0%)"    
# [10] "4\n(0%)"     "412\n(5%)"   "743\n(8%)"   "284\n(3%)"   "8\n(0%)"     "625\n(7%)"

Plot_with_soil[["layers"]][[4]][["data"]][["text"]] <- c("1\n(<1%)", "14\n(<1%)", "3105\n(35%)", "18\n(<1%)", "140\n(2%)", "616\n(7%)",
                                                         "2506\n(28%)", "353\n(4%)", "9\n(<1%)", "4\n(<1%)", "412\n(5%)", "743\n(8%)",
                                                         "284\n(3%)", "8\n(<1%)", "625\n(7%)")

rm(soil_samples, soil_asvs, root_samples, root_asvs, leaf_samples, leaf_asvs, berry_samples, berry_asvs)

## Since the data set used for all the main analyses conducted filtering of ASVs without
## soil samples included in the data set, the number of ASVs is slightly lower as less
## ASVs did met the minimum of being present in 5 or more samples. I will make another
## venn diagram to represent this dataset (just berry, leaf, and root).
phy_without_soil <- readRDS("../Data_files/phyloseq_16s_no_soil_filtered_dataset.rds")
# Root ASVs
root_samples <- subset_samples(phy_without_soil, plant_body_site == "root")
root_samples <- prune_taxa(taxa_sums(root_samples) > 0, root_samples)
root_asvs <- row.names(root_samples@tax_table)
# Leaf ASVs
leaf_samples <- subset_samples(phy_without_soil, plant_body_site == "leaf")
leaf_samples <- prune_taxa(taxa_sums(leaf_samples) > 0, leaf_samples)
leaf_asvs <- row.names(leaf_samples@tax_table)
# Berry ASVs
berry_samples <- subset_samples(phy_without_soil, plant_body_site == "berry")
berry_samples <- prune_taxa(taxa_sums(berry_samples) > 0, berry_samples)
berry_asvs <- row.names(berry_samples@tax_table)
# add lists together and give capitalized names for plot labels
without_soil_list <- c(list(Berry = berry_asvs), list(Leaf = leaf_asvs), list(Root = root_asvs))
Plot_without_soil <- ggvenn(without_soil_list, fill_alpha = 0.5, fill_color = c("#5a1991", "#139d08", "#5c3c0d"),
                         text_color = 'white', show_percentage = TRUE, digits = 0, stroke_color = "black", text_size = 5)
# Change percentages less than zero to be <1% (otherwise they are displayed as 0%)
Plot_without_soil[["layers"]][[4]][["data"]][["text"]]
# [1] "1\n(0%)"     "16\n(0%)"    "4909\n(62%)" "142\n(2%)"   "602\n(8%)"   "1282\n(16%)" "1029\n(13%)"
Plot_without_soil[["layers"]][[4]][["data"]][["text"]] <- c("1\n(<1%)", "16\n(<1%)", "4909\n(62%)", "142\n(2%)",
                                                            "602\n(8%)", "1282\n(16%)", "1029\n(13%)")
# Arrange and save
Figure_S5 <- ggarrange(Plot_without_soil, Plot_with_soil, labels = "AUTO")

ggsave("Figure_S5_venn_diagrams.svg", Figure_S5, height = 6, width = 10)
ggsave("Figure_S5_venn_diagrams.png", Figure_S5, height = 6, width = 10)
