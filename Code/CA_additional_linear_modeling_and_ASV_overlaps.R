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
# load dataset with soil samples included
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

# For Venn Diagram
write.csv(soil_asvs, "Soil_ASVs.csv", row.names = FALSE)
write.csv(root_asvs, "Root_ASVs.csv", row.names = FALSE)
write.csv(leaf_asvs, "Leaf_ASVs.csv", row.names = FALSE)
write.csv(berry_asvs, "Berry_ASVs.csv", row.names = FALSE)

# upset plot
library(UpSetR)

compartment_sets <- list(
  Berry = berry_asvs,
  Leaf = leaf_asvs,
  Soil = soil_asvs,
  Root = root_asvs)

upset_plot <- upset(fromList(compartment_sets),
      sets.x.label = "ASVs")

pdf(file="ASV_overlap_compartments.pdf") # or other device
upset_plot
dev.off()