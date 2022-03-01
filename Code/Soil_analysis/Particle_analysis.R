# Processing of soil hydrometer readings
# Code by: Joel F. Swift
library(envalysis)
library(ggplot2)
library(data.table)
library(ggpubr)
library(emmeans)

#CA2
CA2 <- read.table(file= "CA02.txt", sep = "\t", header = TRUE)
CA2_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA2, plot = T, hydrometer = "152H")
CA2_tex
#CA3
CA3 <- read.table(file= "CA03.txt", sep = "\t", header = TRUE)
CA3_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA3, plot = T, hydrometer = "152H")
CA3_tex
#CA4
CA4 <- read.table(file= "CA04.txt", sep = "\t", header = TRUE)
CA4_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA4, plot = T, hydrometer = "152H")
CA4_tex
#CA5
CA5 <- read.table(file= "CA05.txt", sep = "\t", header = TRUE)
CA5_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA5, plot = T, hydrometer = "152H")
CA5_tex
#CA6
CA6 <- read.table(file= "CA06.txt", sep = "\t", header = TRUE)
CA6_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA6, plot = T, hydrometer = "152H")
CA6_tex
#CA9
CA9 <- read.table(file= "CA09.txt", sep = "\t", header = TRUE)
CA9_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA9, plot = T, hydrometer = "152H")
CA9_tex
#CA10
CA10 <- read.table(file= "CA10.txt", sep = "\t", header = TRUE)
CA10_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA10, plot = T, hydrometer = "152H")
CA10_tex
#CA11
CA11 <- read.table(file= "CA11.txt", sep = "\t", header = TRUE)
CA11_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA11, plot = T, hydrometer = "152H")
CA11_tex
#CA12
CA12 <- read.table(file= "CA12.txt", sep = "\t", header = TRUE)
CA12_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA12, plot = T, hydrometer = "152H")
CA12_tex
#CA13
CA13 <- read.table(file= "CA13.txt", sep = "\t", header = TRUE)
CA13_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA13, plot = T, hydrometer = "152H")
CA13_tex
#CA14
CA14 <- read.table(file= "CA14.txt", sep = "\t", header = TRUE)
CA14_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA14, plot = T, hydrometer = "152H")
CA14_tex
#CA15
CA15 <- read.table(file= "CA15.txt", sep = "\t", header = TRUE)
CA15_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA15, plot = T, hydrometer = "152H")
CA15_tex
#CA16
CA16 <- read.table(file= "CA16.txt", sep = "\t", header = TRUE)
CA16_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA16, plot = T, hydrometer = "152H")
CA16_tex
#CA17
CA17 <- read.table(file= "CA17.txt", sep = "\t", header = TRUE)
CA17_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA17, plot = T, hydrometer = "152H")
CA17_tex
#CA18
CA18 <- read.table(file= "CA18.txt", sep = "\t", header = TRUE)
CA18_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA18, plot = T, hydrometer = "152H")
CA18_tex
#CA19
CA19 <- read.table(file= "CA19.txt", sep = "\t", header = TRUE)
CA19_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA19, plot = T, hydrometer = "152H")
CA19_tex
#CA20
CA20 <- read.table(file= "CA20.txt", sep = "\t", header = TRUE)
CA20_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA20, plot = T, hydrometer = "152H")
CA20_tex
#CA21
CA21 <- read.table(file= "CA21.txt", sep = "\t", header = TRUE)
CA21_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA21, plot = T, hydrometer = "152H")
CA21_tex
#CA22
CA22 <- read.table(file= "CA22.txt", sep = "\t", header = TRUE)
CA22_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA22, plot = T, hydrometer = "152H")
CA22_tex
#CA23
CA23 <- read.table(file= "CA23.txt", sep = "\t", header = TRUE)
CA23_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA23, plot = T, hydrometer = "152H")
CA23_tex
#CA24
CA24 <- read.table(file= "CA24.txt", sep = "\t", header = TRUE)
CA24_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA24, plot = T, hydrometer = "152H")
CA24_tex
#CA25
CA25 <- read.table(file= "CA25.txt", sep = "\t", header = TRUE)
CA25_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA25, plot = T, hydrometer = "152H")
CA25_tex
#CA26
CA26 <- read.table(file= "CA26.txt", sep = "\t", header = TRUE)
CA26_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA26, plot = T, hydrometer = "152H")
CA26_tex
#CA27
CA27 <- read.table(file= "CA27.txt", sep = "\t", header = TRUE)
CA27_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA27, plot = T, hydrometer = "152H")
CA27_tex
#CA28
CA28 <- read.table(file= "CA28.txt", sep = "\t", header = TRUE)
CA28_tex <- texture(Reading ~ Correction_factor + Time + Temperature, CA28, plot = T, hydrometer = "152H")
CA28_tex

#Combine dataframe of % Clay, Silt, Sand
Full_df <- as.data.frame(rbind(CA02 = CA2_tex$usda[1,], 
                               CA03 = CA3_tex$usda[1,],
                               CA04 = CA4_tex$usda[1,],
                               CA05 = CA5_tex$usda[1,],
                               CA06 = CA6_tex$usda[1,],
                               CA09 = CA9_tex$usda[1,],
                               CA10 = CA10_tex$usda[1,],
                               CA11 = CA11_tex$usda[1,],
                               CA12 = CA12_tex$usda[1,],
                               CA13 = CA13_tex$usda[1,],
                               CA14 = CA14_tex$usda[1,],
                               CA15 = CA15_tex$usda[1,],
                               CA16 = CA16_tex$usda[1,],
                               CA17 = CA17_tex$usda[1,],
                               CA18 = CA18_tex$usda[1,],
                               CA19 = CA19_tex$usda[1,],
                               CA20 = CA20_tex$usda[1,],
                               CA21 = CA21_tex$usda[1,],
                               CA22 = CA22_tex$usda[1,],
                               CA23 = CA23_tex$usda[1,],
                               CA24 = CA24_tex$usda[1,],
                               CA25 = CA25_tex$usda[1,],
                               CA26 = CA26_tex$usda[1,],
                               CA27 = CA27_tex$usda[1,],
                               CA28 = CA28_tex$usda[1,]))


Full_df$Site <-c("Ripperdan", "Ripperdan", "Ripperdan", "Livingston", "Livingston", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Liberty", "Livingston", "Livingston", "Livingston", "Livingston", "Ripperdan", "Ripperdan", "Ripperdan", "Ripperdan")
#Save object for other Rscript Ternary_plot_ONLY.R
saveRDS(Full_df, file = 'dataframe_for_analysis.RDS')
