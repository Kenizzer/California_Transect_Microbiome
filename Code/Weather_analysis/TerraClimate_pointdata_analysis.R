# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(ggpubr); packageVersion('ggpubr')
library(lubridate); packageVersion('lubridate')

# Theme set and Color Palettes
theme_set(theme_pubr())
site_palette <- c('#e6ab02', '#281c39', '#12664c')

# Take terraclimate point data and get mean annual precipitation, Evapotranspiration and PDSI
# then calculate an aridity index.
Get_aridity <- function(terraclimdf){
  temp  <- terraclimdf %>% group_by(Year) %>% summarize(Annual_precip_mm = sum(ppt.mm.), Annual_ETo_mm = sum(pet.mm.), Annual_PDSI = mean(PDSI.unitless.), Annual_soilmoisture_mm = sum(soil.mm.))
  # get the last normal AKA 1990 to 2021 
  temp2 <- temp %>% filter(between(Year, 1990, 2021)) %>% summarize(MA_Prec = mean(Annual_precip_mm), MA_ETo = mean(Annual_ETo_mm), M_monthly_PDSI = mean(Annual_PDSI), MA_Soil_moisture = mean(Annual_soilmoisture_mm))
  temp2$Aridity_index <- temp2$MA_Prec / temp2$MA_ETo
  return(temp2)
}

#Aridity Index Value chart https://doi.org/10.1038/s41597-022-01493-1
# Climate Class 
# <0.03 Hyper Arid
# 0.03–0.2 Arid
# 0.2–0.5 Semi-Arid
# 0.5–0.65 Dry sub-humid
# >0.65 Humid

##### Coords #####
#Smoky Valley Ranch prairie [latitude 38.8665, longitude -100.9951]
#Smoky Valley Ranch grain field [38.8791, -100.9828]
#Hayes prairie [38.8355, -99.3033]
#The Land Institute grain field [38.6943, -97.5912]
#The Land Institute prairie [38.9698, -97.4690]
#Konza native prairie [39.1056, -96.6099])

# Load point data
Madera_2018 <- read.csv("Madera_2018.csv", skip = 14, header = TRUE)
Madera_2019 <- read.csv("Madera_2019.csv", skip = 14, header = TRUE)
Madera_data <- rbind(Madera_2018, Madera_2019)
Madera_data$Year <- as.factor(Madera_data$Year)
Madera_data$Month <- as.factor(Madera_data$Month)
Madera_data$Day <- as.factor(Madera_data$Day)
Madera_data$Site <- as.factor("Madera")
Madera_data$Date <- paste(Madera_data$Year, Madera_data$Month, Madera_data$Day, sep = "_")
summary(Madera_data)


Merced_2018 <- read.csv("Merced_2018.csv", skip = 14, header = TRUE)
Merced_2019 <- read.csv("Merced_2019.csv", skip = 14, header = TRUE)
Merced_data <- rbind(Merced_2018, Merced_2019)
Merced_data$Year <- as.factor(Merced_data$Year)
Merced_data$Month <- as.factor(Merced_data$Month)
Merced_data$Day <- as.factor(Merced_data$Day)
Merced_data$Site <- as.factor("Merced")
Merced_data$Date <- paste(Merced_data$Year, Merced_data$Month, Merced_data$Day, sep = "_")
summary(Merced_data)


Sanjoaquin_2018 <- read.csv("Sanjoaquin_2018.csv", skip = 14, header = TRUE)
Sanjoaquin_2019 <- read.csv("Sanjoaquin_2019.csv", skip = 14, header = TRUE)
Sanjoaquin_data <- rbind(Sanjoaquin_2018, Sanjoaquin_2019)
Sanjoaquin_data$Year <- as.factor(Sanjoaquin_data$Year)
Sanjoaquin_data$Month <- as.factor(Sanjoaquin_data$Month)
Sanjoaquin_data$Day <- as.factor(Sanjoaquin_data$Day)
Sanjoaquin_data$Site <- as.factor("San Joaquin")
Sanjoaquin_data$Date <- paste(Sanjoaquin_data$Year, Sanjoaquin_data$Month, Sanjoaquin_data$Day, sep = "_")
summary(Sanjoaquin_data)

full_data <- rbind(Madera_data, Merced_data, Sanjoaquin_data)

full_data$Month_Long <- month(ymd(full_data$Date), label = TRUE, abbr = FALSE)
full_data$Date <- ymd(full_data$Date)
colnames(full_data)

full_data <- full_data[full_data$Month_Long != "May",]
full_data <- full_data[full_data$Month_Long != "September",]


# Plots of weather data
# temperature
full_data$temp_mean <- (full_data$tmmn.degC. + full_data$tmmx.degC.) / 2

a <- ggplot(full_data[full_data$Year == 2018,], aes(x = Date, y = temp_mean, color = Site, fill = Site)) +
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-06-19"), xmax=as.Date("2018-06-21"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-10"), xmax=as.Date("2018-07-13"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-31"), xmax=as.Date("2018-08-02"), ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth() +
  ylab("Mean Temperature (°C)") +
  ggtitle("2018") +
  scale_x_date(limits = as.Date(c("2018-06-15","2018-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'right')

b <- ggplot(full_data[full_data$Year == 2019,], aes(x = Date, y = temp_mean, color = Site, fill = Site)) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-06-25"), xmax=as.Date("2019-06-27"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-07-23"), xmax=as.Date("2019-07-25"), ymin=-Inf, ymax=Inf) + 
  geom_point() +
  geom_smooth() +
  ylab("Mean Temperature (°C)") +
  ggtitle("2019") +
  scale_x_date(limits = as.Date(c("2019-06-15","2019-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'right')

# Evapotranspiration 
c <- ggplot(full_data[full_data$Year == 2018,], aes(x = Date, y = pet.mm., color = Site, fill = Site)) +
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-06-19"), xmax=as.Date("2018-06-21"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-10"), xmax=as.Date("2018-07-13"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-31"), xmax=as.Date("2018-08-02"), ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth() +
  ylab("Evapotranspiration (mm)") +
  scale_x_date(limits = as.Date(c("2018-06-15","2018-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'right')

d <- ggplot(full_data[full_data$Year == 2019,], aes(x = Date, y = pet.mm., color = Site, fill = Site)) +
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-06-25"), xmax=as.Date("2019-06-27"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-07-23"), xmax=as.Date("2019-07-25"), ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth() +
  ylab("Evapotranspiration (mm)") +
  scale_x_date(limits = as.Date(c("2019-06-15","2019-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'right')

# Solar Radiation 
e <- ggplot(full_data[full_data$Year == 2018,], aes(x = Date, y = srad.Wm.2., color = Site, fill = Site)) +
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-06-19"), xmax=as.Date("2018-06-21"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-10"), xmax=as.Date("2018-07-13"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-31"), xmax=as.Date("2018-08-02"), ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth() +
  ylab("Downward Shortwave Rad.") +
  scale_x_date(limits = as.Date(c("2018-06-15","2018-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'right')

f <- ggplot(full_data[full_data$Year == 2019,], aes(x = Date, y = srad.Wm.2., color = Site, fill = Site)) +
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-06-25"), xmax=as.Date("2019-06-27"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-07-23"), xmax=as.Date("2019-07-25"), ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth() +
  ylab("Downward Shortwave Rad.") +
  scale_x_date(limits = as.Date(c("2019-06-15","2019-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'right')

# Relative Humidity 
full_data$rel_humid_mean <- (full_data$rmin... + full_data$rmax...) / 2
g <- ggplot(full_data[full_data$Year == 2018,], aes(x = Date, y = rel_humid_mean, color = Site, fill = Site)) +
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-06-19"), xmax=as.Date("2018-06-21"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-10"), xmax=as.Date("2018-07-13"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2018-07-31"), xmax=as.Date("2018-08-02"), ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth() +
  ylab("Mean Relative Humidity") +
  scale_x_date(limits = as.Date(c("2018-06-15","2018-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.title.x = element_blank(),
        legend.position = 'right')

h <- ggplot(full_data[full_data$Year == 2019,], aes(x = Date, y = rel_humid_mean, color = Site, fill = Site)) +
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-06-25"), xmax=as.Date("2019-06-27"), ymin=-Inf, ymax=Inf) + 
  annotate("rect", fill = "black", alpha= 0.5, xmin=as.Date("2019-07-23"), xmax=as.Date("2019-07-25"), ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth() +
  ylab("Mean Relative Humidity") +
  scale_x_date(limits = as.Date(c("2019-06-15","2019-08-10"))) +
  scale_color_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  theme(legend.text.align = 0,
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'right')


combine_plots <- ggarrange(a,b,c,d,e,f,g,h, nrow = 4, ncol = 2, align = 'hv', common.legend = TRUE, legend = 'right', labels = c("A", "", "B", "", "C", "", "D", ""))

ggsave("Weather_data_plot.svg", combine_plots, height = 10, width = 10)

# Stats
temp_mod <- lm(temp_mean ~ Site + Year, data=full_data)
evap_mod <- lm(pet.mm. ~ Site + Year, data=full_data)
srad_mod <- lm(srad.Wm.2. ~ Site + Year, data=full_data)
rhum_mod <- lm(rel_humid_mean ~ Site + Year, data=full_data)

anova(temp_mod)
anova(evap_mod)
anova(srad_mod)
anova(rhum_mod)

pairs(emmeans::emmeans(temp_mod,~Site))
pairs(emmeans::emmeans(evap_mod,~Site))
pairs(emmeans::emmeans(srad_mod,~Site))
pairs(emmeans::emmeans(rhum_mod,~Site))