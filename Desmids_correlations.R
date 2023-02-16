library("vegan")
library("ggplot2")
library("funrar")
library("dplyr")
library("ggvegan")
library("viridis")
library("wesanderson")
library("ggforce")

# This script performs various data explorations

original_data_set <- read.csv("data/desmid_data_set.csv", row.names = 1)
species_abundances <- original_data_set[ , 1:243]
environmental_variables <- original_data_set [ , 244:253]
# Hurlbert Swamp measurements are approximate
# conductivity is missing for Hedgehog Pond

# first, get species richness and  Shannon index for each site
species_richness <- as.data.frame(specnumber(species_abundances, MARGIN = 1))
shannon_diversity_index <- as.data.frame(diversity(species_abundances, index="shannon"))
#shannon can be replaced by simpson too
names(shannon_diversity_index)[1] <- "Shannon_Diversity_Index"
names(species_richness)[1] <- "Species_Richness"

#genera abundances
Penium <- rowSums(species_abundances%>% 
                    select(starts_with('Penium')))
Netrium <- rowSums(species_abundances%>% 
                     select(starts_with('Netrium')))
Triploceras <- rowSums(species_abundances%>% 
                         select(starts_with('Triploceras')))
Hyalotheca <- rowSums(species_abundances%>% 
                        select(starts_with('Hyalotheca')))
Desmidium <- rowSums(species_abundances%>% 
                       select(starts_with('Desmidium')))
Docidium <- rowSums(species_abundances%>% 
                      select(starts_with('Docidium')))
Bambusina <- rowSums(species_abundances%>% 
                       select(starts_with('Bambusina')))
Spirotaenia <- rowSums(species_abundances%>% 
                         select(starts_with('Spirotaenia')))
Spondylosium <- rowSums(species_abundances%>% 
                          select(starts_with('Spondylosium')))
Tetmemorus <- rowSums(species_abundances%>% 
                        select(starts_with('Tetmemorus')))
Xanthidium <- rowSums(species_abundances%>% 
                        select(starts_with('Xanthidium')))
Teilingia <- rowSums(species_abundances%>% 
                       select(starts_with('Teilingia')))
Mesotaenium <- rowSums(species_abundances%>% 
                         select(starts_with('Mesotaenium')))
Cylindrocystis <- rowSums(species_abundances%>% 
                            select(starts_with('Cylindrocystis')))
Staurodesmus <- rowSums(species_abundances%>% 
                          select(starts_with('Staurodesmus')))
Micrasterias <- rowSums(species_abundances%>% 
                          select(starts_with('Micrasterias')))
Actinotaenium <- rowSums(species_abundances%>% 
                           select(starts_with('Actinotaenium')))
Euastrum <- rowSums(species_abundances%>% 
                      select(starts_with('Euastrum')))
Staurastrum <- rowSums(species_abundances%>% 
                         select(starts_with('Staurastrum')))
Cosmarium <- rowSums(species_abundances%>% 
                       select(starts_with('Cosmarium')))
Pleurotaenium <- rowSums(species_abundances%>% 
                           select(starts_with('Pleurotaenium')))
Plaenotaenium <- rowSums(species_abundances%>% 
                           select(starts_with('Plaenotaenium')))
Closterium <- rowSums(species_abundances%>% 
                        select(starts_with('Closterium')))
#create a data frame with the columns as the abundance of the genera for each site with data.frame(). In the parentheses of the function data.frame(), type the names of the genera, which were assigned to the sums of the rows for each site for the columns corresponding to each genus. This will form a data frame with the sites as rows and the genera as columns. The values will be the abundances of the genera at each site.
genera_abundances <- data.frame(Penium, Netrium, Triploceras, Hyalotheca, Desmidium, Docidium, Bambusina, Spirotaenia, Spondylosium, Tetmemorus, Xanthidium, Teilingia, Mesotaenium, Cylindrocystis, Staurodesmus, Micrasterias, Actinotaenium, Euastrum, Staurastrum, Cosmarium, Pleurotaenium, Plaenotaenium, Closterium)

#genera presence_absence
#convert genera abundances to presence-absence data with the apply() function and the ifelse() function.
genera_presence_absence <- as.data.frame(apply(genera_abundances,2,function(x) ifelse(x>0,1,x)))
genera_presence_absence

#genera richness (number of genera)
#find the genera richness for each site by typing rowSums(yourdata).
rowSums(genera_presence_absence)
#create a data frame with the genera richness by typing as.data.frame(rowSums(genera_presence absence). Assign the data frame to genera_richness.
genera_richness <- as.data.frame(rowSums(genera_presence_absence))
genera_richness
#name the column of genera_richness that contains the number of genera Number_of_Genera
names(genera_richness)[1] <- "Number_of_Genera"
#this step creates the variable genera_richness$Number_of_Genera, which can be used for the number of genera in the correlation and regression analyses.

# make a single matrix
all_data <- environmental_variables
all_data$Species_Richness <- species_richness$Species_Richness
all_data$Shannon_Diversity_Index <- shannon_diversity_index$Shannon_Diversity_Index
all_data$Number_of_Genera <- genera_richness$Number_of_Genera

#make a historgram of each variable to check if the variables are normally distribted. To make a historgram of a certain column in an object, type hist(objectname$columnname)
hist(species_richness$Species_Richness) #does not look normally distributed
hist(genera_richness$Number_of_Genera) #looks somewhat normally distributed
hist(simpson_diversity_index$Simpson_Diversity_Index) #does not look normally distributed
hist(shannon_diversity_index$Shannon_Diversity_Index)#looks somewhat normally distributed
hist(environmental_variables$pH) #looks somewhat normally distributed
hist(environmental_variables$Conductivity) #does not look normally distributed
hist(environmental_variables$Shoreline_Length) #does not look normally distributed
hist(environmental_variables$Total_Area) #does not look normally distributed
hist(environmental_variables$Unvegetated_Open_Water_Area) #does not look normally distributed
hist(environmental_variables$Percent_Vegetation_Cover) #does not look normally distributed
hist(environmental_variables$Elev) #almost normal
#you can also check if the variables are normally distributed with the Shapiro-Wilk test. Type shapiro.test(your.data). To perform the Shapiro-Wilk test on one variable (one column), type shapiro.test(objectname$columnname).
shapiro.test(species_richness$Species_Richness)
#if p is greater than .05, then the result is not significant, and this means that the data are not significantly different from normal distribution. If p is less than .05, this means that the p-value is significant and the sample is significantly different from a normal distribution. 
#the p-value for Species_Richness is 0.00586, so Species_Richness is not normally distributed.
shapiro.test(genera_richness$Number_of_Genera)
#the p-value for Number_of_Genera is 0.3196, so Number_of_Genera is not significantly different from normally distributed.
shapiro.test(shannon_diversity_index$Shannon_Diversity_Index)
#the p-value for Shannon_Diversity_Index is 0.924, so the Simpson_Diversity_Index is not significantly different from normally distributed.
shapiro.test(environmental_variables$pH)
#the p-value for pH is 0.04045, so pH is not normally distributed.
shapiro.test(environmental_variables$Conductivity)
#the p-value for Conductivity is 0.02935, so Conductivity is not normally distributed.
shapiro.test(environmental_variables$Shoreline_Length)
#the p-value for shoreline_length is 0.02637, so shoreline length is not normally distributed.
shapiro.test(environmental_variables$Total_Area)
#the p-value for Total_Area is 0.003563, so Total_Area is not normally distributed.
shapiro.test(environmental_variables$Unvegetated_Open_Water_Area)
#the p-value for Unvegetated_Open_Water_Area is 0.0007423, so Unvegetated_Open_Water_Area is not normally distributed.
shapiro.test(environmental_variables$Percent_Vegetation_Cover)
#the p-value for Percent_Vegetation_Cover is 0.1088, so Percent_Vegetation_Cover is not significantly different from normal.
shapiro.test(environmental_variables$Elev) # 0.3127
shapiro.test(environmental_variables$Lat) #0.00107
shapiro.test(environmental_variables$Long) #0.01033

# vegetation cover, number of genera, elevation, and Shannon index are normally distributed
# no significant correlations there (or almost anywhere)

##rank regression; now using the all_data object for simplicity
#perform rank regression for the dependent variables and pH and create scatterplots of the four measures of diversity vs. pH with the rank regression line displayed.
#the example lines below can be modified to regress and plot other combinations of variables

#pH
pH_model <- lm(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$pH))
summary(pH_model)
plot(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$pH), main = "Rank Regression", xlab = "pH", ylab = "Number of Genera", pch = 16, las = 1)
abline(pH_model)
#Shannon, species richness, and generic richness all have a downward trend with pH
# only generic richness is significant 0.03595
# plain plot
plot(all_data$Number_of_Genera ~ all_data$pH, main = "", xlab = "pH", ylab = "Number of Genera", pch = 16, las = 1)


#conductivity
# conductivity data missing from Hedgehog Pond so we could exclude it
# all_data1 <- all_data[-14,]
# but no need, R will drop missing values automatically
cond_model <- lm(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Conductivity))
summary(cond_model)
plot(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Conductivity), main = "Rank Regression", xlab = "Conductivity", ylab = "Shannon_Diversity_Index", pch = 16, las = 1)
abline(cond_model)
# all diversity metrics downward trend, not significant


#Shoreline Length
shoreline_length_model <- lm(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Shoreline_Length))
summary(shoreline_length_model)
plot(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Shoreline_Length), main = "Rank Regression", xlab = "Shoreline Length", ylab = "Shannon_Diversity_Index", pch = 16, las = 1)
abline(shoreline_length_model)
# all diversity metrics have an upward trend with increasing shoreline length
# none are significant; species richness and Shannon are more pronounced

#Total Wetland Area
area_model <- lm(rank(all_data$Species_Richness) ~ 1 + rank(all_data$Total_Area))
summary(area_model)
plot(rank(all_data$Species_Richness) ~ 1 + rank(all_data$Total_Area), main = "Rank Regression", xlab = "Total Wetland Area", ylab = "Species_Richness", pch = 16, las = 1)
abline(area_model)
#Shannon index, species richness - upward trend, genera less so; none significant

# Open Water Area
open_area_model <- lm(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$Unvegetated_Open_Water_Area))
summary(open_area_model)
plot(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$Unvegetated_Open_Water_Area), main = "Rank Regression", xlab = "Unvegetated Open Water Area", ylab = "Number_of_Genera", pch = 16, las = 1)
abline(open_area_model)
#Shannon, species richness - upward trend, genera less so; none significant

#Percent Vegetation Cover
cover_model <- lm(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Percent_Vegetation_Cover))
summary(cover_model)
plot(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Percent_Vegetation_Cover), main = "Rank Regression", xlab = "Percent Vegetation Cover", ylab = "Shannon_Diversity_Index", pch = 16, las = 1)
abline(cover_model)
# slight downward trend for genera, no trend or slightly upward for richness, Shannon
# none significant

#Elevation
elevation_model <- lm(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$Elev))
summary(elevation_model)
plot(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$Elev), main = "Rank Regression", xlab = "Elevation", ylab = "Number of Genera", pch = 16, las = 1)
abline(elevation_model)
# slight upward trend, only significant for genera (0.01297)

# elevation, number of genera and Shannon are normally distributed, can do Pearson
cor.test(all_data$Elev, all_data$Number_of_Genera, method = "pearson") #0.02373
lm.elev<-lm(Number_of_Genera~Elev, data=all_data)
plot(all_data$Number_of_Genera~all_data$Elev, main = "Linear Regression", xlab = "Elevation", ylab = "Number of Genera", pch = 16, las = 1)
abline(lm.elev)

#plotting the above with ggplot
ggplot(all_data) +
  scale_color_manual(values= wes_palette("GrandBudapest1", n = 4)) +
  geom_point(mapping = aes(x=Elev, y=Number_of_Genera, color=Habitat),size = 5) +
  geom_smooth(aes(x=Elev, y=Number_of_Genera), method='lm') +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
      legend.text = element_text(size=20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title = element_text(size = 20)) + 
  labs(x = "Elevation (m)",
       y = "Number of Genera")


#Latitude
lat_model <- lm(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$Lat))
summary(lat_model)
plot(rank(all_data$Number_of_Genera) ~ 1 + rank(all_data$Lat), main = "Rank Regression", xlab = "Latitude", ylab = "Number_of_Genera", pch = 16, las = 1)
abline(lat_model)
# slight upward trend for Shannon, richness 
# none significant, genera 0.0592... close

#Longitude
long_model <- lm(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Long))
summary(long_model)
plot(rank(all_data$Shannon_Diversity_Index) ~ 1 + rank(all_data$Long), main = "Rank Regression", xlab = "Longitude", ylab = "Shannon_Diversity_Index", pch = 16, las = 1)
abline(long_model)
# slight downward trend for genera, no trend or slightly upward for richness, Shannon
# none significant
