library("vegan")
library("ggplot2")
library("funrar")
library("dplyr")
library("ggvegan")
library("viridis")
library("wesanderson")
library("ggforce")

original_data_set <- read.csv("data/desmid_data_set.csv", row.names = 1)
species_abundances <- original_data_set[ , 1:243]
environmental_variables <- original_data_set [ , 244:253]
# Hurlbert Swamp measurements are approximate
# conductivity is missing for Hedgehog Pond

# shorten names of variables and sites
colnames(environmental_variables) <- c("pH", "cond", "shore", "area", "open", "veg_cover","lat","long","elev","habitat")

# first, get species richness and  Shannon index for each site
specnumber(species_abundances, MARGIN = 1)
diversity(species_abundances, index="shannon")
#shannon can be replaced by simpson too

# species accumulation curve to estimate the completeness of sampling across sites
SA <- specaccum(species_abundances, method = "exact", permutations = 100,
          conditioned =TRUE, gamma = "jack1",  w = NULL)
plot(SA, add = FALSE, random = FALSE)

# rarefaction to to estimate the completeness of sampling within sites
# https://rstudio-pubs-static.s3.amazonaws.com/210845_860b29c5e0a643f987a80179b61bcf16.html
# this may not be the best approach for our data, but it is common
S <- specnumber(species_abundances)
# Number of INDIVIDULAS per site
raremax <- min(rowSums(species_abundances))
# rarefy, w/ raremax as input
Srare <- rarefy(species_abundances, raremax)
#Plot rarefaction results
#par(mfrow = c(1,2))
plot(S, Srare, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species",
     main = " plot(rarefy(desmids, raremax))")
abline(0, 1)
rarecurve(species_abundances, step = 20, 
          sample = raremax, 
          col = "blue", 
          cex = 0.6,
          main = "rarecurve()")
# some sites had very few desmids so we'll exclude them for rarefaction
rarecurve(species_abundances[1:12,], step = 20, sample = raremax, col = "blue", cex = 0.6,
          main = "rarecurve() on subset of data")


#before starting the ordination analyses, we need to standardize the data
# because the sampling effort wasn't equal across sites
#convert the species abundances to species relative abundances with funrar function make_relative(). 
species_rel_abundances <- make_relative(as.matrix(species_abundances))
species_rel_abundances <- as.data.frame(species_rel_abundances)

#calculate the genera abundances from the original data set
#find the abundances of each genera at the sites. 
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

#create a data frame with the columns as the abundance of the genera for each site with data.frame(). 
genera_abundances <- data.frame(Penium, Netrium, Triploceras, Hyalotheca, Desmidium, Docidium, Bambusina, Spirotaenia, Spondylosium, Tetmemorus, Xanthidium, Teilingia, Mesotaenium, Cylindrocystis, Staurodesmus, Micrasterias, Actinotaenium, Euastrum, Staurastrum, Cosmarium, Pleurotaenium, Plaenotaenium, Closterium)

#calculate genera relative abundances similar to species rel. abundances above
genera_rel_abundances <- make_relative(as.matrix(genera_abundances))
genera_rel_abundances <- as.data.frame(genera_rel_abundances)

#first pass at CCA with all variables included
# need to exclude Hedgehog Pond because it is missing conductivity
# also excluding open wetland area; it is redundant
species_rel_abundances_cca <- species_rel_abundances[-14, ]
environmental_variables_cca <- environmental_variables[-14,-5]

# shorten site names for easier plotting
rownames(species_rel_abundances_cca) <- c("Osgood", "Grassy","Casalis", "King","May","Mill","North","Bear","Trail","4thCT","Inlet","Hurlbert","Lime") 
ccamodel_species_all <- cca(species_rel_abundances_cca ~., environmental_variables_cca)

#view the results of the CCA with summary()
summary(ccamodel_species_all)
#test for the significance of the CCA model with anova(ccamodel)
anova(ccamodel_species_all)
#test for the significance of the CCA model axes with anova(ccamodel, by = "axis")
anova(ccamodel_species_all, by = "axis")
# quick plot, do not set xlim or ylim as it causes problems in scaling
ordiplot(ccamodel_species_all, display = c("sp", "bp", "si"), type = "t") 


#attempt to see what changes with absolute abundances
#species_abundances_cca <- species_abundances[-14, ]
#environmental_variables_cca <- environmental_variables[-14,-5]
#ccamodel_species_abs <- cca(species_abundances_cca ~., environmental_variables_cca)
#summary(ccamodel_species_abs)
#anova(ccamodel_species_abs)
#anova(ccamodel_species_abs, by = "axis")
# ordiplot(ccamodel_species_abs, display = c("sp", "bp", "si"), type = "t") 
# poorly sampled sites are driving the analysis, it seems

#cca with genera relative abundances
genera_rel_abundances_cca <- genera_rel_abundances[-14, ]
environmental_variables_cca <- environmental_variables[-14,-5]
rownames(genera_rel_abundances_cca) <- c("Osgood", "Grassy","Casalis", "King","May","Mill","North","Bear","Trail","4thCT","Inlet","Hurlbert","Lime") 
ccamodel_genera_all <- cca(genera_rel_abundances_cca ~., environmental_variables_cca)
summary(ccamodel_genera_all)
anova(ccamodel_genera_all)
anova(ccamodel_genera_all, by = "axis")
# quick plot, do not set xlim or ylim as it causes problems in scaling
ordiplot(ccamodel_genera_all, display = c("sp", "bp", "si"), type = "t") 
# at this level the plot is driven by small genera with few species

# subsetting data based on analyses above
# getting rid of Lime Pond - it only has 5 species, in addition to Hedgehog Pond (no conductivity)
# also excluding open wetland area; it is redundant
# it's weird to have more env. variables than sites 
# also excluding longitude and habitat
# area and shore are correlated, so area is excluded too
species_rel_abundances_cca2 <- species_rel_abundances[-c(13,14), ]
environmental_variables_cca2 <- environmental_variables[-c(13,14),-c(4,5,8)]
#save habitat separately, then exclude
habitat <- environmental_variables_cca2$habitat
environmental_variables_cca2 <- environmental_variables_cca2[,-7]

# shorten site names for easier plotting
rownames(species_rel_abundances_cca2) <- c("Osgood", "Grassy","Casalis", "King","May","Mill","North","Bear","Trail","4thCT","Inlet","Hurlbert") 
ccamodel_species_subset <- cca(species_rel_abundances_cca2 ~., environmental_variables_cca2)

#view the results of the CCA with summary()
summary(ccamodel_species_subset)
#test for the significance of the CCA model with anova(ccamodel)
anova(ccamodel_species_subset)
#test for the significance of the CCA model axes with anova(ccamodel, by = "axis")
anova(ccamodel_species_subset, by = "axis")
# quick plot, do not set xlim or ylim as it causes problems in scaling
ordiplot(ccamodel_species_subset, display = c("sp", "bp", "si"), type = "t") 

# nicer plot with ggvegan
autoplot(ccamodel_species_subset)

pal <- wes_palette(name = "Darjeeling1", type = "continuous")

cca.res<-summary(ccamodel_species_subset)
cca.sites <-data.frame(cca.res$sites)

site.names <- c("Osgood", "Grassy","Casalis", "King","May","Mill","North","Bear","Trail","4thCT","Inlet","Hurlbert") 
ord_df<-data.frame(Site=site.names,CCA1=cca.sites$CCA1,CCA2=cca.sites$CCA2)
ord_df$pH <- environmental_variables_cca2$pH

ord_df$Habitat <-habitat
#row.names(ord_df) <- site.names

cca.species<-data.frame(cca.res$species)
mul <- ggvegan:::arrowMul(cca.res$biplot,
                          rbind(cca.sites, cca.species))
cca.env <- data.frame(cca.res$biplot * mul)
cca.env<-data.frame(cca1=cca.env$CCA1,cca2=cca.env$CCA2, Species = rownames(cca.env))
cca.species<-data.frame(Cca1=cca.species$CCA1,Cca2=cca.species$CCA2,species=rownames(cca.species))

ggplot(ord_df) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  #scale_color_viridis(discrete = FALSE, option = "D") +
  scale_color_manual(values= wes_palette("GrandBudapest1", n = 4)) +
  geom_segment(data = cca.env,
               aes(x = 0, xend = cca1, y = 0, yend = cca2),
               arrow = arrow(length = unit(0.5, "cm")), size=1, colour = "grey") +
  geom_point(mapping = aes(x=CCA1, y=CCA2, color=Habitat),size = 7)+
  geom_text(aes(x=CCA1, y=CCA2, label = Site),
            size = 6,, nudge_x = 0.2, nudge_y = 0.1, colour = "#808080") +
  geom_text(data = cca.env, 
            aes(x = cca1*1.2, y = cca2*1.2, label = Species),
            size = 6) +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
  labs(x = "CCA1",
       y = "CCA2")


#https://rpubs.com/Roeland-KINDT/694021 more plotting options
# https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
#library("BiodiversityR") # has weird dependencies on mac
# but some graphing works well with ggforce

# this theme lacks the grey background which can be nice for visibility
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 12, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

# there is only one HSF site (highland spruce forest), 4th CT Lake. 
# It is not that much higher in elevation than others
# merging it with LSF, which are nearby
ord_df2 <- ord_df
ord_df2$Habitat[ord_df2$Habitat == 'HSF'] <- 'SF'
ord_df2$Habitat[ord_df2$Habitat == 'LSF'] <- 'SF'

plotgg2 <- ggplot(ord_df2) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  scale_color_manual(values= wes_palette("GrandBudapest1", n = 4)) +
  geom_segment(data = cca.env,
               aes(x = 0, xend = cca1, y = 0, yend = cca2),
               arrow = arrow(length = unit(0.5, "cm")), size=1, colour = "grey") +
  geom_point(mapping = aes(x=CCA1, y=CCA2, color=Habitat),size = 7)+
  geom_text(aes(x=CCA1, y=CCA2, label = Site),
            size = 6,, nudge_x = 0.2, nudge_y = 0.1, colour = "#808080") +
  geom_text(data = cca.env, 
            aes(x = cca1*1.2, y = cca2*1.2, label = Species),
            size = 6) +
  labs(x = "CCA1",
       y = "CCA2") +
  geom_mark_ellipse(data=ord_df2, 
                    aes(x=CCA1, y=CCA2, colour=Habitat, 
                        fill=after_scale(alpha(colour, 0.2))), 
                    expand=0, size=0.2, show.legend=FALSE) +
 # geom_mark_hull(data=ord_df2, 
  #               aes(x=CCA1, y=CCA2, colour=Habitat, 
  #                   fill=after_scale(alpha(colour, 0.2))), 
  #               concavity=0.1, size=0.2, show.legend=FALSE) +
  BioR.theme
  
  plotgg2
