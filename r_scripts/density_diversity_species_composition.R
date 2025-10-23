# set the working directory to the r_scripts folder
setwd("~/Desktop/ENVI408/r_scripts")

#                                                                      PACKAGES :::::::::####

  # create vector of required packages
packages_species <- c(
  "magrittr", "tidyr", "dplyr", "tibble", "ggplot2", "stringr",    # tidyverse packages
  "microbiome", "phyloseq", "vegan", "pairwiseAdonis",             # comm. comp. & phyloseqs
  "cowplot",                                                       # making good plots
  "car")                                                           # for levene test

  # install packages (only need to run once)
#lapply(packages, install.packages)

# load them
lapply(packages_species, library, character.only = T)

  # if there's an error message saying "there is no package called lifecycle"
# and you're on windows, install Rtools from:
# https://cran.r-project.org/bin/windows/Rtools/
# and then try and rerun it again.
# per https://stackoverflow.com/questions/70185951/error-using-librarydplyr-in-rstudio-there-is-no-package-called-lifecycle

# set seed
set.seed(2018050519)

#

#                                                                                        SETUP ::::::####
# importing + data manip + creating phyloseq objects                                     ####
# ideally should be able to just run lines 19-363 and it should all work                 ####
######################################################################################## ####
#                                                                   Import data :::::::::####

# import raw data
raw_eucdata <- read.csv("data/eucalypt_data.csv")

#                                                                   Format data :::::::::####

# format the data
raw_eucdata <- raw_eucdata %>% 
  
  # filter for living eucs
  
  filter(alive == T) %>% 
  
  # add column where amp/ter and eug/glob are in stringybark and redgum complexes
  
  mutate(species_complex = case_when(
    (species == "Eucalyptus amplifolia") ~ "Redgum complex",
    (species == "Eucalyptus tereticornis") ~ "Redgum complex",
    (species == "Eucalyptus amplifolia or tereticornis") ~ "Redgum complex",
    
    (species == "Eucalyptus eugenioides") ~ "Stringybark complex",
    (species == "Eucalyptus globoidea") ~ "Stringybark complex",
    (species == "Eucalyptus eugenioides or globoidea") ~ "Stringybark complex",
    
    .default = species)) %>% 
  
  # add column for patch size
  mutate(size = case_when(
    (site_id == "ALR" | site_id == "BBR" | site_id == "CRM") ~ "Large",
    (site_id == "HGA" | site_id == "RGR" | site_id == "WMP") ~ "Medium",
    (site_id == "CSP" | site_id == "MCP" | site_id == "PAP") ~ "Small",
    .default = NA))


# create df for analysing abundance + diversity by patch size

diversities_site <- raw_eucdata %>%
  
  filter(species_complex != "N/A") %>% 
  
  # aggregate by size + site id + species
  
  group_by(site_id, size, species_complex) %>% 
  summarise(n = n()) %>% 
  
  # pivot wide
  pivot_wider(names_from = species_complex, values_from = n) %>% 
  
  # replace nas with 0s 
  replace(is.na(.), 0) %>% 
  
  # add column for the number of sampling points for each site size category
  
  mutate(no_pts = case_when((size == "Large") ~ 5,
                            (size == "Medium") ~ 3,
                            (size == "Small") ~ 1,
                            .default = NA)) %>% 
  
  ungroup() %>% 
  
  # calculate abundance, density/ha, species richness, shannon, simpson
  mutate(abundance = rowSums(.[3:9])/.$no_pts,
         density_ha = (rowSums(.[3:9])/.$no_pts) / (pi * 25 * 2) * 10000,
         richness = rowSums(.[3:9] > 0),
         shannon = diversity(.[3:9], index = "shannon"),
         simpson = diversity(.[3:9], index = "simpson"))

#

#                                                        Create OSU (ASV) files :::::::::####

################################################################################:
#                             SPECIES BY SAMPLING POINT                        #
################################################################################:

# species ASV file
# sum each species for each sampling point, excluding unidentified trees
sum_spp_pt <- subset(raw_eucdata, species_complex != "N/A") %>% 
  group_by(point_id, species_complex) %>% 
  summarise(n = n())

# pivot to format the data properly
pivoted_sum_spp_pt <- pivot_wider(names_from = point_id, values_from = n, data = sum_spp_pt)

# rename "species" to "ASV" since the package we're using only works if the
# species column is named ASV (it's made for microbiomes)
pivoted_sum_spp_pt <- rename(pivoted_sum_spp_pt, "ASV" = "species_complex")

# replace NAs with zeroes
pivoted_sum_spp_pt[is.na(pivoted_sum_spp_pt)] <- 0

# turn column 1 into rownames
ASV_species_pt <- pivoted_sum_spp_pt[,-1]
rownames(ASV_species_pt) <- pivoted_sum_spp_pt$ASV

# turn rownames back into column 1 so that the ASV column name shows up on the table
ASV_species_pt <- ASV_species_pt %>% 
  rownames_to_column('ASV')

# write pivoted_sum_spp as ASV_species.txt
write.table(ASV_species_pt, file = "ASV_metadata_taxonomy/ASV_species_pt.txt",
            col.names = NA, row.names = T, sep = "\t")

################################################################################:
#                                 SPECIES BY SITE                              #
################################################################################:

# species ASV file
# sum each species for each sampling point, excluding unidentified trees
sum_spp_site <- subset(raw_eucdata, species != "N/A") %>% 
  group_by(site_id, species) %>% 
  summarise(n = n())

# pivot to format the data properly
pivoted_sum_spp_site <- pivot_wider(names_from = site_id, values_from = n, data = sum_spp_site)

# rename "species" to "ASV" since the package we're using only works if the
# species column is named ASV (it's made for microbiomes)
pivoted_sum_spp_site <- rename(pivoted_sum_spp_site, "ASV" = "species")

# replace NAs with zeroes
pivoted_sum_spp_site[is.na(pivoted_sum_spp_site)] <- 0

# turn column 1 into rownames
ASV_species_site <- pivoted_sum_spp_site[,-1]
rownames(ASV_species_site) <- pivoted_sum_spp_site$ASV

# turn rownames back into column 1 so that the ASV column name shows up on the table
ASV_species_site <- ASV_species_site %>% 
  rownames_to_column('ASV')

# write pivoted_sum_spp_site as ASV_species_site.txt
write.table(ASV_species_site, file = "ASV_metadata_taxonomy/ASV_species_site.txt",
            col.names = NA, row.names = T, sep = "\t")

# tidying up the environment
rm(pivoted_sum_spp_pt, sum_spp_pt, pivoted_sum_spp_site, sum_spp_site)

#                                                       Create phyloseq objects :::::::::####

# most of this code is from Dr Allison Mertin. thanks Allison

# load useful function
KillZeroRCs <- function(x) {
  x[ which( rowSums(x) != 0) , ] -> x
  x[ , which( colSums(x) != 0) ] -> x
  return(x)
}

################################################################################:
#                         IMPORTING ASV TABLES                                 #
################################################################################:

# removing useless col 1 for these because it's useless
# species by sampling pt
ASV_species_pt <- read.table(
  "ASV_metadata_taxonomy/ASV_species_pt.txt",
  header = TRUE,
  sep = "\t",
  row.names = 2
)
ASV_species_pt <- ASV_species_pt[,-1]

# species by site
ASV_species_site <- read.table(
  "ASV_metadata_taxonomy/ASV_species_site.txt",
  header = TRUE,
  sep = "\t",
  row.names = 2
)
ASV_species_pt <- ASV_species_pt[,-1]

################################################################################:
#                         IMPORTING TAXA TABLES                                #
################################################################################:

# import taxonomy table for ASV_demo and ASV_species
# species
TAXA_species <- read.table(
  "ASV_metadata_taxonomy/Taxonomy_species_noNAs_complex.txt",
  sep = ",", fill = TRUE, row.names = 1)

# add levels of taxonomy to taxa tables
colnames(TAXA_species) <- c("Genus","Complex")

################################################################################:
#                         IMPORTING METADATA TABLES                             #
################################################################################:

# import metadata
META_site <- read.table(
  "ASV_metadata_taxonomy/Metadata_site.txt",
  header = TRUE, sep = ",", row.names = 1)

META_pt <- read.table(
  "ASV_metadata_taxonomy/Metadata_pt.txt",
  header = TRUE, sep = ",", row.names = 1)

################################################################################:
#                      CREATING PHYLOSEQ OBJECTS                               #
################################################################################:

# convert ASVs + taxa tables to matrices in prep for making a phyloseq
ASV_species_pt <- as.matrix(ASV_species_pt)
ASV_species_site <- as.matrix(ASV_species_site)
TAXA_species <- as.matrix(TAXA_species)

# creating phyloseqs
# species
phy_species_pt <- phyloseq(
  otu_table(ASV_species_pt,taxa_are_rows = T),
  tax_table(TAXA_species),
  sample_data(META_pt)
)

phy_species_site <- phyloseq(
  otu_table(ASV_species_site,taxa_are_rows = T),
  tax_table(TAXA_species),
  sample_data(META_site)
)

# clean up environment
rm(ASV_species_pt, ASV_species_site,
   TAXA_species, META_site, META_pt)
#                                                           Create plotting dfs :::::::::####

############################################:
# stacked bar plot of species comp by site #:
############################################:

# merge samples by site
complex_site <- merge_samples(phy_species_pt, "site") %>% 
  
  # turn counts into relative abundance
  
  transform_sample_counts(function(x) 100 * x/sum(x)) %>% 
  
  # psmelt it + turn complex col into character
  
  psmelt() %>% 
  mutate(Complex = as.character(.$Complex)) %>% 
  
  # add cols to order bars + legend + set palette order
  
  mutate(order = case_when((Sample == "ALR" | Sample == "BBR" | Sample == "CRM") ~ 3,
                           (Sample == "HGA" | Sample == "RGR" | Sample == "WMP") ~ 2,
                           (Sample == "CSP" | Sample == "MCP" | Sample == "PAP") ~ 1,
                           .default = NA)) %>% 
  mutate(pattern = case_when((OTU == "Eucalyptus longifolia") ~ 1,
                             (OTU == "Redgum complex") ~ 2,
                             (OTU == "Stringybark complex") ~ 3,
                             (OTU == "Angophora floribunda") ~ 4,
                             (OTU == "Eucalyptus paniculata") ~ 5,
                             (OTU == "Eucalyptus bosistoana") ~ 6,
                             (OTU == "Eucalyptus pilularis") ~ 7,
                             (OTU == "N/A") ~ 8,
                             .default = NA))
##########################################################:
# plots of abundance + density + diversity by patch size #:
##########################################################:

diversities_site_size <- diversities_site %>% 
  
  # group by size and get mean + sd for abundance/density and richness
  group_by(size) %>% 
  summarise(mean_abundance = mean(abundance), 
            sd_abundance = sd(abundance),
            
            mean_density_ha = mean(density_ha), 
            sd_density_ha = sd(density_ha),
            
            mean_richness = mean(richness), 
            sd_richness = sd(richness),
            
            mean_shannon = mean(shannon), 
            sd_shannon = sd(shannon),
            
            mean_simpson = mean(simpson), 
            sd_simpson = sd(simpson)) %>% 
  
  # add col for ordering
  
  mutate(order = case_when((size == "Large") ~ 3,
                           (size == "Medium") ~ 2,
                           (size == "Small") ~ 1,
                           .default = NA))

#####################################:
# plot of raw abundance per species #:
#####################################:

# df for plotting abundance by species
abundance_plot <- raw_eucdata %>%  
  # get confidence of the ID
  mutate(confidence = case_when((species == "N/A" & unconfirmed_id == T) ~ "Low",
                                (unconfirmed_species == "N/A") ~ "Unidentified",
                                (species != "N/A") ~ "High",
                                .default = NA)) %>% 
  
  # group by species and confident; get abundances for each
  
  group_by(unconfirmed_species, confidence) %>% 
  summarise(n = n()) %>% 
  
  # col for ordering spp. on plot
  mutate(order = case_when((unconfirmed_species == "Eucalyptus longifolia") ~ 1,
                           (unconfirmed_species == "Eucalyptus amplifolia or tereticornis") ~ 2,
                           (unconfirmed_species == "Eucalyptus amplifolia") ~ 3,
                           (unconfirmed_species == "Eucalyptus tereticornis") ~ 4,
                           (unconfirmed_species == "Eucalyptus eugenioides or globoidea") ~ 5,
                           (unconfirmed_species == "Eucalyptus eugenioides") ~ 6,
                           (unconfirmed_species == "Eucalyptus globoidea") ~7,
                           (unconfirmed_species == "Angophora floribunda") ~ 8,
                           (unconfirmed_species == "Eucalyptus paniculata") ~ 9,
                           (unconfirmed_species == "Eucalyptus bosistoana") ~ 10,
                           (unconfirmed_species == "Eucalyptus pilularis") ~ 11,
                           (unconfirmed_species == "N/A") ~ 12,
                           .default = NA),
         
  # find and replace generic names with abbreviations + change amp/ter/glob/eug to complex
         
         unconfirmed_species = str_replace(unconfirmed_species,
                                           "Eucalyptus amplifolia or tereticornis",
                                           "Redgum complex"),
         unconfirmed_species = str_replace(unconfirmed_species,
                                           "Eucalyptus eugenioides or globoidea",
                                           "Stringybark complex"),
         unconfirmed_species = str_replace(unconfirmed_species,
                                           "Eucalyptus",
                                           "E."),
         unconfirmed_species = str_replace(unconfirmed_species,
                                           "Angophora",
                                           "A."),
         unconfirmed_species = str_replace(unconfirmed_species,
                                           "N/A",
                                           "Unidentified"))

#

#                                                                                     ANALYSIS ::::::####
######################################################################################## ####
#                                      PERMANOVA for species composition by site :::::::::####

################################################################################:
#                           SPECIES COMPOSITION BY SITE                        #
################################################################################:

# log transform
species_perm <- microbiome::transform(phy_species_pt, transform = "log")

# generate bray distance matrix
species_bray <- phyloseq::distance(species_perm, method = "bray") 
species_jacard <- phyloseq::distance(species_perm, method = "jaccard", binary = TRUE)
species_jacard


# run dispersion test + plot it
disp_species_bray <- vegan::betadisper(species_bray, phyloseq::sample_data(species_perm)$site)
disp_species_bray
plot(disp_species_bray)

# run ADONIS test
vegan::adonis2(species_bray ~ phyloseq::sample_data(species_perm)$site, permutations = 20000)

# run post-hoc comparisons test
pairwise.adonis(species_bray, phyloseq::sample_data(species_perm)$site, sim.method = "bray")

# run SIMPER analysis
simper_sites<- vegan::simper(t(otu_table(phy_species_pt)), sample_data(phy_species_pt)$site, permutations = 999)
summary(simper_sites)

#                                      ANOVA for diversity + pop density by site :::::::::####

# testing assumptions
# normality
par(mfrow = c(2,2))
hist(diversities_site$abundance)
hist(diversities_site$richness)
hist(diversities_site$shannon)
hist(diversities_site$simpson)
par(mfrow = c(1,1))
shapiro.test(diversities_site$abundance)
shapiro.test(diversities_site$richness)
shapiro.test(diversities_site$shannon)
shapiro.test(diversities_site$simpson)

# all normal! 

# test homoscedasticity
leveneTest(abundance ~ size, data = diversities_site)
leveneTest(richness ~ size, data = diversities_site)
leveneTest(shannon ~ size, data = diversities_site)
leveneTest(simpson ~ size, data = diversities_site)

# variance homogeneous! continue 

# run + display anovas for each
abundance_anova_size <- aov(abundance ~ size, data  = diversities_site)
richness_anova_size <- aov(richness ~ size, data  = diversities_site)
shannon_anova_size <- aov(shannon ~ size, data  = diversities_site)
simpson_anova_size <- aov(simpson ~ size, data  = diversities_site)
summary(abundance_anova_size)
summary(richness_anova_size)
summary(shannon_anova_size)
summary(simpson_anova_size)

# only abundance sig. run tukey post-hoc on it
TukeyHSD(abundance_anova_size)

# testing assumptions
# normality
par(mfrow = c(2,2))
hist(diversities_site$abundance)
hist(diversities_site$richness)
hist(diversities_site$shannon)
hist(diversities_site$simpson)
par(mfrow = c(1,1))
shapiro.test(diversities_site$abundance)
shapiro.test(diversities_site$richness)
shapiro.test(diversities_site$shannon)
shapiro.test(diversities_site$simpson)

# all normal! 

# test homoscedasticity
leveneTest(abundance ~ size, data = diversities_site)
leveneTest(richness ~ size, data = diversities_site)
leveneTest(shannon ~ size, data = diversities_site)
leveneTest(simpson ~ size, data = diversities_site)

# variance homogeneous! continue 

# run + display anovas for each
density_anova_size <- aov(density_ha ~ size, data  = diversities_site)
richness_anova_size <- aov(richness ~ size, data  = diversities_site)
shannon_anova_size <- aov(shannon ~ size, data  = diversities_site)
simpson_anova_size <- aov(simpson ~ size, data  = diversities_site)

anova_list <- list(density_anova_size, richness_anova_size,
                   shannon_anova_size, simpson_anova_size)
lapply(anova_list, summary)
# 1 is density, 2 is richness, 3 is shannon, 4 is simpson

#                                                                                     PLOTTING ::::::####
######################################################################################## ####
#                                                   Species composition by site :::::::::####

# plot it
complex_plot_site <- ggplot(complex_site, aes(x = reorder(Sample, order), y = Abundance, fill = reorder(Complex, pattern))) +  
  geom_bar(stat = "identity", width = 0.85, color = "black", linewidth = 0.25) +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Site", y = "Relative abundance (%)") +
  scale_fill_manual(values = c("#94c9b4","#ea9a4f","#698d44","#aa74b1","#5e4c39","#be9f99","#f4d5a2"),
                    name = "Species") + theme_cowplot()
complex_plot_site
#                                                          Abundance by species :::::::::####

# plot it
abundanceplot <- ggplot(abundance_plot, 
                        aes(x = reorder(unconfirmed_species, order), 
                            y = n, fill = confidence)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 315, vjust = 1, hjust=0)) +
  scale_fill_manual(values = c("#79c638","#4d7f24","#2c4814"), name = "Identification 
confidence") + 
  labs(x = "Species",
       y = "Abundance")
abundanceplot


#                                                 Density and diversity by size :::::::::####

# density by size
density_size <- ggplot(diversities_site_size, 
                       aes(x = reorder(size, order), y = mean_density_ha)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_density_ha - sd_density_ha,
                    ymax = mean_density_ha + sd_density_ha),
                width = 0.3) + 
  xlab("Size") + ylab("Mean population density (trees/ha)") + theme_cowplot() 

# richness by size
richness_size <- ggplot(diversities_site_size, 
                        aes(x = reorder(size, order), y = mean_richness)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_richness - sd_richness,
                    ymax = mean_richness + sd_richness),
                width = 0.3) + 
  xlab("Size") + ylab("Mean species richness") + theme_cowplot() 

# shannon by size
shannon_size <- ggplot(diversities_site_size, 
                       aes(x = reorder(size, order), y = mean_shannon)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_shannon - sd_shannon,
                    ymax = mean_shannon + sd_shannon),
                width = 0.3) + 
  xlab("Size") + ylab("Mean Shannon diversity") + theme_cowplot() 

# simpson by size
simpson_size <- ggplot(diversities_site_size, 
                       aes(x = reorder(size, order), y = mean_simpson)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_simpson - sd_simpson,
                    ymax = mean_simpson + sd_simpson),
                width = 0.3) + 
  xlab("Size") + ylab("Mean Simpson diversity") + theme_cowplot() 

# plot in a 4x4 grid
density_diversity_size <- plot_grid(density_size, richness_size, 
                                    shannon_size, simpson_size)

#                                                          Plots used in thesis :::::::::####
  # this is just all the plots made here but i named it this for consistency
# basically just plotting + saving

###############################:
# SPECIES COMPOSITION BY SITE #:
###############################:

complex_plot_site
save_plot("plots/eucalypt_species_stackedbar_site.png", 
          complex_plot_site, 
          base_width = 10,
          base_height = 9)

########################:
# ABUNDANCE BY SPECIES #:
########################:

abundanceplot
save_plot("plots/eucalypt_species_abundance.png",
          abundanceplot, 
          base_height = 9,
          base_width = 10)

#################################:
# DENSITY AND DIVERSITY BY SIZE #:
#################################:

density_diversity_size
save_plot("plots/eucalypt_density_diversity_size.png",
          density_diversity_size,
          base_height = 10, base_width = 10)

