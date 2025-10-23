# set the working directory to the r_scripts folder
setwd("~/Desktop/ENVI408/r_scripts")

#                                                                      PACKAGES :::::::::####

# make vector of packages
packages_insects <- c("magrittr", "tidyr", "dplyr", "ggplot2", # tidyverse packages
              "vegan",                                 # used for diversity indices
              "cowplot",                               # make plots look good
              "car")                                   # levene test

# install them (only need to run once)
#lapply(packages_insects, install.packages)

# load packages
lapply(packages_insects, library, character.only = T)

  # if there's an error message saying "there is no package called lifecycle"
# and you're on windows, install Rtools from:
# https://cran.r-project.org/bin/windows/Rtools/
# and then try and rerun it again.
# per https://stackoverflow.com/questions/70185951/error-using-librarydplyr-in-rstudio-there-is-no-package-called-lifecycle


#                                                                                        SETUP ::::::####
# bunch of data, creation + manipulation to prep for analysis and plots                  ####
# ideally should be able to just run lines 19-128 and it should all work                 ####
######################################################################################## ####
#                                                                   Import data :::::::::####

raw_insectdata <- read.csv("data/insect_data.csv")

#                                                                   Format data :::::::::####
trimmed_insectdata <- raw_insectdata %>%
  
  # select relevant columns
  
  select(c(3:4, 20:28, 31)) %>% 
  
  # add column for patch size category
  
  mutate(size = case_when((site_id == "ALR" | site_id == "BBR" | site_id == "CRM") ~ "Large",
                          (site_id == "HGA" | site_id == "RGR" | site_id == "WMP") ~ "Medium",
                          (site_id == "CSP" | site_id == "MCP" | site_id == "PAP") ~ "Small",
                          .default = NA)) %>% 
  
  # calculate abundance, ordering richness, and shannon and simpson diversity
  
  mutate(insect_abundance = rowSums(.[3:12]),
         insect_richness = rowSums(.[3:12] > 0),
         insect_shannon = diversity(.[3:12], index = "shannon"),
         insect_simpson = diversity(.[3:12], index = "simpson")) %>% 
  
  # column for ordering bars on plots
  
  mutate(ordering = case_when((site_id == "ALR" | site_id == "BBR" | site_id == "CRM") ~ 3,
                           (site_id == "HGA" | site_id == "RGR" | site_id == "WMP") ~ 2,
                           (site_id == "CSP" | site_id == "MCP" | site_id == "PAP") ~ 1,
                           .default = NA))

# make dataframe aggregated by size

trimmed_insect_data_site <- trimmed_insectdata %>% 
  
  # aggregate by site and patch size
  
  group_by(site_id, size, ordering) %>%
  
  # calculate mean + sd for abundance, richness, shannon, simpson by site
  
  dplyr::summarise(
    mean_abundance = mean(insect_abundance), 
    sd_abundance = sd(insect_abundance),
    
    mean_richness = mean(insect_richness), 
    sd_richness = sd(insect_richness),
    
    mean_shannon = mean(insect_shannon), 
    sd_shannon = sd(insect_shannon),
    
    mean_simpson = mean(insect_simpson), 
    sd_simpson = sd(insect_simpson)) %>% 
  
  # replace nas with 0
  
  replace(is.na(.), 0)

# make df for plotting abundance and diversity metrics by patch size

trimmed_insectdata_size <- trimmed_insectdata %>%
  
  # aggregate by patch size
  
  group_by(size, ordering) %>% 
  
  # calculate mean + sd for abundance, richness, shannon, simpson by size
  
  summarise(mean_abundance = mean(insect_abundance), sd_abundance = sd(insect_abundance),
            mean_richness = mean(insect_richness), sd_richness = sd(insect_richness),
            mean_shannon = mean(insect_shannon), sd_shannon = sd(insect_shannon),
            mean_simpson = mean(insect_simpson), sd_simpson = sd(insect_simpson))

# make dfs for plotting stacked barplots

order_table_site <- trimmed_insectdata %>%
  
  # pivot longer + create col of order
  
  pivot_longer(cols = Coleoptera:Trichoptera, 
               names_to = "order", 
               values_to = "n", 
               names_repair = "unique") %>% 
 
  # aggregate by site id and get n of insects in that order per site
  
  group_by(site_id, size, ordering, order) %>% 
  summarise(n = sum(n))

order_table_size <- order_table_site %>% 
  
  # aggregate by size and get n of insects in that order per site
  
  group_by(size, ordering, order) %>% 
  summarise(n = sum(n))

# make df for order abundance plot
order_table <- order_table_site %>% 
  group_by(order) %>% 
  summarise(n = sum(n))


#

#                                                                            STATISTICAL TESTS ::::::####
######################################################################################## ####
#                                                 Diversity & abundance by size :::::::::####

# assumption testing

# test normality with hists + shapiro-wilks

  # plot histograms
par(mfrow = c(2,2))
hist(trimmed_insectdata$insect_abundance)
hist(trimmed_insectdata$insect_richness)
hist(trimmed_insectdata$insect_shannon)
hist(trimmed_insectdata$insect_simpson)
par(mfrow = c(1,1))
# eh. reasonably normal

  # shapiro tests
lapply
shapiro.test(trimmed_insectdata$insect_abundance)
shapiro.test(trimmed_insectdata$insect_richness)
shapiro.test(trimmed_insectdata$insect_shannon)
shapiro.test(trimmed_insectdata$insect_simpson)

# yeah reasonably normal

# test homoscedasticity with levene test
leveneTest(insect_abundance ~ size, data = trimmed_insectdata)
leveneTest(insect_richness ~ size, data = trimmed_insectdata)
leveneTest(insect_shannon ~ size, data = trimmed_insectdata)
leveneTest(insect_simpson ~ size, data = trimmed_insectdata)

# run aovs for abundance and diversity
abundance_aov <- aov(mean_abundance ~ size, data = trimmed_insect_data_site)
richness_aov <- aov(mean_richness ~ size, data = trimmed_insect_data_site)
shannon_aov <- aov(mean_shannon ~ size, data = trimmed_insect_data_site)
simpson_aov <- aov(mean_simpson ~ size, data = trimmed_insect_data_site)
aov_list <- list(abundance_aov, richness_aov, shannon_aov, simpson_aov)


# plot residuals v fitted to check homoscedasticity
plot(abundance_aov, 2)
plot(richness_aov, 2)
plot(shannon_aov, 2)
plot(simpson_aov, 2)
# all fairly ok

# get summaries for all aov
lapply(aov_list, summary)
# 1 is abundance, 2 is richness, 3 is shannon, 4 is simpson


#                                                                                    PLOTTING  ::::::####
######################################################################################## ####
#                                                 Diversity & abundance by size :::::::::####

# plotting mean abundance + diversity by patch size

#############:
# ABUNDANCE #:
#############:

# plot insect abundance by patch size
abundance_size <- ggplot(trimmed_insectdata_size, 
                         aes(x = reorder(size, ordering), 
                             y = mean_abundance)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_abundance - sd_abundance, 
                    ymax = mean_abundance + sd_abundance),
                width = 0.3) + 
  labs(x = "Patch size",
       y = "Mean abundance",
       title = "a)") +
  theme_cowplot()
abundance_size

##################:
# ORDER RICHNESS #:
##################:

# plot insect richness by patch size
richness_size <- ggplot(trimmed_insectdata_size, 
                        aes(x = reorder(size, ordering), 
                            y = mean_richness)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_richness - sd_richness, 
                    ymax = mean_richness + sd_richness),
                width = 0.3) + 
  labs(x = "Patch size",
       y = "Mean order richness",
       title = "b)") +
  theme_cowplot()
richness_size

#####################:
# SHANNON DIVERSITY #:
#####################:

# plot of shannon diversity by patch size
shannon_size <- ggplot(trimmed_insectdata_size, 
                       aes(x = reorder(size, ordering), 
                           y = mean_shannon)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_shannon - sd_shannon, 
                    ymax = mean_shannon + sd_shannon),
                width = 0.3) + 
  labs(x = "Patch size",
       y = "Mean Shannon diversity",
       title = "c)") +
  theme_cowplot()
shannon_size

#####################:
# SIMPSON DIVERSITY #:
#####################:

# plot of simpson diversity by patch size
simpson_size <- ggplot(trimmed_insectdata_size, aes(x = reorder(size, ordering), y = mean_simpson)) +
  geom_bar(stat = "identity", color = "black", fill = "#79c638") + 
  geom_errorbar(aes(ymin = mean_simpson - sd_simpson, ymax = mean_simpson + sd_simpson),
                width = 0.3) + 
  labs(x = "Patch size",
       y = "Mean Simpson diversity",
       title = "d)") +
  theme_cowplot()
simpson_size

####################:
# ALL TOGETHER NOW #:
####################:

# plot all together
diversity_abundance_plot <- plot_grid(abundance_size, richness_size,
                                      shannon_size, simpson_size,
                                      nrow = 2)
#

#                                        Stacked bar of orders by size and site :::::::::####

#######################:
# STACKED BAR BY SITE #:
#######################:

# make palette for barplots
stackedbarpalette <- c("#5a3f15", "#cc3838", "#72c239","#f3ab31","#4bd5df",
                       "#027907","#797971","#c190b9","#728d9f")

#
order_plot_site <- ggplot(order_table_site, 
                          aes(x = reorder(site_id, ordering), 
                              y = n,fill = order)) +  
  geom_bar(stat = "identity", position = "fill",
           width = 0.85, color = "black", linewidth = 0.25) +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Site", y = "Relative abundance",
       title = "a)") +
  scale_fill_manual(values = stackedbarpalette) + theme_cowplot()

order_plot_site

#######################:
# STACKED BAR BY SIZE #:
#######################:

order_plot_size <- ggplot(order_table_size,
                          aes(x = reorder(size, ordering), 
                              y = n, fill = order)) +  
  geom_bar(stat = "identity", position = "fill",
           width = 0.85, color = "black", linewidth = 0.25) +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Patch size", y = "Relative abundance",
       title = "b)") +
  scale_fill_manual(values = stackedbarpalette) + theme_cowplot()
order_plot_size

#####################:
# ALL TOGETHER NOW! #:
#####################:

stackedbar_site_size <- plot_grid(order_plot_site, order_plot_size, nrow = 2)
stackedbar_site_size

# 

#                                                Bar plot of abundance by order :::::::::####
######################:
# ABUNDANCE BY ORDER #:
######################:

order_plot <- ggplot(order_table, aes(x = order, y = n)) +
  geom_bar(stat = "identity", colour = "black", fill = "#79c638") +
  labs(x = "Order",
       y = "Count") + 
  theme_cowplot()

order_plot

#

#                                                          Plots used in thesis :::::::::####

###################################:
# ABUNDANCE AND DIVERSITY BY SIZE #:
###################################:

diversity_abundance_plot

save_plot("plots/insect_diversity_abundance_site.png", diversity_abundance_plot, 
          base_height = 10, base_width = 10)

######################################:
# ORDER STACKED BAR BY SITE AND SIZE #:
######################################:

stackedbar_site_size

save_plot("plots/insect_stackedbar_site_size.png", stackedbar_site_size,
          base_height = 11.7, base_width = 9)

######################:
# ABUNDANCE BY ORDER #:
######################:

order_plot

save_plot("plots/insect_order_plot.png", order_plot,
          base_width = 10, base_height = 9)
