  # set the working directory to the r_scripts folder
setwd("~/Desktop/ENVI408/r_scripts")

#                                                                      PACKAGES :::::::::####
# if you run into problems anywhere this is going to be it

# make vector of packages
packages_angophora <- c("magrittr", "dplyr", "ggplot2", "cowplot")

# install them (only need to run once)
#lapply(packages_angophora, install.packages)

  # if there's an error message saying "there is no package called lifecycle"
# and you're on windows, install Rtools from:
# https://cran.r-project.org/bin/windows/Rtools/
# and then try and rerun it again.
# per https://stackoverflow.com/questions/70185951/error-using-librarydplyr-in-rstudio-there-is-no-package-called-lifecycle

# load packages
lapply(packages_angophora, library, character.only = T)
#
#                                                                                        SETUP ::::::####
# importing + data manip + creating phyloseq objects                                     ####
# ideally should be able to just run lines 19-363 and it should all work                 ####
######################################################################################## ####
#                                                                   Import data :::::::::####
raw_eucdata <- read.csv("data/eucalypt_data.csv")
rainfall_data <- read.csv("data/rainfall_data.csv")
#                                                                   Format data :::::::::####

# create df of epicormic growth by angophora/non-angophora
epic_df <- raw_eucdata %>% 
  
  # filter for living eucs
  
  filter(alive == T) %>% 
  
  # lump non-angophora eucs into one category
  
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ "Angophora floribunda",
    .default = "Non-Angophora")) %>% 
  
  # get rid of nas etc. and just make epicormic growth 0-3
  
  filter(epicormic_growth == 0 | epicormic_growth == 1 | 
           epicormic_growth == 2 | epicormic_growth == 3)

# make into a contingency table
epic_table <- table(epic_df$binary_species, epic_df$epicormic_growth)

# make contingency table within hargraves ave
epic_df_hga <- epic_df %>% 
  filter(site_id == "HGA")

epic_table_hga <- table(epic_df_hga$binary_species, epic_df_hga$epicormic_growth)

# make contingency table outside of hargraves ave
epic_df_no_hga <- epic_df %>% 
  filter(site_id != "HGA")

epic_table_no_hga <- table(epic_df_no_hga$binary_species, epic_df_no_hga$epicormic_growth)

# # #

# create df of epicormic growth by angophora/non angophora
die_df <- raw_eucdata %>% 
  
  # filter for living eucs
  
  filter(alive == T) %>% 
  
  # lump non-angophora eucs into one category
  
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ "Angophora floribunda",
    .default = "Non-Angophora")) %>% 
  
  # get rid of nas and make dieoff just 0-3
  
  filter(leaf_dieoff == 0 | leaf_dieoff == 1 | 
           leaf_dieoff == 2 | leaf_dieoff == 3)

# make into a contingency table
die_table <- table(die_df$binary_species, die_df$leaf_dieoff)

# make contingency table within hargraves ave
die_df_hga <- die_df %>% 
  filter(site_id == "HGA")

die_table_hga <- table(die_df_hga$binary_species, die_df_hga$leaf_dieoff)

# make contingency table outside of hargraves ave
die_df_no_hga <- die_df %>% 
  filter(site_id != "HGA")

die_table_no_hga <- table(die_df_no_hga$binary_species, die_df_no_hga$leaf_dieoff)

# rainfall data

# filter to be since 2000 (and not 2025 since 2025 isn't over), 
# and group by year to get annual rainfall
rainfall_data <- rainfall_data %>%  
  filter(Year > 1999) %>% 
  filter(Year < 2025) %>% 
  group_by(Year) %>% 
  summarise(annual_precip = sum(Monthly.Precipitation.Total..millimetres.))

# get mean annual precipitiation from 2000-2024
mean_precip <- mean(rainfall_data$annual_precip)

#
#                                                       Create dfs for plotting :::::::::####

###########:
# OVERALL #:
###########:
# epicormic growth

epic_plotting_overall <- raw_eucdata %>% 
  
  # filter for alive trees
  
  filter(alive == T) %>% 
  
  # turn into binary species. weird formatting is so it looks good on the plot
  
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ 
      "Angophora 
floribunda",
    .default = "Non-Angophora")) %>% 
  
  # aggregate by species + epicormic growth
  group_by(binary_species, epicormic_growth) %>% 
  summarise(n = n()) %>% 
  
  # filter out nas et al
  filter(epicormic_growth == 0 | epicormic_growth == 1 | 
           epicormic_growth == 2 | epicormic_growth == 3)

# same same but for dieoff
die_plotting_overall <- raw_eucdata %>% 
  filter(alive == T) %>% 
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ 
      "Angophora 
floribunda",
    .default = "Non-Angophora")) %>% 
  group_by(binary_species, leaf_dieoff) %>% 
  summarise(n = n()) %>%
  filter(leaf_dieoff == 0 | leaf_dieoff == 1 | 
           leaf_dieoff == 2 | leaf_dieoff == 3)

#############:
# HARGRAVES #:
#############:
# epicormic growth

epic_plotting_hga <- raw_eucdata %>% 
  
  # filter for alive trees + only for hga
  
  filter(alive == T, site_id == "HGA") %>% 
  
  # turn into binary species. weird formatting is so it looks good on the plot
  
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ 
      "Angophora 
floribunda",
    .default = "Non-Angophora")) %>% 
  
  # aggregate by species + epicormic growth
  group_by(binary_species, epicormic_growth) %>% 
  summarise(n = n()) %>% 
  
  # filter out nas et al
  filter(epicormic_growth == 0 | epicormic_growth == 1 | 
           epicormic_growth == 2 | epicormic_growth == 3)

# same same but for dieoff
die_plotting_hga <- raw_eucdata %>% 
  filter(alive == T,
         site_id == "HGA") %>% 
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ 
      "Angophora 
floribunda",
    .default = "Non-Angophora")) %>% 
  group_by(binary_species, leaf_dieoff) %>% 
  summarise(n = n()) %>%
  filter(leaf_dieoff == 0 | leaf_dieoff == 1 | 
           leaf_dieoff == 2 | leaf_dieoff == 3)

#################:
# NON-HARGRAVES #:
#################:
# epicormic growth

epic_plotting_no_hga <- raw_eucdata %>% 
  
  # filter for alive trees + only for hga
  
  filter(alive == T, site_id != "HGA") %>% 
  
  # turn into binary species. weird formatting is so it looks good on the plot
  
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ 
      "Angophora 
floribunda",
    .default = "Non-Angophora")) %>% 
  
  # aggregate by species + epicormic growth
  group_by(binary_species, epicormic_growth) %>% 
  summarise(n = n()) %>% 
  
  # filter out nas et al
  filter(epicormic_growth == 0 | epicormic_growth == 1 | 
           epicormic_growth == 2 | epicormic_growth == 3)

# same same but for dieoff
die_plotting_no_hga <- raw_eucdata %>% 
  filter(alive == T, site_id != "HGA") %>% 
  mutate(binary_species = case_when(
    (unconfirmed_species == "Angophora floribunda" ) ~ 
      "Angophora 
floribunda",
    .default = "Non-Angophora")) %>% 
  group_by(binary_species, leaf_dieoff) %>% 
  summarise(n = n()) %>%
  filter(leaf_dieoff == 0 | leaf_dieoff == 1 | 
           leaf_dieoff == 2 | leaf_dieoff == 3)

#

#                                                                                     ANALYSIS ::::::####
######################################################################################## ####
#                                                    Dieoff + epicormic overall :::::::::####

# display the contingency tables 
epic_table
die_table

# check that <20% of expected values <5
chisq.test(epic_table)$expected
chisq.test(die_table)$expected

# run chi square tests
chisq.test(epic_table)
chisq.test(die_table)

#                                               Dieoff + epicormic at Hargraves :::::::::####

# display the contingency tables 
epic_table_hga
die_table_hga

# check that <20% of expected values <5
chisq.test(epic_table_hga)$expected
chisq.test(die_table_hga)$expected

# <20% of expected values <5 for dieoff, so run fisher test instead
# run chi square tests
chisq.test(epic_table_hga)
fisher.test(die_table_hga)

#                                       Dieoff + epicormic outside of Hargraves :::::::::####

# display the contingency tables 
epic_table_no_hga
die_table_no_hga

# check that <20% of expected values <5
chisq.test(epic_table_no_hga)$expected
chisq.test(die_table_no_hga)$expected

# run chi square tests
chisq.test(epic_table_no_hga)
chisq.test(die_table_no_hga)

#
#                                                                                     PLOTTING ::::::####
######################################################################################## ####
#                                                                       Overall :::::::::####
# make bar plots of freq for epic_plotting_overall and dieoff
epic_plot_ang <- ggplot(subset(epic_plotting_overall, epic_plotting_overall$binary_species == 
                               "Angophora 
floribunda"), 
                      aes(x = epicormic_growth, y = n)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "#79c638") + 
  xlab("Epicormic growth") + ylab("Frequency") +
  ggtitle("Angophora") +
  theme_cowplot()


die_plot_ang <- ggplot(subset(die_plotting_overall, die_plotting_overall$binary_species == 
                              "Angophora 
floribunda"), 
                     aes(x = leaf_dieoff, y = n, fill = leaf_dieoff)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "grey") + 
  scale_fill_brewer(palette = "Greys", name = "Dieoff") +
  xlab("Dieoff") + ylab("Frequency") +
  ggtitle("Angophora") +
  theme_cowplot()


epic_plot_nang <- ggplot(subset(epic_plotting_overall, epic_plotting_overall$binary_species == "Non-Angophora"), 
                      aes(x = epicormic_growth, y = n)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "#79c638") + 
  xlab("Epicormic growth") + ylab("Frequency") +
  ggtitle("Non-Angophora") +
  theme_cowplot()

die_plot_nang <- ggplot(subset(die_plotting_overall, die_plotting_overall$binary_species == "Non-Angophora"), 
                     aes(x = leaf_dieoff, y = n, fill = leaf_dieoff)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "grey") + 
  xlab("Dieoff") + ylab("Relative frequency (%)") +
  ggtitle("Non-Angophora") + 
  theme_cowplot()


overall_bars <- plot_grid(epic_plot_ang, epic_plot_nang, die_plot_ang, die_plot_nang)

# make bar plots of relative freq. for dieoff
epic_plot_stackedbar <- ggplot(epic_plotting_overall, 
                               aes(x = binary_species, 
                                   y = n, 
                                   fill = epicormic_growth)) +
  geom_bar(stat = "identity", 
           position = "fill",
           colour = "black", 
           linewidth = 0.25) + 
  scale_fill_brewer(palette = "Greens", name = 
                      "Epicormic 
growth") +
  labs(x = "Species", y = "Relative frequency") +
  theme_cowplot()

die_plot_stackedbar <- ggplot(die_plotting_overall, 
                              aes(x = binary_species, 
                                  y = n, 
                                  fill = leaf_dieoff)) +
  geom_bar(stat = "identity", 
           position = "fill",
           colour = "black", 
           linewidth = 0.25) + 
  scale_fill_brewer(palette = "Greys", name = "Dieoff") +
  labs(x = "Species", y = "Relative frequency") +
  theme_cowplot()

overall_stackedbar <- plot_grid(epic_plot_stackedbar, die_plot_stackedbar)
overall_stackedbar
overall_angophora <- plot_grid(overall_bars, overall_stackedbar, 
                               labels = c("a)", "b)"))

# plot all together


#

#                                                                Hargraves only :::::::::####
# make bar plots of freq for epic_plotting_hga and dieoff
epic_plot_ang_hga <- ggplot(subset(epic_plotting_hga, epic_plotting_hga$binary_species == 
                                 "Angophora 
floribunda"), 
                        aes(x = epicormic_growth, y = n)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "#79c638") + 
  xlab("Epicormic growth") + ylab("Frequency") +
  ggtitle("Angophora") +
  theme_cowplot()


die_plot_ang_hga <- ggplot(subset(die_plotting_hga, die_plotting_hga$binary_species == 
                                "Angophora 
floribunda"), 
                       aes(x = leaf_dieoff, y = n, fill = leaf_dieoff)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "grey") + 
  scale_fill_brewer(palette = "Greys", name = "Dieoff") +
  xlab("Dieoff") + ylab("Frequency") +
  ggtitle("Angophora") +
  theme_cowplot()


epic_plot_nang_hga <- ggplot(subset(epic_plotting_hga, epic_plotting_hga$binary_species == "Non-Angophora"), 
                         aes(x = epicormic_growth, y = n)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "#79c638") + 
  xlab("Epicormic growth") + ylab("Frequency") +
  ggtitle("Non-Angophora") +
  theme_cowplot()

die_plot_nang_hga <- ggplot(subset(die_plotting_hga, die_plotting_hga$binary_species == "Non-Angophora"), 
                        aes(x = leaf_dieoff, y = n, fill = leaf_dieoff)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "grey") + 
  xlab("Dieoff") + ylab("Relative frequency (%)") +
  ggtitle("Non-Angophora") + 
  theme_cowplot()


hga_bars <- plot_grid(epic_plot_ang_hga, epic_plot_nang_hga, 
                      die_plot_ang_hga, die_plot_nang_hga)

# make bar plots of relative freq. for dieoff
epic_plot_stackedbar_hga <- ggplot(epic_plotting_hga, 
                               aes(x = binary_species, 
                                   y = n, 
                                   fill = epicormic_growth)) +
  geom_bar(stat = "identity", 
           position = "fill",
           colour = "black", 
           linewidth = 0.25) + 
  scale_fill_brewer(palette = "Greens", name = 
                      "Epicormic 
growth") +
  labs(x = "Species", y = "Relative frequency") +
  theme_cowplot()

die_plot_stackedbar_hga <- ggplot(die_plotting_hga, 
                              aes(x = binary_species, 
                                  y = n, 
                                  fill = leaf_dieoff)) +
  geom_bar(stat = "identity", 
           position = "fill",
           colour = "black", 
           linewidth = 0.25) + 
  scale_fill_brewer(palette = "Greys", name = "Dieoff") +
  labs(x = "Species", y = "Relative frequency") +
  theme_cowplot()

hga_stackedbar <- plot_grid(epic_plot_stackedbar_hga, die_plot_stackedbar_hga)

# plot all together
hga_angophora <- plot_grid(hga_bars, hga_stackedbar, 
                           labels = c("a)", "b)"))

# 

#                                                            Non-Hargraves only :::::::::####
# make bar plots of freq for epic_plotting_no_hga and dieoff
epic_plot_ang_no_hga <- ggplot(subset(epic_plotting_no_hga, epic_plotting_no_hga$binary_species == 
                                        "Angophora 
floribunda"), 
                               aes(x = epicormic_growth, y = n)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "#79c638") + 
  xlab("Epicormic growth") + ylab("Frequency") +
  ggtitle("Angophora") +
  theme_cowplot()


die_plot_ang_no_hga <- ggplot(subset(die_plotting_no_hga, die_plotting_no_hga$binary_species == 
                                       "Angophora 
floribunda"), 
                              aes(x = leaf_dieoff, y = n, fill = leaf_dieoff)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "grey") + 
  scale_fill_brewer(palette = "Greys", name = "Dieoff") +
  xlab("Dieoff") + ylab("Frequency") +
  ggtitle("Angophora") +
  theme_cowplot()


epic_plot_nang_no_hga <- ggplot(subset(epic_plotting_no_hga, epic_plotting_no_hga$binary_species == "Non-Angophora"), 
                                aes(x = epicormic_growth, y = n)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "#79c638") + 
  xlab("Epicormic growth") + ylab("Frequency") +
  ggtitle("Non-Angophora") +
  theme_cowplot()

die_plot_nang_no_hga <- ggplot(subset(die_plotting_no_hga, die_plotting_no_hga$binary_species == "Non-Angophora"), 
                               aes(x = leaf_dieoff, y = n, fill = leaf_dieoff)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25, fill = "grey") + 
  xlab("Dieoff") + ylab("Relative frequency (%)") +
  ggtitle("Non-Angophora") + 
  theme_cowplot()


no_hga_bars <- plot_grid(epic_plot_ang_no_hga, epic_plot_nang_no_hga, 
                         die_plot_ang_no_hga, die_plot_nang_no_hga)

# make bar plots of relative freq. for dieoff
epic_plot_stackedbar_no_hga <- ggplot(epic_plotting_no_hga, 
                                      aes(x = binary_species, 
                                          y = n, 
                                          fill = epicormic_growth)) +
  geom_bar(stat = "identity", 
           position = "fill",
           colour = "black", 
           linewidth = 0.25) + 
  scale_fill_brewer(palette = "Greens", name = 
                      "Epicormic 
growth") +
  labs(x = "Species", y = "Relative frequency") +
  theme_cowplot()

die_plot_stackedbar_no_hga <- ggplot(die_plotting_no_hga, 
                                     aes(x = binary_species, 
                                         y = n, 
                                         fill = leaf_dieoff)) +
  geom_bar(stat = "identity", 
           position = "fill",
           colour = "black", 
           linewidth = 0.25) + 
  scale_fill_brewer(palette = "Greys", name = "Dieoff") +
  labs(x = "Species", y = "Relative frequency") +
  theme_cowplot()

no_hga_stackedbar <- plot_grid(epic_plot_stackedbar_no_hga, die_plot_stackedbar_no_hga)

# plot all together
non_hga_angophora <- plot_grid(no_hga_bars, no_hga_stackedbar, 
                               labels = c("a)", "b)"))

# 
#                                                                      Rainfall :::::::::####

rainfall_plot <- ggplot(rainfall_data, aes(x = Year, 
                                           y = annual_precip)) +
  geom_line(colour = "blue") + 
  lims(x = c(2000,2024), y = c(0,2100)) +
  theme_bw() +
  geom_abline(slope = 0, intercept = mean_precip, colour = "red") +
  labs(x = "Year", y = "Annual precipitation")

#
#                                                          Plots used in thesis :::::::::####
# again. basically just saving plots but it's called this for consistency

###########:
# OVERALL #:
###########:
overall_angophora

save_plot("plots/angophora_dieoff_epicormic_overall.png", 
          overall_angophora,
          base_width = 18,
          base_height = 8)
#############:
# HARGRAVES #:
#############:
hga_angophora

save_plot("plots/angophora_dieoff_epicormic_hargraves.png", 
          hga_angophora,
          base_width = 18,
          base_height = 8)

#################:
# NON-HARGRAVES #:
#################:
non_hga_angophora


save_plot("plots/angophora_dieoff_epicormic_non_hargraves.png", 
          non_hga_angophora,
          base_width = 18,
          base_height = 8)

############:
# RAINFALL #:
############:

rainfall_plot

save_plot("plots/rainfall_plot.png",
          rainfall_plot,
          base_width = 8,
          base_height = 5.2)
