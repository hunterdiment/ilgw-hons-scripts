# set the working directory to the r_scripts folder
setwd("~/Desktop/ENVI408/r_scripts")
#

#                                                                      PACKAGES :::::::::####
# make vector of packages
packages_dbh <- c("magrittr", "tidyr", "dplyr", "ggplot2", # tidyverse packages 
              "cowplot",                         # used to make plots good
              "rstatix",                         # used for dunn_test
              "kSamples")                       # used for ad.test

# install them (only need to run once)
#lapply(packages_dbh, install.packages)

# load them
lapply(packages_dbh, library, character.only = T)

  # if there's an error message saying "there is no package called lifecycle"
# and you're on windows, install Rtools from:
# https://cran.r-project.org/bin/windows/Rtools/
# and then try and rerun it again.
# per https://stackoverflow.com/questions/70185951/error-using-librarydplyr-in-rstudio-there-is-no-package-called-lifecycle

# set seed
set.seed(2018050519)

#
#                                                                                        SETUP ::::::####
# basically a bunch of data, creation + manipulation to prep for analysis + plotting     ####
# ideally should be able to just run lines 17-184 and it should all work                 ####
######################################################################################## ####
#                                                                   Import data :::::::::####

raw_eucdata <- read.csv("data/eucalypt_data.csv")
        ####
#                                                                   Format data :::::::::####

# create a whole load of new cols
eucdata <- raw_eucdata %>%
  
  # select only the columns of interest
  
  dplyr::select(site_id, 
                species, 
                dbh_m, 
                alive, forked_at_ground_level, forked_below_ground_level) %>% 
  
  # create species column where redgum + stringybark species are grouped together
  
  mutate(
    species_complex = 
      case_when(
        (species == "Eucalyptus amplifolia") ~ "Redgum complex",
        (species == "Eucalyptus tereticornis") ~ "Redgum complex",
        (species ==  "Eucalyptus amplifolia or tereticornis") ~ "Redgum complex",
        (species == "Eucalyptus eugenioides") ~ "Stringybark complex",
        (species == "Eucalyptus eugenioides") ~ "Stringybark complex",
        (species == "Eucalyptus globoidea") ~ "Stringybark complex",
        (species ==  "Eucalyptus eugenioides or globoidea") ~ "Stringybark complex",
        .default = species)
    ) %>% 
  
  # column of site size
  
  mutate(size = case_when(
    (site_id == "ALR" | site_id == "BBR" | site_id == "CRM") ~ "Large",
    (site_id == "HGA" | site_id == "RGR" | site_id == "WMP") ~ "Medium",
    (site_id == "CSP" | site_id == "MCP" | site_id == "PAP") ~ "Small",
    .default = NA)) %>% 
  
  # columns for ordering by site when plots come around
  
  mutate(order = case_when(
    (site_id == "ALR" | site_id == "BBR" | site_id == "CRM") ~ 3,
    (site_id == "HGA" | site_id == "RGR" | site_id == "WMP") ~ 2,
    (site_id == "CSP" | site_id == "MCP" | site_id == "PAP") ~ 1,
    .default = NA)) %>% 
  
  mutate(invorder = case_when(
    (site_id == "ALR" | site_id == "BBR" | site_id == "CRM") ~ 1,
    (site_id == "HGA" | site_id == "RGR" | site_id == "WMP") ~ 2,
    (site_id == "CSP" | site_id == "MCP" | site_id == "PAP") ~ 3,
    .default = NA)) 

# create df of eucdata without forks
forkless_eucdata <- filter(eucdata, forked_at_ground_level == F, 
                           forked_below_ground_level == F,
                           alive == T) %>% 
  
  # create column for stacked bar plot
  
  mutate(cat_dbh_m = cut(dbh_m,
                         breaks = c(-Inf, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
                                    0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, Inf),
                         labels = c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2",
                                    "0.2-0.25", "0.25-0.3", "0.3-0.35", "0.35-0.4",
                                    "0.4-0.45", "0.45-0.5", "0.5-0.55", "0.55-0.6",
                                    "0.6-0.65", "0.65-0.7", "0.7-0.75", "0.75-0.8",
                                    "0.8+")),
         
         # create column for colouring histogram (purely aesthetic)
        
         hist_fill = cut(dbh_m, 
           breaks = c(-0.025, 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375,
                      0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825,
                      0.875, 0.925, 0.975, 1.025, 1.075, 1.125)))

#
#                                           DBH vectors for analysis + plotting :::::::::####
  
  # make vectors of dbh for each site
# charles stimson
CSP_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "CSP")
# mcphail reserve
MCP_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "MCP")
# phil adams park
PAP_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "PAP")

# hargraves ave
HGA_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "HGA")
# rosella grove
RGR_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "RGR")
# wiseman park
WMP_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "WMP")

# albion park light rail
ALR_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "ALR")
# blackbutt reserve
BBR_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "BBR")
# croom reserve
CRM_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$site_id == "CRM")

# create a list of site vectors
site_dbh_list <- list(CSP_dbh, MCP_dbh, PAP_dbh, 
                      HGA_dbh, RGR_dbh, WMP_dbh, 
                      ALR_dbh, BBR_dbh, CRM_dbh)

#                                                             Subsampling data  :::::::::####

#########:
# SMALL #:
#########:

# generate a sample from each site w 21 (min(n) for small) observations each 
small_site_sample <- data.frame(CSP = sample(CSP_dbh, 21, replace = F),
                                MCP = sample(MCP_dbh, 21, replace = F),
                                PAP = sample(PAP_dbh, 21, replace = F))%>% 
  pivot_longer(names_to = "site_id", cols = CSP:PAP) %>%
  mutate(size = "Small", order = 1,)

##########:
# MEDIUM #:
##########:

# generate a sample from each site w 57 (min(n) for med) observations each 
med_site_sample <- data.frame(HGA = sample(HGA_dbh, 57, replace = F),
                              RGR = sample(RGR_dbh, 57, replace = F),
                              WMP = sample(WMP_dbh, 57, replace = F))%>% 
  pivot_longer(names_to = "site_id", cols = HGA:WMP) %>% 
  mutate(size = "Medium", order = 2)

#########:
# LARGE #:
#########:

# generate a sample from each site w 228 (min(n) for large) observations each  
large_site_sample <- data.frame(ALR = sample(ALR_dbh, 228, replace = F),
                                BBR = sample(BBR_dbh, 228, replace = F),
                                CRM = sample(CRM_dbh, 228, replace = F))%>% 
  pivot_longer(names_to = "site_id", cols = ALR:CRM) %>% 
  mutate(size = "Large", order = 3)

#######################################:
# COMBINE SUBSAMPLED DATA INTO ONE DF #:
#######################################:
subsampled_data <- bind_rows(small_site_sample, med_site_sample, large_site_sample) %>% 
  mutate(cat_dbh_m = cut(value,
                         breaks = c(-Inf, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
                                    0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, Inf),
                         labels = c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2",
                                    "0.2-0.25", "0.25-0.3", "0.3-0.35", "0.35-0.4",
                                    "0.4-0.45", "0.45-0.5", "0.5-0.55", "0.55-0.6",
                                    "0.6-0.65", "0.65-0.7", "0.7-0.75", "0.75-0.8",
                                    "0.8+")),
         hist_fill = cut(value, 
                         breaks = c(-0.025, 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375,
                                    0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825,
                                    0.875, 0.925, 0.975, 1.025, 1.075, 1.125)))


#
#                                                                            STATISTICAL TESTS ::::::####
######################################################################################## ####
#                                                                   DBH by site :::::::::####

# testing for differences in DBH attributes across sites in this section

########################:
# KRUSKAL TEST BY SITE #:
########################:

  # testing for differences in median DBH across sites
# test
kruskal.test(dbh_m ~ site_id, data = forkless_eucdata)

  # testing for differences in median DBH between site pairs
# dunn's test for post-hoc
dunn_test(dbh_m ~ site_id, data = forkless_eucdata, p.adjust.method = "bonferroni") %>% 
            print(n = 36)

#################################:
# ANDERSON-DARLING TEST BY SITE #:
#################################:

# test dbh distribution among sites
ad.test(site_dbh_list, method = "exact", dist = FALSE, Nsim = 1000)

#                             post-hoc test for distribution by site  :::::::::####

  # this is its own section for tidiness's sake
  # basically 36 pairwise comparisons with bonferroni adjustment
  # specificaly bonferroni because it's quite conservative so as to filter  
  # out less significant differences

# post-hoc pairwise comparison for adtest
# need to do bonferroni adjustment
# 36 hypotheses
# 0.05/36 = 0.013889
# alpha level is 0.013889

#### charles stimson park v everywhere else                                     CSP

ad.test(CSP_dbh, MCP_dbh, method = "exact", dist = FALSE, Nsim = 1000)
ad.test(CSP_dbh, PAP_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

ad.test(CSP_dbh, HGA_dbh, method = "exact", dist = FALSE, Nsim = 1000)
ad.test(CSP_dbh, RGR_dbh, method = "exact", dist = FALSE, Nsim = 1000)
ad.test(CSP_dbh, WMP_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

ad.test(CSP_dbh, ALR_dbh, method = "exact", dist = FALSE, Nsim = 1000)
ad.test(CSP_dbh, BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(CSP_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)

#### mcphail reserve                                                            MCP

ad.test(MCP_dbh, PAP_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

ad.test(MCP_dbh, HGA_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # *
ad.test(MCP_dbh, RGR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(MCP_dbh, WMP_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

ad.test(MCP_dbh, ALR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # **
ad.test(MCP_dbh, BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(MCP_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

#### phil adams park                                                            PAP

ad.test(PAP_dbh, HGA_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(PAP_dbh, RGR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(PAP_dbh, WMP_dbh, method = "exact", dist = FALSE, Nsim = 1000)

ad.test(PAP_dbh, ALR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(PAP_dbh, BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(PAP_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

#### hargraves ave                                                              HGA

ad.test(HGA_dbh, RGR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # *
ad.test(HGA_dbh, WMP_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

ad.test(HGA_dbh, ALR_dbh, method = "exact", dist = FALSE, Nsim = 1000)
ad.test(HGA_dbh, BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(HGA_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # **

#### rosella grove                                                              RGR

ad.test(RGR_dbh, WMP_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

ad.test(RGR_dbh, ALR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # **
ad.test(RGR_dbh, BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(RGR_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # **

#### wiseman park                                                               WMP
ad.test(WMP_dbh, ALR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(WMP_dbh, BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***
ad.test(WMP_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***

#### albion park light rail                                                     ALR
ad.test(ALR_dbh, BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # *
ad.test(ALR_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)

#### blackbutt reserve                                                          BBR
ad.test(BBR_dbh, CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)          # ***


#               test that subsamples are representative of full data  :::::::::####

# use ad.test to test whether subsamples + full data come from same distribution

# small

# csp
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "CSP"), 
        CSP_dbh, method = "exact", dist = FALSE, Nsim = 1000)
# mcp
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "MCP"), 
        MCP_dbh, method = "exact", dist = FALSE, Nsim = 1000)
# pap
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "PAP"), 
        PAP_dbh, method = "exact", dist = FALSE, Nsim = 1000)

# med
# hga
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "HGA"), 
        HGA_dbh, method = "exact", dist = FALSE, Nsim = 1000)
# rgr
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "RGR"), 
        RGR_dbh, method = "exact", dist = FALSE, Nsim = 1000)
# wmp
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "WMP"), 
        WMP_dbh, method = "exact", dist = FALSE, Nsim = 1000)

# large
# alr
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "ALR"), 
        ALR_dbh, method = "exact", dist = FALSE, Nsim = 1000)
# bbr
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "BBR"), 
        BBR_dbh, method = "exact", dist = FALSE, Nsim = 1000)
# crm
ad.test(subset(subsampled_data$value, subsampled_data$site_id == "CRM"), 
        CRM_dbh, method = "exact", dist = FALSE, Nsim = 1000)

#
#                                                                   DBH by size :::::::::####

# here we're testing using the subsampled data, so that each site contributes
# equally to the test (otherwise the sites with the biggest
# sample size would dominate)

########################:
# KRUSKAL TEST BY SIZE #:
########################:

# testing for differences in median dbh between patch size categories

# test
kruskal.test(value ~ size, data = subsampled_data)

# dunn test for post-hoc
dunn_test(value ~ size, data = subsampled_data, p.adjust.method = "bonferroni")

# RESULT:
# median eucalypt DBH differs among site sizes
# median eucalypt DBH differs between small and large sites
# median eucalypt DBH differs between medium and large sites

##############################################:
# ANDERSON-DARLING TEST BY SIZE: SAMPLE DATA #:
##############################################:

# testing for differences in dbh distribution between patch size categories

# create list with small, medium, large

size_dbh_list_subsampled <- list(small_site_sample$value, med_site_sample$value,
                                 large_site_sample$value)

# run anderson-darling test for distr. among sites
ad.test(size_dbh_list_subsampled, method = "exact", dist = F, Nsim = 1000)

# run pairwise tests 
# alpha level is 0.01667
  # using 50000 here because at 1000 the p-values keep fluctuating around 0.016

# small vs medium
ad.test(small_site_sample$value, med_site_sample$value, method = "exact",
        dist = FALSE, Nsim = 50000)

# small vs large
ad.test(small_site_sample$value, large_site_sample$value, method = "exact",
        dist = FALSE, Nsim = 1000)

# medium vs large
ad.test(med_site_sample$value, large_site_sample$value, method = "exact",
        dist = FALSE, Nsim = 1000)

############################################:
# ANDERSON-DARLING TEST BY SIZE: REAL DATA #:
############################################:

# NOT USED IN PAPER

# make vectors of dbh for small, medium, large sites
small_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$size == "Small")
med_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$size == "Medium")
large_dbh <- subset(forkless_eucdata$dbh_m, forkless_eucdata$size == "Large")
size_dbh_list <- list(small_dbh, med_dbh, large_dbh)

# test dbh distribution among sites
ad.test(size_dbh_list, method = "exact", dist = FALSE, Nsim = 1000)

# post-hoc pairwise comparison for adtest
# need to do bonferroni adjustment
# 3 hypotheses (small vs med, small vs large, med vs large)
# 0.05/3 = 0.01667
# alpha level is 0.01667
ad.test(small_dbh, med_dbh, method = "exact", dist = FALSE, Nsim = 1000)
ad.test(small_dbh, large_dbh, method = "exact", dist = FALSE, Nsim = 1000)
ad.test(med_dbh, large_dbh, method = "exact", dist = FALSE, Nsim = 1000)

# RESULT:
# eucalypt DBH does not follow the same distribution at different site sizes
# eucalypt DBH does not follow the same distribution between
# small and medium,
# small and large,
# or medium and large sites


#                                                                                    PLOTTING  ::::::####
######################################################################################## ####
#                                                                   DBH by site :::::::::####
################:
# FULL DATASET #:
################:

# boxplots of dbh by site
boxplot_dbh_site <- ggplot(forkless_eucdata, aes(x = dbh_m, 
                                                 y = reorder(site_id, order))) +
  geom_boxplot(fill = "#94ef5d") + theme_cowplot() +
  labs(x = "Diameter at breast height (m)", 
       y = "Site",
       title = "(a)") +
  scale_y_discrete(limits=rev) +
  annotate("text", y = c("CSP", "MCP", "PAP", 
                         "HGA", "RGR", "WMP", 
                         "ALR", "BBR", "CRM"),
           x = 1.2,
           label = c("n=105","n=45","n=21",
                     "n=141","n=250","n=57",
                     "n=228","n=441","n=290")) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))

# hists of dbh by site
  # palette for bars
histpal <- scales::seq_gradient_pal("#94ef5d", "#2b100c", "Lab")(seq(0,1,length.out=23))

hist_dbh_site <- ggplot(forkless_eucdata, aes(x = dbh_m, fill = hist_fill)) +
  geom_histogram(binwidth = 0.05, colour = "black", width = 0.25) + 
  facet_wrap(~reorder(site_id, order), ncol = 1, scale = "free_y") + 
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values = histpal) +
  labs(x = "Diameter at breast height (m)", 
       y = "Site",
       title = "(b)") +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))

  # plot them together (do not sell seperately)
full_site <- plot_grid(boxplot_dbh_site, hist_dbh_site)
full_site

#####################:
# SUBSAMPLED DATSET #:
#####################:

# boxplots of dbh by site
boxplot_dbh_site_samp <- ggplot(subsampled_data, aes(x = value, 
                                                 y = reorder(site_id, order))) +
  geom_boxplot(fill = "#94ef5d") + theme_cowplot() +
  labs(x = "Diameter at breast height (m)", 
       y = "Site",
       title = "(a)") +
  scale_y_discrete(limits=rev) +
  annotate("text", 
           y = c("CSP", "MCP", "PAP", 
                 "HGA", "RGR", "WMP", 
                 "ALR", "BBR", "CRM"),
           x = 1.2,
           label = c("n=21","n=21","n=21",
                     "n=57","n=57","n=57",
                     "n=228","n=228","n=228")) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))

  
# hists of dbh by site
hist_dbh_site_samp <- ggplot(subsampled_data, aes(x = value, fill = hist_fill)) +
  geom_histogram(binwidth = 0.05, colour = "black", width = 0.25) + 
  facet_wrap(~reorder(site_id, order), ncol = 1, scale = "free_y") + 
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  labs(x = "Diameter at breast height (m)", 
       y = "Site",
       title = "(b)") +
  scale_fill_manual(values = histpal) +
  xlab("Diameter at breast height (m)") + ylab("Count") +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))


# plot them together (do not sell seperately)
subsampled_site <- plot_grid(boxplot_dbh_site_samp, hist_dbh_site_samp)
subsampled_site

###############################:
# FULL VS SUBSAMPLED DATASETS #:
###############################:

# plot them togetherer
full_subsampled_site <- plot_grid(boxplot_dbh_site, hist_dbh_site, 
                                  boxplot_dbh_site_samp, hist_dbh_site_samp,
                                  nrow = 1, ncol = 4, 
                                  labels = c("Full dataset", "", "Subsampled", ""),
                                  label_x = 0.4, label_y = 0.996)
full_subsampled_site

#
#                                                                   DBH by size :::::::::####

######################:
# SUBSAMPLED DATASET #:
######################:

boxplot_dbh_size_samp <- ggplot(subsampled_data, aes(x = value, 
                                                 y = reorder(size, order))) +
  geom_boxplot(fill = "#94ef5d") + 
  theme_cowplot() +
  labs(x = "Diameter at breast height (m)", 
       y = "Site",
       title = "(a)") +
  scale_y_discrete(limits=rev) +
  annotate("text", y = c("Small", "Medium", "Large"),
           x = 1.2,
           label = c("n=63", "n=171", "n=684")) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))


# histogram of dbh by site size
hist_dbh_size_samp <- ggplot(subsampled_data, aes(x = value, fill = hist_fill)) +
  geom_histogram(binwidth = 0.05, colour = "black", width = 0.25) + 
  facet_wrap(~reorder(size, order), ncol = 1, scale = "free_y") + 
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values = histpal) +
  labs(x = "Diameter at breast height (m)", 
       y = "Site",
       title = "(b)") +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))

dbh_size_subsampled <- plot_grid(boxplot_dbh_size_samp, hist_dbh_size_samp)
dbh_size_subsampled

################:
# FULL DATASET #:
################:
  # NOT USED IN THESIS
boxplot_dbh_size <- ggplot(forkless_eucdata, aes(x = dbh_m, 
                                                 y = reorder(size, order))) +
  geom_boxplot(fill = "#94ef5d") + theme_cowplot() +
  xlab("Diameter at breast height (m)") + ylab("Site size") +
  scale_y_discrete(limits=rev)  +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))

# histogram of dbh by site size
hist_dbh_size <- ggplot(forkless_eucdata, aes(x = dbh_m, fill = hist_fill)) +
  geom_histogram(binwidth = 0.05, colour = "black", width = 0.25) + 
  facet_wrap(~reorder(size, order), ncol = 1, scale = "free_y") + 
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values = histpal) +
  xlab("Diameter at breast height (m)") + ylab("Count") +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20))

plot_grid(boxplot_dbh_size, hist_dbh_size)

#############################:
# SUBSAMPLED + FULL DATASET #:
#############################:
  # NOT USED IN THESIS
plot_grid(boxplot_dbh_size, hist_dbh_size, 
          boxplot_dbh_size_samp, hist_dbh_size_samp,
          nrow = 1, ncol = 4, labels = c("Full dataset", "", 
                                         "Subsampled", ""),
          label_x = 0.13)

sampsize_table <- forkless_eucdata %>% 
  group_by(site_id, size) %>% 
  summarise(n = n())
sampsize_table

######################:
# STACKED BAR OF DBH #:
######################:

  # this isn't in the thesis but i used it in my seminar because i think it
# gives a great look at the data at a very quick glance
# plus it looks nice

dbh_stackedbar_data <- subsampled_data %>% 
  group_by(size, cat_dbh_m) %>% 
  summarise(n = n()) %>% 
  mutate(order = case_when((size == "Small") ~ 1,
                           (size == "Medium") ~ 2,
                           (size == "Large") ~ 3))

barpal <- scales::seq_gradient_pal("#94ef5d", "#2b100c", "Lab")(seq(0,1,length.out=17))
ggplot(dbh_stackedbar_data, aes(x = reorder(size, order), y = n, fill = cat_dbh_m)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_cowplot() + 
  scale_fill_manual(values = barpal, name = "DBH (m)") +
  labs(x = "Patch size",
       y = "Relative abundance") +
  scale_y_reverse()
  

#


#                                                                   Overall DBH :::::::::####
# boxplot of dbh overall
boxplot_dbh_overall <- ggplot(forkless_eucdata, aes(y = dbh_m)) +
  geom_boxplot(fill = "#94ef5d") + theme_cowplot() +
  labs(title = "(a)",
       y = "Diameter at breast height (m)",
       x = "you shoudn't be seeing this (unless you're looking thru the script)") +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(colour = NA),
        axis.text.x = element_text(colour = NA),
        axis.ticks.x = element_line(colour = NA)) # yes i am making the 
                                                        # text white instead of 
                                                        # element_blank() so the
                                                        # axes of the hist + 
                                                        # boxplot line up

hist_dbh_overall <- ggplot(forkless_eucdata, aes(x = dbh_m, fill = hist_fill)) +
  geom_histogram(binwidth = 0.05, colour = "black", width = 0.25) + 
  theme_cowplot() +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values = histpal) +
  labs(x = "Diameter at breast height (m)", 
       y = "Site",
       title = "(b)") +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme(plot.title = element_text(size = 20), legend.position = "none")

overall_dbh <- plot_grid(boxplot_dbh_overall, hist_dbh_overall)


#                                                          Plots used in thesis :::::::::####

########################:
# FULL DATASET BY SITE #:
########################:

# plot boxplot + hist of dbh together
full_site

# save plot as a3 landscape
save_plot("plots/dbh_full_data_site.png", 
          full_site, 
          base_width = 16.5,
          base_height = 9)

#######################################:
# FULL VS SUBSAMPLED DATASETS BY SITE #:
#######################################:

# plot of boxplot + hist of full + subsampled data by site
full_subsampled_site

# save plot as a3 landscape
save_plot("plots/dbh_full_subsampled_data_site.png", 
          full_subsampled_site,
          base_width = 25.5,
          base_height = 15)

##############################:
# SUBSAMPLED DATASET BY SIZE #:
##############################:

# plot of boxplot + hist of subsampled data by site
dbh_size_subsampled

# save plot as a3 landscape
save_plot("plots/dbh_subsampled_size.png", 
          dbh_size_subsampled, 
          base_width = 16.5,
          base_height = 8.5)

#######################:
# PLOT OF DBH OVERALL #:
#######################:

# boxplot + hist of full data overall
overall_dbh

# save plot as a3 portrait
save_plot("plots/dbh_full_overall.png", 
          overall_dbh,
          base_width = 10,
          base_height = 13)

#

#                                                           DBH summary tables  :::::::::####
# summary table of DBH by site
site_summary_table <- forkless_eucdata %>% 
  group_by(site_id, size) %>% 
  summarise(n = n(), mean = mean(dbh_m),
            median = median(dbh_m), q1 = quantile(dbh_m, 0.25),
            q3 = quantile(dbh_m, 0.75), sd = sd(dbh_m))
site_summary_table

# summary table of DBH by size
size_summary_table <- subsampled_data %>% 
  group_by(size) %>% 
  summarise(n = n(), mean = mean(value),
            median = median(value), q1 = quantile(value, 0.25),
            q3 = quantile(value, 0.75), sd = sd(value))
size_summary_table


