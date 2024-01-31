# Load required packages
source("utils/load_packages.R")

# Create "outputs" folder
dir.create("outputs", showWarnings=F) 

# Load environment data
environment <- read.csv(file.path("data", "environment.csv"), 
                        header = TRUE,
                        row.names = 1) 

# Load abundance data 
spiders <- read.csv(file.path("data", "spiders.csv"), 
                    sep = ";",
                    header = TRUE,
                    row.names = 1) 

# Load functional trait data
{
  spiders_traits <- read.csv(file.path("data", "spiders_traits.csv"),
                             sep = ";",
                             row.names = 1, 
                             header= TRUE)
  
  spiders_traits$Hunt <- as.factor(spiders_traits$Hunt)
  spiders_traits$Ballo <- as.factor(spiders_traits$Ballo)
}

# # Load thermal index data
# Tpref_spiders <- source("community_temperature_index_spiders.R")

# Pool abundance data into a table of one sum per site and remove sampling sites
# excluded during the sampling site procedure 
spiders_sums <- spiders %>% 
  group_by(Species_full, Code_site) %>% 
  summarise(sum = sum(Abundance)) %>% 
  pivot_wider(names_from = Code_site, values_from = sum) %>% 
  filter(!Species_full == "na") %>% 
  column_to_rownames(var = "Species_full") %>%
  dplyr::select(-c(Lm_01, Hh_01, Lh_03)) %>% # remove 3 sites to exclude  
  replace(is.na(.), 0) %>%
  dplyr::select(order(colnames(.))) %>% 
  mutate(sum_row = rowSums(across(where(is.numeric)))) %>%  
  filter(!sum_row == 0) %>% # remove non-observed species  
  dplyr::select(-c(sum_row))

# Formatting species data for ordination and Indval analysis with rare species omitted
{       
  n = 3 # set the minimal number of individuals per site to consider a species 
  
  spiders_sums_indval <- apply(spiders_sums, 2, 
                             function(x) ifelse(x < n, 0, x)) %>% 
        as.data.frame() %>%
        rownames_to_column(var = "species") %>%
        t() %>%
        row_to_names(row_number = 1) %>%
        as.data.frame()
    
  spiders_sums_indval[ , 
                       c(1:ncol(spiders_sums_indval))] <- apply(spiders_sums_indval[ ,
                                                                                 c(1:ncol(spiders_sums_indval))], 
                                                              2,
                                                              function(x) {
                                                                  as.numeric((x))})
  spiders_sums_indval = spiders_sums_indval[,colSums(spiders_sums_indval) != 0]
}

# Formatting species data for analysis of functional composition with no species omitted
{  
  spiders_sums_fc <- spiders_sums %>% 
        as.data.frame() %>%
        rownames_to_column(var = "species") %>%
        t() %>%
        row_to_names(row_number = 1) %>%
        as.data.frame()
    
  spiders_sums_fc[,c(1:ncol(spiders_sums_fc))] <- apply(spiders_sums_fc[,
                                                                    c(1:ncol(spiders_sums_fc))], 
                                                      2,
                                                      function(x) {
                                                          as.numeric((x))})
  
  spiders_sums_fc <- spiders_sums_fc[,colSums(spiders_sums_fc) != 0]
}

# Formatting functional trait data
{
  spiders_traits <- rownames_to_column(spiders_traits, var = "species")
  Tpref_spiders <- rownames_to_column(Tpref_spiders, var = "species")
  
  spiders_traits_fc <- as.data.frame(t(spiders_sums_fc)) %>% 
    rownames_to_column(var="species") %>%
    left_join(spiders_traits, by = "species") %>%
    left_join(Tpref_spiders, by = "species") %>%
    column_to_rownames(var = "species") %>% 
    dplyr::select(Body, Hunt, Ballo, Tpref_mean, Tpref_max, Tpref_min) %>% 
    rownames_to_column(var = "Species_full") 
  
  
    spiders_traits <- column_to_rownames(spiders_traits, var = "species")
    Tpref_spiders <- column_to_rownames(Tpref_spiders, var = "species")
}

# Run Hierarchical clustering analysis based on environment variables
source("scripts/HCA_environment.R")

# Run NMDS analyses + Permanova to test effect of clusters + 
# Monte-Carlo to test the effect of environment variables
source("scripts/NMDS_spiders.R")

# Run Indval analysis to identify indicator species 
source("scripts/IndVal_spiders.R")

# Compute community trait metrics, compare among cluster and run GLM 
source("scripts/functional_composition.R")

# Analyse shifts in species richness 
source("scripts/shift_species_richness.R")

