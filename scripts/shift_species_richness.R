# Test if shift in community-wide functional shifts are due to species 
# replacement or to a decrease of species with particular trait values

# Classify species in 3 classes according to their thermal affinity
cutoff_Tmean <- bin(Tpref_spiders$Tpref_mean, nbins = 3, method = "content")
cutoff_Tmax <- bin(Tpref_spiders$Tpref_max, nbins = 3, method = "content")
cutoff_Tmin <- bin(Tpref_spiders$Tpref_min, nbins = 3, method = "content")

Tpref_spiders <- rownames_to_column(Tpref_spiders, var = "species")
Tpref_spiders_class <- data.frame(Tpref_spiders$species, 
                                  cutoff_Tmean, 
                                  cutoff_Tmax, 
                                  cutoff_Tmin)
Tpref_spiders <- column_to_rownames(Tpref_spiders, var = "species")

names(Tpref_spiders_class) <- c("Species_full", 
                                "cutoff_Tmean", 
                                "cutoff_Tmax", 
                                "cutoff_Tmin")

# Calculate specific richness on sites considering  
# separately species from each Tmean class  

new_names <- c(class_low = "(7.19,8.99]", 
               class_mid = "(8.99,9.48]",
               class_high = "(9.48,11.4]")
environment_scaled_clus <- rownames_to_column(environment_scaled_clus, var = "Code_site")
spiders_Tmean_richness <- 
    spiders_ab.traits %>% 
    group_by(Code_site, Species_full) %>%
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(n = 1) %>%
    left_join(Tpref_spiders_class, by = "Species_full") %>% 
    pivot_wider(names_from = cutoff_Tmean, 
                values_from = n, 
                values_fill = list(n = 0)) %>% 
    select(-c("NA")) %>% 
    rename(all_of(new_names)) %>% 
    na.omit() %>% 
    mutate(class_tot = 1) %>% #### 
    group_by(Code_site) %>% 
    summarise(Ric_low_sp = sum(class_low),
              Ric_mid_sp = sum(class_mid),
              Ric_high_sp = sum(class_high),
              Ric_tot_sp = sum(class_tot)) %>% ####
    filter(!Code_site == "Lm_01",
           !Code_site == "Hh_01",
           !Code_site == "Lh_03") %>% 
    full_join(environment_scaled_clus, by = "Code_site") %>% 
    column_to_rownames(var = "Code_site") %>% 
    relocate(Ric_tot_sp)
rm(new_names)
environment_scaled_clus <- column_to_rownames(environment_scaled_clus, var = "Code_site")


# Classify species in 3 classes according to their body size
cutoff_size <- bin(spiders_traits_fc$Body, nbins = 3, method = "content")

Size_spiders_class <- spiders_traits_fc %>% 
    rownames_to_column(var="Species") %>% 
    select(Species) 
Size_spiders_class <- data.frame(Size_spiders_class, cutoff_size)

names(Size_spiders_class) <- c("Species_full", 
                               "cutoff_size")

# Calculate specific richness on sites considering  
# separately species from each size class  

new_names <- c(class_small = "(0.983,3]", 
               class_mid = "(3,5.9]",
               class_large = "(5.9,18]")
environment_scaled_clus <- rownames_to_column(environment_scaled_clus, var = "Code_site")
spiders_size_richness <- 
    spiders_ab.traits %>% 
    group_by(Code_site, Species_full) %>%
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(n = 1) %>%
    left_join(Size_spiders_class, by = "Species_full") %>% 
    pivot_wider(names_from = cutoff_size, 
                values_from = n, 
                values_fill = list(n = 0)) %>% 
    select(-c("NA")) %>% 
    rename(all_of(new_names)) %>% 
    na.omit() %>% 
    mutate(class_tot = 1) %>% #### 
    group_by(Code_site) %>% 
    summarise(Ric_small_sp = sum(class_small),
              Ric_mid_sp = sum(class_mid),
              Ric_large_sp = sum(class_large),
              Ric_tot_sp = sum(class_tot)) %>% ####
    filter(!Code_site == "Lm_01",
           !Code_site == "Hh_01",
           !Code_site == "Lh_03") %>% 
    full_join(environment_scaled_clus, by = "Code_site") %>% 
    column_to_rownames(var = "Code_site") %>% 
    relocate(Ric_tot_sp)
rm(new_names)
environment_scaled_clus <- column_to_rownames(environment_scaled_clus, var = "Code_site")


### Thermal affinity ###

# Test differences in Richness of species with high, mid and low  
# temperature affinity among clusters
richness_Tmean_aov_res <- list()
richness_Tmean_hsd_res <- list()

for (i in names(spiders_Tmean_richness)[c(1:4)]) 
{

resAOV.Tmean.rich <- summary(res.aov.Tmean.rich <- aov(get(i) ~ clus, 
                                           data = spiders_Tmean_richness))
richness_Tmean_aov_res[[i]] <- resAOV.Tmean.rich
richness_Tmean_hsd_res[[i]] <- TukeyHSD(res.aov.Tmean.rich)
}

# Plot differences in richness among clusters and save plot in a list 

richness_Tmean_boxplot <- list()
for (i in names(spiders_Tmean_richness)[c(2:4)]) 
{
    richness_Tmean_boxplot[[i]] <- ggboxplot(spiders_Tmean_richness, 
                                            x = "clus", 
                                            y = paste0(i), 
                                            order = c("1", "2", "3"),
                                            ylab = paste0(i), 
                                            xlab = "cluster")
}

NOMpng=paste0("outputs/richness_comparison/richness_Tmean_boxplot.png")
png(file = NOMpng, width = 740, height = 300) 
ggarrange(richness_Tmean_boxplot$Ric_low_sp,
          richness_Tmean_boxplot$Ric_mid_sp,
          richness_Tmean_boxplot$Ric_high_sp,
          labels = c("a)", "b)", "c)"),
          ncol = 3, nrow = 1)
dev.off()

# Plot differences in richness among clusters in a single plot and save
# width=570, height=400

NOMpng=paste0("outputs/richness_comparison/richness_Tmean_boxplot_2.png")
png(file = NOMpng, width = 740, height = 300) 
ggplot(melt(spiders_Tmean_richness[,c(1,2,3,4,10)], 
            id = "clus"), 
       aes(x = variable, 
           y = value, 
           fill = clus)) +
    geom_boxplot() + 
    theme_classic() +
    scale_fill_manual(values=c("lightgrey", "grey55", "grey35"))
dev.off()


#  Averaged low-temp species richness by cluster 
spiders_Tmean_richness %>% 
    group_by(clus) %>% 
    summarise(mean_Ric_low = mean(Ric_low_sp))

#  Averaged mid-temp species richness by cluster 
spiders_Tmean_richness %>% 
    group_by(clus) %>% 
    summarise(mean_Ric_mid = mean(Ric_mid_sp))

#  Averaged high-temp species richness by cluster 
spiders_Tmean_richness %>% 
    group_by(clus) %>% 
    summarise(mean_Ric_high = mean(Ric_high_sp))


# Run glm on species richness according to Tmean classes to identify determinant variables

glm.rich.temp <- list()
rsq.rich.temp <- list()
glm.rich.temp.df <- data.frame()

for (i in names(spiders_Tmean_richness[,c(2:4)])) {
    
    glm.rich.temp[[i]] <- summary(glm(sqrt(get(i)) ~ share_impervious_100 +
                                          mean_UHI_100 +
                                          Mean_local_UHI_total +
                                          Mean_veg_height +
                                          Mean_veg_herb,
                                      data = spiders_Tmean_richness))
    
    rsq.rich.temp[[i]] <- rsq(glm(sqrt(get(i)) ~ share_impervious_100 +
                                          mean_UHI_100 +
                                          Mean_local_UHI_total +
                                          Mean_veg_height +
                                          Mean_veg_herb,
                                      data = spiders_Tmean_richness),
                              adj = TRUE)
    
    glm.rich.temp[[i]] <- glm.rich.temp[[i]]$coefficients
    
    glm.rich.temp.df <- as.data.frame(rbind(glm.rich.temp.df, 
                                            glm.rich.temp[[i]]))
}    

### Body size ###

# Test differences in Richness of species with large, medium and small  
# body size among clusters

richness_size_aov_res <- list()
richness_size_hsd_res <- list()

for (i in names(spiders_size_richness)[c(1:4)]) 
{
    
resAOV.size.rich <- summary(res.aov.size.rich <- aov(get(i) ~ clus, 
                                           data = spiders_size_richness))

richness_size_aov_res[[i]] <- resAOV.size.rich
richness_size_hsd_res[[i]] <- TukeyHSD(res.aov.size.rich)
}

# Plot differences in richness among clusters and save plot in a list 

richness_size_boxplot <- list()
for (i in names(spiders_size_richness)[c(2:4)]) 
    {
    richness_size_boxplot[[i]] <- ggboxplot(spiders_size_richness, 
                                       x = "clus", 
                                       y = paste0(i), 
                                       order = c("1", "2", "3"),
                                       ylab = paste0(i), 
                                       xlab = "cluster")
}

NOMpng=paste0("outputs/richness_comparison/richness_size_boxplot.png")
png(file = NOMpng, width = 740, height = 300) 
ggarrange(richness_size_boxplot$Ric_small_sp,
          richness_size_boxplot$Ric_mid_sp,
          richness_size_boxplot$Ric_large_sp,
          labels = c("a)", "b)", "c)"),
          ncol = 3, nrow = 1)
dev.off()

# Plot differences in richness among clusters in a single plot and save
# width=570, height=400
ggplot(melt(spiders_size_richness[,c(1,2,3,4,10)], 
            id = "clus"), 
       aes(x = variable, 
           y = value, 
           fill = clus)) +
    geom_boxplot() + 
    theme_classic() +
    scale_fill_manual(values=c("lightgrey", "grey55", "grey35"))
    

#  Averaged small species richness by cluster 
spiders_size_richness %>% 
    group_by(clus) %>% 
    summarise(mean_Ric_small = mean(Ric_small_sp))

#  Averaged mid species richness by cluster 
spiders_size_richness %>% 
    group_by(clus) %>% 
    summarise(mean_Ric_mid = mean(Ric_mid_sp))

#  Averaged large species richness by cluster 
spiders_size_richness %>% 
    group_by(clus) %>% 
    summarise(mean_Ric_large = mean(Ric_large_sp))


# Run glm on species richness according to size classes to identify determinant variables

glm.rich.size <- list()
rsq.rich.size <- list()
glm.rich.size.df <- data.frame()

for (i in names(spiders_size_richness[,c(2:4)])) {
    
    glm.rich.size[[i]] <- summary(glm(sqrt(get(i)) ~ share_impervious_100 +
                                          mean_UHI_100 +
                                          Mean_local_UHI_total +
                                          Mean_veg_height +
                                          Mean_veg_herb,
                                      data = spiders_size_richness))
    
    rsq.rich.size[[i]] <- rsq(glm(sqrt(get(i)) ~ share_impervious_100 +
                                          mean_UHI_100 +
                                          Mean_local_UHI_total +
                                          Mean_veg_height +
                                          Mean_veg_herb,
                                      data = spiders_size_richness),
                              adj = TRUE)
    
    glm.rich.size[[i]] <- glm.rich.size[[i]]$coefficients
    
    glm.rich.size.df <- as.data.frame(rbind(glm.rich.size.df, 
                                            glm.rich.size[[i]]))
}    

# Plot differences in richness among clusters in a single plot and save
spiders_Tmean_richness <- spiders_Tmean_richness %>% 
    rownames_to_column(var = "Site_ID")
spiders_richness_full <- spiders_size_richness %>% 
    rownames_to_column(var = "Site_ID") %>% 
    full_join(spiders_Tmean_richness, by = "Site_ID") %>% 
    column_to_rownames(var = "Site_ID")
spiders_Tmean_richness <- spiders_Tmean_richness %>% 
    column_to_rownames(var = "Site_ID")

# Update plot changes in absolute species richness 
# width=1000, height=400
plot_changes_species_richness <- ggplot(melt(spiders_richness_full[,c(11,2,3,4,12,13,14,10)], 
            id = "clus.x"), 
       aes(x = variable, 
           y = value, 
           fill = clus.x)) +
    geom_boxplot() + 
    theme_classic() +
    scale_fill_manual(values=c("#137C8B", "#B8CBD0", "#595959"))

# Update plot changes in relative species richness 
# width=1000, height=400
spiders_richness_full_percent <- read.csv("data/spiders_richness_full_percent.csv", 
                                          sep = ",", 
                                          header = TRUE, 
                                          row.names = 1)
spiders_richness_full_percent$clus.x <- as.factor(spiders_richness_full_percent$clus.x)

plot_changes_species_richness_percent <- ggplot(melt(spiders_richness_full_percent,
                                                     id = "clus.x"), 
                                                aes(x = variable, 
                                                    y = value, 
                                                    fill = clus.x)) +
    geom_boxplot() + 
    theme_classic() +
    scale_fill_manual(values=c("#137C8B", "#B8CBD0", "#595959"))

# Test differences in relative species richness in functional classes among clusters

richness_size_aov_res_rel <- list()
richness_size_hsd_res_rel <- list()

for (i in names(spiders_richness_full_percent)[c(1:6)]) 
{
    
    resAOV.size.rich.rel <- summary(res.aov.size.rich.rel <- aov(get(i) ~ clus.x, 
                                                         data = spiders_richness_full_percent))
    
    richness_size_aov_res_rel[[i]] <- resAOV.size.rich.rel
    richness_size_hsd_res_rel[[i]] <- TukeyHSD(res.aov.size.rich.rel)
}
