# Create trait community-mean data set

# Compute Community Mean indices (non-weighted)
spiders_ab.traits <- spiders %>% 
    full_join(spiders_traits.fc, by = "Species_full")

spiders_traits.cm <-   
    spiders_ab.traits %>%  
    mutate_at("Ballo", as.numeric) %>% 
    group_by(Code_site) %>%   
    distinct(Species_full, .keep_all = TRUE) %>% 
    dplyr::summarize(          
        Body = mean(Body, na.rm = TRUE),
        Ballo = mean(Ballo, na.rm = TRUE),
        Tmean_mean = mean(Tpref_mean, na.rm = TRUE),
        Tmean_sd = sd(Tpref_mean, na.rm = TRUE),
        Tmean_min = min(Tpref_mean, na.rm = TRUE),
        Tmean_max = max(Tpref_mean, na.rm = TRUE),
        Tmin_mean = mean(Tpref_min, na.rm = TRUE),
        Tmin_sd = sd(Tpref_min, na.rm = TRUE),
        Tmin_min = min(Tpref_min, na.rm = TRUE),
        Tmin_max = max(Tpref_min, na.rm = TRUE),
        Tmax_mean = mean(Tpref_max, na.rm = TRUE),
        Tmax_sd = sd(Tpref_max, na.rm = TRUE),
        Tmax_min = min(Tpref_max, na.rm = TRUE),
        Tmax_max = max(Tpref_max, na.rm = TRUE)) %>% 
    filter(!Code_site == "Lm_01", !Code_site == "Hh_01", !Code_site == "Lh_03") %>% # remove 3 sites to exclude 
    column_to_rownames(var = "Code_site")

# Run glm on traits' community mean metrics to identify determinant variables

# formatting dataset 
{
  environment_scaled_clus <- rownames_to_column(environment_scaled_clus, 
                                              var = "Site_ID")

spiders_traits.cm.env <- spiders_traits.cm %>% 
    rownames_to_column(var = "Site_ID") %>% 
    full_join(environment_scaled_clus, by = "Site_ID") %>% 
    column_to_rownames(var = "Site_ID")

environment_scaled_clus <- column_to_rownames(environment_scaled_clus, 
                                              var = "Site_ID")
    }

# save each models' summary

glm.cm <- list()
rsq.cm <- list()
glm.cm.df <- data.frame()

for (i in names(spiders_traits.cm)) {
    
    glm.cm[[i]] <- summary(glm(get(i) ~ share_impervious_100 +
                                   mean_UHI_100 +
                                   Mean_local_UHI_total +
                                   Mean_veg_height +
                                   Mean_veg_herb,
                               data = spiders_traits.cm.env))
    
    rsq.cm[[i]]<- rsq(glm(get(i) ~ share_impervious_100 +
                                  mean_UHI_100 +
                                  Mean_local_UHI_total +
                                  Mean_veg_height +
                                  Mean_veg_herb,
                              data = spiders_traits.cm.env),
                      adj = TRUE)
    
    glm.cm[[i]] <- glm.cm[[i]]$coefficients
    
    glm.cm.df <- as.data.frame(rbind(glm.cm.df, glm.cm[[i]]))
    
}

write.csv(glm.cm.df, "outputs/functional_comp/glm_cm_df.csv")

# Test differences in CM among clusters and save results in lists 
traits_aov_res <- list()
traits_hsd_res <- list()

for (i in names(spiders_traits.cm.env)[c(1:14)]) 
{
resAOV <- summary(res.aov <- aov(get(i) ~ clus, 
                                 data = spiders_traits.cm.env))

traits_aov_res[[i]] <- resAOV
traits_hsd_res[[i]] <- TukeyHSD(res.aov)

}

#  Averaged community mean body size by cluster 
spiders_traits.cm.env %>% 
    group_by(clus) %>% 
    summarise(mean_body = mean(Body))

#  Averaged community mean Tmean by cluster 
spiders_traits.cm.env %>% 
    group_by(clus) %>% 
    summarise(mean_Tmean = mean(Tmean_mean))

#  Averaged community mean Tmax by cluster 
spiders_traits.cm.env %>% 
    group_by(clus) %>% 
    summarise(mean_Tmax = mean(Tmax_mean))

#  Averaged community mean Tmin by cluster 
spiders_traits.cm.env %>% 
    group_by(clus) %>% 
    summarise(mean_Tmin = mean(Tmin_mean))

# Plot differences in CM among clusters and save plot in a list 
# export size = 250*300
traits_boxplot <- list()
for (i in names(spiders_traits.cm.env)[c(1:14)]) 
{
    traits_boxplot[[i]] <- ggboxplot(spiders_traits.cm.env, 
                                     x = "clus", 
                                     y = paste0(i), 
                                     fill = "clus",
                                     palette = c("#137C8B", "#B8CBD0", "#595959"),
                                     add = "jitter",
                                     order = c("1", "2", "3"),
                                     ylab = paste0(i), 
                                     xlab = "cluster")
}

traits_boxplot$Tmean_mean + traits_boxplot$Tmin_mean + traits_boxplot$Tmax_mean +
    plot_layout(guides = 'collect')
