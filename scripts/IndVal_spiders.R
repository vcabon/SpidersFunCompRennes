#  Create Indval folder and set exporting directory 
dir.create("outputs/Indval",showWarnings=F) 
save_path <- paste0("outputs/Indval/")

# Run Indval analysis to identify indicator species 
ind_values <- indval(spiders_sums_indval,df_clus[,-1],numitr=1000)
summary(ind_values)

ind_val_summ <- filter(data.frame(ind_values$maxcls,
                                  ind_values$indcls,
                                  ind_values$pval),
                       ind_values.pval <= 0.05)

ind_val_summ <- arrange(ind_val_summ, ind_values.maxcls, desc(ind_values.indcls))

write.csv(ind_val_summ,paste0(save_path, "species_indicator_values_clusters.csv"))

environment_scaled_clus <- rownames_to_column(environment_scaled_clus, var = "Site_ID")
ab.spiders.indval.glm <- spiders_sums_indval %>% 
    rownames_to_column(var = "Site_ID") %>% 
    full_join(environment_scaled_clus, by = "Site_ID") %>% 
    column_to_rownames(var = "Site_ID")
environment_scaled_clus <- column_to_rownames(environment_scaled_clus, var = "Site_ID")

# save each models' summary
glm.indval <- list()
rsq.indval <- list()
glm.indval.df <- data.frame()

for (i in colnames(spiders_sums_indval)) {
    
    glm.indval[[i]] <- summary(glm(get(i) ~ share_impervious_100 +
                                       mean_UHI_100 +
                                       Mean_local_UHI_total +
                                       Mean_veg_height +
                                       Mean_veg_herb,
                                   data = ab.spiders.indval.glm,
                                   family = "quasipoisson"))
    
    rsq.indval[[i]] <- rsq(glm(get(i) ~ share_impervious_100 +
                                       mean_UHI_100 +
                                       Mean_local_UHI_total +
                                       Mean_veg_height +
                                       Mean_veg_herb,
                                   data = ab.spiders.indval.glm,
                                   family = "quasipoisson"), 
                           adj = TRUE)
    

    glm.indval[[i]] <- glm.indval[[i]]$coefficients
    
    glm.indval.df <- as.data.frame(rbind(glm.indval.df, glm.indval[[i]]))
    
}

write.csv(glm.indval.df, paste0(save_path, "glm_indval_df.csv"))
