# Run NMDS analyses + Permanova to test effect of clusters + 
# Monte-Carlo to test the effect of environment variables

# prepare data for NMDS computing
spiders_env <- environment_scaled_clus
names(spiders_env) <- c("Urbanization", 
                        "UHI intensity", 
                        "Local temperature", 
                        "Vegetation height", 
                        "Vegetation proportion",
                        "Cluster")

# compute NMDS - method Bray
nmds_spiders <- metaMDS(sqrt(spiders_sums_indval), 
                        distance = "bray", autotransform = F)

spiders.env_all.fit <- envfit(nmds_spiders, 
                              environment_scaled) # test all predictors

spiders.env.fit <- envfit(nmds_spiders, 
                          environment_scaled, permutations = 9999) # select predictors

spiders.sp.fit <- envfit(nmds_spiders, 
                         spiders_sums_indval)

source("./utils/function_plot_nmds.R")

plot.nmds.spiders.env <- plot.nmds(nmds_spiders, 
                               spiders.sp.fit, 
                               pval.sp = 0.05,
                               spiders.env.fit,
                               pval.env = 0.05,
                               env.var = "Cluster",
                               ssp.fit = FALSE) 

# test effect of predictors on beta-diversity by PERMANOVA
adonis2(sqrt(spiders_sums_indval) ~ spiders_env$Cluster, permutations = 9999)

spiders_sums_bray <- vegdist(spiders_sums_indval,method="bray")

pairwise.perm.manova(spiders_sums_bray, 
                     spiders_env$Cluster,
                     nperm = 9999, 
                     progress = TRUE, 
                     p.method = "fdr")
