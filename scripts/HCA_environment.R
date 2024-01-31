# Run Hierarchical clustering analysis based on environment variables

#  Create HCA folder and set exporting directory 
dir.create("outputs/HCA",showWarnings=F) 
save_path <- paste0("outputs/HCA/")

# Define global variables
M <- list("manhattan","euclidian") # metrics distance
m <- 2 #choose the correct metrics distance

# Build PCA
res_pca <- PCA(environment,
               quali.sup=1,
               scale.unit = TRUE,
               ncp=5,graph=F)  

# Build HCA

# Computing agglomeration coefficients (AC) to select method:
# AC describes the strength of the clustering structure that has been 
# obtained by group average linkage

# Preparing data for CAH / Scaling for finding best method
environment <- column_to_rownames(environment, var = "Site_ID")
environment_scaled=scale(environment, center = TRUE, scale = TRUE)

METHOD <- c( "average", "single", "complete", "ward") # assessment methods
names(METHOD) <- c("average", "single", "complete", "ward") 

# Function to compute agglomeration coefficients (AC)
agglo_coeff <- function(x) {
  agnes(environment_scaled, method = x)$ac # function from agnes() 
} 

# compute AC with 4 methods:

# Apply a  dbl function to each element of a list or atomic vector 
# thanks to the map_dbl() from the purr package
ac_tab <- purrr::map_dbl(METHOD, agglo_coeff) 
print(ac_tab)
ac_tab <- as.data.frame(ac_tab)
print(ac_tab)
names(ac_tab)=c("ac")
ac_max=max(ac_tab$ac) # Values close to 1 indicate a balanced clustering
print(ac_max)

ac_tab <- cbind(newColName = rownames(ac_tab), ac_tab)
names(ac_tab)[1] ="ac_method"

Best_table=subset(ac_tab,ac==ac_max)
Best_method=Best_table[1,1]
print(paste0("The best agglomeration method is: ", Best_method))

#build the best HCA with the optimum number of clusters
best_hcpc<- HCPC(res_pca, 
                 nb.clust=-1,
                 metric = M[m],
                 method = Best_method, 
                 graph=F) 

plot(res_pca, choix="var", axes= c(1,2))
        
# Inertia histogram
NOMpng <- "HCA_inertia.png"
png(file = paste0(save_path, NOMpng), width = 700, height = 700)
arbre2 <- agnes(environment_scaled, method = Best_method,metric=M[[m]])
as.hclust(arbre2)
inertie <- sort(arbre2$height, decreasing = TRUE)
p=plot(inertie[1:20], type = "s", 
       xlab = "Nombre de classes", ylab = "Inertie")
points(c(3, 4, 5, 6, 7, 8, 9), inertie[c(3, 4, 5, 6, 7, 8, 9)], 
       col = c("green3", "red3", "blue3","orange", "grey", "purple", "pink"), 
       cex = 2, lwd = 3)
legend(x = "topright", # Position
       legend = c("3 clusters", 
                  "4 clusters", 
                  "5 clusters",
                  "6 clusters", 
                  "7 clusters", 
                  "8 clusters",
                  "9 clusters"),  # Legend texts
       fill =c("green3", "red3", "blue3","orange", "grey", "purple", "pink"))
print(p)
dev.off()

# Graph of sites (individuals), colored by group
NOMpng <- "HCA_clusters.png"
png(file = paste0(save_path, NOMpng), width = 700, height = 700)
p=fviz_cluster(best_hcpc,
               repel = TRUE,
               geom="point",# Evite le chevauchement des textes
               show.clust.cent = TRUE, # Montre le centre des clusters
               palette = "jco", # Palette de couleurs, voir eggpubr::ggpar
               ggtheme = theme_minimal(),
               main = paste0(M[m]," distance")) 
p+geom_text(aes(label=""))
print(p)
dev.off()

# Principal components + tree
NOMpng="HCA_Tree_PCA.png"
png(file = paste0(save_path, NOMpng), width = 700, height = 700)
p=plot(best_hcpc, choice = "3D.map",sub=paste0(M[m]," distance"))
print(p)
dev.off()

# Analysing cluster thanks to quantitative variables
print(best_hcpc$desc.var$quanti) # quantitative variables describing each cluster

# determine nb of cluster k :
levels=unique(best_hcpc[["call"]][["X"]][["clust"]])
levels=as.numeric(levels)
K=max(levels)
print(paste0("Le nombre de clusters est de ", K))

# silhouette 
clus <- as.integer(best_hcpc$data.clust$clust)
names(clus) <- best_hcpc$data.clust$Site_ID
NOMpng <- "HCA_silhouette.png"
png(file = paste0(save_path, NOMpng), width = 700, height = 700)
sil <- silhouette(clus, dist(res_pca$ind$coord, method = "euclidian"))
rownames(sil) <- names(clus)
p=plot(sil)
dev.off()

# Clusters' characterization
## Mean indices value
df_clus <- as.data.frame(clus) %>% 
  rownames_to_column(var ="Site_ID") 

environment_scaled_clus <- environment %>%
  mutate(across(1:ncol(.), scale)) %>%
  rownames_to_column(var ="Site_ID") %>%
  inner_join(df_clus, by = "Site_ID") %>%
  mutate(clus = as.character(clus)) %>%
  column_to_rownames(var ="Site_ID")

environment_scaled_clus_avrgd <- aggregate(environment_scaled_clus[,1:length(environment_scaled_clus)-1],
                                      by = list(environment_scaled_clus$clus),
                                      FUN=mean,na.rm=T)
environment_scaled_clus_std <- aggregate(environment_scaled_clus[,1:length(environment_scaled_clus)-1],
                                    by = list(environment_scaled_clus$clus),
                                    FUN=sd,na.rm=T)

colnames(environment_scaled_clus_avrgd)[ 1] <- "cluster" ## rename Cluster
colnames(environment_scaled_clus_std)[ 1] <- "cluster" ## rename Cluster

environment_scaled_clus_avrgd <- melt(environment_scaled_clus_avrgd, id.vars = "cluster")
environment_scaled_clus_std <- melt(environment_scaled_clus_std, id.vars = "cluster")
environment_scaled_clus_avrgd[["std"]] <- environment_scaled_clus_std$value

### Plot clusters' characterization according to environment variables 
clust.labs <- c("Cluster 1", "Cluster 2", "Cluster 3") 
names(clust.labs) <- c("1", "2", "3")

NOMpng="parameters_environment.png"
png(file = paste0(save_path, NOMpng), width = 700, height = 700)
p=ggplot(transform(environment_scaled_clus_avrgd,
                   clus=factor(cluster,levels=c("1","2","3")))) + ## graph en bar
  geom_bar(aes(x = variable, 
               y = value),
           stat = "identity") +
  facet_wrap(~ cluster,labeller = labeller(cluster=clust.labs)) +
  coord_flip() +theme(legend.position = "none") +
  geom_errorbar(aes(ymin=value-std, ymax=value+std,x=variable), width=.2)
plot(p)
dev.off()

