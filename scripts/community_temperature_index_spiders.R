# import spider distribution data per country obtained from: 
# World Spider Catalog https://wsc.nmbe.ch/ (accessed on 1 December 2022)

spiders_country <- read.csv("data/spiders_country.csv",
                                   sep = ";",
                                   header = TRUE)

spiders_country_sorted <- spiders_traits %>% 
  tibble::rownames_to_column(var = "Species_full") %>% 
  dplyr::left_join(spiders_country, by = "Species_full") %>% 
  dplyr::select(-c(Body,Hunt,Ballo)) %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "countries.name")

# Extract countries' temperatures 

countries <- st_read("data/spatial_data/Europe-administrative-boundaries.shp")

Tmean_temp_europe <- raster::stack("data/spatial_data/dataset-insitu-gridded-observations-europe/mean_0.25deg_2011-2022.nc")
Tmean_temp_europe <- calc(Tmean_temp_europe, fun = mean)
Av_mean_Temp_countries <- exact_extract(Tmean_temp_europe,
                                       countries,
                                       'mean')

Tmax_temp_europe <-  raster::stack("data/spatial_data/dataset-insitu-gridded-observations-europe/max_0.25deg_2011-2022.nc")
Tmax_temp_europe <- calc(Tmax_temp_europe, fun = mean)
Av_max_Temp_countries <- exact_extract(Tmax_temp_europe,
                                       countries,
                                       'mean')

Tmin_temp_europe <-  raster::stack("data/spatial_data/dataset-insitu-gridded-observations-europe/min_0.25deg_2011-2022.nc")
Tmin_temp_europe <- calc(Tmin_temp_europe, fun = mean)
Av_min_Temp_countries <- exact_extract(Tmin_temp_europe,
                                       countries,
                                       'mean')

Av_Temp_countries <- data.frame(countries$name, 
                               Av_mean_Temp_countries, 
                               Av_max_Temp_countries, 
                               Av_min_Temp_countries)
Av_Temp_countries <- Av_Temp_countries[order(Av_Temp_countries$countries.name),] 

# Compute thermal index per spider species 

spiders_country_sorted_Tpref <- Av_Temp_countries %>% 
  right_join(spiders_country_sorted, by = "countries.name")

Tpref_spiders_list <- list()

for (i in c("Av_mean_Temp_countries",
            "Av_max_Temp_countries", 
            "Av_min_Temp_countries")) {

Tpref_spiders_list[[i]] <- spiders_country_sorted_Tpref
  
Tpref_spiders_list[[i]][,c(5:length(spiders_country_sorted_Tpref))] <- apply(Tpref_spiders_list[[i]][,c(5:length(spiders_country_sorted_Tpref))],
                                      2,
                                      function(x) {(x = ifelse(x == "1",
                                                                 x,
                                                               Tpref_spiders_list[[i]][[i]]))
                                      })

Tpref_spiders_list[[i]] <- Tpref_spiders_list[[i]][,c(5:length(spiders_country_sorted_Tpref))] %>% 
  mutate_if(is.character, as.numeric)

Tpref_spiders_list[[i]] <- as.data.frame(colMeans(Tpref_spiders_list[[i]],
                                          na.rm = TRUE))
colnames(Tpref_spiders_list[[i]]) <- i

}

Tpref_spiders <- bind_cols(Tpref_spiders_list)
colnames(Tpref_spiders) <- c("Tpref_mean", "Tpref_max", "Tpref_min")



