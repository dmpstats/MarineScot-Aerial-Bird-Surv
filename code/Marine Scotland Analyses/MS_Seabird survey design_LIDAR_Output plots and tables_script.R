### ------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------- ###
### ---                                                                         --- ###
### ---     Script for summary tables and plots for LiDAR results               --- ###
### ---                                                                         --- ###
### ---                                                                         --- ###
### ------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------- ###


## --------------------- ###
## ----   Preamble     ----
## --------------------- ###

require(tidyverse)
require(sf)
require(tmap)
require(tmaptools)
require(units)
require(progress)
require(data.table)
require(scales)
require(furrr)
require(magrittr)



outFolder <- "outputs/Report Plots/LiDAR FHD/"


speciesKey <- tribble(~species_name, ~species_code,
                      "Atlantic Puffin", "FRARC",
                      "Black Legged Kittiwake", "RITRI",
                      "Common Guillemot", "URAAL",
                      "Herring Gull", "LAARG",
                      "Northern Fulmar", "FUGLA",
                      "Northern Gannet", "MOBAS",
                      "Razorbill", "ALTOR", 
                      "Great black-backed gull", "LAMAR")




## ---------------------------------------------- ###
## ----   Upload results and auxiliary data      ----
## ---------------------------------------------- ###

# read in lidar results
sim_lidar <- list.files("outputs/simulations/", pattern = "LiDAR", full.names = TRUE) %>%
  map(function(x){
    read_rds(x)
  })


names(sim_lidar) <- list.files("outputs/simulations/", pattern = "LiDAR", full.names = FALSE)


# upload survey region polygons
Moray_Firth <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Moray_Firth_surevey_area/Wider_Moray_Firth_survey_area.shp")
Herman_Berwick <- st_read(dsn = "data/East coast surveys Marine Scotland/East herman to Berwick/East herman to Berwick.shp")
east_Scotland <- st_read(dsn = "data/East coast surveys Marine Scotland/East_coast_Scotland/East_coast_Scotland.shp")
Forth_Tay <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Forth_Tay_survey_area/Wider_Forth_Tay_survey_area.shp")


# re-project survey region to density surface projections
Forth_Tay_utm <- st_transform(Forth_Tay, st_crs(sim_lidar[[1]][[1]]))
Herman_Berwick_utm <- st_transform(Herman_Berwick, st_crs(sim_lidar[[1]][[1]]))
east_Scotland_utm <- st_transform(east_Scotland, st_crs(sim_lidar[[1]][[1]]))
Moray_Firth_utm <- st_transform(Moray_Firth, st_crs(sim_lidar[[1]][[1]]))



# Upload Coastline & re-project 
coastline <- read_sf("data/GBR_adm/GBR_adm0.shp") %>% select(NAME_LOCAL)
coastline
coastline_utm <- st_transform(coastline, st_crs(sim_lidar[[1]][[1]]))

coastline_utm_cropped <- st_crop(coastline_utm, bb(Forth_Tay_utm, ext = 3))
plot(coastline_utm_cropped)








## --------------------------------------------------------------------- ###
## ----   LiDAR Binning plots - 3 species, all regions, by season       ----
## ------------------------------------------------------------------=-- ###


# Gather region polygons into a grid
surveyRegions_ls <- list('EastHerman Berwick' = Herman_Berwick_utm, 'Moray Firth' = Moray_Firth_utm, 
                         'East Scotland' = east_Scotland_utm, 'Forth Tay' = Forth_Tay_utm)



library(wesanderson)
Lidar_pal <- wes_palette("Zissou1", 15, type = "continuous")


lidar_binMaps <- surveyRegions_ls %>%
  map2(., names(.), function(region_Poly, region_name){
    
    suppressWarnings({
      coastline_utm_cropped <- st_crop(coastline_utm, bb(region_Poly, ext = 3))
    })
    
    #browser()
    
    # subset for region
    cReg_Lidar <- sim_lidar[[str_which(names(sim_lidar), region_name)]]
  
    sp_code <- c("URAAL", "FRARC", "RITRI")
    
    #browser()
    
    sp_plots <- map(sp_code, function(sp){
      
      #browser()
      
      # subset for species
      cReg_cspec_Lidar <- cReg_Lidar[str_which(names(cReg_Lidar), paste0(" ", sp))]
      
      season <- c("winter", "summer", "autumn")
      
      n_breaks <- pretty_breaks(8)(c(10, 50, 100, unlist(map(cReg_cspec_Lidar, "n")), 1500))
      
      pcounter <- 1
      sp_seas_plots <- map(season, function(seas){
        
        dat <- cReg_cspec_Lidar[[str_which(names(cReg_cspec_Lidar), paste0(" ", seas))]] %<>% 
          mutate(n = round(n, 0), 
                 surv_cov_dist = round(set_units(surv_cov_dist, "km"), 0),
                 n_text = paste0("n = ", n),
                 effort_text = paste0(surv_cov_dist ,"km")
          )
        
        pTitle <- str_to_title(seas)
        pMainTitle <- filter(speciesKey, species_code == sp) %>% pull(species_name)
        
        if(pcounter != 1){
          pMainTitle <- " "
        }
        
        p <- tm_shape(region_Poly) +
          tm_borders(lwd = 2) +
          tm_shape(dat) +
          tm_polygons("n", breaks = n_breaks, palette = Lidar_pal, alpha = 0.6, midpoint = 200) +
          tm_text("effort_text", size = 0.6, ymod = -0.5, col = "black") +
          tm_shape(dat) +
          tm_text("n", size = 0.6, ymod = 0, fontface = "bold") +
          tm_shape(coastline_utm_cropped) +
          tm_polygons(col = "tan", border.col = "burlywood4") +
          tm_layout(title = pTitle, title.size = 0.9, title.position = c("LEFT", "TOP"),
                    title.bg.color = "white", title.bg.alpha = 0.8,
                    main.title = pMainTitle, main.title.size = 1.1, 
                    legend.show = FALSE)
        
        pcounter <<- pcounter + 1
        
        p
        
      })
      
      #browser()
      
      sp_seas_plots_arranged <- tmap_arrange(sp_seas_plots, ncol = 3)
      sp_name <- filter(speciesKey, species_code == sp) %>% pull(species_name)
      tmap_save(sp_seas_plots_arranged, 
                filename = paste0(outFolder, "Lidar seasonal maps_", region_name, "_", sp_name, "_mean.png"),
                dpi = 300, asp = 0, width = 13, height = 7, units = "in"
      )
      
      sp_seas_plots_arranged
      
    })
    
    names(sp_plots) <- sp_code
    sp_plots
    
  })








## ------------------------------------------------ ###
## ----          LiDAR summary Table               ----
## ------------------------------------------------ ###


lidarRes_df <- sim_lidar %>%
  map2(., names(.), function(x, y){
    
    #browser()
    
    c_region <- str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+")
    
    x %>%
      map(st_drop_geometry) %>%
      rbindlist(use.names = TRUE, idcol = "id") %>%
      mutate(region = c_region,
             sp_code = str_extract(id, "(?<=Species = )\\w+"),
             season = str_to_title(str_extract(id, "(?<=Season = )\\w+")), 
             predLayer = case_when(
               str_detect(sp_code, "ci025") ~ "ci025" ,
               str_detect(sp_code, "ci975") ~ "ci975",
               TRUE ~ "Mean"
             ), 
             sp_code = str_replace_all(sp_code, "ci025|ci975", replacement = "")) %>%
      left_join(., speciesKey, by = c("sp_code" = "species_code"))
  }) %>%
  rbindlist(use.names = TRUE)  %>%
  group_by(species_name, region, season, predLayer) %>%
  filter(n >= 200) %>%
  summarise(min_dist = min(surv_cov_dist), n = n[which.min(surv_cov_dist)], surv_cov_area = surv_cov_area[which.min(surv_cov_dist)]) %>%
  ungroup() %>%
  complete(species_name, region, season, predLayer) %>%
  mutate(min_dist = round(set_units(min_dist, "km"), 0), 
         surv_cov_area = round(set_units(surv_cov_area, "km^2"), 0), 
         appr_day_surv = round(min_dist/set_units(1600, "km/d"), 2)
         )




write_csv(lidarRes_df, path = "outputs/Lidar_summary Table.csv") 
  











# lidarPlots <- list()
# for(sp in speciesKey$species_code){
#   
#   pcounter <- 1
#   lidarPlots[[sp]] <- Lidar_results[str_which(names(Lidar_results), paste0(" ", sp))] %>%
#     map2(., names(.), function(x , y, regPoly = Forth_Tay_utm, spKey = speciesKey){
#       
#       x %<>% 
#         mutate(n = round(n, 0), 
#                surv_cov_dist = round(set_units(surv_cov_dist, "km"), 0),
#                n_text = paste0("n = ", n),
#                effort_text = paste0(surv_cov_dist ,"km")
#         )
#       
#       pTitle <- str_extract(y, pattern = "(?<= Season = )\\w+") %>%
#         str_to_title()
#       
#       sp_code <- str_extract(y, pattern = "(?<=Species = )\\w+(?=;)")
#       pMainTitle <- filter(spKey, species_code == sp_code) %>% pull(species_name)
#       
#       if(pcounter != 1){
#         pMainTitle <- " "
#       }
#       
#       p <- tm_shape(regPoly) +
#         tm_polygons(lwd = 2, col = "sienna1") +
#         tm_shape(x) +
#         tm_polygons(alpha = 0.6) +
#         tm_text("effort_text", size = 0.6, ymod = -0.4) +
#         tm_shape(x) +
#         #tm_text("subRegion_Id", size = 0.7, ymod = 0.5, fontface = "bold") +
#         tm_text("n", size = 0.6, ymod = 0, fontface = "bold") +
#         tm_shape(coastline_utm_cropped) +
#         tm_polygons(col = "navajowhite1", border.col = "burlywood4") +
#         tm_layout(title = pTitle, title.size = 1, title.position = c("right", "bottom"),
#                   main.title = pMainTitle, main.title.size = 1.2)
#       
#       pcounter <<- pcounter + 1
#       p
#     })
# }
# 
# 
# 
# 
# tmap_arrange(lidarPlots$RITRI, ncol = 3)
# tmap_arrange(lidarPlots$URAAL, ncol = 3)
# 



