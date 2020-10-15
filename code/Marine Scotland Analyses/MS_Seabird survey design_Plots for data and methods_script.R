### -------------------------------------------------------------------- ###
### -------------------------------------------------------------------- ###
### ---                                                              --- ###
### ---              Script for data and methods plots               --- ###
### ---                                                              --- ###
### ---                                                              --- ###
### -------------------------------------------------------------------- ###
### -------------------------------------------------------------------- ###


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


outFolder <- "outputs/Plots and Tables/data maps/"

speciesKey <- tribble(~species_name, ~species_code,
                      "Atlantic Puffin", "FRARC",
                      "Black Legged Kittiwake", "RITRI",
                      "Common Guillemot", "URAAL",
                      "Herring Gull", "LAARG",
                      "Northern Fulmar", "FUGLA",
                      "Northern Gannet", "MOBAS",
                      "Razorbill", "ALTOR")


regionFullNames <- list(Herman_Berwick = "Herman. to Berwick", #"East of Hermaness to Berwick", 
                        east_Scotland = "East Scotland", #"East Coast Scotland", 
                        Forth_Tay = "Forth and Tay", #"Wider Forth and Tay", 
                        Moray_Firth = "Moray Firth")# "Wider Moray Firth")




## -------------------------- ###
## ----   Upload data        ----
## -------------------------- ###

# Upload SeaPop density surfaces (for species of interest) for each season

# Upload density surfaces (for species of interest) for each season
pops_dens_summer <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Summer/Summer.shp") %>%
  select_at(.vars = vars(ID:Ocean, matches("MOBAS"), matches("RITRI"), matches("URAAL"), matches("ALTOR"), 
                         matches("FRARC"), matches("LAARG"), matches("LAMAR")))

pops_dens_autumn <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Autumn/Autumn.shp") %>%
  select_at(.vars = vars(ID:Ocean, matches("MOBAS"), matches("RITRI"), matches("URAAL"), matches("ALTOR"), 
                         matches("FRARC"), matches("LAARG"), matches("LAMAR")))

pops_dens_winter <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Winter/Winter.shp") %>%
  select_at(.vars = vars(ID:Ocean, matches("MOBAS"), matches("RITRI"), matches("URAAL"), matches("ALTOR"), 
                         matches("FRARC"), matches("LAARG"), matches("LAMAR")))

SeaPop_dens_List <- list(summer = pops_dens_summer, autumn = pops_dens_autumn, winter = pops_dens_winter)
rm(pops_dens_summer, pops_dens_autumn, pops_dens_winter)


# Upload MERP density surfaces (for species of interest) for each season
MERP_dens_List <- read_rds("data/MERP_Monthly/MERP_fitted_allData_polys_seasonal_sf.rda")


# Upload Harbour Porpoise density surfaces
hp_dens <- read_rds("data/Porpoise_DensSurf_mean and ciBounds_sf.rds")



# upload survey region polygons
Moray_Firth <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Moray_Firth_surevey_area/Wider_Moray_Firth_survey_area.shp")
Herman_Berwick <- st_read(dsn = "data/East coast surveys Marine Scotland/East herman to Berwick/East herman to Berwick.shp")
east_Scotland <- st_read(dsn = "data/East coast surveys Marine Scotland/East_coast_Scotland/East_coast_Scotland.shp")
Forth_Tay <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Forth_Tay_survey_area/Wider_Forth_Tay_survey_area.shp")




# Upload Coastline
coastline <- read_sf("data/GBR_adm/GBR_adm0.shp") %>% select(NAME_LOCAL)
coastline



st_crs(SeaPop_dens_List$summer)
st_crs(MERP_dens_List$summer)
st_crs(hp_dens)






## --------------------------------------- ###
## ----     SeaPop Maps - by season       ----
## --------------------------------------- ###

#re-project survey region to SeaPop projections
seaPop_CRS <- st_crs(SeaPop_dens_List$summer)

Moray_Firth_utm33 <- st_transform(Moray_Firth, crs = seaPop_CRS)
Herman_Berwick_utm33 <- st_transform(Herman_Berwick, crs = seaPop_CRS)
east_Scotland_utm33 <- st_transform(east_Scotland, crs = seaPop_CRS)
Forth_Tay_utm33 <- st_transform(Forth_Tay, crs = seaPop_CRS)
coastline_utm33 <- st_transform(coastline, crs = seaPop_CRS)

# Gather region polygons into a grid
surveyRegions_ls <- list(Herman_Berwick = Herman_Berwick_utm33, Moray_Firth = Moray_Firth_utm33, 
     east_Scotland = east_Scotland_utm33, Forth_Tay = Forth_Tay_utm33)


# regionCols <- list(Herman_Berwick = "mediumorchid1", east_Scotland = "seagreen3", Forth_Tay = "red3", 
#                    Moray_Firth = "royalblue3")



seaPop_maps_seas <- map2(surveyRegions_ls, names(surveyRegions_ls), function(region_Poly, region_name){
  
  sp_code <- "RITRI"
  
  suppressWarnings({
    coastline_utm_cropped <- st_crop(coastline_utm33, bb(region_Poly, ext = 3))
    SeaPop_seas_cropped <- map(SeaPop_dens_List, ~st_crop(., y = bb(region_Poly, ext = 2.5)))
  })
  
  #browser()
  
  # get breaks for counts scale across all seasons - using the sqrt scale due to the large range of values between the seasons
  dens_brks <- trans_breaks("sqrt", function(x){x ^ 2})(c(SeaPop_seas_cropped$summer[[sp_code]], 
                                                          SeaPop_seas_cropped$autumn[[sp_code]], 
                                                          SeaPop_seas_cropped$winter[[sp_code]]))
  #dens_brks <- pretty_breaks(8)(c(SeaPop_dens_List$summer$RITRI, SeaPop_dens_List$autumn$RITRI, SeaPop_dens_List$winter$RITRI))
  
  seasons <- c("winter", "summer", "autumn")
  
  pcounter <- 1
  seaPop_p <- map(seasons, function(season, survRegion = region_Poly, nameRegion = region_name){
    
    #browser()
    # c_spName <- filter(speciesKey, species_code == sp_code) %>% pull(species_name)
    # cRegName <- regionFullNames[[nameRegion]]
    # pMainTitle <- paste0(cRegName," - ", c_spName)
    pMainTitle <- filter(speciesKey, species_code == sp_code) %>% pull(species_name)
    
    if(pcounter != 1){
      pMainTitle <- " "
      pLegend <- FALSE
    }else{
      pLegend <- TRUE
    }
    
    p <- SeaPop_seas_cropped[[season]] %>%
      tm_shape() +
      tm_polygons(sp_code, breaks = dens_brks, title = "Counts", palette = "YlGnBu") +
      tm_shape(survRegion, is.master = TRUE, bbox = bb(survRegion, ext = 1.4)) +
      #tm_polygons(col = regionCols[[nameRegion]], alpha = 0.2, lwd = 2) +
      tm_borders(lwd = 2, col = "black") +
      tm_shape(coastline_utm_cropped) +
      tm_polygons(col = "tan", border.col = "burlywood4") +
      tm_layout(title = str_to_title(season), title.size = 0.9, title.position = c("LEFT", "TOP"), #title.fontface = "bold",
                title.bg.color = "white", title.bg.alpha = 0.8,
                main.title = pMainTitle, main.title.size = 1.1,# main.title.fontface = "bold",
                legend.show =  pLegend, legend.bg.color = "white", legend.bg.alpha = 0.7, 
                legend.position = c("left", "bottom"))
    
    pcounter <<- pcounter + 1
    p
  })
  
  tmap_arrange(seaPop_p, ncol = 3)
  
})

seaPop_maps_seas$Forth_Tay

# save maps out
seaPop_maps_seas %>%
  walk2(., names(.), function(x, y){
    flPath <- paste0(outFolder, "SeaPop Maps_", y, "_by season_Black Legged Kittiwake.png")
    tmap_save(x, filename = flPath, dpi = 300, asp = 0, width = 12, height = 6, units = "in")
  })











## --------------------------------------------------------- ###
## ----     SeaPop Maps - by prediction layer (summer)      ----
## --------------------------------------------------------- ###


seaPop_maps_predLayer <- map2(surveyRegions_ls, names(surveyRegions_ls), function(region_Poly, region_name){
  
  sp_code <- "RITRI"
  
  suppressWarnings({
    coastline_utm_cropped <- st_crop(coastline_utm33, bb(region_Poly, ext = 3))
    SeaPop_summer_cropped <- st_crop(SeaPop_dens_List$summer, y = bb(region_Poly, ext = 2.5))
  })
  
#  browser()
  
  predLayers <- c(paste0("ci025", sp_code), sp_code, paste0("ci975", sp_code))
  
  
  # get breaks for counts scale across all seasons - using the sqrt scale due to the large range of values between the seasons
  dens_brks <- trans_breaks("sqrt", function(x){x ^ 2})(c(SeaPop_summer_cropped[[predLayers[1]]],
                                                          SeaPop_summer_cropped[[predLayers[2]]],
                                                          SeaPop_summer_cropped[[predLayers[3]]]))

  
  pcounter <- 1
  seaPop_p <- map(predLayers, function(predLayer, survRegion = region_Poly, nameRegion = region_name){
    
    pMainTitle <- filter(speciesKey, species_code == sp_code) %>% pull(species_name)
    
    if(pcounter != 1){
      pMainTitle <- " "
      pLegend <- FALSE
    }else{
      pLegend <- TRUE
    }
    
    subtitle <- case_when(
      str_detect(predLayer, "ci025") ~ "2.5%tile" ,
      str_detect(predLayer, "ci975") ~ "97.5%tile",
      str_detect(predLayer, paste0("^", sp_code, "$")) ~ "Mean"
    )

    p <- SeaPop_summer_cropped %>%
      tm_shape() +
      tm_polygons(predLayer, breaks = dens_brks, title = "Counts", palette = "YlGnBu") +
      tm_shape(survRegion, is.master = TRUE, bbox = bb(survRegion, ext = 1.4)) +
      #tm_polygons(col = regionCols[[nameRegion]], alpha = 0.2, lwd = 2) +
      tm_borders(lwd = 2, col = "black") +
      tm_shape(coastline_utm_cropped) +
      tm_polygons(col = "tan", border.col = "burlywood4") +
      tm_layout(title = subtitle, title.position = c("LEFT", "TOP"), title.size = 0.9, #title.fontface = "bold",
                title.bg.color = "white", title.bg.alpha = 0.8,
                main.title = pMainTitle, main.title.size = 1.1,# main.title.fontface = "bold",
                legend.show =  pLegend, legend.bg.color = "white", legend.bg.alpha = 0.7, 
                legend.position = c("left", "bottom"))
    
    pcounter <<- pcounter + 1
    p
  })
  
  tmap_arrange(seaPop_p, ncol = 3)
  
})


seaPop_maps_predLayer$Moray_Firth

# save maps out
seaPop_maps_predLayer %>%
  walk2(., names(.), function(x, y){
    flPath <- paste0(outFolder, "SeaPop Maps_", y, "_mean and ciBounds_Black Legged Kittiwake.png")
    tmap_save(x, filename = flPath, dpi = 300, asp = 0, width = 12, height = 6, units = "in")
  })








## --------------------------------------------------------- ###
## ----              MERP Maps - by season                 ----
## --------------------------------------------------------- ###

#re-project survey region to MERP projections
MERP_CRS <- st_crs(MERP_dens_List$summer)

Moray_Firth_utm30 <- st_transform(Moray_Firth, crs = MERP_CRS)
Herman_Berwick_utm30 <- st_transform(Herman_Berwick, crs = MERP_CRS)
east_Scotland_utm30 <- st_transform(east_Scotland, crs = MERP_CRS)
Forth_Tay_utm30 <- st_transform(Forth_Tay, crs = MERP_CRS)
coastline_utm30 <- st_transform(coastline, crs = MERP_CRS)

# Gather region polygons into a grid
surveyRegions_utm30_ls <- list(Herman_Berwick = Herman_Berwick_utm30, Moray_Firth = Moray_Firth_utm30, 
                         east_Scotland = east_Scotland_utm30, Forth_Tay = Forth_Tay_utm30)



MERP_maps_seas <- surveyRegions_utm30_ls %>%
  map2(., names(.), function(region_Poly, region_name){
  
  sp_code <- "RITRI"
  
  suppressWarnings({
    coastline_utm_cropped <- st_crop(coastline_utm30, bb(region_Poly, ext = 3))
    MERP_seas_cropped <- map(MERP_dens_List, ~st_crop(., y = bb(region_Poly, ext = 2.5)))
  })
  
  #browser()
  
  # get breaks for counts scale across all seasons - using the sqrt scale due to the large range of values between the seasons
  # dens_brks <- trans_breaks("sqrt", function(x){x ^ 2})(c(MERP_seas_cropped$summer[[sp_code]], 
  #                                                         MERP_seas_cropped$autumn[[sp_code]], 
  #                                                         MERP_seas_cropped$winter[[sp_code]]))
  
  dens_brks <- pretty_breaks(6)(c(MERP_seas_cropped$summer[[sp_code]], 
                                  MERP_seas_cropped$autumn[[sp_code]], 
                                  MERP_seas_cropped$winter[[sp_code]]))
  
  seasons <- c("winter", "summer", "autumn")
  
  pcounter <- 1
  MERP_p <- map(seasons, function(season, survRegion = region_Poly, nameRegion = region_name){
    
    pMainTitle <- filter(speciesKey, species_code == sp_code) %>% pull(species_name)
    
    if(pcounter != 1){
      pMainTitle <- " "
      pLegend <- FALSE
    }else{
      pLegend <- TRUE
    }
    
    p <- MERP_seas_cropped[[season]] %>%
      tm_shape() +
      tm_polygons(sp_code, breaks = dens_brks, title = "Counts", palette = "YlGnBu") +
      tm_shape(survRegion, is.master = TRUE, bbox = bb(survRegion, ext = 1.4)) +
      tm_borders(lwd = 2, col = "black") +
      tm_shape(coastline_utm_cropped) +
      tm_polygons(col = "tan", border.col = "burlywood4") +
      tm_layout(title = str_to_title(season), title.size = 0.9, title.position = c("LEFT", "TOP"), 
                title.bg.color = "white", title.bg.alpha = 0.8, #title.fontface = "bold",
                main.title = pMainTitle, main.title.size = 1.1,# main.title.fontface = "bold",
                legend.show =  pLegend, legend.bg.color = "white", legend.bg.alpha = 0.7, 
                legend.position = c("left", "bottom"))
    
    pcounter <<- pcounter + 1
    p
  })
  
  tmap_arrange(MERP_p, ncol = 3)
})


MERP_maps_seas$east_Scotland

# save maps out
MERP_maps_seas %>%
  walk2(., names(.), function(x, y){
    flPath <- paste0(outFolder, "MERP Maps_", y, "_by season_Black Legged Kittiwake.png")
    tmap_save(x, filename = flPath, dpi = 300, asp = 0, width = 12, height = 6, units = "in")
  })





## -------------------------------------------------------- ###
## ----    Harbour Porpoise Maps - by prediction layer     ----
## -------------------------------------------------------- ###

#re-project survey region to harbour purpose surface
hp_CRS <- st_crs(hp_dens)

Moray_Firth_utm31 <- st_transform(Moray_Firth, crs = hp_CRS)
Herman_Berwick_utm31 <- st_transform(Herman_Berwick, crs = hp_CRS)
east_Scotland_utm31 <- st_transform(east_Scotland, crs = hp_CRS)
Forth_Tay_utm31 <- st_transform(Forth_Tay, crs = hp_CRS)
coastline_utm31 <- st_transform(coastline, crs = hp_CRS)

# Gather region polygons into a grid
surveyRegions_utm31_ls <- list(Herman_Berwick = Herman_Berwick_utm31, Moray_Firth = Moray_Firth_utm31, 
                               east_Scotland = east_Scotland_utm31, Forth_Tay = Forth_Tay_utm31)




hp_maps_predLayer <- surveyRegions_utm31_ls %>%
  map2(., names(.), function(region_Poly, region_name){
  
  suppressWarnings({
    coastline_utm_cropped <- st_crop(coastline_utm31, bb(region_Poly, ext = 3))
    hp_cropped <- st_crop(hp_dens, y = bb(region_Poly, ext = 1.7))
  })
  
  #  browser()
  
  predLayers <- c("lower", "mean", "upper")
  
  
  # get breaks for counts scale across all seasons - using the sqrt scale due to the large range of values between the seasons
  # dens_brks <- trans_breaks("sqrt", function(x){x ^ 2}, n = 6)(c(hp_cropped$lower,
  #                                                         hp_cropped$upper,
  #                                                         hp_cropped$upper))

  dens_brks <- pretty_breaks(6)(c(hp_cropped$lower,
                                  hp_cropped$upper,
                                  hp_cropped$upper))
  
  
  pcounter <- 1
  hp_p <- map(predLayers, function(predLayer, survRegion = region_Poly, nameRegion = region_name){
    
    pMainTitle <- "Harbour Porpoise"
    
    if(pcounter != 1){
      pMainTitle <- " "
      pLegend <- FALSE
    }else{
      pLegend <- TRUE
    }
    
    subtitle <- case_when(
      str_detect(predLayer, "lower") ~ "2.5%tile" ,
      str_detect(predLayer, "upper") ~ "97.5%tile",
      str_detect(predLayer, "mean") ~ "Mean"
    )
    
    #browser()
    
    p <- hp_cropped %>%
      tm_shape() +
      tm_polygons(predLayer, breaks = dens_brks, title = "Counts", palette = "YlGnBu") +
      tm_shape(survRegion, is.master = TRUE, bbox = bb(survRegion, ext = 1.4)) +
      tm_borders(lwd = 2, col = "black") +
      tm_shape(coastline_utm_cropped) +
      tm_polygons(col = "tan", border.col = "burlywood4") +
      tm_layout(title = subtitle, title.position = c("LEFT", "TOP"), title.size = 0.9, #title.fontface = "bold",
                title.bg.color = "white", title.bg.alpha = 0.8,
                main.title = pMainTitle, main.title.size = 1.1,# main.title.fontface = "bold",
                legend.show = pLegend, legend.bg.color = "white", legend.bg.alpha = 0.7, 
                legend.position = c("left", "bottom"))
    
    pcounter <<- pcounter + 1
    p
  })
  
  tmap_arrange(hp_p, ncol = 3)
  
})


hp_maps_predLayer$Moray_Firth

# save maps out
hp_maps_predLayer %>%
  walk2(., names(.), function(x, y){
    flPath <- paste0(outFolder, "Harbour Porpoise Maps_", y, "_mean and ciBounds.png")
    tmap_save(x, filename = flPath, dpi = 300, asp = 0, width = 12, height = 6, units = "in")
  })








## --------------------------------------------------------------------------------------------- ###
## ----    Example of simulated population, transects and detections  - SeaPop, summer, mean    ----
## --------------------------------------------------------------------------------------------- ###

source("code/Main Analysis/MS_Seabird survey design_Main Analysis Tools.R")

# -- set up
spacing = 5000
swath = 500
spatial_units = "m"
impact = 0

spCode <- "RITRI"
spName <- filter(speciesKey, species_code == spCode)$species_name


# -- Select density for one species in one season and rename the density column as 'cell_N' (to meet formatting requirements)
sp_dens <- SeaPop_dens_List$summer %>%
  select(ID:Ocean, spCode) %>%
  rename(cell_N = spCode)


# -- Moray Firth
check_it <- TRUE
ex_sim_results_MorayFirth <- survey_simulator(reps = 5, 
                                     region_polygon = Moray_Firth_utm33, 
                                     sp_dens_surface = sp_dens,
                                     trans_spacing = spacing, 
                                     platform_swath = swath, 
                                     spatial_units = spatial_units,
                                     check_it = check_it, 
                                     do_power = TRUE, 
                                     impact = impact)


survInfo <- paste0("Swath = ", ex_sim_results_MorayFirth$sim_summaries$sim_specs$platform_swath, "m\n",
                   "Transect gap = ", ex_sim_results_MorayFirth$sim_summaries$sim_specs$transect_spacing, "m\n",
                   "Travel distance = ", round(set_units(ex_sim_results_MorayFirth$surveys_results$Survey_summary$covered_distance, "km"), 0), "km\n",
                   "Area covered = ", round(100*ex_sim_results_MorayFirth$surveys_results$Survey_summary$covered_prop, 1), "%\n",
                   "n = ", ex_sim_results_MorayFirth$surveys_results$Survey_summary$n
)


p1 <- ex_sim_results_MorayFirth$density_plot +
  tm_layout(main.title = paste0(spName, " - Summer"), main.title.size = 1.1, main.title.fontface = "plain")

p2 <- ex_sim_results_MorayFirth$survey_plot +
  tm_layout(main.title = " ", 
            title = survInfo, title.size = 0.8, title.position = c("LEFT", "BOTTOM"))

tmap_save(tmap_arrange(p1, p2), filename = paste0("outputs/Report Plots/Example simulated population and survey_Moray Firth_", spName, ".png"), 
          dpi = 300, asp = 0, width = 9, height = 6, units = "in")


          

# -- Forth Tay
check_it <- TRUE
ex_sim_results_Forth_Tay <- survey_simulator(reps = 5, 
                                              region_polygon = Forth_Tay_utm33, 
                                              sp_dens_surface = sp_dens,
                                              trans_spacing = 3000, 
                                              platform_swath = swath, 
                                              spatial_units = spatial_units,
                                              check_it = check_it, 
                                              do_power = TRUE, 
                                              impact = impact)


survInfo_Forth_Tay <- paste0("Swath = ", ex_sim_results_Forth_Tay$sim_summaries$sim_specs$platform_swath, "m\n",
                   "Transect gap = ", ex_sim_results_Forth_Tay$sim_summaries$sim_specs$transect_spacing, "m\n",
                   "Travel distance = ", round(set_units(ex_sim_results_Forth_Tay$surveys_results$Survey_summary$covered_distance, "km"), 0), "km\n",
                   "Area covered = ", round(100*ex_sim_results_Forth_Tay$surveys_results$Survey_summary$covered_prop, 1), "%\n",
                   "n = ", ex_sim_results_Forth_Tay$surveys_results$Survey_summary$n
)


p1 <- ex_sim_results_Forth_Tay$density_plot +
  tm_layout(main.title = paste0(spName, " - Summer"), main.title.size = 1.1, main.title.fontface = "plain")

p2 <- ex_sim_results_Forth_Tay$survey_plot +
  tm_layout(main.title = " ", 
            title = survInfo_Forth_Tay, title.size = 0.8, title.position = c("LEFT", "BOTTOM"))

tmap_save(tmap_arrange(p1, p2), filename = paste0("outputs/Report Plots/Example simulated population and survey_Forth Tay_", spName, ".png"), 
          dpi = 300, asp = 0, width = 9, height = 6, units = "in")

