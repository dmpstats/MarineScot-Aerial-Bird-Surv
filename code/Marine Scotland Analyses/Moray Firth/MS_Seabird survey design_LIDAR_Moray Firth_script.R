### ------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------- ###
### ---                                                                                     --- ###
### ---     Script to simulate survey in a Moray Firth for ALL species & ALL seasons        --- ###
### ---                                                                                     --- ###
### ---                                                                                     --- ###
### ------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------- ###


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


source("code/Main Analysis/MS_Seabird survey design_Main Analysis Tools.R")

survey_region_label <- "Moray Firth"

dens_surface_name <- "SeaPop_fitted and ciBounds"

outFolder <- "outputs/simulations/MS_Seabird survey design_LiDAR_Simulation outputs"



## ------------------------ ###
## ----   Upload data      ----
## ------------------------ ###

# upload survey region polygon
Moray_Firth <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Moray_Firth_surevey_area/Wider_Moray_Firth_survey_area.shp")


# Upload density surfaces (for species of interest) for each season
pops_dens_summer <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Summer/Summer.shp") %>%
  select_at(.vars = vars(ID:Ocean, contains("MOBAS"), contains("RITRI"), contains("URAAL"), contains("ALTOR"), 
                         contains("FRARC"), contains("LAARG"), contains("LAMAR")))

pops_dens_autumn <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Autumn/Autumn.shp") %>%
  select_at(.vars = vars(ID:Ocean, contains("MOBAS"), contains("RITRI"), contains("URAAL"), contains("ALTOR"), 
                         contains("FRARC"), contains("LAARG"), contains("LAMAR")))

pops_dens_winter <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Winter/Winter.shp") %>%
  select_at(.vars = vars(ID:Ocean, contains("MOBAS"), contains("RITRI"), contains("URAAL"), contains("ALTOR"), 
                         contains("FRARC"), contains("LAARG"), contains("LAMAR")))

pops_dens_List <- list(summer = pops_dens_summer, autumn = pops_dens_autumn, winter = pops_dens_winter)
rm(pops_dens_summer, pops_dens_autumn, pops_dens_winter)


# re-project survey region to density surface projections (important to be this way around for , as explained in 'survey_simulator')
Moray_Firth_utm <- st_transform(Moray_Firth, st_crs(pops_dens_List$summer))


# upload data on species flying behaviour and compute percentage of birds flying 
sp_flyBehav <- read_csv("data/HiDef bird behaviour_digested.csv") %>%
  mutate(flying_pctg = n_flying/n_total )






## ------------------------------------------------- ###
## ----   subset density surfaces for survey area    ----
## ------------------------------------------------- ###

pops_dens_List %<>%
  map(function(x, surveyPoly = Moray_Firth_utm){
    x[lengths(st_intersects(x, surveyPoly)) > 0, ]
  })


plot(pops_dens_List$winter, max.plot = 15)



## ----------------------------------------------------- ###
## ----    Correct surfaces to flying animals only      ----
## ----------------------------------------------------- ###

# Correction based on HiDef estimates (Andy Webb, Pers. Comm., email subject "Birdy question", received on 24/06/2019)

sp_flyPctg <- sp_flyBehav %>% 
  select(species_code, flying_pctg) %>% 
  spread(species_code, flying_pctg)


pops_dens_List <- pops_dens_List %>%
  map(function(x){
    x %>%
      select_at(vars(ID, 
                     matches("^ALTOR$|ci\\d\\d\\dALTOR"),
                     matches("^FRARC$|ci\\d\\d\\dFRARC"), 
                     matches("^LAARG$|ci\\d\\d\\dLAARG"),
                     matches("^LAMAR$|ci\\d\\d\\dLAMAR"),
                     matches("^MOBAS$|ci\\d\\d\\dMOBAS"),
                     matches("^RITRI$|ci\\d\\d\\dRITRI"),
                     matches("^URAAL$|ci\\d\\d\\dURAAL")
      )) %>%
      mutate_at(vars(contains("ALTOR")), ~ .*sp_flyPctg$ALTOR) %>%
      mutate_at(vars(contains("FRARC")), ~ .*sp_flyPctg$FRARC) %>%
      mutate_at(vars(contains("LAARG")), ~ .*sp_flyPctg$LAARG) %>%
      mutate_at(vars(contains("LAMAR")), ~ .*sp_flyPctg$LAMAR) %>%
      mutate_at(vars(contains("MOBAS")), ~ .*sp_flyPctg$MOBAS) %>%
      mutate_at(vars(contains("RITRI")), ~ .*sp_flyPctg$RITRI) %>%
      mutate_at(vars(contains("URAAL")), ~ .*sp_flyPctg$URAAL)
    
  })


map(pops_dens_List, function(x){
  x %>% st_drop_geometry() %>% 
    select(c("ALTOR", "FRARC", "LAARG", "LAMAR", "MOBAS", "RITRI", "URAAL")) %>%
    colSums()
})





## ------------------------------------------------- ###
## ----   Test Run (one species in one season)      ----
## ------------------------------------------------- ###

if(0){
  
  # -- set up
  swath = 250
  units(swath) <- "m"
  
  # -- Select density for one species in one season and rename the density column as 'cell_N' (to meet formatting requirements)
  sp_dens_seas <- pops_dens_List$autumn %>%
    select(ID, RITRI) %>%
    rename(cell_N =  RITRI)
  
  
  # -- test run 
  LIDAR_results <- LIDAR_surveyRegions(swath = swath, track_gap = swath, region_polygon = Moray_Firth_utm, 
                                       sp_dens_surface = sp_dens_seas, target = 300, 
                                       clipToRegion = TRUE)
  
  
  LIDAR_results %<>% mutate(n = round(n, digits = 0))
  sp_dens_seas %<>% mutate(cell_N = round(cell_N, digits = 2))
  
  tm_shape(sp_dens_seas) +
    tm_polygons("cell_N") +
    tm_text("cell_N", size = 0.5, just = "top") +
    tm_shape(LIDAR_results) +
    tm_polygons(alpha = 0.7) +
    tm_text("n", ymod = 0.2)
  
  tm_shape(Moray_Firth_utm) +
    tm_polygons(lwd = 2, col = "burlywood1") +
    #tm_shape(LIDAR_results, is.master = TRUE) +
    tm_shape(LIDAR_results) +
    tm_polygons(alpha = 0.7) +
    tm_text("n", size = 0.7) +
    tm_shape(LIDAR_results) +
    tm_text("subRegion_Id", size = 0.7, ymod = 0.4, fontface = "bold")
  
  
  LIDAR_results
  
  
}  


#options(error = NULL)
#options(error = recover)




## ----------------------------------------- ###
## ----    Set up parameter grid            ----
## ----------------------------------------- ###

species_col_Name <- c("ALTOR", "ci025ALTOR", "ci975ALTOR", 
                      "FRARC", "ci025FRARC", "ci975FRARC", 
                      "LAARG", "ci025LAARG", "ci975LAARG", 
                      "LAMAR", "ci025LAMAR", "ci975LAMAR", 
                      "MOBAS", "ci025MOBAS", "ci975MOBAS", 
                      "RITRI", "ci025RITRI", "ci975RITRI", 
                      "URAAL", "ci025URAAL", "ci975URAAL") 

season <- c("summer", "autumn", "winter")


gridPars <- expand.grid(species_col_Name = species_col_Name, season = season) 
nrow(gridPars)
simName <- gridPars %>% transmute(simName = paste0("Species = ", species_col_Name, "; Season = ", season))


# split data for each species/season combo, for purrr/furrr use
sp_dens_ls <- gridPars %>%
  split(., simName) %>%
  map(function(x, dens_list = pops_dens_List){
    
    dens_list[[x$season]] %>%
      select(ID, matches(paste0("^", x$species_col_Name, "$"))) %>%
      rename(cell_N = matches(paste0("^", x$species_col_Name, "$")))
  })






## ------------------------------------------------------------ ###
## ----             RUN LIDAR simulation                       ----
## ------------------------------------------------------------ ###

swath = 250
units(swath) <- "m"
targetAnimals = 200


# set-up progress bar
pb <- progress_bar$new(
  format = " Running LIDAR simulations [:bar] :current/:total (:percent) in :elapsed - eta: :eta",
  total = length(sp_dens_ls), clear = FALSE)  

# run
Lidar_results <- sp_dens_ls %>%
  map(., function(x, strip = swath, survPoly = Moray_Firth_utm, target = targetAnimals){
    
    out <- LIDAR_surveyRegions(swath = strip, track_gap = strip, region_polygon = survPoly, 
                               sp_dens_surface = x, target = target,  clipToRegion = TRUE)
    
    pb$tick()
    
    return(out)
  })




Lidar_results
pryr::object_size(Lidar_results)

Lidar_results[[60]]
Lidar_results[[60]] %<>% mutate(n = round(n, digits = 0))

tm_shape(Moray_Firth_utm) +
  tm_polygons(lwd = 2, col = "burlywood1") +
  tm_shape(Lidar_results[[60]]) +
  tm_polygons(alpha = 0.7) +
  tm_text("n", size = 0.7) +
  tm_shape(Lidar_results[[60]]) +
  tm_text("subRegion_Id", size = 0.7, ymod = 0.7, fontface = "bold")




## ------------------------------------------------------------ ###
## ----               Write results out                        ----
## ------------------------------------------------------------ ###

write_rds(Lidar_results, path = paste0(outFolder, "_", survey_region_label, "_", dens_surface_name, ".rds"))


