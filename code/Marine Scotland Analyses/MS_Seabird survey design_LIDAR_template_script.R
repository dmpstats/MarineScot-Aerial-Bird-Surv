### ----------------------------------------------------------------------------------------- ###
### ----------------------------------------------------------------------------------------- ###
### ---                                                                                   --- ###
### ---     Script to simulate survey in a Forth Tay for ALL species & ALL seasons        --- ###
### ---                                                                                   --- ###
### ---                                                                                   --- ###
### ----------------------------------------------------------------------------------------- ###
### ----------------------------------------------------------------------------------------- ###



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

survey_region_label <- "Forth Tay"


## ------------------------ ###
## ----   Upload data      ----
## ------------------------ ###

# upload survey region polygon
Forth_Tay <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Forth_Tay_survey_area/Wider_Forth_Tay_survey_area.shp")


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
Forth_Tay_utm <- st_transform(Forth_Tay, st_crs(pops_dens_List$summer))



## ------------------------------------------------- ###
## ----   subset density surfaces for survey area    ----
## ------------------------------------------------- ###

pops_dens_List %<>%
  map(function(x, surveyPoly = Forth_Tay_utm){
    x[lengths(st_intersects(x, surveyPoly)) > 0, ]
  })


plot(pops_dens_List$winter, max.plot = 15)


## ------------------------------------------------- ###
## ----   Test Run (one species in one season)      ----
## ------------------------------------------------- ###

if(0){

  # -- set up
  swath = 200
  units(swath) <- "m"
  
  # -- Select density for one species in one season and rename the density column as 'cell_N' (to meet formatting requirements)
  sp_dens_seas <- pops_dens_List$summer %>%
    select(ID:Ocean, LAARG) %>%
    rename(cell_N = LAARG)
  
  
  # -- test run 
  LIDAR_results <- LIDAR_surveyRegions(swath = swath, track_gap = swath, region_polygon = Forth_Tay_utm, 
                                       sp_dens_surface = sp_dens_seas, target = 250, 
                                       clipToRegion = TRUE)
    
  
  LIDAR_results %<>% mutate(n = round(n, digits = 0))
  sp_dens_seas %<>% mutate(cell_N = round(cell_N, digits = 2))

  tm_shape(sp_dens_seas) +
    tm_polygons("cell_N") +
    tm_text("cell_N", size = 0.5, just = "top") +
    tm_shape(LIDAR_results) +
    tm_polygons(alpha = 0.7) +
    tm_text("n", ymod = 0.2)

  tm_shape(Forth_Tay_utm) +
    tm_polygons(lwd = 2, col = "burlywood1") +
    tm_shape(LIDAR_results) +
    tm_polygons(alpha = 0.7) +
    tm_text("n", size = 0.7) +
    tm_shape(LIDAR_results) +
    tm_text("subRegion_Id", size = 0.7, ymod = 0.8, fontface = "bold")
  
  
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

swath = 200
units(swath) <- "m"
targetAnimals = 250


# set-up progress bar
pb <- progress_bar$new(
  format = " Running LIDAR simulations [:bar] :current/:total (:percent) in :elapsed - eta: :eta",
  total = length(sp_dens_ls), clear = FALSE)  

#i <- 1

# run
Lidar_results <- sp_dens_ls %>%
  map(., function(x, strip = swath, survPoly = Forth_Tay_utm, target = targetAnimals){
    
    #print(i)
    
    out <- LIDAR_surveyRegions(swath = strip, track_gap = strip, region_polygon = survPoly, 
                        sp_dens_surface = x, target = target,  clipToRegion = FALSE)
    
    pb$tick()
    
    #i <<- i + 1
    
    return(out)
  })



Lidar_results
pryr::object_size(Lidar_results)


Lidar_results[[60]]
Lidar_results[[60]] %<>% mutate(n = round(n, digits = 0))

tm_shape(sp_dens_ls[[60]]) +
  tm_polygons("cell_N") +
  tm_text("cell_N", size = 0.5, just = "top") +
  tm_shape(Lidar_results[[60]]) +
  tm_polygons(alpha = 0.7) +
  tm_text("n", ymod = 0.2)




# undebug(LIDAR_surveyRegions)
# options(error = recover)



