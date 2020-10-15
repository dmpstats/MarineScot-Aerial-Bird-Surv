### -------------------------------------------------------------------------------------------------------------- ###
### -------------------------------------------------------------------------------------------------------------- ###
### ---                                                                                                        --- ###
### ---   Script to simulate survey in Wider Forth Tay for ALL species & ALL seasons - SeaPop 97.5% layer      --- ###
### ---                                                                                                        --- ###
### ---                                                                                                        --- ###
### -------------------------------------------------------------------------------------------------------------- ###
### -------------------------------------------------------------------------------------------------------------- ###



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


source("code/Main Analysis/MS_Seabird survey design_Main Analysis Tools.R")

survey_region_label <- "Forth Tay"
dens_surface_name <- "SeaPop_ci975"

machine_name <- "sauron"

outFolder <- "outputs/simulations/MS_Seabird survey design_Simulation outputs"


## ------------------------ ###
## ----   Upload data      ----
## ------------------------ ###

# upload survey region polygon
Forth_Tay <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Forth_Tay_survey_area/Wider_Forth_Tay_survey_area.shp")


# Upload density surfaces (for species of interest) for each season
pops_dens_summer <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Summer/Summer.shp") %>%
  select_at(.vars = vars(ID:Ocean, contains("ci975MOBAS"), contains("ci975RITRI"), contains("ci975URAAL"), contains("ci975ALTOR"), 
                         contains("ci975FRARC"), contains("ci975LAARG"), contains("ci975LAMAR")))

pops_dens_autumn <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Autumn/Autumn.shp") %>%
  select_at(.vars = vars(ID:Ocean, contains("ci975MOBAS"), contains("ci975RITRI"), contains("ci975URAAL"), contains("ci975ALTOR"), 
                         contains("ci975FRARC"), contains("ci975LAARG"), contains("ci975LAMAR")))

pops_dens_winter <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Winter/Winter.shp") %>%
  select_at(.vars = vars(ID:Ocean, contains("ci975MOBAS"), contains("ci975RITRI"), contains("ci975URAAL"), contains("ci975ALTOR"), 
                         contains("ci975FRARC"), contains("ci975LAARG"), contains("ci975LAMAR")))

pops_dens_List <- list(summer = pops_dens_summer, autumn = pops_dens_autumn, winter = pops_dens_winter)
rm(pops_dens_summer, pops_dens_autumn, pops_dens_winter)


# re-project survey region to density surface projections (important to be this way around for , as explained in 'survey_simulator')
Forth_Tay_utm <- st_transform(Forth_Tay, st_crs(pops_dens_List$summer))




## ------------------------------------------------- ###
## ----   Test Run (one species in one season)      ----
## ------------------------------------------------- ###

if(0){
  
  # -- set up
  spacing = 5000
  swath = 500
  spatial_units = "m"
  impact = 0
  
  # -- Select density for one species in one season and rename the density column as 'cell_N' (to meet formatting requirements)
  sp_dens_winter <- pops_dens_List$winter %>%
    select(ID:Ocean, ci975LAMAR) %>%
    rename(cell_N = ci975LAMAR)
  
  
  # -- check run - runs 1 iteration for QA and plots of example population and survey
  check_it <- TRUE
  test_sim_results <- survey_simulator(reps = 5, 
                                  region_polygon = Forth_Tay_utm, 
                                  sp_dens_surface = sp_dens_winter,
                                  trans_spacing = spacing, 
                                  platform_swath = swath, 
                                  spatial_units = spatial_units,
                                  check_it = check_it, 
                                  do_power = TRUE, 
                                  impact = impact)
  
  
  # -- test run
  check_it <- FALSE
  test_sim_results <- survey_simulator(reps = 5, 
                                  region_polygon = Forth_Tay_utm, 
                                  sp_dens_surface = sp_dens_winter,
                                  trans_spacing = spacing, 
                                  platform_swath = swath, 
                                  spatial_units = spatial_units,
                                  check_it = check_it, 
                                  do_power = TRUE, 
                                  impact = impact)
  
  # summaries for each simulated survey 
  test_sim_results$surveys_results
  
  # summaries across simulations
  test_sim_results$sim_summaries
  
}




## ------------------------------------------------------------ ###
## ----    Set up parameter grid and parallelisation           ----
## ------------------------------------------------------------ ###

# overall parameters
region = Forth_Tay_utm
swath = 500
spatial_units <- "m" # spatial units, for all distances provided bellow



# Variable parameters
impact <- seq(-0.5, 0.5, by = 0.1)
transectGap <- seq(500, 5500, by = 1000)
species_col_Name <- c("ci975ALTOR", "ci975FRARC", "ci975LAARG", "ci975LAMAR", "ci975MOBAS", "ci975RITRI", "ci975URAAL") 
season <- c("summer", "autumn", "winter")

# generate parameter grid
gridPars <- expand.grid(species_col_Name = species_col_Name, season = season, impact = impact, 
                        transectGap = transectGap, swath = swath) 
nrow(gridPars)

# generate simulation/scenario ID
simName <- gridPars %>% transmute(simName = paste0("Species = ", species_col_Name, "; Season = ", season, 
                                                   "; Impact = ", impact, "; transectGap = ", transectGap, 
                                                   ": Swath = ", swath))

# split as list, for purrr/furrr use
gridPars_ls <- gridPars %>%
  split(., simName)


# Parallelisation function 
surv_sim_parallel <- function(x, dens_surf, region_poly, spatial_units, survey_reps){

  # worker-specific parameters
  spacing = x$transectGap
  impact = x$impact
  species_col_Name = as.character(x$species_col_Name)
  season = x$season
  swath = x$swath
  
  
  # Select density for one species in one season and rename the density column as 'cell_N'
  sp_dens <- dens_surf[[season]] %>%
    select(ID:Y_COORD, matches(paste0("^", species_col_Name, "$"))) %>%
    rename(cell_N = matches(paste0("^", species_col_Name, "$")))
  
  sim_results <- survey_simulator(reps = survey_reps, 
                                  region_polygon = region_poly, 
                                  sp_dens_surface = sp_dens,
                                  trans_spacing = spacing, 
                                  platform_swath = swath, 
                                  spatial_units = spatial_units,
                                  check_it = FALSE, 
                                  do_power = TRUE, 
                                  impact = impact,
                                  progressBar = FALSE)
  
  # return summaries across simulations
  sim_results
}





## ------------------------------------------------------------ ###
## ----             RUN the Beast!                             ----
## ------------------------------------------------------------ ###

# set number of survey replicates per scenario
nSims = 100

# Set parallelization
plan(multiprocess)


# get starting time
start_time = Sys.time()

# Run all scenarios for all species
all_sim_results <- future_map(.x = gridPars_ls, .f = surv_sim_parallel, 
                              dens_surf = pops_dens_List, region_poly = region, 
                              spatial_units = spatial_units, 
                              survey_reps = nSims, .progress = TRUE)
# get end time
end_time = Sys.time()

  

# get simulation's run time
run_metrics <- tibble(machine_name, survey_region_label, dens_surface_name, start_time, 
                      end_time, elapsed_time = end_time - start_time)




## ------------------------------------------------------------ ###
## ----             Gather and save out results                ----
## ------------------------------------------------------------ ###

runOutputs <- list(all_sim_results = all_sim_results, 
                   run_metrics = run_metrics)


write_rds(runOutputs, path = paste0(outFolder, "_", survey_region_label, "_", dens_surface_name, ".rds"))



# # gather results in flat format
# power_results <- data.table::rbindlist(map(runOutputs$all_sim_results, ~ .$sim_summaries$sim_power_stats), idcol = "simName") ; power_results
# 
# N_results <- data.table::rbindlist(map(runOutputs$all_sim_results, ~ .$sim_summaries$sim_Nhat_stats), idcol = "simName"); N_results
#  
# survey_summaries <- data.table::rbindlist(map(runOutputs$all_sim_results, ~ .$sim_summaries$sim_surveys_stats), idcol = "simName") ; survey_summaries

