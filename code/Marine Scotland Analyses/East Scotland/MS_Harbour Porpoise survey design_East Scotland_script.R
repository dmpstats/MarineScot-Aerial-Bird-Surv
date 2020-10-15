### -------------------------------------------------------------------------------------------------------------- ###
### -------------------------------------------------------------------------------------------------------------- ###
### ---                                                                                                        --- ###
### ---    Script to simulate survey in East Scotland for harbour porpoise - mean, lower and upper surfaces    --- ###
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
require(magrittr)

source("code/Main Analysis/MS_Seabird survey design_Main Analysis Tools.R")

survey_region_label <- "East Scotland"

machine_name <- "c5.18xlarge"

outFolder <- "outputs/simulations/MS_Porpoise survey design_Simulation outputs"


## ------------------------ ###
## ----   Upload data      ----
## ------------------------ ###

# upload survey region polygon
east_Scotland <- st_read(dsn = "data/East coast surveys Marine Scotland/East_coast_Scotland/East_coast_Scotland.shp")


# Upload density surfaces
hp_dens <- read_rds("data/Porpoise_DensSurf_mean and ciBounds_sf.rds")

# re-project survey region to density surface projections (important to be this way around for , as explained in 'survey_simulator')
east_Scotland_utm <- st_transform(east_Scotland, st_crs(hp_dens))





#' ## ---------------------------------------------------------------------------------------- ###
#' ## ----     Apply correction to densities for sighting availability due to diving          ----
#' ## ---------------------------------------------------------------------------------------- #### 
#' 
#' #' Predictions in original data were adjusted for availability, so we need to back-calculate to get animals 
#' #' actually available to the observers. 
#' #' Figures: mean surface time ~ 0.06 mins Vs mean dive time ~ 0.44 mins
#' diveCorrection <- 0.06/(0.06+0.44)
#' 
#' hp_dens %<>% 
#'   mutate(mean = diveCorrection * mean,
#'          lower = diveCorrection * lower,
#'          upper = diveCorrection * upper)



## ------------------------------------------------------------------------------------------- ###
## ----   Subset surfaces for survey region for increasing parallelisation performance       ----
## ------------------------------------------------------------------------------------------ ###

hp_dens <- hp_dens[lengths(st_intersects(hp_dens, east_Scotland_utm)) > 0, ]

plot(hp_dens)
hp_dens %>% 
  st_drop_geometry() %>%
  colSums()



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
  hp_dens_test <- hp_dens %>%
    select(mean) %>%
    rename(cell_N = mean)
  
  
  # -- check run - runs 1 iteration for QA and plots of example population and survey
  check_it <- TRUE
  test_sim_results <- survey_simulator(reps = 5, 
                                       region_polygon = east_Scotland_utm, 
                                       sp_dens_surface = hp_dens_test,
                                       trans_spacing = spacing, 
                                       platform_swath = swath, 
                                       spatial_units = spatial_units,
                                       check_it = check_it, 
                                       do_power = TRUE, 
                                       impact = impact)
  
  
  # -- test run
  check_it <- FALSE
  test_sim_results <- survey_simulator(reps = 5, 
                                       region_polygon = east_Scotland_utm, 
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
region = east_Scotland_utm
swath = 500
spatial_units <- "m" # spatial units, for all distances provided bellow



# Variable parameters
impact <- seq(-0.5, 0.5, by = 0.1)
transectGap <- seq(500, 5500, by = 1000)
predLayer <- c("mean", "lower", "upper") 

# generate parameter grid
gridPars <- expand.grid(predLayer = predLayer, impact = impact, transectGap = transectGap, swath = swath) 
nrow(gridPars)

# generate simulation/scenario ID
simName <- gridPars %>% transmute(simName = paste0("predLayer = ", predLayer, "; Impact = ", impact, 
                                                   "; transectGap = ", transectGap, ": Swath = ", swath))

# split as list, for purrr/furrr use
gridPars_ls <- gridPars %>%
  split(., simName)


# Parallelisation function 
surv_sim_parallel <- function(x, dens_surf, region_poly, spatial_units, survey_reps){
  
  # worker-specific parameters
  spacing = x$transectGap
  impact = x$impact
  predLayer = as.character(x$predLayer)
  swath = x$swath
  
  
  # Select density for one species in one season and rename the density column as 'cell_N'
  sp_dens <- dens_surf %>%
    select(matches(paste0("^", predLayer, "$"))) %>%
    rename(cell_N = matches(paste0("^", predLayer, "$")))
  
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
                              dens_surf = hp_dens, region_poly = region, 
                              spatial_units = spatial_units, 
                              survey_reps = nSims, .progress = TRUE)
# get end time
end_time = Sys.time()



# get simulation's run time
run_metrics <- tibble(machine_name, survey_region_label, start_time, 
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

