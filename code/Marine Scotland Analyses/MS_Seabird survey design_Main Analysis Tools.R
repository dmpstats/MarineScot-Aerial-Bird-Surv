
### ---------------------------------------------------------------------------- ###
### ---------------------------------------------------------------------------- ###
### ---                                                                      --- ###
### ---     Functions for simulating survey design for main analysis         --- ###
### ---                                                                      --- ###
### ---                                                                      --- ###
### ---------------------------------------------------------------------------- ###
### ---------------------------------------------------------------------------- ###



draw_transects <- function(region_poly, spacing, units = "m"){

  #' ---------------------------------------------------------------------------
  #' Args:
  #'   - region_poly : sf object with the polygon for the survey region
  #'   - spacing     : distance between the parallel transect lines
  #'   - units       : length units for spacing
  #'
  #' Value:
  #'   sf object with lines of horizontal transects, generated based on 
  #'   systematic sampling (i.e. location of 1st transect drawn randomly 
  #'   between [0, spacing], and subsequent transects placed at fixed spacing)
  #' --------------------------------------------------------------------------
  
  require(sf)
  require(units)
  
  # unit consistency
  units(spacing) <- units
  region_units <- st_crs(region_poly)$units
  spacing <- set_units(spacing, region_units, mode = "standard")
    
  # get projection of survey region
  region_csr <- st_crs(region_poly)
  
  # get region boundaries
  limits <- st_bbox(region_poly)
  xlims <- c(limits[["xmin"]], limits[["xmax"]])
  ylims <- c(limits[["ymin"]], limits[["ymax"]])
  
  # randomly generate starting y point, between 0 and spacing value
  y_start <- runif(1, 0, spacing)
  
  # generate remaining y-coords, following a fixed spacing, i.e. systematic sampling
  ycoords <- seq(ylims[1] + y_start, ylims[2], by = as.numeric(spacing))
  
  # generate transects, as sf's linestring format
  transects <- list()
  for(i in 1:length(ycoords)){
    pts <- rbind(c(xlims[1], ycoords[i]), c(xlims[2], ycoords[i]))
    transects[[i]] <- st_linestring(pts)
  }
  
  # convert to sfc
  transects_sfc <- st_sfc(transects, crs = region_csr)
  
  # clip to area
  transects_sfc <- st_intersection(transects_sfc, region_poly)
  
  
  # convert to sf, adding transect's ID and length
  transects_sf <- st_sf(transect_id = 1:length(transects_sfc), #1:length(ycoords), 
                        transect_lt = st_length(transects_sfc),
                        geometry = transects_sfc)
  
  return(transects_sf)
}

# 
# #library(tmap)
# #debug(draw_transects)
# transects <- draw_transects(region_poly = Forth_Tay, spacing = 3, units = "km")
# tm_shape(Forth_Tay) +
#   tm_polygons(lwd = 3, border.col = "black", col = "grey92") +
#   tm_shape(transects) +
#   tm_lines(lwd = 1)







#'  
#' draw_population <- function(dens_surf, region_poly){
#' 
#'   #' ----------------------------------------------------------------------------------------
#'   #' Args:
#'   #'   - dens_surf   : sf object of a population's density surface. Function expects
#'   #'                  grid cells as polygons. MUST CONTAIN a column (field) named 'cell_N'
#'   #'                  giving the number of animals in each grid cell. For better performance,
#'   #'                  using only the cells intersecting the survey region is highly
#'   #'                  recommended.
#'   #'
#'   #'   - region_poly : sf object of the polygon for the survey region
#'   #'
#'   #' Value:
#'   #'   sf object with point coordinates of the randomly generated population in
#'   #'   the survey region
#'   #' ----------------------------------------------------------------------------------------
#' 
#'   if(!"cell_N" %in% colnames(dens_surf)){
#'     stop("\"dens_surf\" object MUST contain field (column) named \"cell_N\" with number of animals in grid cell" )
#'   }
#' 
#'   # draw numbers in each surface's polygon, based on a Poisson spatial process
#'   dens_surf_sim <- dens_surf %>% mutate(cell_N_sim = rpois(length(cell_N), cell_N))
#' 
#'   # spread animals in each's surface polygon
#'   obj_pts_surface <- st_sample(dens_surf_sim, size = dens_surf_sim$cell_N_sim, exact = TRUE)
#' 
#'   # clip to the survey region
#'   obj_pts_region <- st_intersection(obj_pts_surface, region_poly)
#' 
#'   st_sf(obj_ID = 1:length(obj_pts_region), geometry = obj_pts_region)
#' }
#' 
#' 
#
# pop_dens <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Autumn/Autumn.shp") %>%
#   st_transform(crs = st_crs(Forth_Tay)) %>%
#   select(ID:Ocean, MOBAS) %>%
#   rename(cell_N = MOBAS)
#
# pop_dens_region <- pop_dens[lengths(st_intersects(pop_dens, Forth_Tay)) > 0, ]
#
# pop_locations <- draw_population(dens_surf = pop_dens_region, region_poly = Forth_Tay)
#
# tm_shape(pop_dens_region) +
#   tm_polygons(col = "cell_N", palette = "YlGnBu") +
#   tm_shape(Forth_Tay) +
#   tm_polygons(alpha = 0.5, lwd = 2, border.col = "black") +
#   tm_shape(pop_locations) +
#   tm_dots(size = 0.08)
 



sprinkle_Polys <- function(densSurf_poly, size){
  
  pts_sfc <- list()
  densSurf_crs <- st_crs(densSurf_poly)

  for(j in 1:nrow(densSurf_poly)){

    pts_list <- list()

    poly_flat <- densSurf_poly[j, "geometry"]$geometry[[1]][[1]]
    xlims <- c(min(poly_flat[,1]), max(poly_flat[,1]))
    ylims <- c(min(poly_flat[,2]), max(poly_flat[,2]))

    pts_x <- runif(size[j], xlims[1], xlims[2])
    pts_y <- runif(size[j], ylims[1], ylims[2])
    pts <- matrix(c(pts_x, pts_y), byrow = FALSE, ncol = 2)

    if(nrow(pts)>0){
      for(i in 1:nrow(pts)){ #i <- 1
        pts_list[[i]] <- st_point(pts[i, ])
      }
      pts_sfc[[j]] <- st_sf(a = 1, geometry = st_sfc(pts_list), crs = densSurf_crs)  
    }
    
  }

  if(length(pts_sfc) > 0){
    output <- st_as_sf(data.table::rbindlist(pts_sfc))
  }else{
    output <- st_sf(st_sfc())  # empty sf object
  }
  
  return(output)
}



draw_population <- function(dens_surf, region_poly){

  #' ---------------------------------------------------------------------------------------------
  #' Faster version (~3x faster) from using 'sprinkle_Polys' instead of sf's 'st_sample' to spread 
  #' animals in each grid cell of the density surface. HOWEVER, 'sprinkle_Polys_utm' is accurate 
  #' for only rectangular-shaped polygons for grid cells.
  #' 
  #' Args:
  #'   - dens_surf  : sf object of a population's density surface. Fuction expects grid cells 
  #'                  as polygons. MUST CONTAIN a column named 'cell_N' for the number of animals 
  #'                  in each grid cell. For better performance, using only the cells intersecting the 
  #'                   survey region is highly recommended
  #'                      
  #'   - region_poly : sf object of the polygon for the survey region
  #'
  #' Value:
  #'   sf object of point coordinates of simulated population in the survey region, 
  #'   in UTM projection
  #' ---------------------------------------------------------------------------------------------

  # check if survey region and density surface are in the same projection
  if(st_crs(dens_surf) != st_crs(region_poly)){
    stop("'Ahem... 'dens_surf' and 'region_poly' need to be in the same spatial projection.")
  }
  
  if(!"cell_N" %in% colnames(dens_surf)){
    stop("\"dens_surf\" object MUST contain field (column) named \"cell_N\" with number of animals in grid cell" )
  }

  # draw population based on a Poisson spatial process
  dens_surf_sim <- dens_surf %>% mutate(cell_N_sim = rpois(length(cell_N), cell_N))

  # spread animals in each surface polygon
  obj_pts_surface <- sprinkle_Polys(dens_surf, size = dens_surf_sim$cell_N_sim)

  if(nrow(obj_pts_surface) > 0){
    # select points inside the region
    obj_pts_region <- obj_pts_surface[lengths(st_intersects(obj_pts_surface, region_poly)) > 0, ] %>%
      mutate(obj_ID = seq_along(geometry)) %>%
      select(obj_ID, geometry)
    
    out <- obj_pts_region
  }else{
    out <- st_sf(st_sfc())  # empty sf object
  }

  return(out)
}


# pop_dens <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Autumn/Autumn.shp") %>%
#   select(ID:Ocean, MOBAS) %>%
#   rename(cell_N = MOBAS)
# 
# Forth_Tay_utm <- st_transform(Forth_Tay, st_crs(pop_dens))
# pop_dens_region <- pop_dens[lengths(st_intersects(pop_dens, Forth_Tay_utm)) > 0, ]
# 
# pop_locations<- draw_population(dens_surf = pop_dens_region, region_poly = Forth_Tay_utm)
# 
# tm_shape(pop_dens_region) +
#   tm_polygons(col = "cell_N", palette = "YlGnBu") +
#   tm_shape(Forth_Tay_utm) +
#   tm_polygons(alpha = 0.5, lwd = 2, border.col = "black") +
#   tm_shape(pop_locations) +
#   tm_dots(size = 0.08)
# 
# 
# pop_locations <- st_transform(pop_locations, st_crs(Forth_Tay))
# 
# tm_shape(pop_dens_region) +
#   tm_polygons(col = "cell_N", palette = "YlGnBu") +
#   tm_shape(Forth_Tay_utm) +
#   tm_polygons(alpha = 0.5, lwd = 2, border.col = "black") +
#   tm_shape(pop_locations) +
#   tm_dots(size = 0.08)



get_surveyed_pop <- function(pop_coords, transects, swath, swath_unit = "m"){
  
  #' ----------------------------------------------------------------------------
  #' Args:
  #'   - pop_coords   : sf object of point coordinates of animals
  #'   - transects    : sf object of transect lines
  #'   - swath        : numeric, search strip width of the surveying platform
  #'   - swath_unit   : character, length units for swath
  #'
  #' Value:
  #'   sf object with point coordinates of animals covered by the transects
  #' ---------------------------------------------------------------------------
  
  require(sf)
  require(units)
  
  # units consistency
  units(swath) <- swath_unit
  trans_projUnits <- st_crs(transects)$units
  swath <- set_units(swath, trans_projUnits, mode = "standard")

  # Generate polygons for transects coverage based on swath
  swath_buffer <- as.numeric(swath)/2
  transects_buffer <- st_buffer(transects, dist = swath_buffer, endCapStyle="FLAT")
  
  # get animals covered by the transects
  pop_trans_intc <- st_intersects(pop_coords, transects_buffer) %>% as_tibble()
  pop_surveyed <- pop_coords[pop_trans_intc$row.id, ]
  
  # get perpendicular distance between each animal and closest transect
  obs_trans_dist <- st_distance(pop_surveyed, transects) %>% apply(. , 1, min)
  
  pop_surveyed <- pop_surveyed %>%
    mutate(transect_id = pop_trans_intc$col.id,
           dist_to_transect = obs_trans_dist)
  
  return(pop_surveyed)
}


# pop_surveyed <- get_surveyed_pop(pop_coords = pop_locations, transects = transects, swath = 300, swath_unit = "m")
# 
# tm_shape(transects) +
#   tm_lines(lwd = 1, col = "black") +
#   tm_shape(Forth_Tay) +
#   tm_polygons(alpha = 0.5, lwd = 2, border.col = "black") +
#   tm_shape(pop_locations) +
#   tm_dots(size = 0.07, col = "firebrick1", alpha = 0.2) +
#   tm_shape(pop_surveyed) +
#   tm_dots(size = 0.07, col = "firebrick1")
  




N_hat_bootstrap <- function(nBootRep, transect_obs, region_area, swath, swath_unit){
  
  #' --------------------------------------------------------------------------------
  #' Args:
  #'   - nBootRep        : number of bootstrap samples
  #'   - transect_obs    : sf object of transect lines with number of animals covered.
  #'                       Transects with no animals detected must be included. Must 
  #'                       contain two collumns: 
  #'                          (i)  'transect_lt' (transect lengths)
  #'                          (ii) 'transect_n' (number of animals)
  #'   - region_area     : numeric object of class "units", Area of the surveyed region
  #'   - swath           : numeric, search strip width of the surveying platform
  #'   - swath_unit      : character, length units for swath
  #'
  #' Value:
  #'   boostrap distribution of N_hat = estimator of population abundance in survey region
  #' ------------------------------------------------------------------------------
  
  if(!"transect_lt" %in% colnames(transect_obs)){
    stop("\"transect_obs\" object MUST contain column named \"transect_lt\" for length of each transect")
  }
  
  if(!"transect_n" %in% colnames(transect_obs)){
    stop("\"transect_obs\" object MUST contain column named \"transect_n\" for number of animals 
         in each transect")
  }
  
  # units consistency
  units(swath) <- swath_unit
  trans_lt_projUnits <- units(transect_obs$transect_lt)
  swath <- set_units(swath, trans_lt_projUnits, mode = "standard") # conversion, if non-matching units
  
  # drop units, for performance (massive improvement!)
  transects_lt <- as.numeric(transect_obs$transect_lt)
  region_area <- as.numeric(region_area)
  swath <- as.numeric(swath)
  
  # initialize bootstrap holder object
  N_hat_boot <- rep(0, nBootRep)
  
  # run bootstrap for N_hat
  for(i in 1:nBootRep){
    boot_idx <- sample.int(length(transects_lt), replace = TRUE)
    n_boot <- sum(transect_obs$transect_n[boot_idx])
    N_hat_boot[i] <- region_area/sum(transects_lt[boot_idx] * swath) * n_boot
  }
  
  return(N_hat_boot)
}



# transects_obs <- full_join(data.frame(pop_surveyed), data.frame(transects), by = "transect_id") %>%
#   group_by(transect_id) %>%
#   summarise(transect_n = sum(!is.na(obj_ID)), transect_lt = first(transect_lt))
# 
# 
# N_hat_dst <- N_hat_bootstrap(nBootRep = 10000, transect_obs = transects_obs, region_area = st_area(Forth_Tay), 
#                         swath = 300, swath_unit = "m")
# 
# mean(N_hat_dst)
# var(N_hat_dst)
# quantile(N_hat_dst, probs = c(0.025, 0.975))
# 
# N <- nrow(pop_locations); N




survey_simulator <- function(reps = 100, region_polygon, sp_dens_surface, trans_spacing, 
                             platform_swath, spatial_units, 
                             check_it = FALSE, do_power = FALSE, impact = 0, progressBar = TRUE){
  
  #' -----------------------------------------------------------------------------------------
  #' Args:
  #'   - reps                 : number of bootstrap samples
  #'   - region_polygon       : sf object of the polygon for the survey region. Must be in the same spatial 
  #'                            projection as sp_dens_surface
  #'   - sp_dens_surface      : sf object of geometry type POLYGON, for the species' density surface. 
  #'                            Fuctions expect polygons provided for each grid cell. Polygons ideally 
  #'                            rectangular-shaped, otherwise the spread of spatial points is only an 
  #'                            approximation of the density surface. MUST CONTAIN a column named 'cell_N' for 
  #'                            the number of animals in each grid cell
  #'   - trans_spacing        : numeric, the distance between the parallel survey transect lines
  #'   - platform_swath       : numeric, search strip width of the surveying platform
  #'   - spatial_units        : character, overal spatial units (e.g. "m")
  #'   - check_it             : logical, QA - show plots and outputs for one simulated population and survey
  #'   - do_power             : logical, should power analysis be performed?
  #'   - impact               : numeric, impact to apply to the reference population. Expressed in terms of
  #'                            the percentual change (decimal) in the population. Negative values indicate a
  #'                            decline, and vice-versa.
  #'
  #'
  #' Value:
  #'   List with elements:
  #'   (i) 'surveys_results' : for each simulated survey, surveying summaries and statistics for 
  #'                           population estimates (N = abundance, D = density) and power metrics (if 
  #'                           'do_power = TRUE')
  #'   (ii) 'sim_summaries'  : summaries and statistics for surveys, population estimates and power metrics across 
  #'                           all simulations
  #'   (iii) '*_plot'        : QA plots if 'check_it = TRUE'
  #' ------------------------------------------------------------------------------------------
  
  suppressPackageStartupMessages({
    require(sf)
    require(tmap)
    require(units)
    require(tidyverse)
    require(progress)
    require(data.table)
    require(scales)  
  })
  
  
  #browser()
  
  if(check_it) reps <- 1 
  
  # check if survey region and density surface are in the same projection
  if(st_crs(region_polygon) != st_crs(sp_dens_surface)){
    stop("'Ahem... region_polygon' and 'sp_dens_surface' need to be in the same spatial projection.")
  }
   
  # check if density surface data has the key column 'cell_N'
  if(!"cell_N" %in% colnames(sp_dens_surface)){
    stop("Missing \'cell_N\' column - \'sp_dens_surface\' object MUST contain field (column) named \'cell_N\' providing number of animals in each grid cell")
  }
  
  #' subset the density surface for the surveying area (i.e. returns density polygons overlapped by the 
  #' regions - no clipping involved)
  dens_region_ref <- sp_dens_surface[lengths(st_intersects(sp_dens_surface, region_polygon)) > 0, ]
  
  # Check if density surface and survey region do overlap
  if(nrow(dens_region_ref) == 0){
    stop("OOps! 'sp_dens_surface' and 'region_polygon' do not overlap spatially")
  }
  
  
  # get the area of the survey region
  regionArea <- st_area(region_polygon)
  
  
  # when performing power analysis, get abundance of reference population **in the survey region** & apply impact
  if(do_power){
    
    #' assuming uniform spatial distn inside cell, numbers in each cell adjusted for the fraction of the cell 
    #' overlapped by the region
    suppressWarnings(
      dens_region_ref_clipped <- st_intersection(dens_region_ref, region_polygon) %>%
        mutate(area_poly = st_area(geometry), 
               #prop_area_poly = area_poly/set_units(as_units(10*10, "km^2"), "m^2"),    # WARNING!! resolution specific to SeaPop - need to be changed to be defined upfront by user
               prop_area_poly = area_poly/max(area_poly), 
               cell_N = cell_N * prop_area_poly)
    )
    
    N_ref <- sum(dens_region_ref_clipped$cell_N)
    
    # apply impact
    dens_region <- dens_region_ref %>%
      mutate(cell_N = cell_N + cell_N * impact)

    # impact contrast plots  
    if(check_it){
      
      dens_brks <- pretty_breaks(6)(c(dens_region$cell_N, dens_region_ref$cell_N))
      
      refDens_plot <- tm_shape(dens_region_ref) +
        tm_polygons(col = "cell_N", palette = "YlGnBu", title = "Counts", breaks = dens_brks) +
        tm_shape(region_polygon) +
        tm_borders(lwd = 2, col = "black") +
        tm_layout(asp = 1) +
        tm_layout(main.title = "Density surface: Reference",
                         main.title.size = 0.9, main.title.fontface = "bold")
        
      impDens_plot <- tm_shape(dens_region) +
        tm_polygons(col = "cell_N", palette = "YlGnBu", title = "Counts", breaks = dens_brks) +
        tm_shape(region_polygon) +
        tm_borders(lwd = 2, col = "black") +
        tm_layout(asp = 1) +
        tm_layout(main.title = paste0("Density surface: impact = ", impact),
                  main.title.size = 0.9, main.title.fontface = "bold")
        
    }
  }else{  # if not performing power analysis
    dens_region <- dens_region_ref
  }
  
  
  if(progressBar){
    # set-up progress bar
    pb <- progress_bar$new(
      format = " Running survey simulations [:bar] :current/:total (:percent) in :elapsed - eta: :eta",
      total = reps, clear = FALSE)  
  }
  
  
  
  # initialize lists holding simulations
  N <- D <- Survey_summary <- power_metrics <- list()
  
  #browser()
  
  for(i in 1:reps){  # start of loop over survey replicates
    
    if(progressBar){
      pb$tick()
    }
    
    # draw the true population, with animal positions inside the survey region
    pop_points <- draw_population(dens_surf = dens_region, region_poly = region_polygon)
    
    # skip current iteration if no animals present in the region
    if(nrow(pop_points) == 0){ 
      next 
    }
    
    
    # draw survey transects
    transects <- draw_transects(region_poly = region_polygon, spacing = trans_spacing, 
                                units = spatial_units)
    
    # Get surveyed population
    pop_surveyed <- get_surveyed_pop(pop_coords = pop_points, transects = transects, swath = platform_swath, 
                                     swath_unit = spatial_units)
    
    # Get abundance bootstrap distribution
    transects_obs <- full_join(data.frame(pop_surveyed), data.frame(transects), by = "transect_id") %>%
      group_by(transect_id) %>%
      summarise(transect_n = sum(!is.na(obj_ID)), transect_lt = first(transect_lt))
    
    N_hat_dst <- N_hat_bootstrap(nBootRep = 10000, transect_obs = transects_obs, 
                                 region_area = regionArea, swath = platform_swath, 
                                 swath_unit = spatial_units)
    
    # calculate & store abundance estimates for current survey
    N[[i]] <- tibble(
      survey_sim     = i,
      estimate       = mean(N_hat_dst),
      se             = sd(N_hat_dst),
      lcl            = quantile(N_hat_dst, 0.025),
      ucl            = quantile(N_hat_dst, 0.975),
      truth          = nrow(pop_points),
      bias           = estimate - truth,
      CI_cover_truth = ifelse(truth >= lcl & truth <= ucl, 1, 0) )
    
    
    # Convert to density estimates & store
    D_hat_dist <- N_hat_dst/regionArea
    D[[i]] <- tibble(
      survey_sim      = i,
      estimate        = mean(D_hat_dist),
      se              = sd(D_hat_dist),
      lcl             = quantile(D_hat_dist, 0.025),
      ucl             = quantile(D_hat_dist, 0.975),
      truth           = nrow(pop_points)/regionArea,
      bias            = estimate - truth,
      CI_cover_truth  = ifelse(truth >= lcl & truth <= ucl, 1, 0) )
    
    # get & store current survey summary
    Survey_summary[[i]] <- tibble(
      survey_sim       = i,
      covered_distance = sum(transects$transect_lt),
      covered_area     = covered_distance * set_units(platform_swath, spatial_units, mode = "standard"),
      covered_prop     = covered_area/regionArea,
      n                = nrow(pop_surveyed) )
    
    
    if(do_power){
      # get & store power metrics of current survey
      power_metrics[[i]] <- N[[i]] %>%
        select(estimate, lcl, ucl, truth) %>%
        rename(N_hat = estimate, N_lcl = lcl, N_ucl = ucl, N_truth = truth) %>%
        mutate(N_reference = as.numeric(N_ref),
               impact = impact,
               impactDetected = if_else(N_reference>=N_lcl & N_reference<=N_ucl, FALSE, TRUE))
    }
    
    
    if(check_it){
      
      density_plot <- tm_shape(dens_region) +
        tm_polygons(col = "cell_N", palette = "YlGnBu", title = "Counts", breaks = dens_brks) +
        tm_shape(region_polygon) +
        tm_polygons(alpha = 0.5, lwd = 2, border.col = "black") +
        tm_shape(pop_points) +
        tm_dots(size = 0.05, alpha = 0.6) +
        tm_layout(main.title = "(True) Density surface & simulated animal locations", legend.title.size = 1,
                  main.title.size = 0.9, main.title.fontface = "bold")
        
      survey_plot <- tm_shape(region_polygon) +
        tm_polygons(alpha = 0.5, lwd = 2, border.col = "black") +
        tm_shape(transects) +
        tm_lines(lwd = 2, col = "dodgerblue2") +
        tm_shape(pop_points) +
        tm_dots(size = 0.05, col = "firebrick1", alpha = 0.2) +
        tm_shape(pop_surveyed) +
        tm_dots(size = 0.07, col = "firebrick1") +
        tm_layout(main.title = "Survey region, simulated transects & surveyed animals",
                  main.title.size = 0.9, main.title.fontface = "bold")
      
      #browser()
      
      if(do_power){
        check_plot <- tmap_arrange(refDens_plot, impDens_plot, density_plot, survey_plot, ncol = 2, asp = 1)
                                   
      }else{
        check_plot <- tmap_arrange(density_plot, survey_plot, asp = 1)  
      }
      
      #windows(width = 10, height = 10)
      print(check_plot)
      #check_plot
    }
    
  } # end of loop over survey replicates
  
  #browser()
  
  # Gather survey-level simulation results (using data.table::rbindlist to keep units after binding)
  N <- rbindlist(N) %>% as_tibble()
  D <- rbindlist(D) %>% as_tibble()
  Survey_summary <- rbindlist(Survey_summary) %>% as_tibble()
  
  
  # Summarise results 
  sim_surveys_stats <- Survey_summary %>% 
    summarise(
      region_area     = regionArea,
      mean_cvrg_area  = mean(covered_area), 
      mean_cvrg_dist  = mean(covered_distance),
      mean_cvrg_prop  = mean(covered_prop),
      mean_n          = mean(n),
      No_zero_n       = sum(n == 0),
      mean_ER         = mean_n/mean_cvrg_dist)
  
  sim_Nhat_stats <- N %>% 
    summarise(
      mean_trueN    = mean(truth), 
      mean_Nhat     = mean(estimate),
      mean_Nhat_se  = mean(se),
      percent_bias  = mean(bias/truth * 100),
      CI_cvrg_prop  = sum(CI_cover_truth)/nrow(.))
  
  sim_Dhat_stats <- D %>% 
    summarise(
      mean_trueD    = mean(truth), 
      mean_Dhat     = mean(estimate),
      mean_Dhat_se  = mean(se),
      percent_bias  = mean(bias/truth * 100),
      CI_cvrg_prop  = sum(CI_cover_truth)/nrow(.))
  
  sim_specs <- tibble(
    platform_swath   = set_units(platform_swath, spatial_units, mode = "standard"),
    transect_spacing = set_units(trans_spacing, spatial_units, mode = "standard"),
    No_repetitions   = reps
  )
  
  #browser()
  
  if(do_power){
    power_metrics <- rbindlist(power_metrics) %>% as_tibble()
    
    sim_power_stats <- power_metrics %>% 
      summarise(
        impact         = first(impact),
        imp_mean_Nhat  = mean(N_hat),
        imp_mean_N_lcl = mean(N_lcl),
        imp_mean_N_ucl = mean(N_ucl),
        imp_mean_N_truth = mean(N_truth),
        N_reference      = first(N_reference),
        mean_cvrg_prop   = sim_surveys_stats$mean_cvrg_prop,
        power            = sum(impactDetected)/nrow(.))
  }
  
  
  # Arrange outputs
  outputs <- list(
    surveys_results = list(N = N, D = D, Survey_summary = Survey_summary),
    sim_summaries = list(sim_specs = sim_specs, 
                         sim_surveys_stats = sim_surveys_stats, 
                         sim_Nhat_stats = sim_Nhat_stats,
                         sim_Dhat_stats = sim_Dhat_stats))
  
  if(check_it){
    outputs[["density_plot"]] <- density_plot  
    outputs[["survey_plot"]] <- survey_plot
    
    if(do_power){
      outputs[["refDens_plot"]] <- refDens_plot
      outputs[["impDens_plot"]] <- impDens_plot
    }}
  
  if(do_power){
    outputs[["surveys_results"]][["power_metrics"]] <- power_metrics
    outputs[["sim_summaries"]][["sim_power_stats"]] <- sim_power_stats
  }
  
  return(outputs)
}







break_region_target <- function(region_polygon, sp_dens_surface, target, clipToRegion = FALSE){
  
  #' -----------------------------------------------------------------------------------------
  #' Purpose : 
  #'    Breaks down the survey region into sub-regions comprising a chosen target number of animals
  #' 
  #' Args:
  #'   - region_polygon       : sf object of the polygon for the survey region. Must be in the same spatial 
  #'                            projection as sp_dens_surface
  #'   - sp_dens_surface      : sf object of geometry type POLYGON, for the species' density surface.
  #'                            Fuction expect polygons provided for each grid cell. Polygons ideally 
  #'                            rectangular-shaped. MUST CONTAIN a column named 'cell_N' for 
  #'                            the number of animals in each grid cell
  #'   - target               : integer, the target number of animals in each sub-region
  #'   
  #'   - clipToRegion         : logical, should the density surface be clipped to survey polygon and number 
  #'                          of animals in cells adjusted for proportion of cell covered by survey region?
  #'                          Note: shape of output polygons might look a bit funky when clipping is chosen
  #'
  #' Value:
  #'   sf object with polygons of sub-regions inside the survey region, each containing (at minimum) the 
  #'   number of target animals
  #' ------------------------------------------------------------------------------------------
  
  # set the density surface for the cells intersecting the survey region
  sp_dens_region <- sp_dens_surface[lengths(st_intersects(sp_dens_surface, region_polygon)) > 0, ] %>%
    mutate(cell_Id = 1:n())
  
  #browser()
  
  if(clipToRegion){
    #' assuming uniform spatial distn inside cell, numbers in each cell adjusted for the proportion of the cell 
    #' overlapped by the region
    suppressWarnings(
      sp_dens_init <- st_intersection(sp_dens_region, region_polygon) %>%
        mutate(cell_area = st_area(geometry), 
               cell_prop_cov =  as.numeric(cell_area/max(cell_area)),
               cell_N = cell_N * cell_prop_cov)
    )  
  }else{
    sp_dens_init <- sp_dens_region
  }
  
  # compute total number of animals in the survey area
  regionTotal <- sum(sp_dens_init$cell_N)
  
  
  if(regionTotal < target){
    #warning("There are less than the target number of animals in the whole surveyed region. Returning the whole region")
    subRegions <- st_sf(n = sum(regionTotal), geom = st_union(sp_dens_init)) %>%
      mutate(subRegion_Id = "SR1")
    
  }else{
    
    # set vector with range of lengths from maximum cell as the basis to construct buffers
    bound <- st_bbox(sp_dens_region)
    res <- as.numeric(mean(sqrt(st_area(sp_dens_region[sample.int(nrow(sp_dens_region), size = 50), ]))))  # rough calculation of grid resolution
    maxBuffDistStep <- max(bound$xmax - bound$xmin, bound$ymax - bound$ymin)/res + 5
    buffersDistances <- res * seq(0.5, maxBuffDistStep, by = 1)     # res * seq(0.5, 15, by = 0.5) * res
    
    # initialise objects and counters
    sp_dens_current <- sp_dens_init
    cont = TRUE
    i = 1
    subRegions <- list()
    
    while(cont == TRUE){

      # get the cell with maximum abundance
      maxCell <- sp_dens_current %>% arrange(desc(cell_N)) %>% slice(1)
      
      if(maxCell$cell_N > target){
        cellsBuffInt <- maxCell
        bufferTotal <- maxCell$cell_N
        
      }else{
        
        # print(i)
        # if(i == 74) browser()
        
        # cycle over buffer distances from current maximum cell until target is reached
        prev_numCells = 1
        for(buffDist in buffersDistances){
          
          maxCell_buff <- st_buffer(maxCell, dist = buffDist, endCapStyle = "SQUARE")
          cellsBuffInt <- sp_dens_current[lengths(st_intersects(sp_dens_current, maxCell_buff)) > 0, ]
          
          # drop isolated cell(s) caught in the buffer
          if(nrow(cellsBuffInt) > 1){
            cellsBuffInt <- cellsBuffInt[lengths(st_intersects(cellsBuffInt)) > 1, ]  
            
            #' if above intersection returns an empty geometry (meaning none of the cells in the buffer 
            #' are "connected), make current max cell as the only valid cell in the buffer
            if(nrow(cellsBuffInt) == 0){
              cellsBuffInt <- maxCell
            }
          }
          
          # get total numbers in cells intersecting buffer
          bufferTotal <- sum(cellsBuffInt$cell_N)
          
          #' interrupting loop if number of cells remained constant with increase in buffer size, which 
          #' indicates isolated cell(s) with N < target
          curr_nCells <- nrow(cellsBuffInt)
          if(curr_nCells == prev_numCells){
            break
          }else{
            prev_numCells <- curr_nCells
          }
          
          if(bufferTotal > target) break
        }  
      }
      
      # save current subregion features
      subRegions[[paste0("SR", i)]] <- st_sf(n = bufferTotal, geom = st_union(cellsBuffInt))
      
      # drop current subregion from remaining main survey area
      sp_dens_current <- sp_dens_current %>% filter(!cell_Id %in% cellsBuffInt$cell_Id)
      
      # Check if there are enough animals remaining
      if(sum(sp_dens_current$cell_N) < target){
        cont <- FALSE
      }else{
        i = i + 1
      }
    }
    
    #browser()
     
    # bind subregions and drop isolated cells with less than target number of animals
    subRegions <- do.call(what = sf:::rbind.sf, args = subRegions) %>% 
      rownames_to_column(var = "subRegion_Id") %>%
      filter(n > target)
  }
  
  return(subRegions)
}



LIDAR_surveyRegions <-function(swath, track_gap, ...){

  #' -----------------------------------------------------------------------------------------
  #' Purpose :
  #'    Breaks down the survey region into sub-regions comprising a chosen target number of animals
  #'
  #' Args:
  #'   - swath      : object of class units, giving the width of the swath of the LIDAR equipment
  #'   - track_gap  : object of class units, the gap (distance) between transect tracks
  #'   - ...        : Arguments to be passed to break_region_target()
  #'
  #' Value:
  #'   sf object with polygons of sub-regions inside the survey region, each containing (at minimum) the
  #'   number of target animals
  #' ------------------------------------------------------------------------------------------

  subRegions <- break_region_target(...)
  
  subRegions_effort <- subRegions %>%
    split(.$subRegion_Id) %>%
    map(function(x, track_width = swath, spacing = track_gap){

      # get survey transects for the current sub-region, for a given strip (i.e swath)
      transects <- draw_transects(region_poly = x, spacing = as.numeric(spacing), units = deparse_unit(spacing))

      # compute survey metrics
      x %>%
        mutate(surv_cov_dist = sum(transects$transect_lt),
               surv_cov_area = surv_cov_dist * track_width,
               subRegion_area = st_area(.),
               subRegion_covg_prop = surv_cov_area/subRegion_area
        )
    })
  
  # bind sub-regions and associated survey metrics into a single sf object
  subRegions_effort <- do.call(what = sf:::rbind.sf, args = subRegions_effort)
  rownames(subRegions_effort) <- NULL
  
  return(subRegions_effort)
}


 
# swath <- 200
# units(swath) <- "m"
# result <- LIDAR_surveyRegions(swath = swath, region_polygon = Forth_Tay_utm, sp_dens_surface = Forth_Tay_autumn_sp,
#                               target = 250, clipToRegion = FALSE)
# 
# result %<>% mutate(n = round(n, digits = 0))
# Forth_Tay_autumn_sp %<>% mutate(cell_N = round(cell_N, digits = 2))
# 
# tm_shape(Forth_Tay_autumn_sp) +
#   tm_polygons("cell_N") +
#   tm_text("cell_N", size = 0.5, just = "top") +
#   tm_shape(result) +
#   tm_polygons(alpha = 0.7) +
#   tm_text("n", ymod = 0.2)
# 
# tm_shape(Forth_Tay_utm) +
#   tm_polygons(lwd = 2, col = "orange") +
#   tm_shape(result) +
#   tm_polygons(alpha = 0.7) +
#   tm_text("n", size = 0.7) +
#   tm_shape(result) +
#   tm_text("subRegion_Id", size = 0.7, ymod = 0.3)
#
# 
# st_intersection(st_centroid(Forth_Tay_autumn_sp), result[40, ]) %>%
#   summarise(sum(cell_N), first(n))






# pop_dens <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Winter/Winter.shp") %>%
#   #   st_transform(crs = st_crs(Forth_Tay)) %>%
#   select(ID:Ocean, LACAN) %>%
#   rename(cell_N = LACAN)
# 
# 
# pop_dens <- st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Winter/Winter.shp") %>%
#   #   st_transform(crs = st_crs(Forth_Tay)) %>%
#   select(ID:Ocean, LAARG) %>%
#   rename(cell_N = LAARG)
# 
# 
# # Forth_Tay <- st_read(dsn = "data/East coast surveys Marine Scotland/Wider_Forth_Tay_survey_area/Wider_Forth_Tay_survey_area.shp")
# # Forth_Tay_utm <- st_transform(Forth_Tay, st_crs(pop_dens))
# eastScot <- st_read(dsn = "data/East coast surveys Marine Scotland/East_coast_Scotland/East_coast_Scotland.shp")
# eastScot_utm <- st_transform(eastScot, st_crs(pop_dens))
# 
# 
# set.seed(3)
# sim_results <- survey_simulator(reps = 100, region_polygon = eastScot_utm, sp_dens_surface = pop_dens, 
#                                 trans_spacing = 3000, trans_spacing_units = "m", 
#                                 platform_swath = 1000, platform_swath_units = "m", 
#                                 check_it = FALSE, do_power = TRUE, impact = -0.2)
# 
# 
# 
# 
# #options(error = recover)
# options(error = NULL)
# 
# sim_results$surveys_results
# sim_results$sim_summaries
# 
# 
# 
# #  
# # st_read(dsn = "data/SeaPop/SEAPOP_OpenSea_Distr_Winter/Winter.shp") %>%
# #    st_intersection(., eastScot_utm) %>%
# #    summarise_at(vars(ALALL:MOBAS), sum)
# 
# 
# 
# sim_results$surveys_results

