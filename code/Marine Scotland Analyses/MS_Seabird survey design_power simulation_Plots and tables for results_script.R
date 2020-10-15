### -------------------------------------------------------------------------- ###
### ---                                                                    --- ###
### ---     Script for plotting and tabulate power simulation results      --- ###
### ---                                                                    --- ###
### ---                                                                    --- ###
### -------------------------------------------------------------------------- ###
### -------------------------------------------------------------------------- ###


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
library(ggthemes)
library(ggsci)

outFolder <- "outputs/Plots and Tables/"

speciesKey <- tribble(~species_name, ~species_code,
                      "Atlantic Puffin", "FRARC",
                      "Black Legged Kittiwake", "RITRI",
                      "Common Guillemot", "URAAL",
                      "Herring Gull", "LAARG",
                      "Northern Fulmar", "FUGLA",
                      "Northern Gannet", "MOBAS",
                      "Razorbill", "ALTOR", 
                      "Great black-backed gull", "LAMAR")


regionFullNames <- list(Herman_Berwick = "East of Hermaness to Berwick", #"Herman. to Berwick"
                        east_Scotland = "East Coast Scotland", #"East Scotland",
                        Forth_Tay = "Wider Forth and Tay",  #"Forth and Tay",
                        Moray_Firth = "Wider Moray Firth") #"Moray Firth")



options(dplyr.width = Inf)


## -------------------------------------- ###
## ----   Import simulation results      ----
## -------------------------------------- ###

sim_results_all <- list.files("outputs/simulations/", full.names = TRUE) %>%
  map(function(x){
    read_rds(x)
  })


names(sim_results_all) <- list.files("outputs/simulations/", full.names = FALSE)




## --------------------------------------------- ###
## ----   Re-structure power results data       ----
## --------------------------------------------- ###

### ~~~  SeaPop Data  ~~~ ###

SeaPop_PowerResults_names <- str_subset(names(sim_results_all), "SeaPop") %>%
  str_subset(., "LiDAR", negate = TRUE)

# re-structure to get Seapop power outputs in dataframe format for each survey region and prediction layer
SeaPop_res_power <- sim_results_all[SeaPop_PowerResults_names] %>% 
  map(function(x){
    
    x$all_sim_results %>%
      map(function(x){
        bind_cols(x$sim_summaries$sim_power_stats,
                  x$sim_summaries$sim_surveys_stats,
                  x$sim_summaries$sim_specs)
        
      }) %>%
      data.table::rbindlist(idcol = "Scenario") %>%
      mutate(sp_code = str_extract(Scenario, "(?<=Species = )\\w+"), 
             season = str_to_title(str_extract(Scenario, "(?<=Season = )\\w+")))
    
  })



### ~~~  MERP Data  ~~~ ###

MERP_PowerResults_names <- str_subset(names(sim_results_all), "MERP")

# re-structure to get Seapop power outputs in dataframe format for each survey region and prediction layer
MERP_res_power <- sim_results_all[MERP_PowerResults_names] %>% 
  map(function(x){
    
    x$all_sim_results %>%
      map(function(x){
        bind_cols(x$sim_summaries$sim_power_stats,
                  x$sim_summaries$sim_surveys_stats,
                  x$sim_summaries$sim_specs)
        
      }) %>%
      data.table::rbindlist(idcol = "Scenario") %>%
      mutate(sp_code = str_extract(Scenario, "(?<=Species = )\\w+"),
             season = str_to_title(str_extract(Scenario, "(?<=Season = )\\w+")))
    
  })





### ~~~  Harbour Porpoise Data  ~~~ ###

HP_PowerResults_names <- str_subset(names(sim_results_all), "Porpoise")

# re-structure to get Seapop power outputs in dataframe format for each survey region and prediction layer
HP_res_power <- sim_results_all[HP_PowerResults_names] %>% 
  map(function(x){
    
    x$all_sim_results %>%
      map(function(x){
        bind_cols(x$sim_summaries$sim_power_stats,
                  x$sim_summaries$sim_surveys_stats,
                  x$sim_summaries$sim_specs)
        
      }) %>%
      data.table::rbindlist(idcol = "Scenario") %>%
      mutate(sp_code = str_extract(Scenario, "(?<=predLayer = )\\w+"))
    
  })




## ------------------------------------------------- ###
## ----    SeaPop Power plots - fitted Layer       ----
## ------------------------------------------------- ###

# check species overall abundance per season in each region
SeaPop_res_power[str_which(names(SeaPop_res_power), "fitted")] %>%
  map(function(x){
    x %>%
      group_by(sp_code) %>%
      summarise(N_reference = min(N_reference))
  })



# ----- Most abundant (URAAL) and rarer (FRARC) species - all seasons and all areas ------ #

SeaPop_res_power[str_which(names(SeaPop_res_power), "fitted")] %>%
  map2(., names(.), function(x, y){
    
    #browser()
    
    c_region <- str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+")
    region_name <- regionFullNames[[str_subset(names(regionFullNames), str_split(c_region, pattern = " ")[[1]][2])]]

    # c_region <- tolower(str_replace(str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+"), pattern = " ", "_"))
    # region_name <- regionFullNames[[which(tolower(names(regionFullNames)) == c_region)]]
    
    x %>%
      filter(sp_code %in% c("URAAL", "FRARC")) %>%
      mutate(region_name = region_name) %>%
      group_by(region_name, sp_code, season) %>%
      nest() %>%
      mutate(powPlots = pmap(list(data, sp_code, season, region_name), function(dat, sp, seas, regName){
               
        #browser()
        
        spName <- filter(speciesKey, species_code == sp)$species_name
        
        dat %<>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(Coverage = round(mean_cvrg_prop*100, 0),
                 flight_dist = set_units(mean_cvrg_dist, "km")) 
        
        meanFlightDistPerCov <- dat %>%
          group_by(Coverage) %>%
          summarise(flight_dist_mean = round(mean(flight_dist), 0), 
                    flight_dist_mean = floor(flight_dist_mean/10)*10)
        
        
        p <- dat %>%
          left_join(., meanFlightDistPerCov, by = "Coverage") %>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(panel_title = fct_inorder(as.factor(paste0(Coverage, "% (~", flight_dist_mean, " km)")))) %>%
          ggplot() +
          geom_line(aes(impact, power, col = panel_title), size = 1.2, alpha = 0.9, show.legend = FALSE) +
          geom_hline(yintercept = 0.8, linetype = 'dashed') +
          geom_vline(xintercept = c(-0.3, 0.3), col = 'blue', alpha = 0.5) +
          #scale_color_calc() +
          scale_colour_d3("category10") +
          facet_wrap(.~panel_title, nrow = 2) +
          ggtitle(label = spName, paste0(regName, " - ", seas)) +
          xlab("Impact - proportional change") + #ylab("Power to detect change") +
          theme_minimal()
        
        ggsave(paste0("outputs/Report Plots/power analysis/SeaPop_Power V plots_",
                      regName, "_", spName, "_", seas, ".png"),
               plot = p,
               device = 'png', width = 11, height = 5, units = 'in')
      }))
  })




# ----- all species in strongest season (Autumn) ------ #

SeaPop_res_power[str_which(names(SeaPop_res_power), "fitted")] %>%
  map2(., names(.), function(x, y){
    
    #browser()
    
    c_region <- str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+")
    region_name <- regionFullNames[[str_subset(names(regionFullNames), str_split(c_region, pattern = " ")[[1]][2])]]
    
    # c_region <- tolower(str_replace(str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+"), pattern = " ", "_"))
    # region_name <- regionFullNames[[which(tolower(names(regionFullNames)) == c_region)]]
    
    x %>%
      filter(season == "Autumn") %>%
      mutate(region_name = region_name) %>%
      group_by(region_name, sp_code, season) %>%
      nest() %>%
      mutate(powPlots = pmap(list(data, sp_code, season, region_name), function(dat, sp, seas, regName){
        
        browser()
        
        spName <- filter(speciesKey, species_code == sp)$species_name
        
        dat %<>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(Coverage = round(mean_cvrg_prop*100, 0),
                 flight_dist = set_units(mean_cvrg_dist, "km")) 
        
        meanFlightDistPerCov <- dat %>%
          group_by(Coverage) %>%
          summarise(flight_dist_mean = round(mean(flight_dist), 0), 
                    flight_dist_mean = floor(flight_dist_mean/10)*10)
        
        
        p <- dat %>%
          left_join(., meanFlightDistPerCov, by = "Coverage") %>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(panel_title = fct_inorder(as.factor(paste0(Coverage, "% (~", flight_dist_mean, " km)")))) %>%
          ggplot() +
          geom_line(aes(impact, power, col = panel_title), size = 1.2, alpha = 0.9, show.legend = FALSE) +
          geom_hline(yintercept = 0.8, linetype = 'dashed') +
          geom_vline(xintercept = c(-0.3, 0.3), col = 'blue', alpha = 0.5) +
          #scale_color_calc() +
          scale_colour_d3("category10") +
          facet_wrap(.~panel_title, nrow = 2) +
          ggtitle(label = spName, paste0(regName, " - ", seas)) +
          xlab("Impact - proportional change") + #ylab("Power to detect change") +
          theme_minimal()
        
        ggsave(paste0("outputs/Report Plots/power analysis/SeaPop_Power V plots_",
                      regName, "_", spName, "_", seas, ".png"),
               plot = p,
               device = 'png', width = 11, height = 5, units = 'in')
      }))
  })





## -------------------------------------------------------------------------------------------------- ###
## ----    SeaPop Power plots - EVERYTHING: all regions, all species, all seasons, all surfaces    ----
## -------------------------------------------------------------------------------------------------- ###

SeaPop_res_power %>%
  map2(., names(.), function(x, y){
    
    #browser()
    
    c_region <- str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+")
    region_name <- regionFullNames[[str_subset(names(regionFullNames), str_split(c_region, pattern = " ")[[1]][2])]]
    
    x %>%
      #filter(season == "Autumn") %>%
      mutate(region_name = region_name, 
             layer = case_when(
               str_detect(sp_code, "ci025") ~ "2.5%tile" ,
               str_detect(sp_code, "ci975") ~ "97.5%tile",
               NA ~ "Mean"
             ), 
             sp_code = str_replace_all(sp_code, "ci025|ci975", replacement = "")) %>%
      group_by(region_name, sp_code, season, layer) %>%
      nest() %>%
      mutate(powPlots = pmap(list(data, sp_code, season, region_name, layer), function(dat, sp, seas, regName, layerName){
        
        #browser()
        
        spName <- filter(speciesKey, species_code == sp)$species_name
        
        dat %<>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(Coverage = round(mean_cvrg_prop*100, 0),
                 flight_dist = set_units(mean_cvrg_dist, "km")) 
        
        meanFlightDistPerCov <- dat %>%
          group_by(Coverage) %>%
          summarise(flight_dist_mean = round(mean(flight_dist), 0), 
                    flight_dist_mean = floor(flight_dist_mean/10)*10)
        
        
        p <- dat %>%
          left_join(., meanFlightDistPerCov, by = "Coverage") %>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(panel_title = fct_inorder(as.factor(paste0(Coverage, "% (~", flight_dist_mean, " km)")))) %>%
          ggplot() +
          geom_line(aes(impact, power, col = panel_title), size = 1.2, alpha = 0.9, show.legend = FALSE) +
          geom_hline(yintercept = 0.8, linetype = 'dashed') +
          geom_vline(xintercept = c(-0.3, 0.3), col = 'blue', alpha = 0.5) +
          #scale_color_calc() +
          scale_colour_d3("category10") +
          facet_wrap(.~panel_title, nrow = 2) +
          ggtitle(label = paste0(spName, " - ", layerName), paste0(regName, " - ", seas)) +
          xlab("Impact - proportional change") + #ylab("Power to detect change") +
          theme_minimal()
        
        layerName <- str_replace(layerName, "\\%", "pc")
        
        ggsave(paste0("outputs/Report Plots/power analysis/all SeaPop/SeaPop_Power V plots_",
                      regName, "_", spName, "_", seas, "_", layerName, ".png"), 
               plot = p,
               device = 'png', width = 11, height = 5, units = 'in')
      }))
  })








## ------------------------------------------------- ###
## ----    SeaPop Power Table - EVERYTHING       ----
## ------------------------------------------------- ###

powerTable <- rbindlist(SeaPop_res_power, idcol = "Region_Layer_ID") %>%
  mutate(region_short = str_extract(Region_Layer_ID, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+"),
         predLayer = str_extract(Region_Layer_ID, "(?<=_)[:alnum:]+(?=.rds)")) %>%
  filter(power >= 0.8) %>% 
  select(region_short, sp_code, season, predLayer, impact, mean_cvrg_prop, mean_cvrg_dist, power) %>%
  group_by(region_short, sp_code, season, predLayer, impact) %>%
  arrange(mean_cvrg_dist, power, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(impact = round(impact, 1)) %>%
  filter(impact %in% c(-0.3, 0.3)) %>%
  group_by(region_short, sp_code, season, predLayer) %>%
  summarise(impact = abs(first(impact)), 
            mean_cvrg_prop = mean(mean_cvrg_prop),
            mean_cvrg_dist = round(set_units(mean(mean_cvrg_dist), "km"), 0), 
            power = mean(power)) %>%
  ungroup() %>%
  mutate(sp_code = str_replace(sp_code, "ci025", ""), 
         sp_code = str_replace(sp_code, "ci975", ""), 
         impact = paste0("+/-", impact), 
         app_surv_days = round(mean_cvrg_prop/set_units(2000, "km/d"), 2)) %>%
  left_join(., speciesKey, by = c("sp_code" ="species_code")) %>%
  complete(sp_code, region_short, season, predLayer)


powerTable

write_csv(powerTable, "outputs/SeaPop_Power table_30pct impact.csv")





# powerTable %>%
#   select(species_name, sp_code, region_short, season, predLayer, impact, mean_cvrg_prop) %>%
#   spread(predLayer, mean_cvrg_prop)
# 
# 
# powerTable
# 
# 
# table(powerTable$region_short, powerTable$sp_code, powerTable$season,  powerTable$predLayer)











## ------------------------------------------------- ###
## ----    MERP Power plots - fitted Layer       ----
## ------------------------------------------------- ###

# run for URAAL species, for contrast with SeaPop
MERP_res_power %>%
  map2(., names(.), function(x, y){
    
    #browser()
    
    c_region <- str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+")
    region_name <- regionFullNames[[str_subset(names(regionFullNames), str_split(c_region, pattern = " ")[[1]][2])]]
    
    x %>%
      filter(sp_code == "URAAL", 
             season == "Autumn") %>%
      mutate(region_name = region_name) %>%
      group_by(region_name, sp_code, season) %>%
      nest() %>%
      mutate(powPlots = pmap(list(data, sp_code, season, region_name), function(dat, sp, seas, regName){
        
        #browser()
        
        spName <- filter(speciesKey, species_code == sp)$species_name
        
        dat %<>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(Coverage = round(mean_cvrg_prop*100, 0),
                 flight_dist = set_units(mean_cvrg_dist, "km")) 
        
        meanFlightDistPerCov <- dat %>%
          group_by(Coverage) %>%
          summarise(flight_dist_mean = round(mean(flight_dist), 0), 
                    flight_dist_mean = floor(flight_dist_mean/10)*10)
        
        
        p <- dat %>%
          left_join(., meanFlightDistPerCov, by = "Coverage") %>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(panel_title = fct_inorder(as.factor(paste0(Coverage, "% (~", flight_dist_mean, " km)")))) %>%
          ggplot() +
          geom_line(aes(impact, power, col = panel_title), size = 1.2, alpha = 0.9, show.legend = FALSE) +
          geom_hline(yintercept = 0.8, linetype = 'dashed') +
          geom_vline(xintercept = c(-0.3, 0.3), col = 'blue', alpha = 0.5) +
          #scale_color_calc() +
          scale_colour_d3("category10") +
          facet_wrap(.~panel_title, nrow = 2) +
          ggtitle(label = spName, paste0(regName, " - ", seas)) +
          xlab("Impact - proportional change") + #ylab("Power to detect change") +
          theme_minimal()
        
        ggsave(paste0("outputs/Report Plots/power analysis/MERP_Power V plots_",
                      regName, "_", spName, "_", seas, ".png"),
               plot = p,
               device = 'png', width = 11, height = 5, units = 'in')
      }))
  })






## ---------------------------------------------------------- ###
## ----    Harbour Porpoise Power plots - fitted Layer       ----
## ---------------------------------------------------------- ###

HP_res_power %>%
  map2(., names(.), function(x, y){
    
    #browser()
    
    c_region <- str_extract(y, "(?<=outputs_)[:alpha:]+\\s[:alpha:]+")
    region_name <- regionFullNames[[str_subset(names(regionFullNames), str_split(c_region, pattern = " ")[[1]][2])]]
    
    x %>%
      filter(sp_code == "mean") %>%
      mutate(region_name = region_name) %>%
      group_by(region_name) %>%
      nest() %>%
      mutate(powPlots = pmap(list(data, region_name), function(dat, regName){
        
        #browser()
        
        dat %<>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(Coverage = round(mean_cvrg_prop*100, 0),
                 flight_dist = set_units(mean_cvrg_dist, "km")) 
        
        meanFlightDistPerCov <- dat %>%
          group_by(Coverage) %>%
          summarise(flight_dist_mean = round(mean(flight_dist), 0), 
                    flight_dist_mean = floor(flight_dist_mean/10)*10)
        
        
        p <- dat %>%
          left_join(., meanFlightDistPerCov, by = "Coverage") %>%
          arrange(desc(mean_cvrg_prop)) %>%
          mutate(panel_title = fct_inorder(as.factor(paste0(Coverage, "% (~", flight_dist_mean, " km)")))) %>%
          ggplot() +
          geom_line(aes(impact, power, col = panel_title), size = 1.2, alpha = 0.9, show.legend = FALSE) +
          geom_hline(yintercept = 0.8, linetype = 'dashed') +
          geom_vline(xintercept = c(-0.3, 0.3), col = 'blue', alpha = 0.5) +
          #scale_color_calc() +
          scale_colour_d3("category10") +
          facet_wrap(.~panel_title, nrow = 2) +
          ggtitle(label = "Harbour Porpoise", paste0(regName)) +
          xlab("Impact - proportional change") + #ylab("Power to detect change") +
          theme_minimal()
        
        ggsave(paste0("outputs/Report Plots/power analysis/Harbour Porpoise_Power V plots_",
                      regName, ".png"),
               plot = p,
               device = 'png', width = 11, height = 5, units = 'in')
      }))
  })















