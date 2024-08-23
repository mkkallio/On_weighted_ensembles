# ------------------------------------------------------------------------------
# RUN HBV FOR CARAVAN 

library(sf)
library(dplyr)
library(readr)
library(terra)
library(foreach)
library(doParallel)

source("R/HBV.R")
source("R/kge_non_parametric.R")




# ------------------------------------------------------------------------------
# RUN HBV AT 0.05min PARAMETERS from Beck et al. 2020


library(sf)
library(dplyr)
library(readr)
library(terra)
library(foreach)
library(doParallel)

source("R/HBV.R")
source("R/kge_non_parametric.R")

datasets <- c("camels") #, 
              # "camelsaus",
              # "camelsbr",
              # "camelscl",
              # "camelsgb",
              # "camelsdk",
              # "lamah",
              # "hysets")

cl <- parallel::makeCluster(10, outfile = "")
doParallel::registerDoParallel(cl)

logfile <- "log.txt"
write(paste0(Sys.time(), " -- ", "Starting.."),
      logfile,
      append = TRUE)

for(dataset in datasets) {
    
      message <- paste0(Sys.time(), " -- ",
                        "starting ", dataset)
      write(message, logfile, append=TRUE)
      
      
      # ------------------
      # read caravan data
      attr_file <- paste0("../caravan/attributes/",
                          dataset,
                          "/attributes_other_",
                          dataset,
                          ".csv")
      
      if(dataset == "camelsdk") {
          caravan <- readr::read_delim(attr_file, delim = ";")
      } else {
          caravan <- readr::read_csv(attr_file)
      }
      
      
      gauges <- sf::st_as_sf(caravan, coords = c("gauge_lon", "gauge_lat"))
      
      
      catchment_file <- paste0("../caravan/shapefiles/",
                          dataset,
                          "/",
                          dataset,
                          "_basin_shapes.shp")
      catchments <- sf::read_sf(catchment_file)
      
      message <- paste0("reading caravan for ", dataset, " finished")
      write(message, logfile, append=TRUE)
      
      
      # -------------------
      # HBV
      n <- nrow(caravan)
      
      # do each fold
      for(i in 1:10) {
          
          message <- paste0(Sys.time(), " -- ",
                            "simulations for ", dataset, 
                            " fold ", i, " starting")
          write(message, logfile, append=TRUE)
          
          
          #----------------
          # process parameters
          tif <- paste0("../HBV SIMREG/beck_et_al_2020/parameter_maps_fold_",i,".tif")
          parameter_maps <- terra::rast(tif)
          
          # extract parameters at gauge locations
          params <- terra::extract(parameter_maps, terra::vect(catchments))
          params <- split(params, params$ID)
          
          max_sims <- max(sapply(params, nrow))
          
          #--------------------------
          # run HBV for all locations
          # parallelise loop through each gauge
          result <- foreach(ii = 1:nrow(caravan),
                            .errorhandling = "stop") %dopar% {
              
              gauge <- caravan$gauge_id[ii]
              
              filename <- paste0("Data/HBV05_simulations/",
                                 dataset, 
                                 "/fold_",
                                 i,
                                 "/",
                                 gauge,
                                 ".fst")
              test <- file.exists(filename)
              if(test) {
                  message <- paste(gauge, "already exists")
                  write(message, logfile, append=TRUE)
                  return(NA)
              }
              
              
              file <- paste0("../caravan/timeseries/csv/",
                             dataset,"/", gauge, ".csv")
              test <- file.exists(file)
              if(!test) {
                  message <- paste("missing timeseries file for",
                                    dataset, gauge)
                  write(message, logfile, append=TRUE)
                  return(NA)
              }
              
              # --------------------
              # read and process input timeseries
              timeseries <- data.table::fread(file,
                                              select = c("date", 
                                                         "temperature_2m_mean",
                                                         "total_precipitation_sum",
                                                         "potential_evaporation_sum"))
              timeseries <- tibble::as_tibble(timeseries)
              names(timeseries)[1] <- "Date"
              timeseries$Date <- lubridate::as_date(timeseries$Date)
              rowind <- complete.cases(timeseries)
              ts <- timeseries[rowind,]
              
              spinup_length <- 5*365 # 5 years
              
              test <- nrow(ts) >= spinup_length # min 5 years data for spinup needed
              if(test) {
                  ts_with_spinup <- rbind(ts[1:spinup_length,], ts)
              } else {
                  message <- paste("observation timeseries is too short;",
                                   dataset, gauge)
                  write(message, logfile, append=TRUE)
                  return(NA)
              } 
              
              catchment_params <- params[[ii]]
              parameter_sets <- nrow(catchment_params)
              
              sims <- matrix(NA,
                             nrow = nrow(ts), 
                             ncol = nrow(catchment_params))
              
              # -------------------------------------------
              # run HBV for all parameters within catchment
              # and then average the output
              
              for(iii in 1:parameter_sets) {
                  parameters <- catchment_params[iii,-1]
                  
                  # for beck et al implementation
                  # parameters <- unlist(c(parameters, 1, 1, 1, 1))
                  # 
                  # sim <- HBV(ts_with_spinup$total_precipitation_sum,
                  #            ts_with_spinup$potential_evaporation_sum,
                  #            ts_with_spinup$temperature_2m_mean,
                  #            parameters)$Q
                  # 
                  # sims[,iii] <- sim[-c(1:spinup_length)]

                  
                  # for gannon et al hbv implementation
                  parameters <- unlist(parameters)
                  parameters <- c(parameters[2],
                                  parameters[1],
                                  parameters[6],
                                  SFCF = 1,
                                  parameters[9],
                                  parameters[10],
                                  parameters[11],
                                  parameters[12],
                                  parameters[3],
                                  parameters[4],
                                  parameters[5],
                                  parameters[8],
                                  parameters[7],
                                  MAXBAS = 1)
                  
                  sim <- HBV2(parameters,
                              ts_with_spinup$total_precipitation_sum,
                              ts_with_spinup$temperature_2m_mean,
                              ts_with_spinup$potential_evaporation_sum,
                              0)$q

                  sims[,iii] <- sim[-c(1:spinup_length)]
              }
              
              # --------------
              # prepare output
              simulations <- tibble::tibble(Date = ts$Date,
                                            tibble::as_tibble(sims))
              
              names(simulations)[-1] <- paste0("set_", 
                                               1:parameter_sets)
              
              
              simulations <- dplyr::left_join(timeseries[,1],
                                              simulations,
                                              by = "Date")
              
              
              # ------------------------
              # save simulations to file
              fst::write_fst(simulations, filename)
              
              message <- paste("simulations done for ", gauge)
              write(message, logfile, append=TRUE)
              
              return(NA)
          }
          
          # simulations <- tibble::tibble(Date = timeseries$Date,
          #                               tibble::as_tibble(simulations))
          
          filename <- paste0("Data/HBV05_simulations/",
                             dataset,
                             "_fold_",
                             i,
                             ".fst")
          fst::write_fst(simulations, filename)
          
          message <- paste0(Sys.time(), " -- ",
                            "simulations for ", dataset, 
                            " fold ", i, " finished")
          write(message, logfile, append=TRUE)
      }
  }
parallel::stopCluster(cl)


