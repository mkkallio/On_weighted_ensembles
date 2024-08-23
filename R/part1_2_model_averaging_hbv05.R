# ------------------------------------------------------------------------------
# DO MODEL AVERAGING AND ASSESSMENT FOR HBV 0.05 SIMULATIONS

library(sf)
library(readr)
library(stringr)
library(lubridate)
library(dplyr)
library(foreach)
library(doParallel)
sf::sf_use_s2(FALSE)
source("R/kge_non_parametric.R")
source("R/model_averaging_functions.R")
source("R/HBV.R")



# ------------------------------------------------------------------------------
# Load data


# target <- readRDS("data/isimip/camels_target.RDS")
methods <- c("CLS", "NNLS", "NNLS2", "GRA", "GRB", "GRC", 
             "best", 
             "ensemble_mean")

attr_file <- paste0("../caravan/attributes/camels/attributes_other_camels.csv")
caravan <- readr::read_csv(attr_file)

gauges <- sf::st_as_sf(caravan, coords = c("gauge_lon", "gauge_lat"))




# ------------------------------------------------------------------------------
# compute ensembles and record weights, gof, and hydrological signatures
cl <- parallel::makeCluster(10, outfile = "")
doParallel::registerDoParallel(cl)

logfile <- "log_assessment.txt"
write(paste0(Sys.time(), " -- ", "Starting.."),
      logfile,
      append = TRUE)


ind <- 1:nrow(gauges)
result <- foreach(i = ind,
                  .errorhandling = "stop",
                  .packages = c("dplyr", "lubridate", "pracma",
                                "data.table")) %dopar% {
                      

    # is this station already processed?
    station <- gauges$gauge_id[i]
                      
    message <- paste0(Sys.time()," -- i ", i, " starting ", station)
    write(message, logfile, append=TRUE)
                      
    filename <- paste0("Data/HBV05_simulations/",
                       # dataset, 
                       "camels",
                       "/assessment_final/",
                       station,
                       "_gof.fst")
    test <- file.exists(filename)
    if(test) {
        message(Sys.time(), " -- ", filename, " already exist. skipping.")
        return(NA)
    }

    
    # read_all_hbv simulations
    simulations <- lapply(1:10, \(ii) {
        
        file <- paste0("data/hbv05_simulations/camels/fold_", ii,
                       "/", station, ".fst")
        test <- file.exists(file)
        if(!test) next
        ts <- fst::read_fst(file)
        
        names(ts)[-1] <- paste0("fold", ii, "_", names(ts)[-1])
        
        if(ii > 1) ts <- ts[,-1]
        return(ts)
    }) %>% do.call(cbind, .) %>% as_tibble()
  
    
    # read observations
    file <- paste0("../CARAVAN/timeseries/csv/camels/", station, ".csv")
    test <- file.exists(file)
    if(!test) next
    
    observations <- readr::read_csv(file,
                                    col_types = c("D", rep("d", 39)),
                                    col_select = c(Date = "date", 
                                                   "streamflow"))
    
    # prepare matrix for model averaging and evaluation
    predmat <- simulations
    predmat$observations <- observations$streamflow
    predmat <- predmat[complete.cases(predmat),]
    
    if(nrow(predmat) < 365) {
        message(Sys.time(), " -- ", station, " predmat contains less than 365 ",
                "timesteps. Skipping.")
        next
    } 
 
    

    # Do model averaging
    ind <- bind_rows(evaluate_individual(predmat, FALSE),
                     evaluate_individual(predmat, TRUE))
    
    egof <- list(run_comb(methods, predmat, "none", FALSE),
                 run_comb(methods, predmat, "log", FALSE),
                 run_comb(methods, predmat, "sqrt", FALSE),
                 run_comb(methods, predmat, "none", TRUE),
                 run_comb(methods, predmat, "log", TRUE),
                 run_comb(methods, predmat, "sqrt", TRUE))
    
    
    # extract gof, weights and hydrological signatures from the egof object
    gofs <- bind_rows(ind,
                      lapply(egof, \(x) return(x$gof)) %>%
                          do.call(bind_rows, .))

    weights <- lapply(egof, \(x) return(x$weights)) %>%
        do.call(bind_rows, .)
    
    hydsig <- lapply(egof, \(x) return(x$hydsig)) %>%
      do.call(bind_rows, .)
    
    
    # Save the outputs in files
    filename <- paste0("Data/HBV05_simulations/",
                       # dataset, 
                       "camels",
                       "/assessment_final/",
                       station,
                       "_gof.fst")
    fst::write_fst(gofs, filename)
    
    filename <- paste0("Data/HBV05_simulations/",
                       # dataset,
                       "camels",
                       "/assessment_final/",
                       station,
                       "_weights.fst")
    fst::write_fst(weights, filename)
    
    filename <- paste0("Data/HBV05_simulations/",
                       # dataset,
                       "camels",
                       "/hydsig_final/",
                       station,
                       "_full.fst")
    fst::write_fst(hydsig, filename)
    
    
    # some memory management
    rm(ind, egof, gofs, weights, predmat, simulations, observations)
    gc()
    
    message <- paste0(Sys.time()," -- i ", i, " finished ", station)
    write(message, logfile, append=TRUE)
}
parallel::stopCluster(cl)



