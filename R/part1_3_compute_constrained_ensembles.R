# ------------------------------------------------------------------------------
# Compute gof, weights and hydrological signatures of constrained ensembles
# with 2-200 members

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



attr_file <- paste0("../caravan/attributes/camels/attributes_other_camels.csv")
caravan <- readr::read_csv(attr_file)

gauges <- sf::st_as_sf(caravan, coords = c("gauge_lon", "gauge_lat"))




# ------------------------------------------------------------------------------
# prep
cl <- parallel::makeCluster(10, outfile = "")
doParallel::registerDoParallel(cl)

logfile <- "log_weights.txt"
write(paste0(Sys.time(), " -- ", "Starting.."),
      logfile,
      append = TRUE)

# constrained ensembles computed for only these methods and preprocessing
# options
methods <- c("CLS", "NNLS2", "GRC", "ensemble_mean")
biascor <- c(TRUE)
preprocess <- c("none", "log","sqrt")



result <- foreach(i = 1:nrow(gauges),
                  .errorhandling = "stop",
                  .packages = c("dplyr", "lubridate")) %dopar% {
                      
      
      # is this station already processed?
      station <- gauges$gauge_id[i]
      
      message <- paste0(Sys.time()," -- starting ", station)
      write(message, logfile, append=TRUE)
      
      filename <- paste0("Data/HBV05_simulations/",
                         # dataset, 
                         "camels",
                         "/const_final/",
                         station,
                         "_weights.fst")
      test <- file.exists(filename)
      if(test) return(NA)
      
      
      # read all simulations
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
      
      n <- floor(nrow(predmat)/2)
      train <- 1:n
      test <- (n+1):nrow(predmat)
      
      n <- ncol(predmat)
      nmodels <- c(2, 5, 10, 15, 25, 50, 100, 150, 200)
      nmodels <- nmodels[nmodels < (n-2)]
      
      
      # run the loop
      weights_template <- dplyr::slice(predmat[,-c(1, n)], 0)
      weightlist <- list()
      goflist <- list()
      siglist <- list()
      
      for(bc in biascor) {
          for(pp in preprocess) {
              
              combmat <- do_preprocess(predmat, pp, bc)
              
              
               
              for(nm in nmodels) {
                  for(ii in 1:10) {
                      
                      name <- paste("nmod", nm, "set", ii, sep = "_")
                      # mname <- paste(method, bc, pp, name, sep = "_")
                      
                      
                      mods <- sample(2:(n-1), nm)
                      
                      cmat <- combmat[, c(1,mods,n )]
                      # tbl <- run_comb(methods, cmat, pp, bc)
                      
                      for(method in methods) {

                          mname <- paste(method, bc, pp, name, sep = "_")
                          id <- paste(station, "c", mname, ii, sep = "_")

                          # tbl <- run_comb(methods, predmat, pp, bc)

                          if(method == "ensemble_mean") {
                              models <- names(cmat)[-c(1,ncol(cmat))]  
                              
                              combs <- list(Method = "Ensemble mean",
                                            Models = models,
                                            Weights = rep(1, nm) / nm,
                                            Intercept = 0)
                              
                          } else {
                              combs <- try(R.utils::withTimeout(combinations(cmat[train,], 
                                                                  type = method),
                                                                timeout = 10), 
                                           TRUE)
                              
                              testi <- class(combs) == "try-error"
                              if(testi) next
                          }
                          
                              # process and evaluate full timeseries
                              combs <- process_combs(cmat, 
                                                     combs, 
                                                     pp, 
                                                     train, 
                                                     test)
                              
                              rnp <- assess_comb(combs, 
                                                 predmat$observations, 
                                                 train, 
                                                 test)
                              
                              
                              if(method != "ensemble_mean") {
                                  w <- combs$weights
                                  w <- bind_rows(weights_template,
                                                 as_tibble(t(w))) %>%
                                      mutate(gauge_id = station,
                                             id = id,
                                             n_sim = nm,
                                             model = method,
                                             bias_correct = bc,
                                             preprocess = pp,
                                             id = id) %>%
                                      select(gauge_id, id, n_sim, model, 
                                             bias_correct, preprocess,
                                             intercept, everything())
    
                                  weightlist[[mname]] <- w
                              }

                          # create output
                              tbl_train <- bind_cols(tibble(gauge_id = station,
                                                            id = id,
                                                            model = method,
                                                            method = "ensemble",
                                                            preprocess = pp,
                                                            bias_correct = bc,
                                                            period = "train",
                                                            sig = "full",
                                                            n_sim = nm),
                                                     rnp$train)
                          
                              tbl_test <- bind_cols(tibble(gauge_id = station,
                                                           id = id,
                                                           model = method,
                                                           method = "ensemble",
                                                           preprocess = pp,
                                                           bias_correct = bc,
                                                           period = "test",
                                                           sig = "full",
                                                           n_sim = nm),
                                                    rnp$test)

                          hydsig_train <- do_hydsig_dt(cmat$Date[train], 
                                                    predmat$observations[train], 
                                                    combs$opt_train, 
                                                    method, 
                                                    pp, 
                                                    bc,
                                                    "train",
                                                    id,
                                                    nm)
                          hydsig_test <- do_hydsig_dt(cmat$Date[test], 
                                                   predmat$observations[test], 
                                                   combs$opt_test, 
                                                   method, 
                                                   pp, 
                                                   bc,
                                                   "test",
                                                   id,
                                                   nm)
                          
                          tbl <- bind_rows(tbl_train, tbl_test,
                                           hydsig_train, hydsig_test) %>% 
                              mutate(n_sim = nm) %>% 
                              mutate(set = ii, .before = 2)
                          
                          # hydrological signatures
                          train_sig <- hydrological_signatures(combs$opt_train) %>% 
                            tidyr::pivot_wider(names_from = sig, values_from = value)
                          test_sig <- hydrological_signatures(combs$opt_test) %>% 
                            tidyr::pivot_wider(names_from = sig, values_from = value)
                          
                          sig_train_tbl <- tibble(gauge_id = station,
                                                  id = id,
                                                  model = method,
                                                  n_sim = nm,
                                                  bias_correct = bc,
                                                  preprocess = pp,
                                                  period = "train",
                                                  train_sig)
                          
                          sig_test_tbl <- tibble(gauge_id = station,
                                                 id = id,
                                                 model = method,
                                                 n_sim = nm,
                                                 bias_correct = bc,
                                                 preprocess = pp,
                                                 period = "test",
                                                 test_sig)
                          sig_tbl <- rbind(sig_train_tbl, sig_test_tbl) %>% 
                            mutate(n_sim = nm) %>% 
                            mutate(set = ii, .before = 3)
                          
                          goflist[[mname]] <- tbl
                          siglist[[mname]] <- sig_tbl
                          
                      }
                      message(ii)
                  }
                  message(method, " = ", ii)
              }
              message(bc, " ", pp)
          }
      }

      weightlist <- do.call(rbind, weightlist)
      goflist <- do.call(rbind, goflist)
      siglist <- do.call(rbind, siglist)
      
      # save the output to files
      filename <- paste0("Data/HBV05_simulations/",
                         # dataset, 
                         "camels",
                         "/const_final/",
                         station,
                         "_weights.fst")
      
      fst::write_fst(weightlist, filename)
      
      filename <- paste0("Data/HBV05_simulations/",
                         # dataset, 
                         "camels",
                         "/const_final/",
                         station,
                         "_gof.fst")
      
      fst::write_fst(goflist, filename)
      
      filename <- paste0("Data/HBV05_simulations/",
                         # dataset, 
                         "camels",
                         "/hydsig_final/",
                         station,
                         "_const.fst")
      
      fst::write_fst(siglist, filename)
      
      
      message <- paste0(Sys.time()," -- finished ", station)
      write(message, logfile, append=TRUE)
}
parallel::stopCluster(cl)

