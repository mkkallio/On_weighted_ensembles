# ------------------------------------------------------------------------------
# Compute hydrological signatures if they have not been computed prior

library(sf)
library(readr)
library(stringr)
library(lubridate)
library(dplyr)
# library(units)
library(foreach)
library(doParallel)
library(lfstat)
library(pracma)
sf::sf_use_s2(FALSE)
source("R/kge_non_parametric.R")
source("R/model_averaging_functions.R")
source("R/HBV.R")



# ------------------------------------------------------------------------------
# Load data

attr_file <- paste0("../caravan/attributes/camels/attributes_other_camels.csv")
caravan <- readr::read_csv(attr_file)

gauges <- sf::st_as_sf(caravan, coords = c("gauge_lon", "gauge_lat"))




# ------------------------------------------------------------------------------
# compute ensembles and record weights and gof
cl <- parallel::makeCluster(10, outfile = "")
doParallel::registerDoParallel(cl)

logfile <- "log_hydsig.txt"
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
          
          
          
          filename <- paste0("Data/HBV05_simulations/",
                             # dataset, 
                             "camels",
                             "/hydsig_final/",
                             station,
                             ".fst")
          test <- file.exists(filename)
          if(test) {
            message <- paste0(Sys.time(), " -- ", filename, " already exist. skipping.")
            write(message, logfile, append=TRUE)
            return(NA)
          }
          
          message <- paste0(Sys.time()," -- i ", i, " starting ", station)
          write(message, logfile, append=TRUE)
          
          # read_all_sims at 5arcmin
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
          
          
          # read_obs
          file <- paste0("../CARAVAN/timeseries/csv/camels/", station, ".csv")
          test <- file.exists(file)
          if(!test) next
          
          observations <- readr::read_csv(file,
                                          col_types = c("D", rep("d", 39)),
                                          col_select = c(Date = "date", 
                                                         "streamflow"))#, 
         
          
          
          
          
          # prepare matrix for model averaging and evaluation
          predmat <- simulations
          predmat$observations <- observations$streamflow
          predmat <- predmat[complete.cases(predmat),]
          
          # training and testing periods
          n <- floor(nrow(predmat)/2)
          train <- 1:n
          test <- (n+1):nrow(predmat)
          
          
          # compute hydsig for HBV and observations WITH BIAS CORRECTION
          hydsig_nobc <- lapply(2:ncol(predmat), \(ii) {
            
            x <- pull(predmat,ii)
            name <- names(predmat)[ii]
            
            train_sig <- hydrological_signatures(x[train]) %>% 
              tidyr::pivot_wider(names_from = sig, values_from = value)
            
            train_tbl <- tibble(gauge_id = station,
                                model = name,
                                n_sim = 1,
                                bias_correct = FALSE,
                                preprocess = "none",
                                period = "train",
                                train_sig)
            
            test_sig <- hydrological_signatures(x[test]) %>% 
              tidyr::pivot_wider(names_from = sig, values_from = value)
            
            test_tbl <- tibble(gauge_id = station,
                               model = name,
                               n_sim = 1,
                               bias_correct = FALSE,
                               preprocess = "none",
                               period = "test",
                               test_sig)
            tbl <- rbind(train_tbl, test_tbl)
            
            return(tbl)
            
          }) %>% do.call(rbind, .)
          
          
          # compute hydsig for HBV and observations WITH BIAS CORRECTION
          hydsig_bc <- lapply(2:ncol(predmat), \(ii) {
            
            x <- predmat[,c(1,ii, ncol(predmat))]
            x <- do_preprocess(x, "none", TRUE)
            
            x <- pull(predmat,2)
            name <- names(predmat)[2]
            
            train_sig <- hydrological_signatures(x[train]) %>% 
              tidyr::pivot_wider(names_from = sig, values_from = value)
            
            train_tbl <- tibble(gauge_id = station,
                                model = name,
                                n_sim = 1,
                                bias_correct = FALSE,
                                preprocess = "none",
                                period = "train",
                                train_sig)
            
            test_sig <- hydrological_signatures(x[test]) %>% 
              tidyr::pivot_wider(names_from = sig, values_from = value)
            
            test_tbl <- tibble(gauge_id = station,
                               model = name,
                               n_sim = 1,
                               bias_correct = FALSE,
                               preprocess = "none",
                               period = "test",
                               test_sig)
            tbl <- rbind(train_tbl, test_tbl)
            
            return(tbl)
            
          }) %>% do.call(rbind, .)
          
          
          
          tbl_simobs <- rbind(hydsig_nobc, hydsig_bc)
          rm(hydsig_nobc, hydsig_bc)
          
          
          # --------------------------------------------------------------------
          # hyd sig of ensemble mean
          
          n_sim <- ncol(simulations)-1
            
          em <- rowMeans(simulations[,-1], na.rm=TRUE)
          
          train_sig <- hydrological_signatures(em[train]) %>% 
            tidyr::pivot_wider(names_from = sig, values_from = value)
          test_sig <- hydrological_signatures(em[test]) %>% 
            tidyr::pivot_wider(names_from = sig, values_from = value)
          
          train_tbl <- tibble(gauge_id = station,
                              model = "ensemble_mean",
                              n_sim = n_sim,
                              bias_correct = FALSE,
                              preprocess = "none",
                              period = "train",
                              train_sig)
          
          test_tbl <- tibble(gauge_id = station,
                             model = "ensemble_mean",
                             n_sim = n_sim,
                             bias_correct = FALSE,
                             preprocess = "none",
                             period = "test",
                             test_sig)
          tbl <- rbind(train_tbl, test_tbl)
          
          # with bias correction
          em_bc <- do_preprocess(predmat, "none", TRUE)
          em_bc <- rowMeans(em_bc[,-c(1, ncol(em_bc))])
          
          train_sig <- hydrological_signatures(em_bc[train]) %>% 
            tidyr::pivot_wider(names_from = sig, values_from = value)
          test_sig <- hydrological_signatures(em_bc[test]) %>% 
            tidyr::pivot_wider(names_from = sig, values_from = value)
          
          train_tbl <- tibble(gauge_id = station,
                              model = "ensemble_mean",
                              n_sim = n_sim,
                              bias_correct = TRUE,
                              preprocess = "none",
                              period = "train",
                              train_sig)
          
          test_tbl <- tibble(gauge_id = station,
                             model = "ensemble_mean",
                             n_sim = n_sim,
                             bias_correct = TRUE,
                             preprocess = "none",
                             period = "test",
                             test_sig)
          em_tbl <- rbind(tbl, train_tbl, test_tbl)
          
          
          
          # --------------------------------------------------------------------
          # hyd sig of weighted ensembles
          
          w_file <- paste0("Data/HBV05_simulations/",
                             # dataset, 
                             "camels",
                             "/assessment_final/",
                             station,
                             "_weights.fst")
          test <- file.exists(w_file)
          if(!test) {
            message(Sys.time(), " -- ", filename, " weight file doesnt exist. skipping.")
            return(NA)
          }
          
          weights_all <- fst::read_fst(w_file) %>% 
            as_tibble()
          
          tst <- nrow(weights_all) > 0
          if(tst) {
            # process full ensembles
            
            biascor <- c(FALSE, TRUE)
            prep <- c("none", "log", "sqrt")
            hydsig_full <- list()
            system.time({
              for(bc in biascor) {
                for(pp in prep) {
                  
                  weights <- weights_all %>% 
                    filter(bias_correct == bc,
                           preprocess == pp)
                  
                  combmat <- do_preprocess(predmat, pp, bc)
                  combmat <- t(cbind(1, as.matrix(combmat[,-c(1, ncol(combmat))])))
                  modelnames <- c("intercept", names(predmat)[-c(1, ncol(predmat))])
                  
                  for(iii in 1:nrow(weights)) {
                    
                    # Create the optimized timeseries
                    w <- unlist(weights[iii, -c(1:6)])
                    w[is.na(w)] <- 0
                    models <- names(w)
                    modelind <- which(modelnames %in% models, arr.ind=TRUE)
                    Opt <- as.vector(w %*% combmat[modelind,])
                    
                    
                    # back transform
                    if(pp == "log") {
                      opt_train <- exp(Opt[train])
                      opt_test <- exp(Opt[test])
                    } else if(pp == "sqrt") {
                      opt_train <- Opt[train]^2
                      opt_test <- Opt[test]^2
                    } else if(pp == "none") {
                      opt_train <- Opt[train]
                      opt_test <- Opt[test]
                    }
                    
                    train_sig <- hydrological_signatures(opt_train) %>% 
                      tidyr::pivot_wider(names_from = sig, values_from = value)
                    
                    test_sig <- hydrological_signatures(opt_test) %>% 
                      tidyr::pivot_wider(names_from = sig, values_from = value)
                    
                    
                    train_tbl <- tibble(gauge_id = station,
                                        model = weights$model[iii],
                                        n_sim = weights$n_sim[iii],
                                        bias_correct = bc,
                                        preprocess = pp,
                                        period = "train",
                                        train_sig)
                    
                    test_tbl <- tibble(gauge_id = station,
                                       model = weights$model[iii],
                                       n_sim = weights$n_sim[iii],
                                       bias_correct = bc,
                                       preprocess = pp,
                                       period = "test",
                                       test_sig)
                    tbl <- rbind(train_tbl, test_tbl)
                    
                    hydsig_full <- append(hydsig_full, list(tbl))
                  }
                }
              }
            })
            
          } else {
            hydsig_full <- list(NULL)
          }
          
          
          
          
          
          # --------------------------------------------------------------------
          # process constrained ensembles
          
          
          w_file <- paste0("Data/HBV05_simulations/",
                           # dataset, 
                           "camels",
                           "/const_final/",
                           station,
                           "_weights.fst")
          
         
          weights_all <- fst::read_fst(w_file) %>% 
            as_tibble()
          
          
          biascor <- c(TRUE)
          prep <- c("none", "log", "sqrt")
          hydsig_const <- list()
          system.time({
            for(bc in biascor) {
              for(pp in prep) {
                
                weights <- weights_all %>% 
                  filter(bias_correct == bc,
                         preprocess == pp)
                
                combmat <- do_preprocess(predmat, pp, bc)
                combmat <- t(cbind(1, as.matrix(combmat[,-c(1, ncol(combmat))])))
                
                for(iii in 1:nrow(weights)) {
                  
                  # Create the optimized timeseries
                  w <- unlist(weights[iii, -c(1:6)])
                  w[is.na(w)] <- 0
                  Opt <- as.vector(w %*% combmat)
                  
                  # back transform
                  if(pp == "log") {
                    opt_train <- exp(Opt[train])
                    opt_test <- exp(Opt[test])
                  } else if(pp == "sqrt") {
                    opt_train <- Opt[train]^2
                    opt_test <- Opt[test]^2
                  } else if(pp == "none") {
                    opt_train <- Opt[train]
                    opt_test <- Opt[test]
                  }
                  
                  train_sig <- hydrological_signatures(opt_train) %>% 
                    tidyr::pivot_wider(names_from = sig, values_from = value)
                  
                  test_sig <- hydrological_signatures(opt_test) %>% 
                    tidyr::pivot_wider(names_from = sig, values_from = value)
                  
                  
                  train_tbl <- tibble(gauge_id = station,
                                      model = weights$model[iii],
                                      n_sim = weights$n_sim[iii],
                                      bias_correct = bc,
                                      preprocess = pp,
                                      period = "train",
                                      train_sig)
                  
                  test_tbl <- tibble(gauge_id = station,
                                     model = weights$model[iii],
                                     n_sim = weights$n_sim[iii],
                                     bias_correct = bc,
                                     preprocess = pp,
                                     period = "test",
                                     test_sig)
                  tbl <- rbind(train_tbl, test_tbl)
                  
                  hydsig_const <- append(hydsig_const, list(tbl))
                  
                  # ensemble mean of constrained ensembles
                  
                  # Create the optimized timeseries
                  w <- unlist(weights[iii, -c(1:6)])
                  w[1] <- NA
                  na_ind <- is.na(w)
                  w[na_ind] <- 0
                  w[!na_ind] <- 1 / sum(!na_ind)
                  Opt <- as.vector(w %*% combmat)
                  
                  # back transform
                  if(pp == "log") {
                    opt_train <- exp(Opt[train])
                    opt_test <- exp(Opt[test])
                  } else if(pp == "sqrt") {
                    opt_train <- Opt[train]^2
                    opt_test <- Opt[test]^2
                  } else if(pp == "none") {
                    opt_train <- Opt[train]
                    opt_test <- Opt[test]
                  }
                  
                  train_sig <- hydrological_signatures(opt_train) %>% 
                    tidyr::pivot_wider(names_from = sig, values_from = value)
                  
                  test_sig <- hydrological_signatures(opt_test) %>% 
                    tidyr::pivot_wider(names_from = sig, values_from = value)
                  
                  
                  train_tbl <- tibble(gauge_id = station,
                                      model = "ensemble_mean",
                                      n_sim = weights$n_sim[iii],
                                      bias_correct = bc,
                                      preprocess = pp,
                                      period = "train",
                                      train_sig)
                  
                  test_tbl <- tibble(gauge_id = station,
                                     model = "ensemble_mean",
                                     n_sim = weights$n_sim[iii],
                                     bias_correct = bc,
                                     preprocess = pp,
                                     period = "test",
                                     test_sig)
                  tbl <- rbind(train_tbl, test_tbl)
                  
                  hydsig_const <- append(hydsig_const, list(tbl))
                  
                  
                }
              }
            }
          })
          
          
          # --------------------------------------------------------------------
          # combine and write hydrological sigs
          
          hydsig_full <- do.call(rbind, hydsig_full)
          hydsig_const <- do.call(rbind, hydsig_const)
          
          hydsig <- rbind(hydsig_full, em_tbl, hydsig_const,
                          tbl_simobs)
          
          
          # save output
          filename <- paste0("Data/HBV05_simulations/",
                             # dataset, 
                             "camels",
                             "/hydsig_final/",
                             station,
                             ".fst")
          fst::write_fst(hydsig, filename)
          

          
          message <- paste0(Sys.time()," -- i ", i, " finished ", station)
          write(message, logfile, append=TRUE)
        }
parallel::stopCluster(cl)
