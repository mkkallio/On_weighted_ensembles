# ------------------------------------------------------------------------------
# ASSESS THE VARIABILITY OF GOODNESS OF FIT

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

# Function to compute Gini Mean Difference
gmd <- function (x, alpha, na.rm=TRUE) {
  
  if(na.rm) x <- x[!is.na(x)]
  
  if (missing(alpha)) 
    alpha <- 1
  if (alpha == 1) {
    if (is.null(dim(x))) {
      n = length(x)
      r = 2 * seq(1, n) - n - 1
      covgxy = (2 * sum(r * x[order(x)])/(n * (n - 1)))
      return(abs(covgxy))
    }
    else {
      x <- as.matrix(x)
      n <- nrow(x)
      x <- as.matrix(dist(x))
      return(ifelse(n == 1, 0, abs(sum(x)/(n * (n - 1)))))
    }
  }
  else if (0 < alpha & 2 >= alpha) {
    x <- as.matrix(x)
    n <- nrow(x)
    x <- as.matrix(dist(x))
    return(ifelse(n == 1, 0, abs(sum(x^alpha)/(n * (n - 
                                                      1)))))
  }
  else {
    stop("alpha must be the a number between 0 and 2")
  }
}

# ------------------------------------------------------------------------------
# Load data

attr_file <- paste0("../caravan/attributes/camels/attributes_other_camels.csv")
caravan <- readr::read_csv(attr_file)

gauges <- sf::st_as_sf(caravan, coords = c("gauge_lon", "gauge_lat"))




# ------------------------------------------------------------------------------
# compute ensembles and record weights and gof
cl <- parallel::makeCluster(5, outfile = "")
doParallel::registerDoParallel(cl)


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
      # predmat$observations <- units::set_units(observations$streamflow, "mm/d")
      predmat$observations <- observations$streamflow
      predmat <- predmat[complete.cases(predmat),]
      
      # training and testing
      n <- floor(nrow(predmat)/2)
      train <- 1:n
      test <- (n+1):nrow(predmat)
      
      
      
      
      w_file <- paste0("Data/HBV05_simulations/",
                       # dataset, 
                       "camels",
                       "/assessment_final/",
                       station,
                       "_weights.fst")
      
      w_const <- paste0("Data/HBV05_simulations/",
                       # dataset, 
                       "camels",
                       "/const_final/",
                       station,
                       "_weights.fst")
      test <- file.exists(w_file)
      if(!test) {
        return(NA)
      }
      
      weights_all <- fst::read_fst(w_file) %>% 
        bind_rows(fst::read_fst(w_const)) %>% 
        as_tibble()
      
      tst <- nrow(weights_all) > 0
      if(tst) {
        # process full ensembles
        
        biascor <- c(FALSE,TRUE)
        prep <- c("none", "log", "sqrt")
        metrics <- tibble()
        system.time({
          for(bc in biascor) {
            for(pp in prep) {
              
              weights <- weights_all %>% 
                filter(bias_correct == bc,
                       preprocess == pp)
              
              tst <- nrow(weights) == 0
              if(tst) next
              
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
                
                
                temp <- tibble(Date = predmat$Date[test],
                               obs = predmat$observations[test], 
                               sim = opt_test) %>% 
                  mutate(y = year(Date),
                         m = month(Date)) %>% 
                  group_by(y,m) %>% 
                  summarise(rnp = 1 / (2 - RNP(sim, obs)[4]),
                            cv = sd(sim) / mean(sim)) %>% 
                  ungroup() %>%
                  # group_by(method, train) %>% 
                  summarise(mean = mean(rnp, na.rm=TRUE),
                            # cv = sd(rnp, na.rm=TRUE) / mean(rnp, na.rm=TRUE),
                            gmd = gmd(rnp, na.rm=TRUE),
                            # V = (med + (1-cv)),
                            # V_gmd = (med + (1-gmd)),
                            # V_weighted_cv = weighted.mean(c(med, (1-cv)), w = weights),
                            V = weighted.mean(c(mean, (1-gmd)), w = c(0.5,0.5))) 
                
                out <- temp %>% 
                  mutate(gauge_id = station,
                         model = weights$model[iii],
                         bias_correct = bc,
                         preprocess = pp,
                         n_sim = weights$n_sim[iii],
                         .before = 1)
                metrics <- rbind(metrics, out)
                
                gc()
              }
            }
          }
        })
        
        
        gc()
        # metrics <- do.call(rbind, metrics)
      } 
      
      filename <- paste0("Data/HBV05_simulations/",
                         # dataset, 
                         "camels",
                         "/variability_analysis/",
                         station,
                         ".fst")
      fst::write_fst(metrics, filename)
      
      return(NULL)
}
parallel::stopCluster(cl)


# load the data
file_list <- list.files("data/hbv05_simulations/camels/variability_analysis/",
                        full.names = TRUE)

var <- lapply(file_list, \(f2) {
  x <- fst::read_fst(f2)
  x <- as_tibble(x)
  return(x)
}) %>% do.call(bind_rows, .)


library(ggplot2)
var %>% 
  filter(bias_correct == TRUE,
         preprocess == "none") %>% 
  ggplot() +
  geom_point(aes(mean, gmd, color = V)) +
  facet_wrap(~model)

var %>% 
  filter(bias_correct == TRUE,
         preprocess == "none") %>% 
  ggplot() +
  geom_point(aes(n_sim, gmd, color = mean)) +
  geom_smooth(aes(n_sim, gmd)) +
  facet_wrap(~model)

var %>% 
  filter(bias_correct == TRUE,
         preprocess == "none") %>% 
  ggplot() +
  geom_smooth(aes(n_sim, gmd, color = model)) 



var %>% 
  group_by(model, bias_correct, preprocess, n_sim) %>% 
  select(-gauge_id) %>% 
  summarise_all(mean) %>% 
  arrange(n_sim, -V) %>% 
  print(n=Inf)




