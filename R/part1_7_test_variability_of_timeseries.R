# ------------------------------------------------------------------------------
# LOOK AT THE VARIABILITY OF TIMESERIES AS FUNCTION OF ENSEMBLE COMPOSITION

library(sf)
library(readr)
library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)
library(units)
library(foreach)
library(doParallel)
library(ggplot2)
library(patchwork)
library(ggridges)
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
# load data


# ------------------------------------------------------------------------------
# Assess variability of timeseries

attr_file <- paste0("../caravan/attributes/",
                    "camels", 
                    "/attributes_caravan_",
                    "camels",
                    ".csv")
gauges <- readr::read_csv(attr_file)

n_sims <- c(2, 5, 10, 15, 25, 50, 100, 150, 200)
stations <- unique(gauges$gauge_id)
gauges <- stations

cl <- parallel::makeCluster(10, outfile = "")
doParallel::registerDoParallel(cl)

result <- foreach(stat = seq_along(gauges),
                  .errorhandling = "stop",
                  .packages = c("dplyr", "lubridate")) %dopar% {
        
                    station <- stations[stat]
        # variability <- lapply(stations, \(station) {
        timeseries_var <- lapply(n_sims, \(ns) {
          
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
          
          # training and testing
          n <- floor(nrow(predmat)/2)
          train <- 1:n
          test <- (n+1):nrow(predmat)
          
          # get testing period
          predmat <- predmat[test,]
          
          
          # apply weights
          
          filename <- paste0("data/hbv05_simulations/camels/const_final/",
                              station, 
                             "_weights.fst")
          
          weights <- fst::read_fst(filename) %>% 
            as_tibble()
          
          cw <- weights %>% 
            filter(gauge_id == station,
                   n_sim == ns, 
                   preprocess == "none",
                   bias_correct == TRUE) %>% 
            arrange(model)
          
          
          cv <- function(x) {sd(x / mean(x))}
          pred_var <- lapply(unique(cw$model), \(m) {
            
            cwt <- cw %>% 
              filter(model == m)
            
            opt_frames <- lapply(1:nrow(cwt), \(i) {
              
              w <- unlist(cwt[i,8:ncol(cw)]) 
              w <- round(w, 3)
              w <- w[!is.na(w)]
              
              pm <- predmat %>% 
                select(names(w)) 
              
              p <- pm %>% 
                as.matrix() %>% 
                cbind(1, .) %>% 
                t()
              
              w <- c(cwt$intercept[i], w)
              
              pred <- as.vector(w %*% p)
              
              frames <- tibble(opt = pred)
              
            }) %>% do.call(cbind, .)
            
            model <- cwt$model[1]
            
            out <- tibble(min = apply(opt_frames, 1, min),
                          median = apply(opt_frames, 1, median),
                          max = apply(opt_frames, 1, max),
                          cv = apply(opt_frames, 1, cv),
                          gmd = apply(opt_frames, 1, gmd)) %>% 
              summarise_all(.funs = median) %>% 
              mutate(station = station,
                     model = model,
                     n_sim = cw$n_sim[1],
                     .before = 1)
            
            
          }) %>% do.call(rbind, .)
        }) %>% do.call(rbind, .)
        return(timeseries_var)
      }

variability <- do.call(rbind, result)

variability %>% 
  ggplot() +
  geom_smooth(aes(n_sim, gmd, color = model),
              method = "loess") +
  theme_minimal() 
  
  
  variability %>% 
    ggplot() +
    geom_line(aes(n_sim, gmd, group = station), alpha = 0.1) +
    theme_minimal() +
    coord_cartesian(ylim = c(0,3)) +
    facet_wrap(~model)

variability %>% 
  group_by(model, n_sim) %>% 
  summarise(mean_cv = mean(cv, na.rm=TRUE)) %>% 
  pivot_wider(names_from = model, values_from = mean_cv)



