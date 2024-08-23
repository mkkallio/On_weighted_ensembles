#### Draws heavily from ForecastComb package
combinations <- function(flowmat, type="CLS", ...) {
    
    # remove collinear timeseries
    flowmat <- remove_collinear(flowmat)
    
    ncols <- ncol(flowmat)
    obs <- flowmat$observations %>% unname #%>% units::drop_units()
    predmat <- flowmat[,-c(1, ncols)] %>% as.matrix
    models <- colnames(predmat)
    intercept <- 0
    dates <- dplyr::pull(flowmat, "Date")
    
    test <- is.function(type)
    if(test) {
        method <- "Custom function optimisation"
        args <- list(...)
        
        test <- is.null(args$control)
        if(test) args$control <- list()
        
        test <- is.null(args$hessian)
        if(test) args$hessian <- FALSE
        
        test <- is.null(args$method)
        if(test) args$method <- "Nelder-Mead"
        
        test <- is.null(args$lower)
        if(test) args$lower <- -Inf
        
        test <- is.null(args$upper)
        if(test) args$upper <- Inf
        
        tst <- optim(par = rep(1/ncol(predmat), ncol(predmat)),
                     fn = type,
                     gr = args$gr,
                     obs, predmat,
                     dates,
                     method = args$method,
                     lower = args$lower,
                     upper = args$upper,
                     control = args$control,
                     hessian = args$hessian)
        
        weights <- tst$par
        
        
    } else if(type == "CLS") {
        # This is borrows from ForecastComb::comb_CLS(), but edited  
        # according to https://stackoverflow.com/a/28388394. 
        method <- c("Constrained Least Squares; 0 < weights < 1; ",
                    "sum(weights) = 1; intercept = 0")
        p <- NCOL(predmat)
        Rinv <- solve(safe_chol(t(predmat) %*% predmat))
        C <- cbind(rep(1, p), diag(p))
        b <- c(1, rep(0, p))
        d <- t(obs) %*% predmat
        nn2 <- sqrt(norm(d,"2"))
        qp1 <- quadprog::solve.QP(Dmat = Rinv*nn2, 
                                  factorized = TRUE, 
                                  dvec = d/(nn2^2), 
                                  Amat = C, 
                                  bvec = b, 
                                  meq = 1)
        weights <- unname(qp1$sol)
        
    } else if (type == "GRA") {
        method <- "Granger-Ramanathan Type A - Least Squares; intercept = 0"
        lin_model <- lm(obs ~ 0 + predmat)
        weights <- lin_model$coef
        
    } else if (type == "GRB") {
        method <- paste0("Granger-Ramanathan Type B - Least Squares; ",
                         "sum(weights) = 1; intercept = 0")
        
        # solved with pracma
        # Equality constraint: We want the sum of the coefficients to be 1.
        # I.e. Aeq x == beq  
        p <- ncol(predmat)
        Aeq <- matrix(rep(1, p), nrow= 1)
        beq <- c(1)
        
        # Lower and upper bounds of the parameters
        lb <- rep(-1e6, p)
        ub <- rep(1e6, p)
        
        # And solve:
        weights <- pracma::lsqlincon(predmat, obs, Aeq= Aeq, beq= beq, lb= lb, ub= ub)
        
    } else if (type %in% c("GRC", "OLS") ) {
        method <- paste0("Granger-Ramanathan Type C - Least Squares; ",
                         "Ordinary Least Squares")
        lin_model <- lm(obs ~ predmat)
        weights <- lin_model$coef[-1]
        intercept <- lin_model$coef[1]
        
    } else if (type == "NNLS") {
        method <- paste0("Non-Negative Least Squares; weights > 0;",
                         "intercept = 0")
        
        mod <- glmnet::glmnet(predmat, 
                              obs, 
                              lambda = 0, 
                              lower.limits = 0, 
                              intercept = FALSE)
        weights <- coef(mod)[-1]
        
    } else if (type == "NNLS2") {
        method <- paste0("Non-Negative Least Squares; weights > 0")
        
        mod <- glmnet::glmnet(predmat, 
                              obs, 
                              lambda = 0, 
                              lower.limits = 0, 
                              intercept = TRUE)
        weights <- coef(mod)[-1]
        intercept <- coef(mod)[1]
        
    }  else if (type == "BG") {
        method <- "Bates-Granger; sum(weights) = 1; intercept = 0"
        errmat <- obs - predmat
        sample_msqu_pred_error <- (t(errmat) %*% errmat)/length(obs)
        weights <- diag(sample_msqu_pred_error)^(-1)/
            sum(diag(sample_msqu_pred_error)^(-1))
        
    } else if (type == "InvW") {
        method <- "Inverse Rank combination; sum(weights) = 1; intercept = 0"
        errmat <- obs - predmat
        sse <- colSums((errmat)^2)
        ranks <- rank(sse)
        inv_ranks <- 1/ranks
        weights <- inv_ranks/sum(inv_ranks)
        
    } else if (type == "EIG1") {
        method <- paste0("Standard Eigenvector combination;",
                         " sum(weights) = 1; intercept = 0")
        errmat <- obs - predmat
        sample_msqu_pred_error <- (t(errmat) %*% errmat)/length(obs)
        eigen_decomp <- eigen(sample_msqu_pred_error)
        ds <- colSums(eigen_decomp$vectors)
        adj_eigen_vals <- eigen_decomp$values/(ds^2)
        min_idx <- which.min(adj_eigen_vals)
        weights <- eigen_decomp$vectors[, min_idx]/ds[min_idx]
        
    } else if (type == "EIG2") {
        method <- paste0("Bias-corrected Eigenvector combination;",
                         " sum(weights) = 1")
        mean_obs <- mean(obs)
        mean_preds <- colMeans(predmat)
        centered_obs <- obs - mean_obs
        centered_preds <- scale(predmat, scale = FALSE)
        omega_matrix <- t(centered_obs - centered_preds) %*% 
            (centered_obs - centered_preds)/length(obs)
        eigen_decomp <- eigen(omega_matrix)
        ds <- colSums(eigen_decomp$vectors)
        adj_eigen_vals <- eigen_decomp$values/(ds^2)
        min_idx <- which.min(adj_eigen_vals)
        weights <- eigen_decomp$vectors[, min_idx]/ds[min_idx]
        intercept <- as.numeric(mean_obs - t(mean_preds) %*% weights)
        
    } else if (type == "best") {
        method <- "Best individual ensemble member"
        
        rmse <- rep(NA, ncol(predmat))
        for(i in 1:ncol(predmat)) {
            rmse[i] <- hydroGOF::rmse(predmat[,i], obs)
        }
        min <- which(rmse == min(rmse))
        weights <- rep(0, ncol(predmat))
        weights[min] <- 1
        names(weights) <- colnames(predmat)
        intercept <- 0
        
    } else if (type == "BMA") {
        
        method <- c("Bayesian Model Averaging; gaussian family")
        
        BMA_fit <- BMA::bic.glm(obs ~ .,
                                data = data.frame(cbind(obs, predmat)),
                                glm.family = "gaussian")
        
        weights <- BMA_fit$postmean[-1]
        intercept <- BMA_fit$postmean[1]
        
    } else {
        stop("Unknown optimization method: ", type)
    }
    
    names(weights) <- models
    output <- list(Method = method,
                   Models = models,
                   Weights = weights,
                   Intercept = intercept)
    
    return(output)
}

# from package 'lme4'
safe_chol <- function(m) {
    if (all(m==0)) return(m)
    if (nrow(m)==1) return(sqrt(m))
    if (all(dmult(m,0)==0)) {  ## diagonal
        return(diag(sqrt(diag(m))))
    }
    ## attempt regular Chol. decomp
    if (!inherits(try(cc <- chol(m),silent=TRUE),"try-error"))
        return(cc)
    ## ... pivot if necessary ...
    cc <- suppressWarnings(chol(m,pivot=TRUE))
    oo <- order(attr(cc,"pivot"))
    cc[,oo]
}

# from package 'lme4'; needed in safe_chol()
dmult <- function(m,s) {
    diag(m) <- diag(m)*s
    m
}
# remove collinear timeseries. modified from ForecastComb::remove_collinear
remove_collinear <- function (flow) {
    
    ncols <- ncol(flow)
    obs <- flow$observations
    predmat <- flow[,-c(1,ncols)] %>% as.matrix
    
    # remove collinear timeseries according to RMSE, until prediction matrix
    # is full rank.
    repeat {
        repeat_this <- Matrix::rankMatrix(predmat)[1] == ncol(predmat)
        if (repeat_this == TRUE) {
            break
        }
        ranks <- rep(NA, ncol(predmat))
        for (i in 1:ncol(predmat)) {
            ranks[i] <- Matrix::rankMatrix(predmat[, -i])[1]
        }
        maxrank <- which(ranks == max(ranks))
        remove <- rep(0, ncol(predmat))
        for (i in maxrank) {
            remove[i] <- sqrt(mean((predmat[,i] - obs)^2, na.rm = TRUE))
        }
        remove <- which(remove == max(remove))[1]
        predmat <- predmat[, -remove]
    }
    
    keep <- colnames(flow) %in% c("Date", "observations", colnames(predmat))
    output <- flow[,which(keep)]
    return(output)
}


bc_trans <- function(x) {
    
    # lambda <- MASS::boxcox(x ~ 1, interp = TRUE, plotit = FALSE)
    # lambda <- lambda$x[which.max(lambda$y)]
    lambda <- 0.2
    if(lambda == 0) {
        trans <- log(x)
    } else {
        trans <- (x^lambda-1)/lambda
    }
    
    return(list(x = trans,
                lambda = lambda))
}

bc_back <- function(x, lambda) {
    
    if(lambda == 0) {
        x <- exp(x)
    } else {
        x <- (lambda*x+1)^(1/lambda)
    }
    
   # x <- exp(log(lambda * x + 1) / lambda)
}


process_combs <- function(combmat, combs, preprocess, train, test) {
    
    
    # gather weights
    models <- combs$Models
    weights <- combs$Weights
    #names(weights) <- models
    weights[is.na(weights)] <- 0
    
    # get intercept
    int <- combs$Intercept
    names(int) <- "intercept"
    
    # Create the optimized timeseries
    w <- c(int, weights)
    modelind <- which(colnames(combmat) %in% models, arr.ind=TRUE)
    p <- t(cbind(1,as.matrix(combmat[,modelind])))
    Opt <- as.vector(w %*% p)
    Opt <- round(Opt, 6)
    
    
    neg_train <- sum(Opt[train] < 0, na.rm = TRUE) / length(train)
    neg_test <- sum(Opt[test] < 0, na.rm = TRUE) / length(test)
    
    # back transform
    if(preprocess == "log") {
        opt_train <- exp(Opt[train])
        opt_test <- exp(Opt[test])
        
        neg_train <- 0
        neg_test <- 0
        
    } else if(preprocess == "sqrt") {
        opt_train <- Opt[train]^2
        opt_test <- Opt[test]^2
    } else if(preprocess == "none") {
        opt_train <- Opt[train]
        opt_test <- Opt[test]
    }
    
    
    out <- list(weights = w,
                opt_train = opt_train,
                opt_test = opt_test,
                neg_train = neg_train,
                neg_test = neg_test)
    return(out)
}

assess_comb <- function(combs, obs, train, test) {
    
    obs_train <- obs[train]
    obs_test <- obs[test]
    
    # TRAIN PERIOD
    rnp <- RNP(combs$opt_train, obs_train) %>% 
        t() %>% 
        as_tibble
    nrmse <- hydroGOF::rmse(combs$opt_train, obs_train) / mean(obs_train, 
                                                               na.rm=TRUE)
    rnp_train <- cbind(combs$neg_train, nrmse, rnp)
    names(rnp_train) <- c("neg", "nrmse", "var", "bias", "dyn", "rnp")
    
    # TEST PERIOD
    rnp <- RNP(combs$opt_test, obs_test) %>% 
        t() %>% 
        as_tibble
    nrmse <- hydroGOF::rmse(combs$opt_test, obs_test) / mean(obs_test, 
                                                               na.rm=TRUE)
    rnp_test <- cbind(combs$neg_test, nrmse, rnp)
    names(rnp_test) <- c("neg", "nrmse", "var", "bias", "dyn", "rnp")
    
    out <- list(train = rnp_train, test = rnp_test)
    return(out)
}



run_comb <- function(methods, 
                     predmat, 
                     preprocess = "none", 
                     biascor = FALSE,
                     traintest = NULL,
                     return_ts = FALSE,
                     timeout = 10) {
    
    # combmat <- predmat
    
    test <- is.null(traintest)
    if(test) {
      n <- floor(nrow(predmat)/2)
      train <- 1:n
      test <- (n+1):nrow(predmat)
    } else {
      train <- traintest$train
      test <- traintest$test
    }
 
    
    n_sim <- length(starts_with("fold", vars = names(predmat)))
    
    #preprocess
    observations <- predmat$observations
    combmat <- do_preprocess(predmat, preprocess, biascor)
    if(preprocess == "boxcox") {
        lambda <- combmat$lambda
        combmat <- combmat$combmat
    }
    
    # --------------------------------------------------------------------------
    # process all combination methods
    tbls <- list()
    tbls <- lapply(methods, \(m) {
        
        if(m == "ensemble_mean") {
            return(run_ens_mean(m, 
                                predmat, 
                                observations, 
                                preprocess, 
                                biascor,
                                traintest = list(train = train,
                                                 test = test)))
            next 
        }
        
        id <- paste(station, "f", m, sep = "_")
        
        # combinations
        combs <- try(R.utils::withTimeout(combinations(combmat[train,], type = m),
                                          timeout = timeout), 
                     TRUE)
        
        if(class(combs) == "try-error") return(NULL)
        

        combs <- process_combs(combmat, combs, preprocess, train, test)
        
        rnp <- assess_comb(combs, observations, train, test)
        
        
        tbl_train <- bind_cols(tibble::tibble(gauge_id = station,
                                              id = id, 
                                              model = m,
                                              n_sim = n_sim,
                                              preprocess = preprocess,
                                              bias_correct = biascor,
                                              period = "train",
                                              sig = "full"),
                               rnp$train)
        tbl_test <- bind_cols(tibble::tibble(gauge_id = station,
                                             id = id,
                                             model = m,
                                             n_sim = n_sim,
                                             preprocess = preprocess,
                                             bias_correct = biascor,
                                             period = "test",
                                             sig = "full"),
                              rnp$test)
        
        
        # ----------
        # HYDROLOGICAL SIGNATURES
        
        hydsig_train <- do_hydsig_dt(predmat$Date[train], 
                                     predmat$observations[train], 
                                     combs$opt_train, 
                                     m, 
                                     preprocess, 
                                     biascor,
                                     "train",
                                     id,
                                     n_sim)
        hydsig_test <- do_hydsig_dt(predmat$Date[test], 
                                    predmat$observations[test], 
                                    combs$opt_test, 
                                    m, 
                                    preprocess, 
                                    biascor,
                                    "test",
                                    id,
                                    n_sim)
        
          train_sig <- hydrological_signatures(combs$opt_train) %>% 
            tidyr::pivot_wider(names_from = sig, values_from = value)
          test_sig <- hydrological_signatures(combs$opt_test) %>% 
            tidyr::pivot_wider(names_from = sig, values_from = value)
          
          sig_train_tbl <- tibble(gauge_id = station,
                                  model = m,
                                  n_sim = n_sim,
                                  bias_correct = biascor,
                                  preprocess = preprocess,
                                  period = "train",
                                  train_sig)
          
          sig_test_tbl <- tibble(gauge_id = station,
                                 model = m,
                                 n_sim = n_sim,
                                 bias_correct = biascor,
                                 preprocess = preprocess,
                                 period = "test",
                                 test_sig)
          sig_tbl <- rbind(sig_train_tbl, sig_test_tbl)
          
          
          w <- as_tibble(t(combs$weights)) %>% 
            mutate(gauge_id = station,
                   id = id,
                   n_sim = n_sim,
                   model = m,
                   bias_correct = biascor,
                   preprocess = preprocess,
                   .before = 1)
        
        output <- bind_rows(tbl_train, tbl_test, hydsig_train, hydsig_test)
        
        output <- list(gof = output,
                       weights = w,
                       hydsig = sig_tbl)
        
        if(return_ts) {
          out_ts <- predmat[,c(1:2)]
          names(out_ts)[2] <- m
          out_ts[,2] <- NA
          out_ts[train,2] <- combs$opt_train
          out_ts[test,2] <- combs$opt_test
          output[["opt_ts"]] <- out_ts
        } 
        
        
        return(output)

    }) 
    
    gofs <- lapply(tbls, \(x) return(x$gof)) %>% do.call(bind_rows, .)
    weights <- lapply(tbls, \(x) return(x$weights)) %>% do.call(bind_rows, .)
    hydsigs <- lapply(tbls, \(x) return(x$hydsig)) %>% do.call(bind_rows, .)
    
    if(return_ts) {
      ts <- lapply(tbls, \(x) return(x$opt_ts[,2])) %>% do.call(bind_cols, .)
      ts <- bind_cols(ts, predmat[, "observations"])
      ts <- bind_cols(predmat[,"Date"], ts)
      # names(ts) <- c("Date", methods, "observations")
    } else {
      ts <- NULL
    }
    
    return(list(gof = gofs,
                weights = weights,
                hydsig = hydsigs,
                opt_ts = ts))
    
}


run_ens_mean <- function(methods, 
                         combmat,
                         observations,
                         preprocess = "none", 
                         biascor = FALSE,
                         traintest = NULL,
                         return_ts = FALSE,
                         timeout = 10) {
  
  # combmat <- predmat
  
  test <- is.null(traintest)
  if(test) {
    n <- floor(nrow(predmat)/2)
    train <- 1:n
    test <- (n+1):nrow(predmat)
  } else {
    train <- traintest$train
    test <- traintest$test
  }
    
    n_sim <- length(starts_with("fold", vars = names(combmat)))
    
    id <- paste(station, "f", "ensemble_mean", sep = "_")
    
    #preprocess
    combmat <- do_preprocess(combmat, preprocess, biascor)
    if(preprocess == "boxcox") {
        lambda <- combmat$lambda
        combmat <- combmat$combmat
    }

    
    Opt <- rowMeans(combmat[, -c(1, ncol(combmat))])
    
    # back transform
    if(preprocess == "log") {
        Opt <- exp(Opt)
    } else if(preprocess == "sqrt") {
        Opt <- Opt^2
    } #else if(preprocess == "boxcox") {
    #     Opt <- bc_back(Opt, lambda)
    # }
    
    
    # ASSESSMENT
    obs_train <- observations[train]
    obs_test <- observations[test]
    
    opt_train <- Opt[train]
    opt_test <- Opt[test]
    
    # TRAIN PERIOD
    rnp <- RNP(opt_train, obs_train) %>% 
        t() %>% 
        as_tibble
    nrmse <- hydroGOF::rmse(opt_train, obs_train) / mean(obs_train, 
                                                               na.rm=TRUE)
    neg_train <- sum(opt_train < 0) / length(opt_train)
    rnp_train <- c(neg_train, nrmse, rnp)
    names(rnp_train) <- c("neg", "nrmse", "var", "bias", "dyn", "rnp")
    
    # TEST PERIOD
    rnp <- RNP(opt_test, obs_test) %>% 
        t() %>% 
        as_tibble
    nrmse <- hydroGOF::rmse(opt_test, obs_test) / mean(obs_test, 
                                                             na.rm=TRUE)
    neg_test <- sum(opt_test < 0) / length(opt_test)
    rnp_test <- c(neg_test, nrmse, rnp)
    names(rnp_test) <- c("neg", "nrmse", "var", "bias", "dyn", "rnp")
    
    
    tbl_train <- bind_cols(tibble(gauge_id = station,
                                  id = id,
                            model = "ensemble_mean",
                            n_sim = n_sim,
                            preprocess = preprocess,
                            bias_correct = biascor,
                            period = "train",
                            sig = "full"),
                     rnp_train)
    
    tbl_test <- bind_cols(tibble(gauge_id = station,
                                 id = id,
                                  model = "ensemble_mean",
                                 n_sim = n_sim,
                                  preprocess = preprocess,
                                  bias_correct = biascor,
                                  period = "test",
                                  sig = "full"),
                           rnp_test)
    
    # ----------
    # HYDROLOGICAL SIGNATURES
    
    
    hydsig_train <- do_hydsig_dt(combmat$Date[train], 
                        obs_train, 
                        opt_train, 
                        "ensemble_mean", 
                        preprocess, 
                        biascor,
                        "train",
                        id,
                        n_sim)
    
    hydsig_test <- do_hydsig_dt(combmat$Date[test], 
                              obs_test, 
                              opt_test, 
                              "ensemble_mean", 
                              preprocess, 
                              biascor,
                              "test",
                              id,
                              n_sim)
    
    gofs <- bind_rows(tbl_train, tbl_test, hydsig_train, hydsig_test)
    
    
    train_sig <- hydrological_signatures(opt_train) %>% 
      tidyr::pivot_wider(names_from = sig, values_from = value)
    test_sig <- hydrological_signatures(opt_test) %>% 
      tidyr::pivot_wider(names_from = sig, values_from = value)
    
    sig_train_tbl <- tibble(gauge_id = station,
                            model = "ensemble_mean",
                            n_sim = n_sim,
                            bias_correct = biascor,
                            preprocess = preprocess,
                            period = "train",
                            train_sig)
    
    sig_test_tbl <- tibble(gauge_id = station,
                           model = "ensemble_mean",
                           n_sim = n_sim,
                           bias_correct = biascor,
                           preprocess = preprocess,
                           period = "test",
                           test_sig)
    sig_tbl <- rbind(sig_train_tbl, sig_test_tbl)
    
    
    if(return_ts) {
      ts <- lapply(tbls, \(x) return(x$opt_ts[,2])) %>% do.call(bind_cols, .)
      ts <- bind_cols(ts, predmat[, "observations"])
      ts <- bind_cols(predmat[,"Date"], ts)
    } else {
      ts <- NULL
    }
    
    return(list(gof = gofs,
                weights = NULL,
                hydsig = sig_tbl,
                opt_ts = ts))
    
  
}



evaluate_individual <- function(predmat, biascor = FALSE) {
    
    
    n <- floor(nrow(predmat)/2)
    train <- 1:n
    test <- (n+1):nrow(predmat)
    
    cols <- 2:(ncol(predmat)-1)
    
    # bias correct with linear transformation
    if(biascor) {
        predmat[,c(-1,-ncol(predmat))] <- lapply(predmat[,c(-1,-ncol(predmat))], 
                                                 function(x){
                                                     x <- x*mean(predmat$observations)/mean(x)
                                                 })
        
    }    
    
    obs_train <- predmat$observations[train]
    obs_test <- predmat$observations[test]
    
    stat_eval <- lapply(cols, \(ii) {
        
        m <- names(predmat)[ii]
        Opt <- pull(predmat, ii)
        
        opt_train <- Opt[train]
        opt_test <- Opt[test]
        
        id <- paste(station, m, sep = "_")
        
        # ASSESSMENT
        # TRAIN PERIOD
        rnp <- RNP(opt_train, obs_train) %>% 
            t() %>% 
            as_tibble
        nrmse <- hydroGOF::rmse(opt_train, obs_train) / mean(obs_train, 
                                                             na.rm=TRUE)
        neg_train <- sum(opt_train < 0) / length(opt_train)
        rnp_train <- c(neg_train, nrmse, rnp)
        names(rnp_train) <- c("neg", "nrmse", "var", "bias", "dyn", "rnp")
        
        # TEST PERIOD
        rnp <- RNP(opt_test, obs_test) %>% 
            t() %>% 
            as_tibble
        nrmse <- hydroGOF::rmse(opt_test, obs_test) / mean(obs_test, 
                                                           na.rm=TRUE)
        neg_test <- sum(opt_test < 0) / length(opt_test)
        rnp_test <- c(neg_test, nrmse, rnp)
        names(rnp_test) <- c("neg", "nrmse", "var", "bias", "dyn", "rnp")
        
        name <- names(simulations)[ii]
        
        tbl_train <- bind_cols(tibble(gauge_id = station,
                                      id = id,
                                      model = name,
                                      n_sim = 1,
                                      preprocess = "none",
                                      bias_correct = biascor,
                                      period = "train",
                                      sig = "full"),
                               rnp_train)
        
        tbl_test <- bind_cols(tibble(gauge_id = station,
                                     id = id,
                                     model = name,
                                     n_sim = 1,
                                     preprocess = "none",
                                     bias_correct = biascor,
                                     period = "test",
                                     sig = "full"), 
                              rnp_test)
        
        
        
        hydsig_train <- do_hydsig_dt(predmat$Date[train], 
                                  obs_train, 
                                  opt_train, 
                                  m, 
                                  preprocess = "none", 
                                  biascor,
                                  "train",
                                  id,
                                  1) 
        
        hydsig_test <- do_hydsig_dt(predmat$Date[test], 
                                 obs_test, 
                                 opt_test, 
                                 m, 
                                 preprocess = "none", 
                                 biascor,
                                 "test",
                                 id,
                                 1) 
        
        output <- bind_rows(tbl_train, tbl_test, hydsig_train, hydsig_test)
        
    }) %>% do.call(bind_rows, .)
    
    gc()
    
    return(stat_eval)
}


run_comb <- compiler::cmpfun(run_comb)
run_ens_mean <- compiler::cmpfun(run_ens_mean)
evaluate_individual <- compiler::cmpfun(evaluate_individual)


hydrological_signatures <- function(x, 
                                    quantiles = c(0.95, 0.75, 0.5, 0.25, 0.05)) {
    
    x <- x[!is.na(x)]
    
    quant <- quantile(x, quantiles, na.rm=TRUE)
    
    std <- sd(x, na.rm = TRUE)
    mn <- mean(x, na.rm=TRUE)
    cv <- std / mn
    
    med <- median(x, na.rm=TRUE)
    
    ac <- acf(x, lag.max = 7, plot = FALSE, na.action = na.omit)
    
    hyd_sigs <- c(quant, cv, mn, med, ac$acf[2], ac$acf[7])
    
    hyd_sig_names <- c(paste0("q", quantiles),
                       "cv",
                       "mean",
                       "median",
                       "autocor_7day",
                       "autocor_1day")
    
    out <- tibble(sig = hyd_sig_names,
                  value = hyd_sigs)
    return(out)
}

hydrological_signatures_dt <- function(x, 
                                    quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    
    x <- x[!is.na(x)]
    
    quant <- quantile(x, quantiles, na.rm=TRUE)
    
    std <- sd(x, na.rm = TRUE)
    mn <- mean(x, na.rm=TRUE)
    cv <- std / mn
    
    med <- median(x, na.rm=TRUE)
    
    ac <- acf(x, lag.max = 7, plot = FALSE, na.action = na.omit)
    
    hyd_sigs <- c(as.list(quant), list(cv, mn, med, ac$acf[2], ac$acf[7]))
    
    names(hyd_sigs) <- c(paste0("q", quantiles),
                       "cv",
                       "mean",
                       "median",
                       "autocor_7day",
                       "autocor_1day")
    
    return(hyd_sigs)
}


do_preprocess <- function(combmat, preprocess, biascor) {
    
    n <- ncol(combmat)
    
    # bias correct with linear transformation
    if(biascor) {
        combmat[,c(-1,-n)] <- lapply(combmat[,c(-1,-n)], 
             function(x){
                 x <- x*mean(combmat$observations, na.rm = TRUE)/mean(x, na.rm=TRUE)
             })
        
    }    
    
    # preprocessing transformation
    if(preprocess == "log") {
        combmat[,c(-1)] <- lapply(combmat[,c(-1)] ,
                                  function(x) {
                                      # x <- units::drop_units(x)
                                      x[x == 0] <- 1e-6
                                      x <- log(x)
                                      # x <- units::set_units(x, "mm/d")
                                  })
        
    } else if(preprocess == "sqrt") {
        combmat[,c(-1)] <- lapply(combmat[,c(-1)] ,
                                  function(x) {
                                      # x <- units::drop_units(x)
                                      x <- sqrt(x)
                                      # x <- units::set_units(x, "mm/d")
                                  })
        
    } else if(preprocess == "boxcox") {
        
        # x <- units::drop_units(combmat$observations)
        x <- combmat$observations
        x[x <= 0] <- 1e-6
        x <- bc_trans(x)
        lambda <- x$lambda
        # x <- units::set_units(x$x, "mm/d")
        combmat$observations <- x$x
        
        return(list(combmat = combmat, lambda = lambda))
    }
    
    return(combmat)
}

do_hydsig <- function(dates, obs, Opt, m, preprocess, biascor, period) {
    # ----------
    # HYDROLOGICAL SIGNATURES
    
    #annual timeseries
    hyd_sig <- tibble(Date = dates,
                      observations = obs,
                      !!m := Opt) %>%
        tidyr::gather(model, value, -Date) %>%
        mutate(value = value,
               year = lubridate::year(Date),
               month = month(Date),
               .before = 2) %>%
        select(-Date, -month) %>%
        group_by(year, model) %>%
        summarise(sigs = hydrological_signatures(value)) %>%
        tidyr::unnest(sigs)
    
    obs_sig <- hyd_sig %>%
        ungroup() %>%
        filter(model == "observations") %>%
        select(year, sig, obs = value)
    
    # model goodness-of-fit on hydrological signatures
    annual_gofs <- hyd_sig %>%
        filter(model != "observations") %>%
        left_join(obs_sig, by = c("year", "sig")) %>%
        group_by(model, sig) %>%
        summarise(value = RNP(value, obs)) %>%
        mutate(gauge_id = station,
               model = m,
               preprocess = preprocess,
               bias_correct = biascor,
               period = period,
               gof = c("var", "bias", "dyn", "rnp"),
               .before = "sig") %>%
        tidyr::pivot_wider(names_from = gof,
                           values_from = value)
    
    
    #full_ts
    hyd_sig <- tibble(Date = dates,
                      observations = obs,
                      !!m := Opt) %>%
        tidyr::gather(model, value, -Date) %>%
        select(-Date) %>%
        group_by(model) %>%
        summarise(sigs = hydrological_signatures(value)) %>%
        tidyr::unnest(sigs) 
    
    obs_sig <- hyd_sig %>%
        ungroup() %>%
        filter(model == "observations") %>%
        select(sig, obs = value)
    
    # model goodness-of-fit on hydrological signatures
    full_gofs <- hyd_sig %>%
        filter(model != "observations") %>%
        left_join(obs_sig, by = c("sig")) %>% 
        # tidyr::pivot_wider(values_from = x, names_from = q) %>%
        group_by(model, sig) %>%
        summarise(value = value / obs) %>%
        mutate(gauge_id = station,
               model = m,
               preprocess = preprocess,
               bias_correct = biascor,
               .before = "sig") %>% 
        rename(relative_value = value)
    
    
    output <- cbind(annual_gofs, full_gofs[,"relative_value"])
    return(output)
}


# using data.table much faster than using dplyr
do_hydsig_dt <- function(dates, obs, Opt, m, 
                         preprocess, biascor, period, 
                         id, n_sim) {
    # ----------
    # HYDROLOGICAL SIGNATURES
    
        hyd_sig <- data.table::data.table(Date = dates,
                              obs = obs,
                              opt = Opt)
        hyd_sig[, year := year(Date)]
        
        # full timeseries relative value
        opt_sig <- hyd_sig[, hydrological_signatures_dt(opt)]
        obs_sig <- hyd_sig[, hydrological_signatures_dt(obs)]
        
        relative <- t(opt_sig/obs_sig)
        
        # annual timeseries
        opt_sig <- hyd_sig[, hydrological_signatures_dt(opt), by = "year"]
        opt_sig <- data.table::melt(opt_sig, id.vars = "year",
                        variable.name = "sig",
                        value.name = "opt") 
        
        obs_sig <- hyd_sig[, hydrological_signatures_dt(obs), by = "year"]
        obs_sig <- data.table::melt(obs_sig, id.vars = "year",
                        variable.name = "sig",
                        value.name = "obs")
        
        gofs <- opt_sig[obs_sig, on = c("year", "sig")]
        gofs <- gofs[, as.list(RNP(opt, obs)), 
                    by = "sig"]
        
        # add columns
        gofs$relative_value <- as.vector(relative)
        gofs$gauge_id <- station
        gofs$id <- id
        gofs$model <- m
        gofs$n_sim <- n_sim
        gofs$preprocess <- preprocess
        gofs$bias_correct <- biascor
        gofs$period <- period
        
        # rename and reorder
        data.table::setnames(gofs, 
                 c("alpha_var", "beta_bias", "dyn_spearman", "RNP"),
                 c("var", "bias", "dyn", "rnp"))
        data.table::setcolorder(gofs, c("gauge_id","id", "model", "n_sim", 
                                        "preprocess", "bias_correct", "period", 
                                        "sig", "var", "bias", "dyn", "rnp", 
                                        "relative_value"))
        gofs <- as_tibble(gofs)
    
    return(gofs)
}
