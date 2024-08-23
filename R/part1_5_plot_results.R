# ------------------------------------------------------------------------------
# Code to plot most of the figures in the paper. 

library(sf)
library(readr)
library(stringr)
library(lubridate)
library(dplyr)
library(units)
library(ggplot2)
library(tidyr)
library(boot)
sf::sf_use_s2(FALSE)
source("R/kge_non_parametric.R")
source("R/model_averaging_functions.R")
source("R/HBV.R")



# ------------------------------------------------------------------------------
# read_data


# READ CARAVAN DATA
attr_file <- paste0("../caravan/attributes/camels/attributes_other_camels.csv")
caravan <- readr::read_csv(attr_file)

gauges <- sf::st_as_sf(caravan, coords = c("gauge_lon", "gauge_lat"))

catchments <- sf::read_sf("../caravan/shapefiles/camels/camels_basin_shapes.shp")
    

# READ GOODNESS-OF-FIT TABLES
files <- list.files("data/hbv05_simulations/camels/assessment_final/",
                    full.names = TRUE,
                    pattern = "_gof")

data <- lapply(files, \(f) {
    tbl <- fst::read_fst(f) # %>% 
    return(tbl)
            
}) %>% do.call(bind_rows, .) %>% 
    as_tibble() 

# count the sizes of full ensembles
n_hbv <- data %>%
  filter(n_sim == 1,
         sig == "full",
         bias_correct == FALSE,
         period == "test") %>%
  group_by(gauge_id) %>%
  summarise(n_sim = n())


# GOODNESS-OF-FIT FROM CONSTRAINED ENSEMBLES
file_list <- list.files("data/hbv05_simulations/camels/const_final/",
                        full.names = TRUE)

const <- file_list[grepl("gof.fst", file_list)]

const <- lapply(const, \(f2) {
  x <- fst::read_fst(f2)
  x <- as_tibble(x)
  return(x)
}) %>% do.call(bind_rows, .)


# read weights
files <- list.files("data/hbv05_simulations/camels/assessment_final/",
                    full.names = TRUE,
                    pattern = "_weights")

weights <- lapply(files, \(f) {
  tbl <- fst::read_fst(f)
  
  return(tbl)
  
}) %>% do.call(bind_rows, .) %>% 
  as_tibble() 


# read hydrological signatures
files <- list.files("data/hbv05_simulations/camels/hydsig_final/",
                    full.names = TRUE)

sigs <- lapply(files, \(f) {
  tbl <- fst::read_fst(f) 
  return(tbl)
  
}) %>% do.call(bind_rows, .) %>% 
  as_tibble() 


sigs <- sigs %>%
  select(1:16) %>%
  left_join(n_hbv %>% rename(n_hbv = n_sim),
            by = "gauge_id")


# setuo for plotting
# colours for plot
cols <- inlmisc::GetColors(10)

# plotting order
model_levels <- c("hbv", "ensemble_mean", "best", " ", "CLS", "NNLS", "NNLS2",
                  "GRB", "GRA", "GRC")
sig_levels <- c("q0.05", "q0.25", "median", "q0.75", "q0.95",
                "cv", "autocor_1day", "autocor_7day", "full")
preprocess_levels <- c("none", "log", "sqrt")#, "boxcox")

# apply plotting order to the data
# change names of hbv simulations to simply "hbv"
# rnp to bounded rnp (rnp here is the non-parametric KGE)
# round numbers to two decimals
data <- data %>%
    filter(!sig %in% c("mean", "q0.5")) %>% # remove redundant
    mutate(model = ifelse(grepl("fold", model), "hbv", model),
           group = paste(model, preprocess, bias_correct, sep = "_"),
           model = factor(model, levels = model_levels),
           sig = factor(sig, levels = sig_levels),
           preprocess = factor(preprocess, levels = preprocess_levels),
           nrmse = round(nrmse, 2),
           rnp = round(rnp / (2-rnp),2),
           var = round(var, 2),
           bias = round(bias, 2),
           dyn = round(dyn, 2)) %>% 
    ungroup() 

# bound KGE to -1 to 1: https://www.researchgate.net/publication/282319260_A_bounded_version_of_the_Nash-Sutcliffe_criterion_for_better_model_assessment_on_large_sets_of_basins

# thresholds for plotting for KGE = 0.75, 0.5, and -0.41 (equivalent to mean 
# flow benchmark)
kgeb75 <- 0.75 / (2-0.75)
kgeb50 <- 0.5 / (2-0.5)
kgeb0 <- -0.41 / (2- -0.41)




# ------------------------------------------------------------------------------
# ANALYSIS

# number of HBV simulations and ensembles: 116 360 sims
data %>% 
  filter(model == "hbv",
         sig == "full",
         period == "test",
         bias_correct == FALSE)

# how many ensemble fits with full ensembles: 18 108 optimised combinations
# (+ 421 721 constrained combinations) = 439829
data %>% 
  filter(model != "hbv",
         sig == "full",
         period == "test")


# ------------------------------------------------------------------------------
# mean model performance vs ensemble performance

# mean performance of hbv
# RNP
data %>% 
  filter(model == "hbv",
         sig == "full",
         period == "test") %>% 
  group_by(bias_correct) %>% 
  summarise(hbv = mean(rnp, na.rm=TRUE),
            iqr = IQR(rnp, na.rm=TRUE),
            q025 = quantile(rnp, probs = 0.25, na.rm=TRUE),
            q075 = quantile(rnp, probs = 0.75, na.rm=TRUE))

# VARIABILITY
data %>% 
  filter(model == "hbv",
         sig == "full",
         period == "test") %>% 
  group_by(bias_correct) %>% 
  summarise(hbv = mean(var, na.rm=TRUE),
            iqr = IQR(var, na.rm=TRUE),
            q025 = quantile(var, probs = 0.25, na.rm=TRUE),
            q075 = quantile(var, probs = 0.75, na.rm=TRUE))

# DYNAMICS
data %>% 
  filter(model == "hbv",
         sig == "full",
         period == "test") %>% 
  group_by(bias_correct) %>% 
  summarise(hbv = mean(dyn, na.rm=TRUE),
            iqr = IQR(dyn, na.rm=TRUE),
            q025 = quantile(dyn, probs = 0.25, na.rm=TRUE),
            q075 = quantile(dyn, probs = 0.75, na.rm=TRUE))


# ------------------------------------------------------------------------------
# PLOT FIGURES


hbv <- data %>% 
    filter(model == "hbv",
           sig == "full",
           period == "test") %>% 
    group_by(gauge_id, bias_correct) %>% 
    summarise(rnp_hbv = mean(rnp))




# ---------
# FIGURE 3
data %>% 
    filter(sig == "full",
           preprocess == "none",
           # bias_correct == TRUE,
           model != "hbv") %>% 
    select(gauge_id, model, bias_correct, preprocess, period, sig, rnp) %>% 
    pivot_wider(names_from = period, values_from = rnp)  %>% 
    mutate(model = factor(model, levels = model_levels[-1])) %>% 
    ggplot() +
    geom_abline(slope = 1, intercept = 0) +
    geom_hline(yintercept = c(kgeb0), linetype = 2) +
    geom_vline(xintercept = c(kgeb0), linetype = 2) +
    geom_point(aes(train, test, color = bias_correct), 
               alpha = 0.25, shape = 16) +
    geom_point(data = data %>% 
                   filter(sig == "full",
                          preprocess == "none",
                          model != "hbv") %>% 
                   group_by(model, bias_correct, period) %>% 
                   summarise(rnp = median(rnp, na.rm=TRUE)) %>% 
                   pivot_wider(names_from = period, values_from = rnp) %>% 
                   mutate(model = factor(model, levels = model_levels[-1])) %>% 
                   arrange(bias_correct, model),
               aes(train, test, shape = bias_correct),
               size = 2) +
    scale_color_manual(values = cols[c(9,4)]) +
    scale_shape_manual( values = 3:4) +
    facet_wrap(~model) +
    theme_minimal() +
    labs(x = "KGEb in training period",
         y = "KGEb in testing period")

# ggsave("performance_in_train_test.pdf", height  = 6, width = 14)


# ------------------------------------------------------------------------------
# FIGURE 4
# Get supplementary figures 1 and 2 by changing the preprocess filter to
# "log" or "sqrt".
sig_obs <- sigs %>% 
  filter(model == "observations")

sigs %>% 
  select(-mean, -q0.5) %>%
  filter(period == "test",
         model != "observations",
         preprocess == "none", # "log" for suppl. fig. 1, and "sqrt" for suppl. fig. 2
         n_hbv == n_sim) %>% 
  gather(sig, value, -gauge_id, -model, -n_sim,
         -bias_correct,-preprocess,-period) %>% 
  mutate(model = ifelse(grepl("fold", model), "hbv", model),
         model = factor(model, levels = c("observations", model_levels)),
         sig = factor(sig, levels = sig_levels)) %>%
  left_join(sig_obs %>% 
              filter(period == "test",
                     preprocess == "none") %>% 
              gather(sig, value, -gauge_id, -model, -n_sim,
                               -bias_correct,-preprocess,-period) %>% 
              select(-model, -n_sim, -preprocess, -bias_correct, -period, 
                     obs = value),
            by = c("gauge_id", "sig")) %>% 
  filter(model %in% c("ensemble_mean", "CLS", "NNLS2", "GRC")) %>% 
  mutate(sig = factor(sig, levels = sig_levels)) %>% 
  ggplot() +
  geom_point(aes(value, obs, color = bias_correct),
             alpha = 0.3, pch = 16, size = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~model + sig, scales = "free", ncol = 9) +
  scale_color_manual(values = cols[c(9,4)]) +
  theme_minimal() +
  theme(strip.text.x = element_text(margin = margin(0,0,0,0, "cm"),
                                    size = 5)) +
  labs(x = "Modelled hydrological signature",
       y = "Observed hydrological signature")

# ggsave("full_ens_hydsig.pdf", height  = 6, width = 14)


# correlations to manually add to the plot
correlations_for_the_plot <- sigs %>% 
  select(-mean, -q0.5) %>%
  filter(period == "test",
         model != "observations",
         preprocess == "none",
         # cv < 60,
         # cv > -30,
         n_hbv == n_sim) %>% 
  gather(sig, value, -gauge_id, -model, -n_sim,
         -bias_correct,-preprocess,-period) %>% 
  mutate(model = ifelse(grepl("fold", model), "hbv", model),
         model = factor(model, levels = c("observations", model_levels)),
         sig = factor(sig, levels = sig_levels)) %>%
  left_join(sig_obs %>% 
              filter(period == "test",
                     preprocess == "none") %>% 
              gather(sig, value, -gauge_id, -model, -n_sim,
                     -bias_correct,-preprocess,-period) %>% 
              select(-model, -n_sim, -preprocess, -bias_correct, -period, 
                     obs = value),
            by = c("gauge_id", "sig")) %>% 
  filter(model %in% c("hbv", "ensemble_mean", "CLS", "NNLS2", "GRC")) %>% 
  mutate(sig = factor(sig, levels = sig_levels)) %>% 
  group_by(model, sig, bias_correct) %>% 
  summarise(r = round(cor(value, obs), 2)) %>% 
  pivot_wider(names_from = sig, values_from = r)

# ------------------------------------------------------------------------------
# FIGURE 5
data2 <- data %>% 
  ungroup() %>% 
  filter(#neg <= 0.05, ## unfair to omit timeseries with many negatives
    !sig %in% c("mean", "q0.5"),
    period == "test") %>% 
  group_by(model, gauge_id, bias_correct, preprocess, sig) %>%
  summarise(mean_rnp = mean(rnp, na.rm=TRUE)) %>%
  pivot_wider(names_from = model, values_from = mean_rnp) %>% 
  gather(model, mean_rnp, -gauge_id, -bias_correct, 
         -preprocess, -sig,-ensemble_mean) %>% 
  mutate(ss = (mean_rnp - ensemble_mean) / (1-ensemble_mean),
         model = factor(model, levels = model_levels)) %>% 
  ungroup()

data2 %>% 
  filter(bias_correct ==  TRUE,
         model %in% c("CLS", "NNLS2", "GRC"),
         sig == "full") %>% 
  ggplot() +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  geom_vline(xintercept = kgeb0, linetype = 2) +
  geom_hline(yintercept = kgeb0, linetype = 2) +
  geom_point(aes(ensemble_mean, mean_rnp, color = preprocess),
             size = 0.8, alpha = 0.3, pch = 16) +
  geom_smooth(aes(ensemble_mean, mean_rnp, color = preprocess),
              size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 1) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  facet_wrap(~model) +
  scale_color_manual(values = cols[c(4,8,10)]) +
  theme_minimal() +
  labs(x = "Ensemble mean KGEb",
       y = "KGEb")

# ggsave("preprocess_vs_ensemble_mean.pdf", height  = 4, width = 12)


# ------------------------------------------------------------------------------
# FIGURE 6
data2 <- data %>% 
  filter(sig == "full",
         period == "test") %>% 
  group_by(model, gauge_id, bias_correct, preprocess) %>% 
  summarise(mean_rnp = mean(rnp, na.rm=TRUE)) %>% 
  pivot_wider(names_from = model, values_from = mean_rnp) %>% 
  mutate(ss = (ensemble_mean - best) / (1-best)) %>% 
  ungroup()


data %>% 
  filter(!model %in% c("hbv", "best", "ensemble_mean", "NNLS", "GRB", "GRA"),
         sig == "full",
         period == "test") %>% 
  left_join(select(data2, gauge_id, bias_correct, preprocess, ensemble_mean), 
            by = c("gauge_id", "bias_correct", "preprocess")) %>% 
  group_by(gauge_id, model, preprocess, bias_correct) %>% 
  summarise(em_better = ensemble_mean >= rnp,
            em_ss = (rnp - ensemble_mean)/(1-ensemble_mean),
            n_sim = mean(n_sim),
            rnp = mean(rnp, na.rm=TRUE)) %>% 
  ggplot() +
  geom_point(aes(n_sim, em_ss, color = bias_correct), alpha = 0.1) +
  geom_smooth(aes(n_sim, em_ss, color = bias_correct)) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-1,1), xlim = c(0,1000), expand = FALSE) +
  # ylim(c(-1,1)) +
  facet_wrap(~model, scales = "free_y") +
  scale_color_manual(values = cols[c(9,4)]) +
  theme_minimal()

# ggsave("skill_score_vs_em.pdf", height  = 4, width = 12)




# ------------------------------------------------------------------------------
# SUPPLEMENTARY FIGURE S3
pdata <- data %>%
  filter(model != "hbv",
         preprocess == "none",
         sig == "full",
         period == "test") %>%
  select(gauge_id, model, bias_correct, rnp) %>%
  # left_join(hbv, by = c("gauge_id", "bias_correct")) %>%
  left_join(data %>%
              filter(model == "hbv", sig == "full", period == "test") %>%
              select(gauge_id, bias_correct, rnp_hbv = rnp),
            by = c("gauge_id", "bias_correct"))
p <- pdata %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = c(kgeb0), linetype = 2) +
  geom_vline(xintercept = c(kgeb0), linetype = 2) +
  geom_point(aes(rnp, rnp_hbv, color = bias_correct),
             alpha = 0.01, pch = 16, size = 0.5) +
  # geom_boxplot(aes(rnp, rnp_hbv, color = bias_correct,
  #                  group = gauge_id), fill = "transparent") +
  geom_point(data = pdata %>%
               group_by(model, bias_correct) %>%
               summarise(rnp = median(rnp, na.rm=TRUE),
                         rnp_hbv = median(rnp_hbv, na.rm=TRUE)),
             aes(rnp, rnp_hbv, shape = bias_correct),
             size = 2) +
  coord_cartesian(ylim = c(-1, 1),
                  xlim = c(-1, 1)) +
  facet_wrap(~model) +
  scale_color_manual(values = cols[c(9,4)]) +
  scale_shape_manual( values = 3:4) +
  theme_minimal() +
  # theme(legend.position = "none") +
  labs(x = "KGEb of optimal combinations (n = 482)",
       y = "KGEb of individual ensemble members")
p
ggsave("opt_comb_vs_hbv_none.png", p, height  = 6, width = 14)











