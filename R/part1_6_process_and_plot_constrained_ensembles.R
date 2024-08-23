# ------------------------------------------------------------------------------
# PLOT FIGURES USING CONSTRAINED ENSEMBLES

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


# ------------------------------------------------------------------------------
# load data

file_list <- list.files("data/hbv05_simulations/camels/const_final/", 
                        full.names = TRUE)

weights <- file_list[grepl("weights.fst", file_list)]
gof <- file_list[grepl("gof.fst", file_list)]
# gof <- gof[!grepl("NNLS2", gof)]

weights <- lapply(weights, \(f) {
    stat <- stringr::word(f, -1, sep = "/")
    stat <- stringr::str_remove(stat, "_weights.fst")
    x <- fst::read_fst(f)
    x <- as_tibble(x)
    x <- mutate(x, gauge_id = stat, .before = 1)
    return(x)
}) %>% do.call(bind_rows, .) #%>%
    # rename(n_sim = nmod)

gof <- lapply(gof, \(f2) {
    x <- fst::read_fst(f2)
    x <- as_tibble(x) %>% 
        filter(period == "test")
    return(x)
}) %>% do.call(bind_rows, .) 

gc()

# read signatures

files <- list.files("data/hbv05_simulations/camels/hydsig_final/",
                    full.names = TRUE)
files <- files[!grepl("_const.fst", files)]
files <- files[!grepl("_full.fst", files)]

sigs <- lapply(files, \(f) {
  tbl <- fst::read_fst(f)  %>% 
  filter(period == "test")
  
  return(tbl)
  
}) %>% do.call(bind_rows, .) %>% 
  as_tibble() %>% 
  filter(period == "test")

sig_obs <- sigs %>% 
  filter(model == "observations")

sigs <- sigs %>%
  filter(n_sim %in% c(2,5,10,15,25,50,100,150,200))


sim_levels <- as.character(sort(unique(gof$n_sim)))
sim_levels <- as.character(sort(unique(weights$n_sim)))

model_levels <- c("hbv", "ensemble_mean", "best", "CLS", "NNLS", "NNLS2",
                  "GRB", "GRA", "GRC")
sig_levels <- c("q0.05", "q0.25", "median", "q0.75", "q0.95",
                "cv", "autocor_1day", "autocor_7day", "full")
preprocess_levels <- c("none", "log", "sqrt")#, "boxcox")

cols <- inlmisc::GetColors(10)

kgeb75 <- 0.75 / (2-0.75)
kgeb50 <- 0.5 / (2-0.5)
kgeb0 <- -0.41 / (2- -0.41)


# apply factors and round 
gof <- mutate(gof,
              rnp = rnp / (2-rnp))
gof <- gof %>% 
    filter(!sig %in% c("mean", "q0.5")) %>% 
    mutate(model = factor(model, model_levels),
           sig = factor(sig, sig_levels),
           preprocess = factor(preprocess, preprocess_levels),
           nrmse = round(nrmse, 2),
           var = round(var, 2),
           bias = round(bias, 2),
           dyn = round(dyn, 2),
           rnp = round(rnp, 2),
           relative_value = round(relative_value, 2))


# ------------------------------------------------------------------------------
# goodness of fit

# median performance as a function of ensemble size
gof %>%
    filter(sig == "full",
           !is.na(model),
           # model %in% c("GRC", "NNLS2"),
           bias_correct == TRUE,
           period == "test") %>% 
    # mutate(n_sim = as.numeric(stringr::word(model,2, sep = "_"))) %>%
    group_by(bias_correct, preprocess, model, n_sim) %>%
    summarise(rnp = median(rnp, na.rm=TRUE)) %>% 
    pivot_wider(names_from = c(n_sim), 
                values_from = rnp) 

gof %>%
  filter(sig == "full",
         !is.na(model),
         # model %in% c("GRC", "NNLS2"),
         bias_correct == TRUE,
         period == "test") %>% 
  # mutate(n_sim = as.numeric(stringr::word(model,2, sep = "_"))) %>%
  group_by(bias_correct, preprocess, model, n_sim) %>%
  summarise(median = median(rnp, na.rm=TRUE),
            q5 = quantile(rnp, 0.05, na.rm=TRUE),
            q95 = quantile(rnp, 0.95, na.rm=TRUE)) %>% 
  ggplot() +
  geom_ribbon(aes(x = n_sim, ymin = q5, ymax = q95, fill = model), 
              alpha = 0.25) +
  geom_line(aes(n_sim, median, color = model), size = 1) +
  facet_wrap(~preprocess)


gof %>%
  filter(sig == "full",
         preprocess == "none",
         !is.na(model)) %>% 
  group_by(n_sim, gauge_id, preprocess, bias_correct, model) %>%
  summarise(rnp = median(rnp)) %>%
  ggplot() +
  geom_line(aes(n_sim, rnp, group = gauge_id), alpha = 0.05, size = 1) +
  geom_line(data = gof %>%
              filter(sig == "full", 
                     preprocess == "none",
                     !is.na(model)) %>% 
              # mutate(n_sim = as.numeric(stringr::word(model,2, sep = "_"))) %>%
              group_by(n_sim, preprocess, bias_correct, model) %>%
              summarise(rnp = median(rnp)),
            aes(n_sim, rnp), color = "red", size = 1) +
  theme_minimal() +
  facet_wrap(~model, nrow = 2, drop = TRUE)



# ------------------------------------------------------------------------------
# FIGURE 7

sigs %>% 
  select(-mean, -q0.5) %>%
  filter(period == "test",
         model != "observations",
         preprocess == "log",
         bias_correct = TRUE,
         # cv < 30,
         # cv > -30,
         # n_sim %in% c(2, 10, 25, 50, 100, 200)
         ) %>% 
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
  mutate(sig = factor(sig, levels = sig_levels),
         relative_value = value/obs) %>% 
  group_by(n_sim, sig, model) %>%
  summarise(relative_value = median(relative_value, na.rm=TRUE)) %>%
  ggplot() +
  geom_abline(slope = 0, intercept = 1, size = 1) +
  # geom_point(aes(n_sim, relative_value, color = model)) +
  # geom_line(aes(n_sim, relative_value, colour = model,fill = model),
  #           size = 1) +
  geom_smooth(aes(n_sim, relative_value, colour = model,fill = model),
            method = "loess") +
  coord_cartesian(ylim = c(0,3), expand = FALSE) +
  facet_grid(~sig) +
  scale_color_manual(values = cols[c(3,5,8,10)]) +
  scale_fill_manual(values = cols[c(3,5,8,10)]) +
  scale_linetype_manual(values = c(1, 3, 5)) +
  theme_minimal() +
  labs(x = "Ensemble size", 
       y = "Relative value")

ggsave("ensemble_size_vs_hydsig_log.pdf", height = 4, width = 12)


# how many stations in each n_sim category
sigs %>% 
  group_by(gauge_id, n_sim) %>% 
  slice(1) %>% 
  group_by(n_sim) %>% 
  summarise(n_station = n())



# ------------------------------------------------------------------------------
# FIGURE 8
{
    p1 <- gof %>% 
        filter(#model == "GRC",
            !is.na(model),
            sig == "full",
            bias_correct == TRUE,
            period == "test") %>% 
        group_by(model, preprocess, n_sim) %>% 
        mutate(TFneg = ifelse(neg < 0.05, FALSE, TRUE),
               neg = ifelse(neg < 0.05, NA, neg)) %>% 
        summarise(prop_neg = sum(TFneg, na.rm=TRUE) / n(),
                  mean_neg = mean(neg, na.rm=TRUE)) %>% 
        ggplot() +
        geom_line(aes(n_sim, prop_neg, color = model, linetype = preprocess)) +
        labs(x = "Ensemble size",
             y = "Proportion of optimized timeseries \nwith >5% negative values") +
        scale_y_continuous(labels = scales::percent, limits = c(0,0.75)) +
        theme_minimal() +
        scale_colour_manual(values = cols[c(3,5,8,10)])
    
    p2 <- gof %>% 
        filter(#model == "GRC",
            !is.na(model),
            sig == "full",
            bias_correct == TRUE,
            period == "test") %>% 
        group_by(model, preprocess, n_sim) %>% 
        mutate(TFneg = ifelse(neg < 0.05, FALSE, TRUE),
               neg = ifelse(neg < 0.05, NA, neg)) %>% 
        summarise(prop_neg = sum(TFneg, na.rm=TRUE) / n(),
                  mean_neg = mean(neg, na.rm=TRUE)) %>% 
        ggplot() +
        geom_line(aes(n_sim, mean_neg, color = model, linetype = preprocess)) +
        geom_hline(yintercept = 0.05) +
        labs(x = "Ensemble size",
             y = "Mean proportion of negative values \nwhen >5% negative values") +
        scale_y_continuous(labels = scales::percent, limits = c(0,0.75)) +
        theme_minimal() +
        scale_colour_manual(values = cols[c(3,5,8,10)]) 
    
    
    p3 <- gof %>%
        filter(#model == "GRC",
            !is.na(model),
            sig == "full",
            bias_correct == TRUE,
            period == "test") %>% 
        # mutate(n_sim = as.numeric(stringr::word(model,2, sep = "_"))) %>%
        group_by(model, bias_correct, preprocess, n_sim) %>%
        summarise(n = sum(neg > 0.05) / n(),
                  neg = mean(neg[neg > 0.05])) %>% 
        ggplot() +
        geom_point(aes(n, neg, color = model, pch = preprocess)) +
        scale_x_continuous(labels = scales::percent, limits = c(0,0.75)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,0.75)) +
        labs(x = "Proportion of optimized timeseries \nwith >5% negative values 5%",
             y = "Mean proportion of negative values \nwhen >5% negative values 5%") +
        theme_minimal() +
        scale_colour_manual(values = cols[c(3,5,8,10)]) +
        scale_shape_manual(values = c(16,17,3))
    
}
p <- p1+p2+p3 + plot_layout(guides = "collect")
p
# ggsave("negative_vals.pdf", p, height  = 4, width = 8)





# ------------------------------------------------------------------------------
# FIGURE 9

weights <- mutate(weights, 
                  nw = rowSums(abs(weights[,-c(1:6)]) > 0.01, na.rm=TRUE),
                  .before =3)
gc()

p1 <- weights %>% 
    filter(model %in% c("GRC", "NNLS2", "CLS"),
           bias_correct == TRUE,
           preprocess == "none") %>% 
    mutate(preprocess = factor(preprocess, preprocess_levels),
           model = factor(model, levels = c("CLS", "NNLS2", "GRC"))) %>% 
    select(1:7) %>% 
    mutate(group = paste(preprocess, n_sim, sep = "_")) %>% 
    ggplot() +
    geom_boxplot(aes(n_sim, nw, group = group, 
                     # color = preprocess, 
                     # fill = preprocess
                     ), 
                 outlier.size = 0.05,
                 width = 3,
                 color = "grey60") +
    geom_abline(intercept = 0, slope = 1) + 
    geom_smooth(data = weights %>% 
                  filter(#preprocess == "none",
                    #bias_correct == FALSE,
                    model %in% c("GRC", "NNLS2", "CLS")) %>% 
                  mutate(preprocess = factor(preprocess, preprocess_levels),
                         model = factor(model, levels = c("CLS", "NNLS2", "GRC"))) %>%
                  group_by(model, bias_correct, preprocess, n_sim) %>% 
                  summarise(mean_nw = mean(nw)),
                aes(n_sim, mean_nw, 
                    # linetype = preprocess, 
                    # color = preprocess
                    ),
                method = "loess",
                size = 0.5,
                color = "grey60") +
    facet_wrap(~model, 
               ncol = 3, scales = "free") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    scale_color_manual(values = cols[c(3,6,9)]) +
    scale_fill_manual(values = cols[c(3,6,9)]) +
    scale_linetype_manual(values = c(1, 3, 5)) +
    labs(x = "Ensemble size",
         y = "Number of members used in \nthe optimal combination") 
p1

plotdata <- weights %>% 
    filter(#preprocess == "none",
           bias_correct == TRUE,
           model %in% c("GRC", "NNLS2", "CLS")) %>% 
    mutate(nw = nw / n_sim,
           preprocess = factor(preprocess, preprocess_levels)) %>%
    group_by(model, bias_correct, preprocess, n_sim) %>% 
    summarise(mean_nw = mean(nw))

p2 <- ggplot() +
    geom_line(data = filter(plotdata, model %in% c("GRC", "NNLS2", "CLS")),
              aes(n_sim, mean_nw, 
                  color = model, 
                  linetype = preprocess),
              size = 0.5) +
    # geom_line(data = filter(plotdata, model == "NNLS2"),
    #           aes(n_sim, mean_nw, 
    #               color = bias_correct, 
    #               linetype = preprocess)) +
    # geom_line(data = filter(plotdata, model == "CLS"),
    #           aes(n_sim, mean_nw, 
    #               color = bias_correct, 
    #               linetype = preprocess)) +
    # facet_wrap(~model) +
    theme_minimal() +
    # scale_color_manual(values = cols[c(9,4)]) +
    scale_color_manual(values = cols[c(5,10,8)]) +
    scale_linetype_manual(values = c(1, 3, 5)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Ensemble size",
         y = "Ensemble members with weight > 1 %")
    
# ggsave("used_weights_p1.pdf", p1)
# ggsave(p2, "used_weights_p2.pdf")
p <- p1 + p2 + plot_layout(guides = "collect", 
                           widths = c(3,1))
p
ggsave("used_weights.pdf", p1, height = 3, width = 8)

gc()

weights %>% 
    # mutate(nw = nw / n_sim) %>% 
    group_by(model, bias_correct, preprocess, n_sim) %>% 
    summarise(mean_nw = mean(nw)) %>% 
    pivot_wider(names_from= c(model, bias_correct, preprocess), 
                values_from = mean_nw)
