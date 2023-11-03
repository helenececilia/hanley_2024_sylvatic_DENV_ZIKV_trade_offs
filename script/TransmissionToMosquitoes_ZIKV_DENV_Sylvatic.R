## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Helene Cecilia
##
## Date Created: 2022-11-21

rm(list=ls())

## Loading Packages  ------------------
library(tidyverse)
library(bbmle) # for mle2
library(emdbook) # for dbetabinom
library(ggplot2)
library(mgcv) # for rmvn, gam
library(scales) # for alpha function
library(funrar) # for matrix_to_stack
library(readxl)
library(patchwork)
library(tmvtnorm) # for truncated multi variate normal rtmvnorm
library(DHARMa)
library(effects)
library(nlstools)
library(glmmTMB)
library(janitor) # for clean_names

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# DENGUE ----
## DENV squirrel data -----
df_tmp <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                   dec = ".", sep = "\t")
df_tmp <- df_tmp %>% clean_names()

df_tmp <- df_tmp[,c("id","sex","final_treatment","day_post_infection","viremia_deduced",
                    "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                    "x_mosquitoes_infected_body_legs_or_saliva")]

df_tmp$k <- round((df_tmp$x_mosquitoes_infected_body_legs_or_saliva/100)*df_tmp$number_of_mosquitoes_that_engorged_and_survived_to_titer)

select_inf <- df_tmp[!is.na(df_tmp$viremia_deduced) & !is.na(df_tmp$k) & df_tmp$final_treatment != "Control",]

# log transform viremia
select_inf$Viremia_log10_corr <- log10(select_inf$viremia_deduced+1)

denv_sq_inf <- select_inf[,c("id","Viremia_log10_corr",
                             "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                             "k")]

colnames(denv_sq_inf) <- c("ID","log_V","N","k")
denv_sq_inf$paper <- "denv"
denv_sq_inf$disease_class <- "unknown"
denv_sq_inf$serotype <- NA


## DENV cyno data -----
df_tmp <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                   dec = ".", sep = "\t")

df_tmp <- df_tmp %>% clean_names()

df_tmp <- df_tmp[,c("id","sex","final_treatment","day_post_infection","viremia_deduced",
                    "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                    "x_mosquitoes_infected_body_or_legs")]

df_tmp$k <- round((df_tmp$x_mosquitoes_infected_body_or_legs/100)*df_tmp$number_of_mosquitoes_that_engorged_and_survived_to_titer)

select_inf <- df_tmp[!is.na(df_tmp$viremia_deduced) & !is.na(df_tmp$k) & df_tmp$final_treatment != "Control",]

# viremia already in log, but back-transform to be log(vir+1) as in squirrels
select_inf$vir_back <- 10**select_inf$viremia_deduced
select_inf$vir_back[select_inf$vir_back == 1] <- 0
select_inf$Viremia_log10_corr <- log10(select_inf$vir_back+1)

denv_cy_inf <- select_inf[,c("id","Viremia_log10_corr",
                             "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                             "k")]

colnames(denv_cy_inf) <- c("ID","log_V","N","k")
denv_cy_inf$paper <- "denv"
denv_cy_inf$disease_class <- "unknown"
denv_cy_inf$serotype <- NA


## Merging datasets ----
denv_sq_inf$NHP <- "Squirrel"
denv_cy_inf$NHP <- "Cyno"

df_den <- rbind(denv_cy_inf,
                denv_sq_inf)

## General additive model ----
# mod_den <- gam(cbind(k,N-k) ~ s(log_V, k = 6),
#                data = df_den,
#                family = binomial)
# # the model was saved as an rds object if you want to load it directly (see below)
# # do not run this unless you want to overwrite the object
# # write_rds(mod_den, file = "../output/result_files/transmission_to_mosquitoes/GAM_zika_dengue_sylvatic/GAM_dengue.rds")
# 
# # quality check of the model
# plot(mod_den, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
#      col = c("black", rep("red", length(mod_den$residuals))))
# gam.check(mod_den)

# Load the model object instead of fitting it again if you want to make sure you're using the same model as the one published
mod_den <- readRDS("../output/result_files/transmission_to_mosquitoes/GAM_zika_dengue_sylvatic/GAM_dengue.rds")

summary(mod_den)
# Family: binomial
# Link function: logit
# 
# Formula:
#   cbind(k, N - k) ~ s(log_V, k = 6)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -17.43      37.04  -0.471    0.638
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value
# s(log_V) 1.841  1.995  4.173   0.118
# 
# R-sq.(adj) =  0.409   Deviance explained = 67.4%
# UBRE = -0.47875  Scale est. = 1         n = 157

df_den$prob <- df_den$k/df_den$N
# for axis limits and point size breaks
# max(df_den$log_V) # 3.9
# max(df_den$N) # 14
pred_nl_den <- predict(mod_den, type = "response",
                       newdata = data.frame(log_V = seq(0,3.9,length.out = 100)))
pred_den <- data.frame(log_V = seq(0,3.9,length.out = 100),
                       pred_proba = pred_nl_den)

breaks <- c(5, 10, 15)
n_breaks <- length(breaks)
labels <- c(breaks, rep("", n_breaks))
shapes <- c(rep(16, n_breaks), rep(6, n_breaks))
breaks2 <- rep(breaks, 2)

df_den$NHP <- factor(df_den$NHP, levels = c("Squirrel","Cyno")) # change order of legend

p_den <- ggplot() + geom_vline(xintercept = log10(21), color = "darkgrey",
                               linewidth = 1.3) +
  geom_point(data = df_den,
             aes(x = log_V, y = prob, pch = NHP,
                 size = N),
             alpha = 0.5, stroke = 1.5) + 
  geom_line(data = pred_den, aes(x = log_V, y = pred_proba),
            color = "#006837", linewidth = 2, alpha = 0.8) +
  annotate(geom = "text", label = "LOD",
           color = "darkgrey", size = 8,
           x = 1.57, y = 0.96) +
  scale_shape_manual(name = bquote(bold("Monkey species")),
                     values = c("Squirrel" = 16,
                                "Cyno" = 6),
                     labels = c("Squirrel monkeys",
                                "Cynomolgus macaques")) +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = 4))) +
  coord_cartesian(ylim = c(0,1), xlim = c(0,6.35)) +
  scale_x_continuous(breaks = seq(0,6.35),
                     expand = expansion(add = 0.06)) +
  scale_y_continuous(expand = expansion(add = c(0.02,0.03))) +
  scale_size_continuous(limits=c(1,15), range = c(0,10),
                         breaks = breaks2, labels = labels,
                         guide = guide_legend(order = 2, ncol = 2, byrow = FALSE,
                                              override.aes = list(shape = shapes),
                                              direction = "vertical", label.hjust = 1)) + # , label.vjust = -.5 trans = "log10
  theme_classic() +
  labs(x = bquote("Dengue virus titer ("*log[10]~"PFU/ml)"),
       y = "Prob mosquito infection",
       size = bquote(bold("N mosq. titered"))) + 
  theme(axis.text = element_text(size = 25),
        axis.title.y = element_text(size = 26,
                                    margin = margin(r=20)),
        axis.title.x = element_text(size = 26,
                                    margin = margin(t=15)),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 26),
        legend.position = c(0.8,0.75),
        legend.background = element_rect(color = NA, fill = NA))

# plot(p_den)


# ZIKA ----
## ZIKV squirrel data -----
df_tmp <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                   dec = ".", sep = "\t")

df_tmp <- df_tmp %>% clean_names()

df_tmp <- df_tmp[,c("id","sex","final_treatment","day_post_infection","viremia_deduced",
                    "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                    "x_mosquitoes_infected_body_legs_or_saliva")]

df_tmp$k <- round((df_tmp$x_mosquitoes_infected_body_legs_or_saliva/100)*df_tmp$number_of_mosquitoes_that_engorged_and_survived_to_titer)

select_inf <- df_tmp[!is.na(df_tmp$viremia_deduced) & !is.na(df_tmp$k) & df_tmp$final_treatment != "Control",]

# log transform viremia
select_inf$Viremia_log10_corr <- log10(select_inf$viremia_deduced+1)

zikv_sq_inf <- select_inf[,c("id","Viremia_log10_corr",
                             "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                             "k")]

colnames(zikv_sq_inf) <- c("ID","log_V","N","k")
zikv_sq_inf$paper <- "zikv"
zikv_sq_inf$disease_class <- "unknown"
zikv_sq_inf$serotype <- NA

## ZIKV cyno data ----
df_tmp <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Cynomolgus_Macaques.csv",
                   dec = ".", sep = "\t")

df_tmp <- df_tmp %>% clean_names()

df_tmp <- df_tmp[,c("id","sex","final_treatment","day_post_infection","viremia_deduced",
                    "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                    "x_mosquitoes_infected_body_legs_or_saliva")]

df_tmp$k <- round((df_tmp$x_mosquitoes_infected_body_legs_or_saliva/100)*df_tmp$number_of_mosquitoes_that_engorged_and_survived_to_titer)

select_inf <- df_tmp[!is.na(df_tmp$viremia_deduced) & !is.na(df_tmp$k) & df_tmp$final_treatment != "Control",]

# viremia already in log(vir+1) / no need to do anything

zikv_cy_inf <- select_inf[,c("id","viremia_deduced",
                             "number_of_mosquitoes_that_engorged_and_survived_to_titer",
                             "k")]

colnames(zikv_cy_inf) <- c("ID","log_V","N","k")
zikv_cy_inf$paper <- "zikv"
zikv_cy_inf$disease_class <- "unknown"
zikv_cy_inf$serotype <- NA

## Merging datasets ----
zikv_sq_inf$NHP <- "Squirrel"
zikv_cy_inf$NHP <- "Cyno"

df_zik <- rbind(zikv_cy_inf,
                zikv_sq_inf)

# General additive model ----

# mod_zik <- gam(cbind(k,N-k) ~ s(log_V, k = 6),
#                data = df_zik,
#                family = binomial)
# # the model was saved as an rds object if you want to load it directly (see below)
# # do not run this unless you want to overwrite the object
# write_rds(mod_zik, file = "../output/result_files/transmission_to_mosquitoes/GAM_zika_dengue_sylvatic/GAM_zika.rds")
# 
# # quality check of the model
# plot(mod_zik, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
#      col = c("black", rep("red", length(mod_zik$residuals))))
# gam.check(mod_zik)

# Load the model object instead of fitting it again if you want to make sure you're using the same model as the one published
mod_zik <- readRDS("../output/result_files/transmission_to_mosquitoes/GAM_zika_dengue_sylvatic/GAM_zika.rds")

summary(mod_zik)
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   cbind(k, N - k) ~ s(log_V, k = 6)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.0964     0.8098  -3.824 0.000131 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(log_V) 4.623  4.888   90.9  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.857   Deviance explained = 85.2%
# UBRE = 0.28301  Scale est. = 1         n = 74

df_zik$prob <- df_zik$k/df_zik$N
# for axis limits and point size breaks
# max(df_zik$log_V) # 6.3
# max(df_zik$N) # 14
pred_nl_zik <- predict(mod_zik, type = "response",
                       newdata = data.frame(log_V = seq(0,6.3,length.out = 75)))
pred_zik <- data.frame(log_V = seq(0,6.3,length.out = 75),
                       pred_proba = pred_nl_zik)

df_zik$NHP <- factor(df_zik$NHP, levels = c("Squirrel","Cyno")) # change order of legend

p_zik <- ggplot() +   geom_vline(xintercept = log10(21), color = "darkgrey",
                                 linewidth = 1.3) +
  geom_point(data = df_zik,
             aes(x = log_V, y = prob, pch = NHP,
                 size = N),
             alpha = 0.5, stroke = 1.5) + 
  geom_line(data = pred_zik, aes(x = log_V, y = pred_proba),
            color = "#2b8cbe", linewidth = 2, alpha = 0.8) +
  annotate(geom = "text", label = "LOD",
           color = "darkgrey", size = 8,
           x = 1.57, y = 0.96) +
  scale_shape_manual(name = bquote(bold("Monkey species")),
                     values = c("Squirrel" = 16,
                                "Cyno" = 6),
                     labels = c("Squirrel monkeys",
                                "Cynomolgus macaques")) +
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  coord_cartesian(ylim = c(0,1), xlim = c(0,6.35)) +
  scale_x_continuous(breaks = seq(0,6.35),
                     expand = expansion(add = 0.06)) +
  scale_y_continuous(expand = expansion(add = c(0.02,0.03))) +
  scale_size(limits=c(1,15),breaks=c(5,10,15),
             range = c(0,10)) +
  theme_classic() +
  labs(x = bquote("Zika virus titer ("*log[10]~"PFU/ml)"),
       y = "Prob mosquito infection",
       size = bquote(bold("N mosq. titered"))) + 
  theme(axis.text = element_text(size = 25),
        axis.title.y = element_text(size = 26,
                                    margin = margin(r = 20,)),
        axis.title.x = element_text(size = 26,
                                    margin = margin(t = 15)),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 26),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = NA))

# plot(p_zik)

# Assemble DENV and ZIKV plots -----
p <- (p_den / p_zik)
p <- p + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 23),
                                                   plot.tag.position = c(0.14,0.98))
png(filename = "../output/figures/main/Figure_5.png",
    width = 900, height = 1200)
print(p)
dev.off()

