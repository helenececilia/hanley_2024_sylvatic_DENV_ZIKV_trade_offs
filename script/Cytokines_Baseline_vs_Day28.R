## ---------------------------
##
## Script name: Cytokines_Baseline_vs_Day28.R
##
## Purpose of script: Executes the function to test differences between cytokine concentrations
## at baseline and at day 28 of the experiments, in infected monkeys
## When differences are detected, plots the data for control monkeys to qualitatively assess
## if the change was similar, and therefore not due to infection
## Detailed results are produced in .txt files in output/result_files/cytokines/baseline_vs_day28,
## but you can find the summary of these results in Table S.5.2
##
## Author: Helene Cecilia
##
## Date Created: 2023-07-05

rm(list=ls())

## Loading Packages  ------------------
library(tidyverse)
library(kableExtra)
library(patchwork)

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
source("Cytokines_Functions.R")

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------


data <- read.csv("../data/Cytokines_All_Exp.csv")
df_denv_cyno <- data[data$virus == "Dengue virus" & data$NHP == "Cyno",]
df_denv_squirrel <- data[data$virus == "Dengue virus" & data$NHP == "Squirrel",]
df_zikv_squirrel <- data[data$virus == "Zika virus" & data$NHP == "Squirrel",]

# The function already deals with problems of LOD and such,
# so disable the warnings for final rendering
defaultW <- getOption("warn")
options(warn = -1) 
cyto_return_baseline_denv_cyno_noLOD <- compare_baseline_day28(df_denv_cyno, exclude_LOD = T,
                                                            verbose = T, plot = F,
                                                            filename = "../output/result_files/cytokines/baseline_vs_day28/baseline_day28_excludeLOD_DENV_cyno")

cyto_return_baseline_denv_cyno_LOD <- compare_baseline_day28(df_denv_cyno, exclude_LOD = F,
                                                            verbose = T, plot = F,
                                                            filename = "../output/result_files/cytokines/baseline_vs_day28/baseline_day28_withLOD_DENV_cyno")

cyto_return_baseline_denv_squirrel_noLOD <- compare_baseline_day28(df_denv_squirrel, exclude_LOD = T,
                                                            verbose = T, plot = F,
                                                            filename = "../output/result_files/cytokines/baseline_vs_day28/baseline_day28_excludeLOD_DENV_squirrel")

cyto_return_baseline_denv_squirrel_LOD <- compare_baseline_day28(df_denv_squirrel, exclude_LOD = F,
                                                            verbose = T, plot = F,
                                                            filename = "../output/result_files/cytokines/baseline_vs_day28/baseline_day28_withLOD_DENV_squirrel")

cyto_return_baseline_zikv_squirrel_noLOD <- compare_baseline_day28(df_zikv_squirrel, exclude_LOD = T,
                                                            verbose = T, plot = F,
                                                            filename = "../output/result_files/cytokines/baseline_vs_day28/baseline_day28_excludeLOD_ZIKV_squirrel")

cyto_return_baseline_zikv_squirrel_LOD <- compare_baseline_day28(df_zikv_squirrel, exclude_LOD = F,
                                                            verbose = T, plot = F,
                                                            filename = "../output/result_files/cytokines/baseline_vs_day28/baseline_day28_withLOD_ZIKV_squirrel")

options(warn = defaultW)

# Significant differences detected for some cytokines in infected cynomolgus macaques
# Qualitative comparisons with control cynomolgus macaques
# We concluded that differences were similar and therefore not due to infection

my_colors <- c("Control" = "darkgrey",
               "1 Mosquito" = "#c2e699",
               "10 Mosquitos" = "#006837")

data <- read.csv("../data/Cytokines_All_Exp.csv")
cyno <- data[data$virus == "Dengue virus" & data$NHP == "Cyno",]

i_tac <- cyno[cyno$day %notin% seq(0,27) & cyno$cytokine == "I.TAC",]

p1 <- ggplot(i_tac) +
  geom_line(aes(x = day, y = log10(value), group = ID, color = group)) +
  geom_point(aes(x = day, y = log10(value), color = group)) +
  geom_hline(aes(yintercept = log10(LOD)), color = "black") + ggtitle("A - I-TAC") +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(limits = c(-7,28),
                     breaks = c(-7,28),
                     labels = c("-7","28")) +
  labs(y = bquote("Concentration ("*log[10]~"pg/\u03BCl)"), x = "", color = "") +
  annotate(geom = "text", label = "LOD",
         color = "black", size = 6,
         x = -5.7, y = 1.54) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 25,
                                  margin = margin(r = 15)),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 23),
        plot.title = element_text(size = 24),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 17))

il12 <- cyno[cyno$day %notin% seq(0,27) & cyno$cytokine == "IL.12",]

p2 <- ggplot(il12) +
  geom_line(aes(x = day, y = log10(value), group = ID, color = group)) +
  geom_point(aes(x = day, y = log10(value), color = group)) +
  geom_hline(aes(yintercept = log10(LOD)), color = "black") + ggtitle("B - IL-12") +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(limits = c(-7,28),
                     breaks = c(-7,28),
                     labels = c("-7","28")) +
  labs(y = "", x = "Day post infection", color = "") +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25,
                                  margin = margin(t = 15)),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 23),
        plot.title = element_text(size = 24),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 17))


ilra <- cyno[cyno$day %notin% seq(0,27) & cyno$cytokine == "IL.RA",]

p3 <- ggplot(ilra) +
  geom_line(aes(x = day, y = log10(value), group = ID, color = group)) +
  geom_point(aes(x = day, y = log10(value), color = group)) +
  geom_hline(aes(yintercept = log10(LOD)), color = "black")  + ggtitle("C - IL-RA") +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(limits = c(-7,28),
                     breaks = c(-7,28),
                     labels = c("-7","28")) +
  labs(y = "", x = "", color = "") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25,
                                  margin = margin(t = 15)),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 23),
        plot.title = element_text(size = 24),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 0))


p <- (p1|p2|p3)
plot(p)

