## ---------------------------
##
## Script name: Cytokines_Sex_Diff_at_Baseline.R
##
## Purpose of script: Executes the function to determine if differences
## exist between sexes for cytokine concentrations before the start of the experiments
## Detailed results are produced in .txt files in output/result_files/cytokines/sex_diff_at_baseline,
## but you can find the summary of these results in Table S.5.1
##
## Author: Helene Cecilia
##
## Date Created: 2023-05-31

rm(list=ls())

## Loading Packages  ------------------
library(tidyverse)
library(kableExtra)
library(patchwork)
library(car) # for leveneTest

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
source("Cytokines_Functions.R")

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

data <- read.csv("../data/Cytokines_All_Exp.csv")
df_cyno <- data[data$NHP == "Cyno",]
df_squirrel <- data[data$NHP == "Squirrel",]


defaultW <- getOption("warn")
options(warn = -1) 
# excluding measures below LOD
cyto_sex_diff_cyno <- sex_diff_at_baseline(df_cyno, exclude_LOD = T,
                                           plot = T, verbose = T,
                                           filename = "../output/result_files/cytokines/sex_diff_at_baseline/sex_diff_cyno_excludeLOD")
cyto_sex_diff_squirrel <- sex_diff_at_baseline(df_squirrel, exclude_LOD = T,
                                               plot = T, verbose = T,
                                               filename = "../output/result_files/cytokines/sex_diff_at_baseline/sex_diff_squirrel_excludeLOD")
# including measures below LOD
cyto_sex_diff_cyno <- sex_diff_at_baseline(df_cyno, exclude_LOD = F,
                                           plot = T, verbose = T,
                                           filename = "../output/result_files/cytokines/sex_diff_at_baseline/sex_diff_cyno_withLOD")
cyto_sex_diff_squirrel <- sex_diff_at_baseline(df_squirrel, exclude_LOD = F,
                                               plot = T, verbose = T,
                                               filename = "../output/result_files/cytokines/sex_diff_at_baseline/sex_diff_squirrel_withLOD")
options(warn = defaultW)


