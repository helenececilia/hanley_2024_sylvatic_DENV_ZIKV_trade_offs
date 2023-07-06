## ---------------------------
##
## Script name: Cytokines_Data_Formatting.R
##
## Purpose of script: Takes the supplementary Tables S1,S2,S3 and turns them into
## tidy dataframes for cytokines concentrations, used for further testing.
##
## Author: Helene Cecilia
##
## Date Created: 2022-05-05

rm(list=ls())

## Loading Packages  ------------------
library(pzfx)
library(readxl)
library(ggplot2)
library(pracma) # for trapz
library(tidyverse)
library(BAS) # for Bayesian regression
library(stringdist) # for amatch
library(ggh4x) # for facet_nested

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location 
getwd()

## Load source files ------------------
# source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

# DENV Sylv Cyno -----
df_complete <- read.csv(file = "../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                        sep = "\t", dec = ".")
df_cyto <- df_complete[,c(1,2,4,5,seq(59,87))] 

# transform columns into rows, to have one observation per row (tidy rule)
tidy_cyto <- df_cyto %>% pivot_longer(-c(`ID`,`Sex`,`Final.Treatment`,`Day.Post.Infection`),
                                      names_to = "description",
                                      values_to = "level") # we specify which column not to pivot (shorter because many to pivot)

df_cyto_LOD <- read.csv("../data/Cytokines_LOD_Sylvatic_DENV-2_Cynomolgus_Macaques.csv")
colnames(df_cyto_LOD) <- c("description","LOD")
tidy_cyto <- merge(x = tidy_cyto, y = df_cyto_LOD, by = "description", all = T) 
tidy_cyto <- tidy_cyto[complete.cases(tidy_cyto$level),]

# Assign values below LOD to LOD

tidy_cyto$level[tidy_cyto$level <= tidy_cyto$LOD] <- tidy_cyto$LOD[tidy_cyto$level <= tidy_cyto$LOD]
test <- tidy_cyto
colnames(test) <- c("cytokine","ID","sex","group","day","value","LOD")
test$NHP <- "Cyno"
test$virus <- "Dengue virus"
write.csv(test, file = "../data/Cytokines_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
          row.names = F)

# Sylv DENV Squirrel -----
df_complete <- read.csv(file = "../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                        sep = "\t", dec = ".")
df_cyto <- df_complete[,c(1,3,5,6,seq(62,90))] # keep bites (not titer results)
# transform columns into rows, to have one observation per row (tidy rule)
tidy_cyto <- df_cyto %>% pivot_longer(-c(`ID`,`Sex`,`Final.Treatment`,`Day.Post.Infection`),
                                      names_to = "description",
                                      values_to = "level") # we specify which column not to pivot (shorter because many to pivot)

df_cyto_LOD <- read.csv("../data/Cytokines_LOD_Sylvatic_DENV-2_Squirrel_Monkeys.csv")
colnames(df_cyto_LOD) <- c("description","LOD")

tidy_cyto <- merge(x = tidy_cyto, y = df_cyto_LOD, by = "description", all = T) 
tidy_cyto <- tidy_cyto[complete.cases(tidy_cyto$level),]

# Assign values below LOD to LOD

tidy_cyto$level[tidy_cyto$level <= tidy_cyto$LOD] <- tidy_cyto$LOD[tidy_cyto$level <= tidy_cyto$LOD]
test <- tidy_cyto
colnames(test) <- c("cytokine","ID","sex","group","day","value","LOD")
test$NHP <- "Squirrel"
test$virus <- "Dengue virus"
write.csv(test, file = "../data/Cytokines_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
          row.names = F)

# Sylv ZIKV Squirrel ----
df_complete <- read.csv(file = "../data/Table_S3_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                        sep = "\t", dec = ".")
df_cyto <- df_complete[,c(1,3,5,6,seq(62,90))]
# transform columns into rows, to have one observation per row (tidy rule)
tidy_cyto <- df_cyto %>% pivot_longer(-c(`ID`,`Sex`,`Final.Treatment`,`Day.Post.Infection`),
                                      names_to = "description",
                                      values_to = "level") # we specify which column not to pivot (shorter because many to pivot)

df_cyto_LOD <- read.csv("../data/Cytokines_LOD_Sylvatic_ZIKV_Squirrel_Monkeys.csv")
colnames(df_cyto_LOD) <- c("description","LOD")

tidy_cyto <- merge(x = tidy_cyto, y = df_cyto_LOD, by = "description", all = T) 
tidy_cyto <- tidy_cyto[complete.cases(tidy_cyto$level),]

# Assign values below LOD to LOD

tidy_cyto$level[tidy_cyto$level <= tidy_cyto$LOD] <- tidy_cyto$LOD[tidy_cyto$level <= tidy_cyto$LOD]
test <- tidy_cyto
colnames(test) <- c("cytokine","ID","sex","group","day","value","LOD")
test$NHP <- "Squirrel"
test$virus <- "Zika virus"
write.csv(test, file = "../data/Cytokines_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
          row.names = F)

df1 <- read.csv("../data/Cytokines_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                sep = ",", dec = ".")
df2 <- read.csv("../data/Cytokines_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                sep = ",", dec = ".")
df3 <- read.csv("../data/Cytokines_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                sep = ",", dec = ".")

df <- rbind(df1,df2,df3)
write.csv(df, file = "../data/Cytokines_All_Exp.csv",
          row.names = F)

