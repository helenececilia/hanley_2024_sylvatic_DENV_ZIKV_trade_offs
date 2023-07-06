## ---------------------------
##
## Script name: Cytokines_GLMM_Infection_Virus_Host_Effects.R
##
## Purpose of script: Test through generalized linear mixed effect models 
## if cytokine concentrations differ by infection status, virus species, or host species
## Detailed results are produced in .txt files in output/result_files/cytokines/mixed_models/effect_infection_virus_monkey_species,
## but you can find the summary of these results in Table S.5.3 and S.5.4
##
## Author: Helene Cecilia
##
## Date Created: 2023-07-05

rm(list=ls())

## Loading Packages  ------------------
library(ggplot2)
library(tidyverse)
library(BAS)
library(patchwork)
library(cowplot)
library(gghalves)
library(car) # for leveneTest
library(MASS) # for glm.nb and multidimensional scaling (isoMDS, sammon)
library(VGAM) # for vglm
library(lme4) # for lmer
library(lmerTest)
library(modelr)
library(performance) # for check_collinearity
library(scales) # for trans_breaks
library(MuMIn) # for r.squaredGLMM
library(lattice) # for qqmath
library(sjPlot) # for plot_model
library(nlme) # for varIdent
library(report)
library(kableExtra)
library(bbmle)
library(clusrank) # wilcoxon test for clustered (non-independent) data
library(DHARMa)
library(effects)
library(glmmTMB)
library(gplots) # for textplot
library(mgcv) # for gam
library(ggtext) # for element_markdown

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
# source()

## Global command
`%notin%` <- Negate(`%in%`)


## -------------------------------

# Data ----
data <- read.csv("../data/Cytokines_All_Exp.csv")
baseline <- data[data$day == -7,c("ID","cytokine","value")]
colnames(baseline)[3] <- "baseline"
data <- merge(data, baseline, by = c("ID","cytokine"))

data$inf_status <- "Control"
data[data$group != "Control",]$inf_status <- "Infected"
data$inf_status <- as.factor(data$inf_status)
data <- within(data, inf_status <- relevel(inf_status, ref = "Infected"))
data$group <- interaction(data$NHP,data$virus)
data$group <- as.factor(data$group)
data <- within(data, group <- relevel(group, ref = "Squirrel.Dengue virus"))
data <- data[data$day != 28 & data$day !=-7,]

# Models excluding measures below LOD ----
pdf(file = "../output/result_files/cytokines/mixed_effect_models/effect_infection_virus_monkey_species/REF_MODELS_full_report_RE_ID_day_no_LOD_no_dispformula.pdf")
for(c in unique(data$cytokine)){
  df <- data[data$cytokine == c,]
  textplot(c, cex = 2)
  #########
  textplot("Infection in DENV-cyno", cex = 2)
  d_cyno <- df[df$value != df$LOD & df$group == "Cyno.Dengue virus",]
  nb_excl <- length(df[df$value == df$LOD & df$group == "Cyno.Dengue virus",]$value)
  nb_inf <- length(d_cyno[d_cyno$inf_status == "Infected",]$value)
  nb_cont <- length(d_cyno[d_cyno$inf_status == "Control",]$value)

  textplot(paste0("Nb excluded (LOD): ", nb_excl,"\n",
                  "Nb obs infection: ", nb_inf, "\n",
                  "Nb obs control: ", nb_cont), cex = 2)
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ inf_status +
                       (1|ID) + (1|day),
                     data = d_cyno)
    simulateResiduals(m_inf, plot = T)
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  textplot("Infection in DENV-squirrel", cex = 2)
  inf_d_sq <- df[df$group == "Squirrel.Dengue virus",]
  # add the controls from the ZIKV exp
  cont_z <- df[df$inf_status == "Control" & df$group == "Squirrel.Zika virus",]
  d_sq <- rbind(inf_d_sq, cont_z)
  nb_excl <- length(d_sq[d_sq$value == d_sq$LOD,]$value)
  d_sq <- d_sq[d_sq$value != d_sq$LOD,]
  nb_inf <- length(d_sq[d_sq$inf_status == "Infected",]$value)
  nb_cont <- length(d_sq[d_sq$inf_status == "Control",]$value)
  
  textplot(paste0("Nb excluded (LOD): ", nb_excl,"\n",
                  "Nb obs infection: ", nb_inf, "\n",
                  "Nb obs control: ", nb_cont), cex = 2)  
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ inf_status +
                       (1|ID) + (1|day),
                     data = d_sq)
    simulateResiduals(m_inf, plot = T)
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})

  textplot("Infection in ZIKV-squirrel", cex = 2)
  inf_z_sq <- df[df$group == "Squirrel.Zika virus",]
  # add the controls from the ZIKV exp
  cont_d_sq <- df[df$inf_status == "Control" & df$group == "Squirrel.Dengue virus",]
  z_sq <- rbind(inf_z_sq, cont_d_sq)
  nb_excl <- length(z_sq[z_sq$value == z_sq$LOD,]$value)
  z_sq <- z_sq[z_sq$value != z_sq$LOD,]
  nb_inf <- length(z_sq[z_sq$inf_status == "Infected",]$value)
  nb_cont <- length(z_sq[z_sq$inf_status == "Control",]$value)
  
  textplot(paste0("Nb excluded (LOD): ", nb_excl,"\n",
                  "Nb obs infection: ", nb_inf, "\n",
                  "Nb obs control: ", nb_cont), cex = 2)   
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ inf_status +
                       (1|ID) + (1|day),
                     data = z_sq)
    simulateResiduals(m_inf, plot = T)
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  #########
  textplot("Effect of experiment", cex = 2)
  my_df <- df[df$value != df$LOD & df$inf_status == "Infected",]
  nb_excl <- length(df[df$value == df$LOD & df$inf_status == "Infected",]$value)
  
  nb_d_cy <- length(my_df[my_df$group == "Cyno.Dengue virus",]$value)
  nb_d_sq <- length(my_df[my_df$group == "Squirrel.Dengue virus",]$value)
  nb_z_sq <- length(my_df[my_df$group == "Squirrel.Zika virus",]$value)
  
  textplot(paste0("Nb excluded (LOD): ", nb_excl,"\n",
                  "Nb obs DENV-squirrel (infected only): ", nb_d_sq, "\n",
                  "Nb obs DENV-cyno (infected only): ", nb_d_cy, "\n",                
                  "Nb obs ZIKV-squirrel (infected only): ", nb_z_sq), cex = 1.2) 
  tryCatch({
    m_group <- glmmTMB(log10(value) ~ group +
                         (1|ID) + (1|day),
                       data = my_df)
    simulateResiduals(m_group, plot = T)
    textplot(capture.output(summary(m_group)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  #########
}
dev.off()


# Including measures below LOD ----
pdf(file = "../output/result_files/cytokines/mixed_effect_models/effect_infection_virus_monkey_species/full_report_RE_ID_day_with_LOD_no_dispformula.pdf")
for(c in unique(data$cytokine)){
  df <- data[data$cytokine == c,]
  textplot(c, cex = 2)
  #########
  textplot("Infection in DENV-cyno", cex = 2)
  d_cyno <- df[df$group == "Cyno.Dengue virus",]
  nb_excl <- length(df[df$value == df$LOD & df$group == "Cyno.Dengue virus",]$value)
  # browser()
  textplot(paste0("Nb obs (total):", length(d_cyno$value),
                  "\nNb LOD (included): ", nb_excl), cex = 2)
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ inf_status +
                       (1|ID) + (1|day),
                     data = d_cyno)
    simulateResiduals(m_inf, plot = T)
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  textplot("Infection in DENV-squirrel", cex = 2)
  inf_d_sq <- df[df$group == "Squirrel.Dengue virus",]
  # add the controls from the ZIKV exp
  cont_z <- df[df$inf_status == "Control" & df$group == "Squirrel.Zika virus",]
  d_sq <- rbind(inf_d_sq, cont_z)
  nb_excl <- length(d_sq[d_sq$value == d_sq$LOD,]$value)
  textplot(paste0("Nb obs (total):", length(d_sq$value),
                  "\nNb LOD (included): ", nb_excl), cex = 2)
  
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ inf_status +
                       (1|ID) + (1|day),
                     data = d_sq)
    simulateResiduals(m_inf, plot = T)
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  textplot("Infection in ZIKV-squirrel", cex = 2)
  inf_z_sq <- df[df$group == "Squirrel.Zika virus",]
  # add the controls from the ZIKV exp
  cont_d_sq <- df[df$inf_status == "Control" & df$group == "Squirrel.Dengue virus",]
  z_sq <- rbind(inf_z_sq, cont_d_sq)
  nb_excl <- length(z_sq[z_sq$value == z_sq$LOD,]$value)
  textplot(paste0("Nb obs (total):", length(z_sq$value),
                  "\nNb LOD (included): ", nb_excl), cex = 2)
  
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ inf_status +
                       (1|ID) + (1|day),
                     data = z_sq)
    simulateResiduals(m_inf, plot = T)
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  #########
  textplot("Effect of experiment", cex = 2)
  my_df <- df[df$inf_status == "Infected",]
  nb_excl <- length(df[df$value == df$LOD & df$inf_status == "Infected",]$value)
  # browser()
  textplot(paste0("Nb obs (total):", length(my_df$value),
                  "\nNb LOD (included): ", nb_excl), cex = 2)
  tryCatch({
    m_group <- glmmTMB(log10(value) ~ group +
                         (1|ID) + (1|day),
                       data = my_df)
    simulateResiduals(m_group, plot = T)
    textplot(capture.output(summary(m_group)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  #########
}
dev.off()



# Plots & inspection/confirmation ----

# Confirmation of effect of infection in ZIKV squirrel ----
c <- "RANTES"
df <- data[data$cytokine == c,]
inf_z_sq <- df[df$group == "Squirrel.Zika virus",]
# add the controls from the other squirrel exp
cont_d_sq <- df[df$inf_status == "Control" & df$group == "Squirrel.Dengue virus",]
z_sq <- rbind(inf_z_sq, cont_d_sq)
lod <- unique(z_sq$LOD)
p_lod <- ggplot() +
  geom_hline(aes(yintercept = lod), color = "black") +
  geom_boxplot(data = z_sq, aes(x = inf_status, y = value, fill = inf_status),
               alpha = 0.5, outlier.shape = NA) +
  geom_jitter(data = z_sq,aes(x = inf_status, y = value, color = inf_status),
              size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Infected" = "#2b8cbe",
                                "Control" = "darkgrey")) +
  scale_fill_manual(values = c("Infected" = "#2b8cbe",
                               "Control" = "darkgrey")) +
  # annotate(geom = "text", label = "LOD", color = "grey", size = 7,
  #          x = 2.4, y = min(1.2*lod,lod+0.5), hjust=0,vjust=0) +
  annotate(geom = "text", label = "    B", color = "black", size = 7,
           x = 0.3, y = 10**1.61, hjust=0,vjust=0) + # for some reason the x position doesn't work properly, adjust with spaces
  annotate("text",
           x = 1:length(table(z_sq$inf_status)),
           y = 1.1*max(z_sq$value) - 5,
           label = table(z_sq$inf_status),
           vjust = -1,
           size = 7) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(.x))) +
  coord_cartesian(ylim = c(lod-1,1.1*max(z_sq$value))) +
  theme_classic() + labs(y = "including measures &le; LOD", x = "", color = "", fill = "") +
  # ggtitle(paste0(c," - effect of Zika infection in squirrels")) +
  theme(axis.text = element_text(size = 22),
        axis.title.y = element_markdown(size = 22,
                                    margin = margin(r = 15)),
        legend.position = "none",
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 22))

z_sq <- z_sq[z_sq$value != z_sq$LOD,]
p_no_lod <- ggplot() +
  geom_hline(aes(yintercept = lod), color = "black") +
  geom_boxplot(data = z_sq, aes(x = inf_status, y = value, fill = inf_status),
               alpha = 0.5, outlier.shape = NA) +
  geom_jitter(data = z_sq,aes(x = inf_status, y = value, color = inf_status),
              size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Infected" = "#2b8cbe",
                                "Control" = "darkgrey")) +
  scale_fill_manual(values = c("Infected" = "#2b8cbe",
                               "Control" = "darkgrey")) +
  annotate(geom = "text", label = "    A", color = "black", size = 7,
           x = 0.3, y = 10**1.61, hjust=0,vjust=0) + # for some reason the x position doesn't work properly, adjust with spaces
  annotate(geom = "text", label = "LOD", color = "black", size = 7,
           x = 2.4, y = 10**0.9, hjust=0,vjust=0) +
  annotate("text",
           x = 1:length(table(z_sq$inf_status)),
           y = 1.1*max(z_sq$value) - 5,
           label = table(z_sq$inf_status),
           vjust = -1,
           size = 7) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(.x))) +
  coord_cartesian(ylim = c(lod-1,1.1*max(z_sq$value))) +
  theme_classic() + labs(x = "", color = "", fill = "",
                         y = "Concentration (log<sub>10</sub> pg/\u03BCl)<br>excluding measures &le; LOD") +
  # ggtitle(paste0(c," - effect of Zika infection in squirrels")) +
  theme(axis.text = element_text(size = 22),
        axis.title.y = element_markdown(size = 22,
                                    margin = margin(r = 15)),
        legend.position = "none",
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 22))

p <- p_no_lod | p_lod

png(filename = paste0("../output/figures/suppl/Figure_S7.png"), width = 1600, height = 500)
print(p)
dev.off()


# Confirmation of effect of species ----
# Effect detected in infected monkeys, models ran for the same cytokines in control monkeys
# and comparing coeff estimates and confidence intervals
select_cyto <- c("I.TAC","MIP.1a","RANTES","Eotaxin")
# exclude measures below LOD
pdf(file = "../output/result_files/cytokines/mixed_effect_models/effect_infection_virus_monkey_species/report_species_effect_in_controls_no_LOD.pdf")
for(c in select_cyto){
  df <- data[data$cytokine == c,]
  textplot(c, cex = 2)
  my_df <- df[df$inf_status == "Control" & df$virus == "Dengue virus",]
  # add the controls from the ZIKV exp
  cont_z <- df[df$inf_status == "Control" & df$group == "Squirrel.Zika virus",]
  my_df <- rbind(my_df, cont_z)
  
  nb_excl <- length(my_df[my_df$value == my_df$LOD,]$value)
  
  my_df$NHP <- as.factor(my_df$NHP)
  my_df <- my_df[my_df$value != my_df$LOD,]
  nb_cy <- length(my_df[my_df$NHP == "Cyno",]$value)
  nb_sq <- length(my_df[my_df$NHP == "Squirrel",]$value)

  textplot(paste0("Nb excluded (LOD): ", nb_excl,"\n",
                  "Nb obs control squirrel: ", nb_sq, "\n",
                  "Nb obs control cyno: ", nb_cy), cex = 1.2)
  tryCatch({
    my_df <- within(my_df, NHP <- relevel(NHP, ref = "Squirrel"))
    m_group <- glmmTMB(log10(value) ~ NHP +
                         (1|ID) + (1|day),
                       data = my_df)
    simulateResiduals(m_group, plot = T)
    textplot(capture.output(summary(m_group)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
}
dev.off()

# include measures below LOD
pdf(file = "../output/result_files/cytokines/mixed_effect_models/effect_infection_virus_monkey_species/report_species_effect_in_controls_with_LOD.pdf")
for(c in select_cyto){
  df <- data[data$cytokine == c,]
  textplot(c, cex = 2)
  my_df <- df[df$inf_status == "Control" & df$virus == "Dengue virus",]
  # add the controls from the ZIKV exp
  cont_z <- df[df$inf_status == "Control" & df$group == "Squirrel.Zika virus",]
  my_df <- rbind(my_df, cont_z)
  
  my_df$NHP <- as.factor(my_df$NHP)
  # my_df <- my_df[my_df$value != my_df$LOD,]
  nb_cy <- length(my_df[my_df$NHP == "Cyno",]$value)
  nb_sq <- length(my_df[my_df$NHP == "Squirrel",]$value)
  
  textplot(paste0("Nb obs control squirrel: ", nb_sq, "\n",
                  "Nb obs control cyno: ", nb_cy), cex = 1.2)
  tryCatch({
    my_df <- within(my_df, NHP <- relevel(NHP, ref = "Squirrel"))
    m_group <- glmmTMB(log10(value) ~ NHP +
                         (1|ID) + (1|day),
                       data = my_df)
    simulateResiduals(m_group, plot = T)
    textplot(capture.output(summary(m_group)), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
}
dev.off()

# Figure
# MIP.1a was removed because too few observations and inconsistent results between models
select_cyto <- c("I.TAC","RANTES","Eotaxin")

tags_no_lod <- c("Eotaxin" = "A",
                 "I.TAC" = "B",
                 "RANTES" = "C")

tags_lod <- c("Eotaxin" = "D",
              "I.TAC" = "E",
              "RANTES" = "F")
# it has to be a list for the key to give both positions at once (2 dimensions)
tags_pos <- list("Eotaxin" = c(0.17,1),
              "I.TAC" = c(0.09,1),
              "RANTES" = c(0.07,1))

fig_title <- c("Eotaxin" = "Figure_S5_A_D",
               "I.TAC" = "Figure_S5_B_E",
               "RANTES" = "Figure_S5_C_F")

for(c in select_cyto){
  df <- data[data$cytokine == c,]
  my_df <- df[df$inf_status == "Control" & df$virus == "Dengue virus",]
  # add the controls from the ZIKV exp
  cont_z <- df[df$inf_status == "Control" & df$group == "Squirrel.Zika virus",]
  my_df <- rbind(my_df, cont_z)
  
  my_df$NHP <- as.factor(my_df$NHP)
  inf_df <- df[df$inf_status == "Infected" & df$virus == "Dengue virus",]
  plot_df <- rbind(my_df,inf_df)
  plot_df$group <- interaction(plot_df$NHP,plot_df$inf_status)
  lod <- min(plot_df$LOD) # in that case several LODs possible
  
  global_title <- c
  lgd_pos <- "none"
  if(c != "Eotaxin"){
    y_title <- ""
  }else{
    y_title <- "Concentration (log<sub>10</sub> pg/\u03BCl)<br>including measures &le; LOD"
  }
  if(c == "I.TAC"){
    x_title <- "Monkey species"
  }else{
    x_title <- ""
  }
  
  p_lod <- ggplot() +
    geom_hline(aes(yintercept = lod), color = "black") +
    geom_boxplot(data = plot_df, aes(x = group, y = value, fill = inf_status,
                                     group = group),
                 outlier.shape = NA, alpha = 0.5) +
    geom_jitter(data = plot_df,aes(x = group, y = value, color = inf_status)) +
    scale_color_manual(values = c("Infected" = "#006837",
                                  "Control" = "darkgrey")) +
    scale_fill_manual(values = c("Infected" = "#006837",
                                 "Control" = "darkgrey")) +
    annotate("text",
             x = 1:length(table(plot_df$group)),
             y = min(max(plot_df$value)-50,10**3.98),
             label = table(plot_df$group),
             vjust = -1,
             size = 8) +
    scale_y_continuous(trans = "log10",
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(.x))) +
    scale_x_discrete(labels = c("Cyno","Squirrel","Cyno","Squirrel")) +
    coord_cartesian(ylim = c(lod-1,1.15*max(plot_df$value))) +
    theme_classic() + labs(y = y_title,
                           x = x_title, color = "", fill = "",
                           tag = tags_lod[[c]]) +
    theme(axis.text = element_text(size = 25),
          axis.title.y = element_markdown(size = 25,
                                          margin = margin(r=15)),
          axis.title.x = element_text(size = 25,
                                      margin = margin(t = 15)),
          legend.position = lgd_pos,
          legend.text = element_text(size = 25),
          plot.title = element_text(size = 25))
  
  plot_df <- plot_df[plot_df$value != plot_df$LOD,]
  
  x_title <- ""
  if(c != "Eotaxin"){
    my_table <- table(plot_df$group)
    my_x_labels <- c("Cyno","Squirrel","Cyno","Squirrel")
    lod_label <- ""
    y_title <- ""
    lgd_pos <- "none"
  }else{
    my_table <- table(plot_df$group)[which(table(plot_df$group)!=0)]
    my_x_labels <- c("Cyno","Squirrel","Cyno")
    y_title <- "Concentration (log<sub>10</sub> pg/\u03BCl)<br>excluding measures &le; LOD"
    lod_label <- "LOD"
    lgd_pos <- c(0.87,0.33)
  }
  
  p_no_lod <- ggplot() +
    geom_hline(aes(yintercept = lod), color = "black") +
    geom_boxplot(data = plot_df, aes(x = group, y = value, fill = inf_status,
                                     group = group),
                 outlier.shape = NA, alpha = 0.5) +
    geom_jitter(data = plot_df,aes(x = group, y = value, color = inf_status)) +
    scale_color_manual(values = c("Infected" = "#006837",
                                  "Control" = "darkgrey")) +
    scale_fill_manual(values = c("Infected" = "#006837",
                                 "Control" = "darkgrey")) +
    annotate(geom = "text", label = lod_label, color = "black", size = 8,
             x = 3.2, y = 10**0.24, hjust=0,vjust=0) +
    annotate("text",
             x = 1:length(my_table),
             y = min(max(plot_df$value)-50,10**3.98),
             label = my_table,
             vjust = -1,
             size = 8) +
    scale_y_continuous(trans = "log10",
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(.x))) +
    scale_x_discrete(labels = my_x_labels) +
    coord_cartesian(ylim = c(lod-1,1.15*max(plot_df$value))) +
    theme_classic() + labs(y = y_title,
                           x = x_title, color = "", fill = "",
                           tag = tags_no_lod[[c]]) +
    # ggtitle(paste0(c," - DENV")) +
    theme(axis.text = element_text(size = 25),
          axis.title.y = element_markdown(size = 25,
                                          margin = margin(r = 15)),
          axis.title.x = element_text(size = 25,
                                          margin = margin(t = 15)),
          legend.position = lgd_pos,
          legend.text = element_text(size = 25),
          plot.title = element_text(size = 25))
  
  p <- p_no_lod / p_lod
  p <- p + plot_annotation(title = global_title) & theme(plot.title = element_text(size = 25),
                                                         plot.tag = element_text(size = 24),
                                                         plot.tag.position = tags_pos[[c]])
  
  png(filename = paste0("../output/figures/suppl/",fig_title[[c]],".png"), width = 800, height = 1000)
  print(p)
  dev.off()
}


