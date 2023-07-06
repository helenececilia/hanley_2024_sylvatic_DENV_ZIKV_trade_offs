## ---------------------------
##
## Script name: Cytokines_NK_Cells.R
##
## Purpose of script: Determine if some cytokines are associated with NK cell activation
##
## Author: Helene Cecilia
##
## Date Created: 2023-04-21

rm(list=ls())

## Loading Packages  ------------------
library(tidyverse)
library(janitor) # for clean_names
library(DHARMa)
library(effects)
library(glmmTMB)
library(mgcv)
library(patchwork)
library(gammit) # for predict_gamm
library(ggtext) # for element_markdown

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# Data ----

df_control_cyno_complete <- read.csv("../data/Cytokines_NK_Cells_Control_Cynomolgus_Macaques.csv")
df_denv_cyno_complete <- read.csv("../data/Cytokines_NK_Cells_Sylvatic_DENV-2_Cynomolgus_Macaques.csv")

# The initial goal was to test associations of IL-12, IL-15 and TNF-alpha with NK cells,
# in both cynomolgus macaques and squirrel monkeys.
# But the number of measures below the limit of detection
# prevented most model fitting (see Supplemental text S.4.2),
# so we focused on IL-12 in cynomolgus macaques

# Model selection : control monkeys ----
df <- df_control_cyno_complete
df <- df[df$cytokine == "IL.12" & df$day %in% seq(1,7),]
df <- df[complete.cases(df$value, df$perc_nk_cells),]
df$log_cy <- log10(df$value)
df$ID <- as.factor(df$ID)
df$day <- as.factor(df$day)

ggplot(df) + geom_point(aes(x = log_cy, y = perc_nk_cells,
                            color = ID))
# no clear relationship within individuals, 
# but across individuals the negative association seems strong

################### LMs #################
# RE = random effect
## 1 : lm without RE ----
m1 <- glmmTMB(perc_nk_cells ~ log_cy,
              REML = FALSE,
              data = df)
simulateResiduals(m1, plot = T)
testDispersion(m1)
MuMIn::AICc(m1) # 119.8581
plot(allEffects(m1, partial.residuals = T)) 

## 2 : lm with RE intercept ID ----
m2 <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID),
              REML = FALSE,
              data = df)
simulateResiduals(m2, plot = T)
testDispersion(m2)
MuMIn::AICc(m2) # 114.8253
plot(allEffects(m2, partial.residuals = T))

anova(m1,m2) # signif pref for m2

## 3 : lm with RE intercept day ----
m3 <- glmmTMB(perc_nk_cells ~ log_cy + (1|day),
              REML = FALSE,
              data = df)
simulateResiduals(m3, plot = T)
testDispersion(m3)
MuMIn::AICc(m3) # 122.5972
plot(allEffects(m3, partial.residuals = T)) 

anova(m1,m3) # no pref for m3 over m1

## 4 : lm with RE intercept ID & day ----
m4 <- glmmTMB(perc_nk_cells ~ log_cy + + (1|ID) + (1|day),
              REML = FALSE,
              data = df)
simulateResiduals(m4, plot = T) 
testDispersion(m4) 
MuMIn::AICc(m4) # 117.8134
plot(allEffects(m4, partial.residuals = T)) 

anova(m2,m4) # no pref for m4 over m2
anova(m3,m4) # pref for m4 over m3
# m2 and m3 are not nested

## 5 : lm with RE intercept and slope uncorrelated (ID) to be comparable with the gam later ----
m5b <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID) + (0 + log_cy|ID),
               REML = FALSE,
               data = df)
simulateResiduals(m5b, plot = T) 
testDispersion(m5b) 
MuMIn::AICc(m5b) # 117.1914
plot(allEffects(m5b, partial.residuals = T)) 

anova(m2,m5b) # not in favour of m5b

## 6 : lm with RE intercept and slope uncorrelated (day) to be comparable with the gam later ----
m6b <- glmmTMB(perc_nk_cells ~ log_cy + (1|day) + (0 + log_cy|day),
               REML = FALSE,
               data = df)
simulateResiduals(m6b, plot = T) 
testDispersion(m6b) 
MuMIn::AICc(m6b) # 125.5853
plot(allEffects(m6b, partial.residuals = T)) 

anova(m3,m6b) # not in favour of m6b

## 7 : lm with RE intercept and slope uncorrelated (ID & day) ----
m7b <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID) + (0 + log_cy|ID) +
                 (1|day) + (0 + log_cy|day),
               REML = FALSE,
               data = df)
simulateResiduals(m7b, plot = T) 
testDispersion(m7b) 
MuMIn::AICc(m7b) # 124.0642
plot(allEffects(m7b, partial.residuals = T)) 

anova(m6b,m7b) # slightly in favour of m7b
anova(m5b,m7b) # not in favour of m7b

## 8 : lm with RE intercept (ID and day) and slope ID, uncorrelated ----
m8b <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID) + (0 + log_cy|ID) +
                 (1|day),
               REML = FALSE,
               data = df)
simulateResiduals(m8b, plot = T) 
testDispersion(m8b) 
MuMIn::AICc(m8b) # 120.4642
plot(allEffects(m8b, partial.residuals = T)) 

anova(m7b,m8b) # in favour of m8b
anova(m5b,m8b) # in favour of m5b
anova(m4,m8b) # in favour of m4
anova(m3,m8b) # slightly in favour of m8b

## 9 : lm with RE intercept (ID and day) and slope day, uncorrelated ----
m9b <- glmmTMB(perc_nk_cells ~ log_cy + (1|day) + (0 + log_cy|day) +
                 (1|ID),
               REML = FALSE,
               data = df)
simulateResiduals(m9b, plot = T) 
testDispersion(m9b) 
MuMIn::AICc(m9b) # 121.0861
plot(allEffects(m9b, partial.residuals = T)) 

anova(m7b,m9b) # in favour of m9b
anova(m6b,m9b) # in favour of m9b
anova(m4,m9b) # in favour of m4
anova(m3,m9b) # slightly in favour of m9b

################## GAMs ###################
## 1 : gam without RE ----
mg1 <- gam(perc_nk_cells ~ s(log_cy, k = 4),
           data = df,
           method = "ML")
plot(mg1, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg1$residuals))))
MuMIn::AICc(mg1) # 116.6261

## 2 : gam with RE intercept ID ----
mg2 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(ID, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg2, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg2$residuals))))
MuMIn::AICc(mg2) # 106.7367

anova(mg1,mg2) # deviance substantial

## 3 : gam with RE intercept day ----
mg3 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg3, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg3$residuals))))
MuMIn::AICc(mg3) # 116.6263

anova(mg1,mg3) # really small deviance

## 4 : gam with RE intercept ID & day ----
mg4 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(ID, bs = "re", k = 2) +
             s(day, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg4, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg4$residuals))))
MuMIn::AICc(mg4) # 106.7368

anova(mg3,mg4) # deviance substantial
anova(mg2,mg4) # really small deviance

## 5 : gam with RE intercept and slope ID ----
mg5 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(ID, bs = "re", k = 2) +
             s(ID, log_cy, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg5, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg5$residuals))))
MuMIn::AICc(mg5) # 106.6315

anova(mg5,mg2) # small deviance
# AICc would indicate mg5 but by a notch, so mg2 prefered, more parsimonious
MuMIn::AICc(mg2, mg5)
MuMIn::AICc(m2,mg2) # mg2 prefered

## 6 : gam with intercept and slope day ----
mg6 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, log_cy, bs = "re", k = 2) +
             s(day, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg6) # 115.5697

anova(mg3,mg6) #  small
MuMIn::AICc(mg6,mg3) # in favour of mg6, not by much

## 7 : gam with intercept and slope day and ID ----
mg7 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, log_cy, bs = "re", k = 2) +
             s(day, bs = "re", k = 2) +
             s(ID, bs = "re", k = 2) +
             s(ID, log_cy, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg7) # 104.0768

anova(mg6, mg7) # deviance substantial
MuMIn::AICc(mg7,mg6) # in favour of mg7
anova(mg5,mg7) # 
MuMIn::AICc(mg7,mg5) # mg7

## 8 : gam with intercept day and ID, slope ID ----
mg8 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, bs = "re", k = 2) +
             s(ID, bs = "re", k = 2) +
             s(ID, log_cy, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg8) # 104.0768

anova(mg7,mg8) # 
MuMIn::AICc(mg8,mg7) # identical, prefer mg8 more parsimonious
anova(mg5,mg8) # 
MuMIn::AICc(mg8,mg5) # mg8
anova(mg4,mg8) # 
MuMIn::AICc(mg8,mg4) # mg8
anova(mg3,mg8) # 
MuMIn::AICc(mg8,mg3) # mg8

## 9 : gam with intercept day and ID, slope day ----
mg9 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, log_cy, bs = "re", k = 2) +
             s(day, bs = "re", k = 2) +
             s(ID, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg9) # 104.0768

anova(mg7,mg9) # 
MuMIn::AICc(mg9,mg7) # identical, mg9 prefered
anova(mg6,mg9) # 
MuMIn::AICc(mg9,mg6) # mg9
anova(mg3,mg9) # 
MuMIn::AICc(mg9,mg3)
anova(mg4,mg9) # 
MuMIn::AICc(mg9,mg4) # mg9

MuMIn::AICc(mg9,mg2)
MuMIn::AICc(mg9,mg8) # identical, equally parsimonious!

########################################################################
# In control cynos, the final model retained would be mg8 or mg9 #######
# mg8 chosen to be consistent with selected model for infected monkeys #
# (see later in the script) ###
# but fits were identical anyway ###
########################################################################
# don't run this unless you want to overwrite our model object
# write_rds(mg8, file = "../output/result_files/cytokines/mixed_effect_models/correlation_NK_cells/Control_cyno_selected_model_IL12.rds")

m_select <- mg8
# if you want to make sure you are using the model we fitted
# m_select <- readRDS("../output/result_files/cytokines/mixed_effect_models/correlation_NK_cells/Control_cyno_selected_model_IL12.rds")
summary(m_select)
min_ <- min(df$log_cy)
max_ <- max(df$log_cy)
newDat <- data.frame(log_cy = seq(min_,max_,0.01))
newDat$ID <- "NV259" # provided to avoid error but won't be used for now
newDat$day <- 1
pred_nl_fix <- predict_gamm(m_select, type = "response",
                            newdata = newDat, se = T,
                            re_form = NA) # NULL = with ranef / NA = no ranef
pred_fix <- data.frame(log_cy = newDat$log_cy,
                       pred_nk = pred_nl_fix$prediction,
                       lwr = pred_nl_fix$prediction + qnorm(0.025)*pred_nl_fix$se,
                       upr = pred_nl_fix$prediction + qnorm(0.975)*pred_nl_fix$se)
# note that the value fixed for day doesn't have an impact
newDat <- data.frame()
for(i in unique(df$ID)){
  test <- df[df$ID == i,]
  min_ <- min(test$log_cy)
  max_ <- max(test$log_cy)
  new <- data.frame(log_cy = seq(min_,max_,0.01))
  new$ID <- i
  new$day <- 1
  newDat <- rbind(newDat,new)
}
pred_nl_RE <- predict_gamm(m_select, type = "response",
                           newdata = newDat, se = F,
                           re_form = NULL) # NULL = with ranef / NA = no ranef
pred_RE <- data.frame(log_cy = newDat$log_cy,
                      pred_nk = pred_nl_RE$prediction,
                      ID = newDat$ID)

p_il12 <- ggplot() + geom_line(data = pred_fix, aes(x = log_cy, y = pred_nk),
                               linewidth = 1.5) +
  geom_line(data = pred_RE, aes(x = log_cy, y = pred_nk, group = ID, col = ID),
            linewidth = 1.2, alpha = 0.5) +
  geom_ribbon(data = pred_fix, aes(x = log_cy, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.5) +
  geom_point(data = df,
             aes(x = log_cy, y = perc_nk_cells,
                 col = ID),
             size = 2.5, alpha = 0.8) +
  labs(x = "IL-12 concentration (log<sub>10</sub> pg/\u03BCl)",
       y = "% NK cells",
       color = "") +
  ggtitle("A - Control") +
  scale_x_continuous(limits = c(2,3.65),
                     expand = expansion(add = 0.1)) +
  scale_y_continuous(expand = expansion(add = 0.1)) +
  coord_cartesian(ylim = c(0,12.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 22),
        axis.title.y = element_markdown(size = 23,
                                        margin = margin(r = 10)),
        axis.title.x = element_markdown(size = 23,
                                        margin = margin(t = 10)),
        legend.position = c(0.15,0.7),
        legend.text = element_text(size = 22),
        plot.title = element_text(size = 23))
png(filename = "../output/figures/suppl/Figure_S3A.png",
    height = 500, width = 800)
plot(p_il12)
dev.off()
##########################################################

# Model selection : infected monkeys -----
df <- df_denv_cyno_complete
df <- df[df$cytokine == "IL.12" & df$day %in% seq(1,7),]
df <- df[complete.cases(df$value, df$perc_nk_cells),]
df$log_cy <- log10(df$value)
df$ID <- as.factor(df$ID)
df$day <- as.factor(df$day)
ggplot(df) + geom_point(aes(x = log_cy, y = perc_nk_cells,
                            color = ID))
# no clear relationship within individuals, 
# but across individuals the positive association seems strong

################### LMs #################
## 1 : lm without RE ----
m1 <- glmmTMB(perc_nk_cells ~ log_cy,
              REML = FALSE,
              data = df)
simulateResiduals(m1, plot = T) 
testDispersion(m1) 
MuMIn::AICc(m1) # 319.9846
plot(allEffects(m1, partial.residuals = T)) 

## 2 : lm with RE intercept ID ----
m2 <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID),
              REML = FALSE,
              data = df)
simulateResiduals(m2, plot = T) 
testDispersion(m2) 
MuMIn::AICc(m2) # 280.3971
plot(allEffects(m2, partial.residuals = T)) 

anova(m1,m2) # signif pref for m2

## 3 : lm with RE intercept day ----
m3 <- glmmTMB(perc_nk_cells ~ log_cy + (1|day),
              REML = FALSE,
              data = df)
simulateResiduals(m3, plot = T) 
testDispersion(m3) 
MuMIn::AICc(m3) # 322.2495
plot(allEffects(m3, partial.residuals = T)) 

anova(m1,m3) # in favour of m1

## 4 : lm with RE intercept ID & day ----
m4 <- glmmTMB(perc_nk_cells ~ log_cy + + (1|ID) + (1|day),
              REML = FALSE,
              data = df)
simulateResiduals(m4, plot = T) 
testDispersion(m4) 
MuMIn::AICc(m4) # 282.2184
plot(allEffects(m4, partial.residuals = T)) 

anova(m2,m4) # no pref for m4 over m2
anova(m3,m4) # pref for m4 over m3

## 5 : lm with RE intercept and slope uncorrelated (ID) to be comparable with the gam later ----
m5b <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID) + (0 + log_cy|ID),
               REML = FALSE,
               data = df)
simulateResiduals(m5b, plot = T) 
testDispersion(m5b) 
MuMIn::AICc(m5b) # 282.7601
plot(allEffects(m5b, partial.residuals = T)) 

anova(m2,m5b) # not in favour of m5b

## 6 : lm with RE intercept and slope uncorrelated (day) to be comparable with the gam later ----
m6b <- glmmTMB(perc_nk_cells ~ log_cy + (1|day) + (0 + log_cy|day),
               REML = FALSE,
               data = df)
# Model convergence problem; non-positive-definite Hessian matrix
diagnose(m6b) # elements related to day unnecessary
simulateResiduals(m6b, plot = T) 
testDispersion(m6b) 
MuMIn::AICc(m6b) # NA
plot(allEffects(m6b, partial.residuals = T)) 

anova(m3,m6b) # not working
anova(m6a,m6b) # not working

## 7 : lm with RE intercept and slope uncorrelated (ID & day) ----
m7b <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID) + (0 + log_cy|ID) +
                 (1|day) + (0 + log_cy|day),
               REML = FALSE,
               data = df)
simulateResiduals(m7b, plot = T) 
testDispersion(m7b) 
MuMIn::AICc(m7b) # 286.6163
plot(allEffects(m7b, partial.residuals = T)) 

anova(m6b,m7b) # not working
anova(m5b,m7b) # not in favour of m7b

## 8 : lm with RE intercept (ID and day) and slope ID, uncorrelated ----
m8b <- glmmTMB(perc_nk_cells ~ log_cy + (1|ID) + (0 + log_cy|ID) +
                 (1|day),
               REML = FALSE,
               data = df)
simulateResiduals(m8b, plot = T) 
testDispersion(m8b) 
MuMIn::AICc(m8b) # 284.6658
plot(allEffects(m8b, partial.residuals = T)) 

anova(m7b,m8b) # in favour of m8b
anova(m5b,m8b) # in favour of m5b
anova(m4,m8b) # in favour of m4
anova(m3,m8b) # in favour of m8b

## 9 : lm with RE intercept (ID and day) and slope day, uncorrelated ----
m9b <- glmmTMB(perc_nk_cells ~ log_cy + (1|day) + (0 + log_cy|day) +
                 (1|ID),
               REML = FALSE,
               data = df)
simulateResiduals(m9b, plot = T) 
testDispersion(m9b) 
MuMIn::AICc(m9b) # 284.0799
plot(allEffects(m9b, partial.residuals = T)) 

anova(m7b,m9b) # in favour of m9b
anova(m6b,m9b) # not working
anova(m4,m9b) # in favour of m4
anova(m3,m9b) # in favour of m9b

################## GAMs ###################
## 1 : gam without RE ----
mg1 <- gam(perc_nk_cells ~ s(log_cy, k = 4),
           data = df,
           method = "ML")
plot(mg1, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg1$residuals))))
MuMIn::AICc(mg1) # 310.8467

## 2 : gam with RE intercept ID ----
mg2 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(ID, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg2, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg2$residuals))))
MuMIn::AICc(mg2) # 259.9148

anova(mg1,mg2) # deviance substantial

## 3 : gam with RE intercept day ----
mg3 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg3, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg3$residuals))))
MuMIn::AICc(mg3) # 309.8193

anova(mg1,mg3) # deviance substantial (AICc quite close)
MuMIn::AICc(mg1,mg3)

## 4 : gam with RE intercept ID & day ----
mg4 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(ID, bs = "re", k = 2) +
             s(day, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg4, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg4$residuals))))
MuMIn::AICc(mg4) # 266.3259

anova(mg3,mg4) # deviance substantial
anova(mg2,mg4) # deviance not negligible but AICc favours mg2

## 5 : gam with RE intercept and slope ID ----
mg5 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(ID, bs = "re", k = 2) +
             s(ID, log_cy, bs = "re", k = 2),
           data = df,
           method = "ML")
plot(mg5, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
     col = c("black", rep("red", length(mg5$residuals))))
MuMIn::AICc(mg5) # 259.907

anova(mg5,mg2) # deviance really small 
MuMIn::AICc(mg2,mg5)
# AICc would indicate mg5 but by a notch, so mg2 prefered

MuMIn::AICc(m2,mg2) # mg2 prefered

## 6 : gam with intercept and slope day ----
mg6 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, log_cy, bs = "re", k = 2) +
             s(day, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg6) # 304.2643

anova(mg3,mg6) # 
MuMIn::AICc(mg6,mg3) # in favour of mg6

## 7 : gam with intercept and slope day and ID ----
mg7 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, log_cy, bs = "re", k = 2) +
             s(day, bs = "re", k = 2) +
             s(ID, bs = "re", k = 2) +
             s(ID, log_cy, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg7) # 257.1102

anova(mg6, mg7) # 
MuMIn::AICc(mg7,mg6) # in favour of mg7
anova(mg5,mg7) # 
MuMIn::AICc(mg7,mg5) # mg7

## 8 : gam with intercept day and ID, slope ID ----
mg8 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, bs = "re", k = 2) +
             s(ID, bs = "re", k = 2) +
             s(ID, log_cy, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg8) # 256.6784

anova(mg8,mg7) # deviance quite small
MuMIn::AICc(mg8,mg7) # mg8 but not by much
anova(mg5,mg8) # 
MuMIn::AICc(mg8,mg5) # mg8
anova(mg4,mg8) # 
MuMIn::AICc(mg8,mg4) # mg8
anova(mg3,mg8) # 
MuMIn::AICc(mg8,mg3) # mg8

## 9 : gam with intercept day and ID, slope day -----
mg9 <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
             s(day, log_cy, bs = "re", k = 2) +
             s(day, bs = "re", k = 2) +
             s(ID, bs = "re", k = 2),
           REML = FALSE,
           data = df)
MuMIn::AICc(mg9) # 259.679

anova(mg7,mg9) # 
MuMIn::AICc(mg9,mg7) # mg7
anova(mg4,mg9) # 
MuMIn::AICc(mg9,mg4) # mg9
anova(mg9,mg3) # 
MuMIn::AICc(mg9,mg3)
anova(mg6,mg9) # 
MuMIn::AICc(mg9,mg6)

###############################################################
# In infected cynos, the final model retained would be mg8 ####
###############################################################
# Trend seems driven by one point with low cytokine,
# see how it looks without this point
# conclusions remain the same / keep initial model
# mg8_no_outlier <- gam(perc_nk_cells ~ s(log_cy, k = 4) +
#              s(day, bs = "re", k = 2) +
#              s(ID, bs = "re", k = 2) +
#              s(ID, log_cy, bs = "re", k = 2),
#            REML = FALSE,
#            data = df[df$log_cy > 2.5,])
# summary(mg8_no_outlier)

# don't run this unless you want to overwrite our model object
# write_rds(mg8, file = "../output/result_files/cytokines/mixed_effect_models/correlation_NK_cells/Infected_cyno_selected_model_IL12.rds")

m_select <- mg8
# if you want to make sure you are using the model we fitted
# m_select <- readRDS("../output/result_files/cytokines/mixed_effect_models/correlation_NK_cells/Infected_cyno_selected_model_IL12.rds")

min_ <- min(df$log_cy)
max_ <- max(df$log_cy)
newDat <- data.frame(log_cy = seq(min_,max_,0.01))
newDat$ID <- "BC116A" # provided to avoid error but won't be used for now
newDat$day <- 1 # provided to avoid error but won't be used for now
pred_nl_fix <- predict_gamm(m_select, type = "response",
                            newdata = newDat, se = T,
                            re_form = NA) # NULL = with ranef / NA = no ranef
pred_fix <- data.frame(log_cy = newDat$log_cy,
                       pred_nk = pred_nl_fix$prediction,
                       lwr = pred_nl_fix$prediction + qnorm(0.025)*pred_nl_fix$se,
                       upr = pred_nl_fix$prediction + qnorm(0.975)*pred_nl_fix$se)

newDat <- data.frame()
# Note : changing the day doesn't change the final plot here
for(i in unique(df$ID)){
  test <- df[df$ID == i,]
  min_ <- min(test$log_cy)
  max_ <- max(test$log_cy)
  new <- data.frame(log_cy = seq(min_,max_,0.01))
  new$ID <- i
  new$day <- 1
  newDat <- rbind(newDat,new)
}
pred_nl_RE <- predict_gamm(m_select, type = "response",
                           newdata = newDat, se = F,
                           re_form = NULL) # NULL = with ranef / NA = no ranef
pred_RE <- data.frame(log_cy = newDat$log_cy,
                      pred_nk = pred_nl_RE$prediction,
                      ID = newDat$ID,
                      day = newDat$day)

p_il12 <- ggplot() + 
  geom_line(data = pred_fix, aes(x = log_cy, y = pred_nk),
            linewidth = 1.5) +
  geom_line(data = pred_RE, aes(x = log_cy, y = pred_nk, group = ID, col = ID),
            linewidth = 1.2, alpha = 0.5) +
  geom_ribbon(data = pred_fix, aes(x = log_cy, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.5) +
  geom_point(data = df,
             aes(x = log_cy, y = perc_nk_cells,
                 color = ID),
             size = 2.5, alpha = 0.8) +
  labs(x = "IL-12 concentration (log<sub>10</sub> pg/\u03BCl)",
       y = "% NK cells",
       color = "") +
  ggtitle("B - DENV-infected") +
  guides(color = guide_legend(ncol = 2)) +
  scale_x_continuous(limits = c(2,3.65),
                     expand = expansion(add = 0.1)) +
  scale_y_continuous(expand = expansion(add = 0.1)) +
  coord_cartesian(ylim = c(0,12.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 22),
        axis.title.y = element_markdown(size = 23,
                                        margin = margin(r = 10)),
        axis.title.x = element_markdown(size = 23,
                                        margin = margin(t = 10)),
        legend.position = c(0.33,0.8),
        legend.text = element_text(size = 22),
        plot.title = element_text(size = 23))

png(filename = "../output/figures/suppl/Figure_S3B.png",
    height = 500, width = 800)
plot(p_il12)
dev.off()


# Plots focusing on temporal trend, no model fit ----

df <- df_control_cyno_complete
df <- df[df$cytokine == "IL.12" & df$day %in% seq(1,7),]
df <- df[complete.cases(df$value, df$perc_nk_cells),]
df$log_cy <- log10(df$value)
df$ID <- as.factor(df$ID)
df$day <- as.factor(df$day)

p_il12 <- ggplot() + 
  geom_point(data = df,
             aes(x = log_cy, y = perc_nk_cells,
                 color = ID),
             size = 3, alpha = 0.8) +
  scale_color_discrete(drop = FALSE) +
  labs(x = "IL-12 concentration (log<sub>10</sub> pg/\u03BCl)",
       y = "% NK cells",
       color = "") +
  ggtitle("A - Control") +
  guides(color = guide_legend(ncol = 3)) +
  facet_wrap(~ day) +
  scale_x_continuous(limits = c(2,3.65),
                     expand = expansion(add = 0.1)) +
  scale_y_continuous(limits = c(0,12.5),
                     expand = expansion(add = 0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 22),
        axis.title.y = element_markdown(size = 24,
                                        margin = margin(r = 10)),
        axis.title.x = element_markdown(size = 24,
                                        margin = margin(t = 10)),
        legend.position = c(0.67,0.13),
        strip.text.x = element_text(size = 19),
        legend.text = element_text(size = 21),
        plot.title = element_text(size = 24))

png(filename = "../output/figures/suppl/Figure_S4A.png",
    height = 750, width = 1200)
plot(p_il12)
dev.off()

df <- df_denv_cyno_complete
df <- df[df$cytokine == "IL.12" & df$day %in% seq(1,7),]
df <- df[complete.cases(df$value, df$perc_nk_cells),]
df$log_cy <- log10(df$value)
df$ID <- as.factor(df$ID)
df$day <- as.factor(df$day)

# TTM : transmission to mosquitoes
df$vir_TTM <- FALSE
df$vir_TTM[df$ID == "FR469A" & df$day == 4] <- TRUE
df$vir_TTM[df$ID == "BC116A" & df$day == 4] <- TRUE
df$vir_TTM[df$ID == "FR423A" & df$day == 4] <- TRUE
df$vir_TTM[df$ID == "FR423A" & df$day == 5] <- TRUE
df$vir_TTM[df$ID == "BC167" & df$day == 5] <- TRUE
df$vir_TTM[df$ID == "FR423A" & df$day == 6] <- TRUE
df$vir_TTM[df$ID == "FR423A" & df$day == 7] <- TRUE

p_il12 <- ggplot() + 
  geom_point(data = df[df$vir_TTM == FALSE,],
             aes(x = log_cy, y = perc_nk_cells,
                 color = ID),
             size = 3, alpha = 0.8) +
  geom_point(data = df[df$vir_TTM == TRUE,],
             aes(x = log_cy, y = perc_nk_cells,
                 fill = ID),
             stroke = 1.5,
             shape = 21,
             size = 3, alpha = 0.8) +
  scale_color_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE, guide = "none") +
  labs(x = "IL-12 concentration (log<sub>10</sub> pg/\u03BCl)",
       y = "% NK cells",
       color = "") +
  ggtitle("B - DENV-infected") +
  guides(color = guide_legend(ncol = 3)) +
  facet_wrap(~ day) +
  scale_x_continuous(limits = c(2,3.65),
                     expand = expansion(add = 0.1)) +
  scale_y_continuous(limits = c(0,12.5),
                     expand = expansion(add = 0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 22),
        axis.title.y = element_markdown(size = 24,
                                        margin = margin(r = 10)),
        axis.title.x = element_markdown(size = 24,
                                        margin = margin(t = 10)),
        legend.position = c(0.67,0.13),
        strip.text.x = element_text(size = 19),
        legend.text = element_text(size = 21),
        plot.title = element_text(size = 24))

png(filename = "../output/figures/suppl/Figure_S4B.png",
    height = 750, width = 1200)
plot(p_il12)
dev.off()

