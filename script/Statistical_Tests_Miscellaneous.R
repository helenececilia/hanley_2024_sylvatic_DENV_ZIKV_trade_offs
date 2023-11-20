## ---------------------------
##
## Script name: Statistical_Tests_Miscellaneous
##
## Purpose of script: Regroups all basic tests performed on our datasets
##
## Author: Helene Cecilia
##
## Date Created: 2023-02-17

rm(list=ls())

## Loading Packages  ------------------
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(effects)
library(janitor)
library(MuMIn)
library(multcomp) # for glht
library(lsr) # for etaSquared

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# Cynomolgus macaques / ZIKV ----

## Just some data plotting : too few individuals to run tests ----
# Really close doses and peak titers between individuals
df <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
df0 <- df[df$day_post_infection == 0,]
df0 <- df0[,c("id",
              "corrected_estimated_titer_delivered_log10")]
colnames(df0) <- c("id","dose_passage_log10")
vir <- df %>% group_by(id) %>% summarise(max_vir_raw = max(viremia_deduced, na.rm = T))
test <- merge(df0,vir,by = "id")

ggplot(test) + geom_point(aes(x = dose_passage_log10, y = max_vir_raw))

# early % NK cells of ZIKV-infected within the range of control cynomolgus macaques ----
control <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                    sep = "\t", dec = ".")
control <- control %>% clean_names()
control <- control[control$final_treatment == "Control",]
control <- control[control$day_post_infection == 1,]
df1 <- df[df$day_post_infection == 1,]
df1 <- df1[,c("id","final_treatment","x_nk_cells")]
control <- control[,c("id","final_treatment","x_nk_cells")]
test <- rbind(df1,control)

ggplot(test) + geom_jitter(aes(x = final_treatment, y = x_nk_cells)) +
  geom_boxplot(aes(x = final_treatment, y = x_nk_cells),
               fill = NA)

# early NK cell mobilization and neutralizing antibody titers (PRNT80) 
df1 <- df[df$day_post_infection == 1,]
df1 <- df1[,c("id","x_nk_cells")]
df28 <- df[df$day_post_infection == 28,]
df28 <- df28[,c("id","prnt80")]
test <- merge(df1,df28, by = "id")

test$prnt80[test$prnt80 == ">640"] <- "640" #underestimation
test$prnt80 <- as.numeric(test$prnt80)
ggplot(test) + geom_point(aes(x = x_nk_cells, y = log10(prnt80)))


# early NK cells and peak virus titer
# later analyzed combined with squirrel monkeys
df_nk <- df[df$day_post_infection == 1,c("id","x_nk_cells")]
df_nk <- df_nk[complete.cases(df_nk),]
df_vir <- df %>% group_by(id) %>% summarise(max_vir = max(viremia_deduced, na.rm = T))
test <- merge(df_nk, df_vir,by = "id")

ggplot(test) + geom_point(aes(x = x_nk_cells, y = max_vir))


# Cynomolgus macaques / DENV ----
## Impact of number of mosquitoes with virus-positive saliva on the likelihood of becoming detectably viremic ----
df <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")
# used for nb_saliva below
df[df$Day.Post.Infection == 0 & df$Final.Treatment != "Control",
   c("ID","Number..IT.injected.mosquitoes.with.saliva.positive.for.DENV.2.after.one.passage.of.saliva.in.C6.36.cells")]

# used for vir_status_detect below
# viremic status defined as detection through serum titer
df1 <- df[df$Viremia..log10pfu.ml. != 0,]
unique(df1$ID)

# used for vir_status_deduced below
# viremic status defined as transmission to mosquitoes (body, legs, or saliva)
df2 <- df[df$Viremia.deduced != 0,]
unique(df2$ID)

denv_cyno <- data.frame(nb_saliva = c(0,1,1,1,6,5,3,4,3),
                        vir_status_detect = c(0,1,0,0,1,1,0,0,0),
                        vir_status_deduced = c(0,1,0,0,1,1,0,1,1))

denv_cyno$vir_status_deduced <- as.factor(denv_cyno$vir_status_deduced)
m <- glmmTMB(nb_saliva ~ vir_status_deduced,
              data = denv_cyno,
              family = poisson)
simulateResiduals(m, plot = T) # ok
testDispersion(m) # ok 
summary(m)
car::Anova(m, type = "II")

## Differences in early NK cells between treatments ----
df <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[,c("id","final_treatment","day_post_infection",
            "x_nk_cells","number_nk_cells")]
df1 <- df[df$day_post_infection == 1,]

df1$final_treatment <- as.factor(df1$final_treatment)
m1 <- lm(x_nk_cells ~ final_treatment,
         data = df1)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
plot(allEffects(m1, partial.residuals = T))
summary(m1) # difference between low/high exposure groups detected (pvalue close to 0.05)

# when we adjust pvalues, individual comparisons between groups are not significant
tuk = glht(m1, linfct = mcp(final_treatment = "Tukey"))
summary(tuk, test = adjusted("fdr"))

car::Anova(m1, type = "II") # overall the contribution of the treatment variable is not significant


## Differences in day 28 PRNT80 values between treatments ----
df <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[,c("id","final_treatment","day_post_infection",
            "prnt80")]
df1 <- df[df$day_post_infection == 28,]
df1 <- df1[df1$final_treatment != "Control",]

m1 <- lm(log10(prnt80) ~ final_treatment,
         data = df1)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
plot(allEffects(m1, partial.residuals = T))

summary(m1)
car::Anova(m1, type = "II")
t.test(log10(prnt80) ~ final_treatment,
       data = df1)

## Association between early NK cell mobilization and neutralizing antibody titers : REMOVED ----
df <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[,c("id","final_treatment","day_post_infection",
            "x_nk_cells","number_nk_cells","prnt80")]
# infected only
df <- df[df$final_treatment != "Control",]
df1 <- df[df$day_post_infection == 1,]
df1 <- df1[,c("id","x_nk_cells")]
df28 <- df[df$day_post_infection == 28,]
df28 <- df28[,c("id","prnt80")]
test <- merge(df1,df28, by = "id")

m0 <- lm(log10(prnt80) ~ x_nk_cells,
         data = test)
simulateResiduals(m0, plot = T) # ok
testDispersion(m0) # ok
plot(allEffects(m0, partial.residuals = T))
summary(m0)
# the two ways of computing the CIs do not give exactly the same results 
# (one includes 0 and not the other one)
# we report confint
confint(m0)
-0.12380 - 1.96*0.05377
-0.12380 + 1.96*0.05377
car::Anova(m0, type = "II")


# Comparison cynos / ZIKV vs DENV ----
## Just some data plotting : too few individuals in the ZIKV group to run tests ----
# PRNT80 values between viruses 
denv <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                 sep = "\t", dec = ".")
zikv <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Cynomolgus_Macaques.csv",
                 sep = "\t", dec = ".")
denv <- denv %>% clean_names()
zikv <- zikv %>% clean_names()

denv <- denv[denv$final_treatment != "Control",]
zikv <- zikv[zikv$final_treatment != "Control",]

denv_df28 <- denv[denv$day_post_infection == 28,]
denv_df28 <- denv_df28[,c("id","prnt80")]
denv_df28$virus <- "DENV"

zikv_df28 <- zikv[zikv$day_post_infection == 28,]
zikv_df28 <- zikv_df28[,c("id","prnt80")]
zikv_df28$virus <- "ZIKV"

df28 <- rbind(denv_df28, zikv_df28)
df28$prnt80[df28$prnt80 == ">640"] <- "640" # LOWER BOUND 
df28$prnt80 <- as.numeric(df28$prnt80)

ggplot(df28) + geom_point(aes(x = virus, y = prnt80))

# Early % NK cells between viruses 
denv_df <- denv[denv$day_post_infection == 1,]
denv_df <- denv_df[,c("id","x_nk_cells")]
denv_df$virus <- "DENV"

zikv_df <- zikv[zikv$day_post_infection == 1,]
zikv_df <- zikv_df[,c("id","x_nk_cells")]
zikv_df$virus <- "ZIKV"

df <- rbind(denv_df, zikv_df)

ggplot(df) + geom_point(aes(x = virus, y = x_nk_cells))

# Squirrel monkeys / DENV ----
## Impact of the dose of DENV delivered on the likelihood of squirrel monkeys becoming viremic ----

# here viremic = detectable virus in raw or passaged serum or infection of mosquitoes 
df <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df1 <- df[df$viremia_pfu_ml != 0,]
unique(df1$id) # initially 3 monkeys viremic (detected in serum)
df <- df[,c("id","estimated_titer_delivered","corrected_estimated_titer_delivered_log10","viremia_deduced")]
df2 <- df[df$viremia_deduced != 0,]
unique(df2$id) # 7 viremic monkeys once we account for transmission to mosquitoes
df <- df[!is.na(df$estimated_titer_delivered),]
df$vir_status_deduced <- c(1,0,1,0,1,1,1,1,1,0)

m0 <- glmmTMB(vir_status_deduced ~ corrected_estimated_titer_delivered_log10,
              data = df,
              family = binomial)
simulateResiduals(m0, plot = T) # issues
testDispersion(m0) # ok
summary(m0) # non signif 
car::Anova(m0, type = "II")

## Differences in early % NK cells between DENV-2 infected and control squirrel monkeys ----
df <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
control <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                    sep = "\t", dec = ".")
control <- control %>% clean_names()
control <- control[control$final_treatment == "Control",]
control <- control[control$day_post_infection %in% c(1,2),]
df <- df[df$day_post_infection %in% c(1,2),]
df <- df[,c("id","final_treatment","x_nk_cells")]
control <- control[,c("id","final_treatment","x_nk_cells")]
test <- rbind(df,control)

m0 <- lm(x_nk_cells ~ final_treatment,
         data = test)
simulateResiduals(m0, plot = T) # ok
testDispersion(m0) # ok
summary(m0)
car::Anova(m0, type = "II")

t.test(x_nk_cells ~ final_treatment,
       data = test)

## Relationship between early % NK cells and PRNT80 values on Day 28 : REMOVED ----
df <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[,c("id","final_treatment","day_post_infection",
            "x_nk_cells","number_nk_cells","prnt80")]
# infected only
df <- df[df$final_treatment != "Control",]
df1 <- df[df$day_post_infection %in% c(1,2),]
df1 <- df1[,c("id","x_nk_cells")]
df28 <- df[df$day_post_infection == 28,]
df28 <- df28[,c("id","prnt80")]
test <- merge(df1,df28, by = "id")
test <- test[!is.na(test$x_nk_cells),]

m2 <- lm(log10(prnt80) ~ x_nk_cells,
         data = test)
# warning
simulateResiduals(m2, plot = T) # ok
testDispersion(m2) # ok
summary(m2) 
# not exactly the same CI / we report confint
confint(m2)
-0.11437 - 1.96*0.05057
-0.11437 + 1.96*0.05057

car::Anova(m2, type = "II")

# Comparison squirrel - cyno / DENV ----
df_sq <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                  sep = "\t", dec = ".")
df_sq <- df_sq %>% clean_names()
df_sq <- df_sq[,c("id","final_treatment","day_post_infection",
                  "x_nk_cells","number_nk_cells","prnt80")]
df_sq <- df_sq[df_sq$final_treatment != "Control",]
df_sq1 <- df_sq[df_sq$day_post_infection %in% c(1,2),]
df_sq1 <- df_sq1[,c("id","x_nk_cells")]
df_sq28 <- df_sq[df_sq$day_post_infection == 28,]
df_sq28 <- df_sq28[,c("id","prnt80")]
test_sq <- merge(df_sq1,df_sq28, by = "id")
test_sq <- test_sq[!is.na(test_sq$x_nk_cells),]
test_sq$species <- "squirrel"

df_cy <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                  sep = "\t", dec = ".")
df_cy <- df_cy %>% clean_names()
df_cy <- df_cy[,c("id","final_treatment","day_post_infection",
                  "x_nk_cells","number_nk_cells","prnt80")]
df_cy <- df_cy[df_cy$final_treatment != "Control",]
df_cy1 <- df_cy[df_cy$day_post_infection == 1,]
df_cy1 <- df_cy1[,c("id","x_nk_cells")]
df_cy28 <- df_cy[df_cy$day_post_infection == 28,]
df_cy28 <- df_cy28[,c("id","prnt80")]
test_cy <- merge(df_cy1,df_cy28, by = "id")
test_cy <- test_cy[!is.na(test_cy$x_nk_cells),]
test_cy$species <- "cyno"

test <- rbind(test_cy, test_sq)
test$log_prnt <- log10(test$prnt80)

## Differences between species, for both NK cells and PRNT80 ----
t.test(log_prnt ~ species, data = test, var.equal = TRUE)
t.test(x_nk_cells ~ species, data = test, var.equal = TRUE)

## Possible species effect in relationship between early NK cells and PRNT80 during DENV infection ----
m1 <- lm(log_prnt ~ x_nk_cells*species,
         data = test)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
summary(m1) 
confint(m1)

M <- aov(m1)
etaSquared(M) # interaction < 0.01 so type II preferred
car::Anova(m1, type = "II")


# ZIKV - squirrel ----
## Relationship between number of mosquitoes salivating virus and total dose delivered ----
df <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
df <- df[df$day_post_infection == 0,]
df <- df[,c("id",
            "number_of_mosquitoes_with_detectable_denv_2_in_collected_saliva",
            "number_it_injected_mosquitoes_with_saliva_positive_for_denv_2_after_one_passage_of_saliva_in_c6_36_cells",
            "raw_estimated_titer_delivered_ffu",
            "corrected_estimated_titer_delivered_log10")]
colnames(df) <- c("id","nb_saliva_raw","nb_saliva_passage",
                  "dose_raw","dose_passage_log10")

m1 <- lm(dose_passage_log10 ~ nb_saliva_passage,
         data = df)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
summary(m1)
confint(m1)
car::Anova(m1, type = "II")

## Effect of total dose delivered peak titer of ZIKV ----
df <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
df0 <- df[df$day_post_infection == 0,]
df0 <- df0[,c("id",
              "corrected_estimated_titer_delivered_log10")]
colnames(df0) <- c("id","dose_passage_log10")
vir <- df %>% group_by(id) %>% summarise(max_vir_raw = max(viremia_pfu_ml, na.rm = T),
                                         max_vir_corr = max(viremia_deduced, na.rm = T))
test <- merge(df0,vir,by = "id")

# without monkey 4683
m0 <- lm(log10(max_vir_raw) ~ dose_passage_log10,
         data = test[test$id != 4683,])
simulateResiduals(m0, plot = T) # ok
testDispersion(m0) # ok
summary(m0)
car::Anova(m0, type = "II")

# with monkey 4683
m1 <- lm(log10(max_vir_raw) ~ dose_passage_log10,
         data = test)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
summary(m1)
car::Anova(m1, type = "II")

## Differences in early % NK cells between ZIKV-infected and control squirrel monkeys ----
df <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
control <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                    sep = "\t", dec = ".")
df <- df %>% clean_names()
control <- control %>% clean_names()
control <- control[control$final_treatment == "Control",]
control <- control[control$day_post_infection %in% c(1,2),]
df <- df[df$day_post_infection %in% c(1,2),]
df <- df[,c("id","final_treatment","x_nk_cells")]
control <- control[,c("id","final_treatment","x_nk_cells")]
test <- rbind(df,control)

test <- test[complete.cases(test),]
test$final_treatment <- as.factor(test$final_treatment)
test <- within(test, final_treatment <- relevel(final_treatment, ref = "Control"))

# with monkey 4683
m0 <- lm(x_nk_cells ~ final_treatment,
         data = test)
simulateResiduals(m0, plot = T) # ok
testDispersion(m0) # ok
summary(m0) 
car::Anova(m0, type = "II")

# without monkey 4683 (results reported in manuscript)
m1 <- lm(x_nk_cells ~ final_treatment,
         data = test[test$id != 4683,])
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
summary(m1)
car::Anova(m1, type = "II")
t.test(x_nk_cells ~ final_treatment,
       data = test[test$id != 4683,])

## Relationship between early NK cells and peak virus titer ----
df <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
df_nk <- df[df$day_post_infection %in% c(1,2),c("id","x_nk_cells")]
df_nk <- df_nk[complete.cases(df_nk),]
df_vir <- df %>% group_by(id) %>% summarise(max_vir = max(log10(viremia_deduced), na.rm = T))
test <- merge(df_nk, df_vir,by = "id")

# without monkey 4683 : significant
m0 <- lm(max_vir ~ x_nk_cells,
         data = test[test$id != 4683,])
simulateResiduals(m0, plot = T) # ok
testDispersion(m0) # ok
summary(m0) # signif (0.009)
car::Anova(m0, type = "II")

# with monkey 4683 : not significant
m1 <- lm(max_vir ~ x_nk_cells,
         data = test)
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
summary(m1) # not signif
car::Anova(m1, type = "II")

# Comparison squirrels / ZIKV vs DENV ----
## Differences in PRNT80 values between viruses ----
denv <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                 sep = "\t", dec = ".")
zikv <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                 sep = "\t", dec = ".")
denv <- denv %>% clean_names()
zikv <- zikv %>% clean_names()

denv <- denv[denv$final_treatment != "Control",]
zikv <- zikv[zikv$final_treatment != "Control",]

denv_df28 <- denv[denv$day_post_infection == 28,]
denv_df28 <- denv_df28[,c("id","prnt80")]
denv_df28$virus <- "DENV"

zikv_df28 <- zikv[zikv$day_post_infection == 28,]
zikv_df28 <- zikv_df28[,c("id","prnt80")]
zikv_df28$virus <- "ZIKV"

df28 <- rbind(denv_df28, zikv_df28)

m0 <- lm(log10(prnt80) ~ virus,
         data = df28)
simulateResiduals(m0, plot = T) # ok
testDispersion(m0) # ok
summary(m0)
car::Anova(m0, type = "II")
t.test(log10(prnt80) ~ virus,
       data = df28)

## Differences in early % NK cells between viruses ----
denv <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                 sep = "\t", dec = ".")
zikv <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                 sep = "\t", dec = ".")
denv <- denv %>% clean_names()
zikv <- zikv %>% clean_names()

denv <- denv[denv$final_treatment != "Control",]
zikv <- zikv[zikv$final_treatment != "Control",]

denv_df <- denv[denv$day_post_infection %in% c(1,2),]
denv_df <- denv_df[,c("id","x_nk_cells")]
denv_df$virus <- "DENV"

zikv_df <- zikv[zikv$day_post_infection %in% c(1,2),]
zikv_df <- zikv_df[,c("id","x_nk_cells")]
zikv_df$virus <- "ZIKV"

df <- rbind(denv_df, zikv_df)

# with monkey 4683
m0 <- lm(x_nk_cells ~ virus,
         data = df)
simulateResiduals(m0, plot= T) # ok
testDispersion(m0) # ok
summary(m0)
car::Anova(m0, type = "II")
t.test(x_nk_cells ~ virus,
       data = df)

# without monkey 4683
m1 <- lm(x_nk_cells ~ virus,
         data = df[df$id != 4683,])
simulateResiduals(m1, plot= T) # ok
testDispersion(m1) # ok
summary(m1)
car::Anova(m1, type = "II")
t.test(x_nk_cells ~ virus,
       data = df[df$id != 4683,])

## Differences between viruses in effect of dose on peak viremia ----
denv <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                 sep = "\t", dec = ".")
zikv <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                 sep = "\t", dec = ".")
denv <- denv %>% clean_names()
zikv <- zikv %>% clean_names()

denv <- denv[denv$final_treatment != "Control",]
zikv <- zikv[zikv$final_treatment != "Control",]

denv$virus <- "DENV"
zikv$virus <- "ZIKV"
denv <- denv[,c("id",
                "virus",
                "day_post_infection",
                "viremia_deduced",
                "corrected_estimated_titer_delivered_log10")]
zikv <- zikv[,c("id",
                "virus",
                "day_post_infection",
                "viremia_deduced",
                "corrected_estimated_titer_delivered_log10")]
df <- rbind(denv,zikv)
df0 <- df[df$day_post_infection == 0,]
df0 <- df0[,c("id",
              "virus",
              "corrected_estimated_titer_delivered_log10")]
colnames(df0) <- c("id","virus","dose_passage_log10")
vir <- df %>% group_by(id,virus) %>% summarise(max_vir = max(viremia_deduced, na.rm = T)) %>% ungroup()
test <- merge(df0,vir,by = c("id","virus"))
test$log_vir_max <- log10(test$max_vir + 1)

# with monkey 4683
m0 <- lm(log_vir_max ~ dose_passage_log10*virus,
         data = test)
simulateResiduals(m0, plot = T) # ok
testDispersion(m0) # ok
plot(allEffects(m0, partial.residuals = T))
M <- aov(m0)
etaSquared(M) # interaction < 0.01, type II prefered
car::Anova(m0, type = "II") # virus signif

# without monkey 4683 : results reported in manuscript
m1 <- lm(log_vir_max ~ dose_passage_log10*virus,
         data = test[test$id != 4683,])
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
plot(allEffects(m1, partial.residuals = T))
M <- aov(m1)
etaSquared(M) # interaction < 0.01, type II prefered
car::Anova(m1, type = "II") # virus signif

# if we focus on a subset of doses common to DENV and ZIKV monkeys
select <- test[test$dose_passage_log10 <= 3.3 & test$dose_passage_log10 >=2.3,]
# note : this does not include monkey 4683

ms0 <- lm(log_vir_max ~ virus, 
          data = select)
simulateResiduals(ms0, plot = T) # ok
testDispersion(ms0) # ok
plot(allEffects(ms0, partial.residuals = T))
summary(ms0) 
car::Anova(ms0, type = "II") 


# Grouped analysis squirrel - cyno / ZIKV ----
## Relationship between peak viremia and duration ----
df <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
vir <- df[!is.na(df$viremia_deduced) & df$viremia_deduced != 0,]
vir <- vir[,c("id","day_post_infection","viremia_deduced","viremia_pfu_ml")]
vir$log_V <- log10(vir$viremia_deduced)
test1 <- vir %>% group_by(id) %>% summarise(max_vir = max(log_V, na.rm = T),
                                           min_day = min(day_post_infection),
                                           max_day = max(day_post_infection)) %>%
  ungroup()
test1$duration <- (test1$max_day - test1$min_day) +1

df <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
vir <- df[!is.na(df$viremia_deduced) & df$viremia_deduced != 0,]
vir <- vir[,c("id","day_post_infection","viremia_deduced","viremia_pfu_ml")]
vir$log_V <- vir$viremia_deduced # already in log
test2 <- vir %>% group_by(id) %>% summarise(max_vir = max(log_V, na.rm = T),
                                           min_day = min(day_post_infection),
                                           max_day = max(day_post_infection)) %>%
  ungroup()
test2$duration <- (test2$max_day - test2$min_day) +1

my_df <- rbind(test1,test2)
my_df$duration_factor <- as.factor(my_df$duration)
with_Minerva <- my_df
no_Minerva <- my_df[my_df$id != 4683,]

# Plot : 4683 is at the top right
ggplot(my_df) + geom_point(aes(x = max_vir, y = duration, color = as.factor(min_day)))

# results of ordinal logistic regression reported in manuscript, not run in R

## Effect of early % NK cells on peak viremia ----
df <- read.csv("../data/Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
df_nk <- df[df$day_post_infection %in% c(1,2),c("id","x_nk_cells")]
df_nk <- df_nk[complete.cases(df_nk),]
df_vir <- df %>% group_by(id) %>% summarise(max_vir = max(log10(viremia_deduced), na.rm = T))
df_sq <- merge(df_nk, df_vir,by = "id")
df_sq$species <- "Squirrel"

df <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")
df <- df %>% clean_names()
df <- df[df$final_treatment != "Control",]
df_nk <- df[df$day_post_infection == 1,c("id","x_nk_cells")]
df_nk <- df_nk[complete.cases(df_nk),]
df_vir <- df %>% group_by(id) %>% summarise(max_vir = max(viremia_deduced, na.rm = T))
df_cy <- merge(df_nk, df_vir,by = "id")
df_cy$species <- "Cyno"

test <- rbind(df_sq, df_cy)
ggplot(test) + geom_point(aes(x = x_nk_cells, y = max_vir,
                              color = species))

m <- lm(max_vir ~ x_nk_cells,
         data = test)
simulateResiduals(m, plot = T) # no issues 
testDispersion(m) # ok 
summary(m) # not signif negative 


m1 <- lm(max_vir ~ x_nk_cells,
         data = test[test$id != "4683",])
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
summary(m1) # signif negative
