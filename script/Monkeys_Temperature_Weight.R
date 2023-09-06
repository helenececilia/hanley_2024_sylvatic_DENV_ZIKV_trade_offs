## ---------------------------
##
## Script name: Monkeys_Temperature_Weight.R
##
## Purpose of script: Plot supplementary figures for temperature and weight of monkeys
## over the course of the experiments. Also tests is weight change differs across treatment groups
##
## Author: Helene Cecilia
##
## Date Created: 2022-12-16

rm(list=ls())

## Loading Packages  ------------------
library(ggplot2)
library(patchwork)
library(lubridate)
library(chron)
library(glmmTMB) # inconsistency between Matrix versions / install.packages('TMB', type = 'source')
library(DHARMa)

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# Function that takes a value (or array) in hms format 
# (hour-minute-second, from the package hms or lubridate)
# and turns it into a fraction of a day
hms_to_decimal_day <-function(time_hms){
  h <- hours(time_hms)
  m <- minutes(time_hms)
  s <- seconds(time_hms)
  sum_minutes <- h*60 + m + s/60
  dec <- sum_minutes / (24*60)
  return(dec)
}

# Temperature ----

df1 <- read.csv("../data/Temperatures_Sylvatic_DENV-2_Cynomolgus_Macaques.csv")
df2 <- read.csv("../data/Temperatures_Sylvatic_DENV-2_Squirrel_Monkeys.csv")
df3 <- read.csv("../data/Temperatures_Sylvatic_ZIKV_Squirrel_Monkeys.csv")
df4 <- read.csv("../data/Temperatures_Sylvatic_ZIKV_Cynomolgus_Macaques.csv")

treat1 <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                   dec = ".", sep  ="\t")
treat2 <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                   dec = ".", sep  ="\t")
treat3 <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                   dec = ".", sep  ="\t")
treat4 <- read.csv("~/Documents/POSTDOC/trade_offs_NMSU/data/ZIKV_Sylv_Cyno/master_ZIKV_Cyno.csv",
                   dec = ".", sep  ="\t")

treat1 <- unique(treat1[,c("ID","Final.Treatment")])
treat2 <- unique(treat2[,c("ID","Final.Treatment")])
treat3 <- unique(treat3[,c("ID","Final.Treatment")])
treat4 <- unique(treat4[,c("ID","Final.Treatment")])

df1 <- merge(df1,treat1, by = "ID")
df2 <- merge(df2,treat2, by = "ID")
df3 <- merge(df3,treat3, by = "ID")
df4 <- merge(df4,treat4, by = "ID")

# merge controls for same species experiments
cont1 <- df1[df1$Final.Treatment == "Control",]
cont2 <- df2[df2$Final.Treatment == "Control",]
cont3 <- df3[df3$Final.Treatment == "Control",]

df2 <- rbind(df2,cont3)
df3 <- rbind(df3,cont2)
df4 <- rbind(df4,cont1)

df1$chron_time <- chron(times. = as.character(df1$time),
                        format = "h:m:s")
df2$chron_time <- chron(times. = as.character(df2$time),
                        format = "h:m:s")
df3$chron_time <- chron(times. = as.character(df3$time),
                        format = "h:m:s")
df4$chron_time <- chron(times. = as.character(df4$time),
                        format = "h:m:s")

df1$Study.Day <- df1$Day + hms_to_decimal_day(df1$chron_time)
df2$Study.Day <- df2$Day + hms_to_decimal_day(df2$chron_time)
df3$Study.Day <- df3$Day + hms_to_decimal_day(df3$chron_time)
df4$Study.Day <- df4$Day + hms_to_decimal_day(df4$chron_time)

# for axis limit
# max(c(df1$temp,df2$temp,df3$temp,df4$temp)) # 40.65

# change order for the legends
df1$Final.Treatment <- factor(df1$Final.Treatment, levels = c("Control",
                                                              "1 Mosquito",
                                                              "10 Mosquitos"))
df2$Final.Treatment <- factor(df2$Final.Treatment, levels = c("Control",
                                                              "15 Mosquitos"))
df3$Final.Treatment <- factor(df3$Final.Treatment, levels = c("Control",
                                                              "15 Mosquitos"))
df4$Final.Treatment <- factor(df4$Final.Treatment, levels = c("Control",
                                                              "15 Mosquitos"))

pt1 <- ggplot(df1) + geom_line(aes(x = Study.Day, y = temp,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.6) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "1 Mosquito" = "#c2e699",
                                "10 Mosquitos" = "#006837")) +
  geom_vline(xintercept = 0) +
  labs(x = "", y = expression("Temperature " ( degree*C)),
       color = "") +
  coord_cartesian(xlim = c(-7,28), ylim = c(34,41)) +
  scale_x_continuous(breaks = c(-7,0,5,10,
                                15,20,25,28),
                     labels = c("-7","0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
#pt1

pt2 <- ggplot(df2) + geom_line(aes(x = Study.Day, y = temp,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.6) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "15 Mosquitos" = "#006837")) +
  geom_vline(xintercept = 0) +
  labs(x = "", y = "",
       color = "") +
  coord_cartesian(xlim = c(-7,28), ylim = c(34,41)) +
  scale_x_continuous(breaks = c(-7,0,5,10,
                                15,20,25,28),
                     labels = c("-7","0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
#pt2

pt3 <- ggplot(df3) + geom_line(aes(x = Study.Day, y = temp,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.6) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "15 Mosquitos" = "#2b8cbe")) +
  geom_vline(xintercept = 0) +
  labs(x = "Day post infection", y = "",
       color = "") +
  coord_cartesian(xlim = c(-7,28), ylim = c(34,41)) +
  scale_x_continuous(breaks = c(-7,0,5,10,
                                15,20,25,28),
                     labels = c("-7","0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30,
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
#pt3

pt4 <- ggplot(df4) + geom_line(aes(x = Study.Day, y = temp,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.6) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "15 Mosquitos" = "#2b8cbe")) +
  geom_vline(xintercept = 0) +
  labs(x = "Day post infection", y = expression("Temperature " ( degree*C)),
       color = "") +
  coord_cartesian(xlim = c(-7,28), ylim = c(34,41)) +
  scale_x_continuous(breaks = c(-7,0,5,10,
                                15,20,25,28),
                     labels = c("-7","0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30,
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
# pt4

pt <- (pt1 | pt2) / (pt4 | pt3)
pt <- pt + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 27))
png(filename = "../output/figures/suppl/Figure_S1.png",
    width = 1950, height = 1300)
print(pt)
dev.off()

# Weight ----

df1 <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                dec = ".", sep  ="\t")
df2 <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                dec = ".", sep  ="\t")
df3 <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                dec = ".", sep  ="\t")
df4 <- read.csv("~/Documents/POSTDOC/trade_offs_NMSU/data/ZIKV_Sylv_Cyno/master_ZIKV_Cyno.csv",
                   dec = ".", sep  ="\t")

df1 <- df1[,c("ID","Final.Treatment","Day.Post.Infection","Weight..kg.")]
df2 <- df2[,c("ID","Final.Treatment","Day.Post.Infection","Weight..g.")]
df3 <- df3[,c("ID","Final.Treatment","Day.Post.Infection","Weight..g.")]
df4 <- df4[,c("ID","Final.Treatment","Day.Post.Infection","Weight..g.")]

colnames(df1)[4] <- "weight"
colnames(df2)[4] <- "weight"
colnames(df3)[4] <- "weight"
colnames(df4)[4] <- "weight"

# Everything in grams
df1$weight <- 1000 * df1$weight
df2$weight <- as.numeric(df2$weight)
df3$weight <- as.numeric(df3$weight)
df4$weight <- as.numeric(df4$weight)

df1 <- df1[complete.cases(df1),]
df2 <- df2[complete.cases(df2),]
df3 <- df3[complete.cases(df3),]
df4 <- df4[complete.cases(df4),]

cont1 <- df1[df1$Final.Treatment == "Control",]
cont2 <- df2[df2$Final.Treatment == "Control",]
cont3 <- df3[df3$Final.Treatment == "Control",]

df2 <- rbind(df2,cont3)
df3 <- rbind(df3,cont2)
df4 <- rbind(df4,cont1)

# df1 also has measurements at day -20, but we'll take all baseline at day 0
baseline1 <- df1[df1$Day.Post.Infection == 0,]
baseline2 <- df2[df2$Day.Post.Infection == 0,]
baseline3 <- df3[df3$Day.Post.Infection == 0,]
baseline4 <- df4[df4$Day.Post.Infection == 0,]

baseline1 <- baseline1[,c("ID","weight")]
baseline2 <- baseline2[,c("ID","weight")]
baseline3 <- baseline3[,c("ID","weight")]
baseline4 <- baseline4[,c("ID","weight")]

colnames(baseline1)[2] <- "baseline"
colnames(baseline2)[2] <- "baseline"
colnames(baseline3)[2] <- "baseline"
colnames(baseline4)[2] <- "baseline"

df1 <- merge(df1,baseline1, by = "ID")
df2 <- merge(df2,baseline2, by = "ID")
df3 <- merge(df3,baseline3, by = "ID")
df4 <- merge(df4,baseline4, by = "ID")

# change order for the legends
df1$Final.Treatment <- factor(df1$Final.Treatment, levels = c("Control",
                                                              "1 Mosquito",
                                                              "10 Mosquitos"))
df2$Final.Treatment <- factor(df2$Final.Treatment, levels = c("Control",
                                                              "15 Mosquitos"))
df3$Final.Treatment <- factor(df3$Final.Treatment, levels = c("Control",
                                                              "15 Mosquitos"))
df4$Final.Treatment <- factor(df4$Final.Treatment, levels = c("Control",
                                                              "15 Mosquitos"))
df1$species <- "Cyno"
df2$species <- "Squirrel"
df3$species <- "Squirrel"
df4$species <- "Cyno"
df <- rbind(df1,df2,df3,df4)

# Check species difference
m0 <- glmmTMB(weight ~ species,
              data = df)
simulateResiduals(m0, plot = T) # issues
m0 <- glmmTMB(weight ~ species + (1|ID) + (1|Day.Post.Infection),
              data = df)
simulateResiduals(m0, plot = T) # issues
summary(m0)

# Checking that weight change is not different between treatments
m1 <- glmmTMB(weight/baseline ~ Final.Treatment,
              data = df1)
simulateResiduals(m1, plot = T) # ok
summary(m1)

m2 <- glmmTMB(weight/baseline ~ Final.Treatment,
              data = df2)
simulateResiduals(m2, plot = T) # ok
summary(m2)

m3 <- glmmTMB(weight/baseline ~ Final.Treatment,
              data = df3)
simulateResiduals(m3, plot = T) # ok
summary(m3)

m4 <- glmmTMB(weight/baseline ~ Final.Treatment,
              data = df4)
simulateResiduals(m4, plot = T) # Levene test signif
summary(m4)

# Figure
pw1 <- ggplot(df1) + geom_line(aes(x = Day.Post.Infection, y = weight/baseline,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.8) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "1 Mosquito" = "#c2e699",
                                "10 Mosquitos" = "#006837")) +
  labs(x = "", y = "Relative weight",
       color = "") +
  scale_x_continuous(limits = c(0,28),
                     breaks = c(0,5,10,
                                15,20,25,28),
                     labels = c("0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0.8,1.2)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.title.x = element_text(size = 30),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
#pw1

pw2 <- ggplot(df2) + geom_line(aes(x = Day.Post.Infection, y = weight/baseline,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.8) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "15 Mosquitos" = "#006837")) +
  labs(x = "", y = "",
       color = "") +
  scale_x_continuous(limits = c(0,28),
                     breaks = c(0,5,10,
                                15,20,25,28),
                     labels = c("0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0.8,1.2)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.title.x = element_text(size = 30),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
#pw2

pw3 <- ggplot(df3) + geom_line(aes(x = Day.Post.Infection, y = weight/baseline,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.8) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "15 Mosquitos" = "#2b8cbe")) +
  labs(x = "Day post infection", y = "",
       color = "") +
  scale_x_continuous(limits = c(0,28),
                     breaks = c(0,5,10,
                                15,20,25,28),
                     labels = c("0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0.8,1.2)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.title.x = element_text(size = 30,
                                    margin = margin(t = 15)),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
#pw3

pw4 <- ggplot(df4) + geom_line(aes(x = Day.Post.Infection, y = weight/baseline,
                                   color = Final.Treatment, group = ID),
                               alpha = 0.8) +
  scale_color_manual(values = c("Control" = "darkgrey",
                                "15 Mosquitos" = "#2b8cbe")) +
  labs(x = "Day post infection", y = "Relative weight",
       color = "") +
  scale_x_continuous(limits = c(0,28),
                     breaks = c(0,5,10,
                                15,20,25,28),
                     labels = c("0","","10",
                                "","20","","28"),
                     expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0.8,1.2)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 30,
                                    margin = margin(r = 15)),
        axis.title.x = element_text(size = 30,
                                    margin = margin(t = 15)),
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom")
#pw4

pw <- (pw1 | pw2) / (pw4 | pw3)
pw <- pw + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 27))
png(filename = "../output/figures/suppl/Figure_S2.png",
    width = 1950, height = 1300)
print(pw)
dev.off()
