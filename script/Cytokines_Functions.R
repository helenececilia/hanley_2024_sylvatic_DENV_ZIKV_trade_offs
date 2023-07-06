## ---------------------------
##
## Script name:
##
## Purpose of script: Functions to ease the analysis of cytokine data
##
## Author: Helene Cecilia
##
## Date Created: 2023-05-31

sex_diff_at_baseline <- function(data, exclude_LOD = T,
                                 verbose = F, plot = F, filename = NULL){
  exclude <- NULL
  var_het <- NULL
  sex_diff <- NULL
  no_sex_diff <- NULL
  pvalues <- NULL
  
  if(exclude_LOD == TRUE){
    data <- data[data$value != data$LOD,]
  }
  
  resname <- paste0(filename,".txt")
  for(c in unique(data$cytokine)){
    df <- data[data$cytokine == c & data$day == -7,]
    df_m <- df[df$sex == "M",]
    df_f <- df[df$sex == "F",]
    all_obs_m <- dim(df_m)[1]
    all_obs_f <- dim(df_f)[1]
    LOD_female <- sum(df_f$value == df_f$LOD)
    LOD_male <- sum(df_m$value == df_m$LOD)

    if(verbose & !is.null(filename)){
      write(c, file = resname, append = T)
      write(paste0(all_obs_m," male obs"), file = resname, append = T)
      write(paste0(all_obs_f," female obs"), file = resname, append = T)
      write(paste0(LOD_male," LOD male"), file = resname, append = T)
      write(paste0(LOD_female," LOD female"), file = resname, append = T)
    }
    tryCatch({
      # Wilcoxon test even if Levene's test significant 
      lev <- leveneTest(log10(value) ~ sex, data = df)
      if(verbose & !is.null(filename)){
        res1 <- capture.output(lev)
        write(res1, file = resname, append = T)
      }
      wil <- wilcox.test(log10(df_m$value), log10(df_f$value))
      if(verbose & !is.null(filename)){
        res2 <- capture.output(wil)
        write(res2, file = resname, append = T)
      }
      if(lev$`Pr(>F)`[1] <= 0.05){
        var_het <- c(var_het,c)
      }
      if(wil$p.value <= 0.05){
        sex_diff <- c(sex_diff,c)
        pvalues <- c(pvalues, wil$p.value)
      }else{
        no_sex_diff <- c(no_sex_diff,c)
      }
    }, error=function(e){cat(c, "ERROR :",conditionMessage(e), "\n")}) #
    if(plot){
      if(c %in% sex_diff){
        plotname <- paste0(filename,"_",c,".png")
        p <- ggplot(df, aes(x = sex, y = log10(value), color = sex)) +
          geom_boxplot() + geom_jitter() +
          geom_hline(aes(yintercept = log10(LOD))) +
          scale_color_viridis_d() +
          ggtitle(c) +
          labs(y = "Concentration (log10 pg/\u03BCl)", x = "", color = "") +
          theme_bw() +
          theme(axis.title = element_text(size = 25),
                axis.text = element_text(size = 24),
                legend.text = element_text(size = 23),
                plot.title = element_text(size = 24),
                legend.position = "none",
                panel.grid.minor = element_blank()) 
        
        png(filename = plotname, width = 450, height = 700)
        plot(p)
        dev.off()
        
      }
    }
  }
  return(list("exclude" = exclude, "var_het" = var_het,
              "sex_diff" = sex_diff, "no_sex_diff" = no_sex_diff, "pvalues" = pvalues))
}

compare_baseline_day28 <- function(data, exclude_LOD = T,
                                   verbose = F, plot = F, filename = NULL){
  keep_day28 <- NULL
  pvalues <- NULL
  remove_day28 <- NULL
  exclude <- NULL
  
  if(exclude_LOD == TRUE){
    data <- data[data$value != data$LOD,]
  }
  data_c <- data
  data <- data[data$group != "Control",]
  
  day28 <- data[data$day == 28,]
  baseline <- data[data$day == -7,]
  resname <- paste0(filename,".txt")
  for(c in unique(baseline$cytokine)){
    tryCatch({
      df1 <- baseline[which(baseline$cytokine == c),]
      df2 <- day28[which(day28$cytokine == c),]
      
      # For squirrels, not all individuals appear on day 28, so keep the right subset for baseline
      # Also if an ID is not present at one day (because LOD) then should not be present the other day
      IDs <- intersect(unique(df1$ID), unique(df2$ID))
      df1 <- df1[which(df1$ID %in% IDs),]
      df2 <- df2[which(df2$ID %in% IDs),]
      
      df1 <- df1[with(df1, order(ID)),]
      df2 <- df2[with(df2, order(ID)),]
      df <- rbind(df1,df2)
      
      all_obs1 <- dim(df1)[1]
      all_obs2 <- dim(df2)[1]
      LOD_baseline <- sum(df1$value == df1$LOD)
      LOD_day28 <- sum(df2$value == df2$LOD)

      if(verbose & !is.null(filename)){
        write(c, file = resname, append = T)
        write(paste0(all_obs1," obs in baseline"), file = resname, append = T)
        write(paste0(all_obs2," obs in day 28"), file = resname, append = T)
        write(paste0(LOD_baseline," LOD baseline"), file = resname, append = T)
        write(paste0(LOD_day28," LOD day 28"), file = resname, append = T)
      }
      
      test_value1 <- wilcox.test(log10(value) ~ day, data = df, paired = T, alternative = "two.sided")

      if(verbose & !is.null(filename)){
        res1 <- capture.output(test_value1)
        write(res1, file = resname, append = T)
      }
      
      if(test_value1$p.value <= 0.05){ 
        keep_day28 <- c(keep_day28,c)
        pvalues <- c(pvalues, test_value1$p.value)
      }else if(test_value1$p.value > 0.05){ 
        remove_day28 <- c(remove_day28,c)
      }else if(is.na(test_value1$p.value)){ 
        # I think this one is useless, even though the summary of the model shows an NA p-value, it cannot be extracted
        print(c)
        print("NA p-value, put in exclude")
        exclude <- c(exclude,c)
      }
      
    }, error=function(e){}) 
    if(plot){
      if(c %in% keep_day28){
        if(unique(data$virus) == "Dengue virus"){
          if(unique(data$NHP) == "Cyno"){
            my_colors <- c("Control" = "darkgrey",
                           "1 Mosquito" = "#c2e699",
                           "10 Mosquitos" = "#006837")
          }else if(unique(data$NHP) == "Squirrel"){
            my_colors <-c("Control" = "darkgrey",
                          "15 Mosquitos" = "#006837")
          }
        }else if(unique(data$virus) == "Zika virus"){
          my_colors <-c("Control" = "darkgrey",
                        "15 Mosquitos" = "#2b8cbe")
        }
        plotname <- paste0(filename,"_",c,".png")
        test <- data_c[data_c$day %notin% seq(0,27) & data_c$cytokine == c,]
        # browser()
        p1 <- ggplot(test) + geom_line(aes(x = day, y = log10(value), group = ID, color = group)) +
          geom_point(aes(x = day, y = log10(value), color = group)) +
          geom_hline(aes(yintercept = log10(LOD))) + ggtitle(c) +
          scale_color_manual(values = my_colors) +
          scale_x_continuous(limits = c(-7,28),
                             breaks = c(-7,28),
                             labels = c("-7","28")) +
          labs(y = "Concentration (log10 pg/\u03BCl)", x = "Day", color = "") +
          theme_bw() +
          theme(axis.title = element_text(size = 25),
                axis.text = element_text(size = 24),
                legend.text = element_text(size = 23),
                plot.title = element_text(size = 24),
                legend.position = "bottom",
                panel.grid.minor = element_blank(),
                plot.margin = unit(c(1,2,1,1),"cm"))
        
        png(filename = plotname, width = 450, height = 700)
        plot(p1)
        dev.off()
      }
    }
  }
  return(list("signif_diff" = keep_day28, "no_diff" = remove_day28, "exclude" = exclude, "pvalues" = pvalues))
}


