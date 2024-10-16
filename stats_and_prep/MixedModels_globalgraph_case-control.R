## C. Vriend - Amsterdam UMC - July '24
## perform multivariate mixed model analysis on CORE data


# clear variables
rm(list = ls())

#library(tidyverse)
library(dplyr)
library(readr)
library(moderndive)
library(reshape2)
library(lme4)
library(readxl)
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
library(tidyr)
library(stringr)
library(openxlsx)
library(ggplot2)
library(ggdist)

# functions #
clean_dataframe <- function(df) {
  # multiply by 100 and round to two decimals
  df[c('B',
       'B.1',
       '95CI-',
       '95CI+',
       '95CI-.1',
       '95CI+.1',
       'SE',
       'SE.1')] <- df[c('B',
                        'B.1',
                        '95CI-',
                        '95CI+',
                        '95CI-.1',
                        '95CI+.1',
                        'SE',
                        'SE.1')] * 100
  df[c('95CI-', '95CI+', '95CI-.1', '95CI+.1', 'SE', 'SE.1')] <- sapply(df[c('95CI-', '95CI+', '95CI-.1', '95CI+.1', 'SE', 'SE.1')], round, 2)
  # unite columns for table
  df <- df %>%
    unite("CI_crude",
          c("95CI-", "95CI+"),
          sep = " | ",
          remove = TRUE) %>% unite("CI_adj",
                                   c("95CI-.1", "95CI+.1"),
                                   sep = " | ",
                                   remove = TRUE)
  df$BSE_crude <- sprintf("%2.3f [%2.3f]", df$B, df$SE)
  df$BSE_adj <- sprintf("%2.3f [%2.3f]", df$B.1, df$SE.1)
  
  df <- select(df, -c(B.1, B, SE, SE.1))
  
  # calculate FDR correction ##
  
  #df['Pfdr_crude']<-p.adjust(df[['P-value']], method = "BH")
  #df['Pfdr_adj']<-p.adjust(df[['P-value.1']], method = "BH")
  
  # round to 3 decs
  #df[c('Pfdr_crude','Pfdr_adj')] <- sapply(df[c('Pfdr_crude','Pfdr_adj')],round,3)
  
  #df$Pvalues_crude <- sprintf("%2.3f (%2.3f)", df$Pfdr_crude, df[['P-value']])
  #df$Pvalues_adj <- sprintf("%2.3f (%2.3f)", df$Pfdr_adj, df[['P-value.1']])
  df$Pvalues_crude <- sprintf("%2.3f", df[['P-value']])
  df$Pvalues_adj <- sprintf("%2.3f", df[['P-value.1']])
  # remove cols
  df <- select(df, -c('P-value', 'P-value.1'))
  
  
  df <- df %>% relocate("Pvalues_crude", .after = "BSE_crude") %>%
    relocate("CI_crude", .before = "Pvalues_crude") %>%
    relocate("Pvalues_adj", .after = "BSE_adj") %>%
    relocate("CI_adj", .before = "Pvalues_adj")
  return(df)
}

# input variables 
date <- Sys.Date()
basedir = '/data/anw/anw-work/NP/projects/data_chris/CORE'
session = 'ses-T0'
dptvars <- c('Dx')


for (acq in c('multi', 'func', 'dwi')) {
  print(acq)
  if (acq == 'multi') {
    
    graphfile<-file.path(basedir,acq,session,'output',paste0('CORE_ML_',session,'_global.csv'))
    graphfileHC<-file.path(basedir,acq,'HC','output',paste0('CORE_ML_','HC','_global.csv'))
    graphmeasures<-c('eigenvector','eccentricity')
    
  }  else  {
   
  graphfile <- file.path(basedir,
                         acq,
                         session,
                         'graph',
                         paste0('graph_global_acq-', acq, '_CORE.csv'))
  graphfileHC <- file.path(basedir,
                           acq,
                           'HC',
                           'graph',
                           paste0('graph_global_acq-', acq, '_CORE.csv'))
  
  graphmeasures = c('GE', 'Qyeo', 'PC_glob_Schaef', 'Sigma','Q','PC_glob_ind')
  # only to check c('Q','PC_glob_ind')
  }
  
  clinfile <- file.path(basedir, 'clin', 'OCDPTSD_trials4analyses.xlsx')
  clinfileHC <- file.path(basedir, 'clin', 'HC4analyses.xlsx')
  
  outputfile <- file.path(
    basedir,
    'stats',
    paste0(
      "CORE_MixedModel_global_acq-",
      acq,
      "_case-control_",
      date,
      ".xlsx"
    )
  )
  outputfile2 <- file.path(
    basedir,
    'stats',
    paste0(
      "CORE_LOSO_global_acq-",
      acq,
      "_",
      date,
      "_case-control.xlsx"
    )
  )
  summaryfile <- file.path(
    basedir,
    'stats',
    paste0("summaryTables_acq-", acq, "_", date, "_case-control.xlsx")
  )
  

  # Import datasets # 
  dfwide <- read_csv(graphfile)
  dfclin <- read_excel(clinfile)
  dfwideHC <- read_csv(graphfileHC)
  dfclinHC <- read_excel(clinfileHC)
  

  dfclinHC['subjID'] <- dfclinHC$Subj
  
  dfclin <- dfclin %>%
    mutate(subjID = paste0(Subj, '_', session))
  dfwideHC <- dfwideHC %>%
    mutate(subjID = str_replace(subjID, "_ses-T0$", ""))
  
  # drop column T1MRI
  dfclin <- select(dfclin, -c(T1MRI))
  
  # merge clinical and connectome df 
  dfwide <- merge(dfwide, dfclin, by = 'subjID')
  dfwideHC <- merge(dfwideHC, dfclinHC, by = 'subjID')
  
  
  # additional exclusions >65 year
  #dfwideHC <-dfwideHC %>% filter(age<=65)
  
  
  # select columns
  dfwide <- select(
    dfwide,
    -c(
      Subj,
      baseline_sev,
      posttrial_sev,
      dYBOCS_dCAPS,
      responder,
      perc_improv
    )
  )
  # add diagnostic variable 
  dfwide['Dx'] = 'pat'
  dfwideHC['Dx'] = 'HC'
  
  dfclinsubset <- dfclin[c('Subj', 'age', 'seks')]
  dfclinsubset['Dx'] = 'pat'
  dfclinHCsubset <- dfclinHC[c('Subj', 'age', 'seks')]
  dfclinHCsubset['Dx'] <- 'HC'
  
  dfclin2 <- rbind(dfclinsubset, dfclinHCsubset)
  dfclin2[sapply(dfclin2, is.character)] <- lapply(dfclin2[sapply(dfclin2, is.character)], as.factor)
  t.test(dfclin2$age ~ dfclin2$Dx)
  chisq.test(dfclin2$Dx, dfclin2$seks)
  
  
  ## perform match-it
  
  # Perform matching based on age and sex
  library(MatchIt)
  match_it <- matchit(
    Dx ~ age + seks,
    data = dfclin2,
    method = "nearest",
    ratio = 4,
    replace = TRUE
  )
  
  # Get the matched data
  matched_data <- match.data(match_it)
  
  # Split the matched data back into patients and controls
  #matched_patients <- matched_data[matched_data$Dx == 'pat', ]
  matched_controls <- matched_data[matched_data$Dx == 'HC', ]
  
  matched_controls['Dx'] = 'HC'
  matched_controls <- matched_controls %>%
    select(-c(distance, weights))
  
  dfclin3 <- rbind(dfclinsubset, matched_controls)
  dfclin3[sapply(dfclin3, is.character)] <- lapply(dfclin3[sapply(dfclin3, is.character)], as.factor)
  
  #t.test(dfclin3$age ~ dfclin3$Dx)
  #chisq.test(dfclin3$Dx, dfclin3$seks)
  ########
  
  Subjbefore <- unique(dfclin2$Subj)
  Subjafter <- unique(dfclin3$Subj)
  Subjdiff <- setdiff(Subjbefore, Subjafter)
  
  # exclude HC's that are not a match on age/sex with patients
  dfwideHC_filt <- dfwideHC[!dfwideHC$subjID %in% Subjdiff, ]
  
  dfwide <- bind_rows(dfwide, dfwideHC_filt)
  dfwide <- select(dfwide, -c(Subj))
  
  # convert characters to factors
  dfwide <-dfwide %>% mutate(across(c(trial,treatment), ~as.character(.x))) %>%
    mutate(across(c(trial,treatment), ~replace_na(.x,'HC'))) %>%
    mutate(across(c(trial,treatment), ~as.factor(.x)))
  
  if (acq == 'dwi') {
    dfwide <- dfwide %>% mutate_at(vars('GE'), ~ (. * 100))
    
  } else if (acq == 'func') {
    dfwide <- dfwide %>% mutate_at(vars('GE'), ~ (. * 10))
    
  }
  
  
  
  
  ##################
  ### create plot ##
  ##################
  
  
  
  dfwide4plot<-merge(dfwide, dfclin[, c('subjID','perc_improv','responder')], by = 'subjID',all.x=TRUE)
  dfwide4plot$responder <- factor(dfwide4plot$responder, levels = c(0, 1), labels = c("non-resp", "resp"))
  
    for (i in 1:length(graphmeasures)) {
      graph = graphmeasures[i]
  
      yintercept_min <- 0.95 * min(dfwide4plot[[graph]], na.rm = TRUE)
      yintercept_max <- 1.05 * max(dfwide4plot[[graph]], na.rm = TRUE)
  # Create the scatterplot and boxplot
  ggplot() +
    # Scatterplot for Group A (filtering out NA values for X)
    geom_point(data = filter(dfwide4plot, Dx == "pat"), aes(x = perc_improv, y = !!sym(graph), color = responder), alpha = 0.7) +
    geom_smooth(data = filter(dfwide4plot, Dx == "pat"), aes(x = perc_improv, y = !!sym(graph)), method = "lm", color = "black", linetype = "dashed") +
    
    # Raincloud plot for Group B (positioned at a fixed x position)
    geom_violin(data = filter(dfwide4plot, Dx == "HC"), aes(x = -1.5, y = !!sym(graph)), 
                     width = 0.3, fill = "orange", alpha = 0.5) +
    geom_jitter(data = filter(dfwide4plot, Dx == "HC"), aes(x = -1.5, y = !!sym(graph)), 
                width = 0.05, alpha = 0.5) +
    stat_summary(data = filter(dfwide4plot, Dx == "HC"), aes(x = -1.5, y = !!sym(graph)), 
                 fun = mean, geom = "point", shape = 95, size = 6, color = "black") +
        # Boxplot for Group B (position it at a fixed x position)
    geom_boxplot(data = filter(dfwide4plot, Dx == "HC"), aes(x = -1.5, y = !!sym(graph)), width = 0.2, fill = "orange4", color = "black", alpha = 0.5) +
    
    
    # Raincloud plot for Group B (positioned at a fixed x position)
    geom_violin(data = filter(dfwide4plot, Dx == "pat"), aes(x = -1, y = !!sym(graph)), 
                width = 0.3, fill = "green", alpha = 0.5) +
    geom_jitter(data = filter(dfwide4plot, Dx == "pat"), aes(x = -1, y = !!sym(graph)), 
                width = 0.05, alpha = 0.5) +
    stat_summary(data = filter(dfwide4plot, Dx == "pat"), aes(x = -1, y = !!sym(graph)), 
                 fun = mean, geom = "point", shape = 95, size = 6, color = "black") +
    # Boxplot for Group B (position it at a fixed x position)
    geom_boxplot(data = filter(dfwide4plot, Dx == "pat"), aes(x = -1, y = !!sym(graph)), width = 0.2, fill = "green4", color = "black", alpha = 0.5) +
    
    
   # geom_hline(yintercept = yintercept_min, linetype = "dotted", color = "grey") +
  #  annotate("text", x = 2, y = yintercept_min, label = "Break", color = "grey", size = 4, angle = 0) +
    
    
   # geom_segment(aes(x = -2, xend = -1.8, y = yintercept_min, yend = yintercept_min), 
   #              linetype = "solid", color = "black", size = 0.5) +  
    # Adjust size as needed
    geom_vline(xintercept = -0.75, linetype = "solid", color ="black", size = 0.5) +
    # Add 0 where the x and y axes intersect (doesn't work)
    #annotate("text", x = 0, y = 0, label = "0", color = "black", size = 4, vjust = 1, hjust = 1)

    # Adjust the x-axis limits to make room for the boxplot and the scatterplot
    
    scale_x_continuous(limits = c(-1.75, 1.5), 
                       breaks = c(-1.5, -1, -0.5, 0, 0.5, 1), 
                       labels = c("healthy controls", "patients", "-0.5", "0", "0.5", "1")) +
  #  theme(axis.text.x = element_text(angle = 45,hjust =1)) +
    #scale_x_continuous(limits = c(-2, 2), breaks = c(seq(-1.5, 1, by = 0.5), -1.5, 1.5), labels = c(seq(-1.5, 1, by = 0.5), "healthy controls","patients")) +
    # Set the y-axis limits to 0 to 1
    scale_y_continuous(limits = c(yintercept_min, yintercept_max)) +
    
    # Labels and theme adjustments
    labs(x = "% improvement", y = graph) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 10))
  ggsave(file.path(basedir,"stats",paste0('plot_',acq,'_',graph,'_Dx-vs-percimprov.pdf')),width=8,height=6)
  

    }
  ##################
  ### run models ###
  ##################
  
  # Lists to store resulting dataframes
  results_list <- list()
  results_list_treat <- list()
  results_list_trial <- list()
  
  
  for (j in 1:length(dptvars)) {
    covinterest <- dptvars[j]
    print(paste("dependent variable = ", covinterest))
    
    # pre-allocation
    temps_graph_coeff <- matrix(
      0,
      nrow = length(graphmeasures),
      ncol = 10,
      dimnames = list(graphmeasures, rep(
        c("B", "SE", "P-value", "95CI-", "95CI+"), 2
      ))
    )
    
    temps_graph_coeff_trial <- matrix(
      0,
      nrow = length(graphmeasures),
      ncol = 10,
      dimnames = list(graphmeasures, rep(
        c("B", "SE", "P-value", "95CI-", "95CI+"), 2
      ))
    )
    temps_graph_coeff_treat <- matrix(
      0,
      nrow = length(graphmeasures),
      ncol = 10,
      dimnames = list(graphmeasures, rep(
        c("B", "SE", "P-value", "95CI-", "95CI+"), 2
      ))
    )
    
    
    
    
    for (i in 1:length(graphmeasures)) {
      graph = graphmeasures[i]
      
      print(paste("global graph measure = ", graph))
      
      # crude model
      model_graph_crude <- lm(formula = as.formula(paste(graph, "~ ", covinterest)), data =
                                dfwide)
      # adj model
      model_graph_adj <- lm(formula = as.formula(paste(
        graph, "~ ", covinterest, " + age + seks"
      )), data = dfwide)
      
      # rand intercept - trial
      model_graph_crude_trial <- lmer(formula = as.formula(paste(
        graph, "~ ", covinterest, " + (1|trial)"
      )), data = dfwide)
      model_graph_adj_trial <- lmer(formula = as.formula(paste(
        graph, "~ ", covinterest, " + age + seks + (1|trial)"
      )), data = dfwide)
      
      # rand intercept - treatment
      model_graph_crude_treat <- lmer(formula = as.formula(paste(
        graph, "~ ", covinterest, " + (1|treatment)"
      )), data = dfwide)
      model_graph_adj_treat <- lmer(formula = as.formula(
        paste(graph, "~ ", covinterest, " + age + seks + (1|treatment)")
      ), data = dfwide)
      
      
      
      # extract stat parameters
      
      # standard models
      rows_crude <- grep(paste0("^", covinterest),
                         rownames(summary(model_graph_crude)$coefficients),
                         value = TRUE)
      rows_adj <- grep(paste0("^", covinterest),
                       rownames(summary(model_graph_adj)$coefficients),
                       value = TRUE)
      
      coeff_graph_crude <- cbind(summary(model_graph_crude)$coefficients[rows_crude, c(1, 2, 4), drop =
                                                                           FALSE],
                                 confint(model_graph_crude)[rows_crude, , drop =
                                                              FALSE])
      coeff_graph_adj <- cbind(summary(model_graph_adj)$coefficients[rows_adj, c(1, 2, 4), drop =
                                                                       FALSE],
                               confint(model_graph_adj)[rows_adj, , drop = FALSE])
      
      # random intercept -trial models
      rows_crude <- grep(paste0("^", covinterest),
                         rownames(summary(model_graph_crude_trial)$coefficients),
                         value = TRUE)
      rows_adj <- grep(paste0("^", covinterest),
                       rownames(summary(model_graph_adj_trial)$coefficients),
                       value = TRUE)
      
      # in case of random intercept the P-value is in column 5 instead of 4
      if (colnames(summary(model_graph_crude_trial)$coefficients)[3] == 'df') {
        Pcol <- 5
      } else {
        Pcol <- 4
      }
      
      coeff_graph_crude_trial <- cbind(
        summary(model_graph_crude_trial)$coefficients[rows_crude, c(1, 2, Pcol), drop =
                                                        FALSE],
        confint(model_graph_crude_trial)[rows_crude, , drop =
                                           FALSE]
      )
      coeff_graph_adj_trial <- cbind(
        summary(model_graph_adj_trial)$coefficients[rows_adj, c(1, 2, Pcol), drop =
                                                      FALSE],
        confint(model_graph_adj_trial)[rows_adj, , drop =
                                         FALSE]
      )
      rm(Pcol)
      
      # random intercept -treatment models
      rows_crude <- grep(paste0("^", covinterest),
                         rownames(summary(model_graph_crude_treat)$coefficients),
                         value = TRUE)
      rows_adj <- grep(paste0("^", covinterest),
                       rownames(summary(model_graph_adj_treat)$coefficients),
                       value = TRUE)
      
      # in case of random intercept the P-value is in column 5 instead of 4
      if (colnames(summary(model_graph_crude_treat)$coefficients)[3] == 'df') {
        Pcol <- 5
      } else {
        Pcol <- 4
      }
      
      coeff_graph_crude_treat <- cbind(
        summary(model_graph_crude_treat)$coefficients[rows_crude, c(1, 2, Pcol), drop =
                                                        FALSE],
        confint(model_graph_crude_treat)[rows_crude, , drop =
                                           FALSE]
      )
      coeff_graph_adj_treat <- cbind(
        summary(model_graph_adj_treat)$coefficients[rows_adj, c(1, 2, Pcol), drop =
                                                      FALSE],
        confint(model_graph_adj_treat)[rows_adj, , drop =
                                         FALSE]
      )
      rm(Pcol)
      
      
      # round to 3 decimals
      coeff_graph_crude <- sapply(coeff_graph_crude, round, 5)
      coeff_graph_adj <- sapply(coeff_graph_adj, round, 5)
      coeff_graph_crude_trial <- sapply(coeff_graph_crude_trial, round, 3)
      coeff_graph_adj_trial <- sapply(coeff_graph_adj_trial, round, 3)
      coeff_graph_crude_treat <- sapply(coeff_graph_crude_treat, round, 3)
      coeff_graph_adj_treat <- sapply(coeff_graph_adj_treat, round, 3)
      
      
      temps_graph_coeff[i, 1:5] = coeff_graph_crude
      temps_graph_coeff[i, 6:10] = coeff_graph_adj
      
      rm(coeff_graph_crude, coeff_graph_adj)
      temps_graph_coeff_trial[i, 1:5] = coeff_graph_crude_trial
      temps_graph_coeff_trial[i, 6:10] = coeff_graph_adj_trial
      
      temps_graph_coeff_treat[i, 1:5] = coeff_graph_crude_treat
      temps_graph_coeff_treat[i, 6:10] = coeff_graph_adj_treat
      
      rm(
        coeff_graph_crude,
        coeff_graph_crude_trial,
        coeff_graph_crude_treat,
        coeff_graph_adj,
        coeff_graph_adj_trial,
        coeff_graph_adj_treat
      )
      
      
    }
    # save to dataframes
    all_coeff <- as.data.frame(temps_graph_coeff)
    all_coeff_trial <- as.data.frame(temps_graph_coeff_trial)
    all_coeff_treat <- as.data.frame(temps_graph_coeff_treat)
    names(all_coeff) <- make.unique(names(all_coeff))
    names(all_coeff_trial) <- make.unique(names(all_coeff_trial))
    names(all_coeff_treat) <- make.unique(names(all_coeff_treat))
    
    df_globalall <- tibble::rownames_to_column(all_coeff, "globalmeasure")
    df_globalall_trial <- tibble::rownames_to_column(all_coeff_trial, "globalmeasure")
    df_globalall_treat <- tibble::rownames_to_column(all_coeff_treat, "globalmeasure")
    
    
    results_list[[paste0("df_", j)]] <- clean_dataframe(df_globalall)
    results_list_trial[[paste0("df_", j)]] <- clean_dataframe(df_globalall_trial)
    results_list_treat[[paste0("df_", j)]] <- clean_dataframe(df_globalall_treat)
    
  
  }
  
  # Rename and save the dataframes outside the for loop
  df_dx <- results_list$df_1
  df_dx_trial <- results_list_trial$df_1
  df_dx_treat <- results_list_treat$df_1
  ### save and format output ###
  
  write.list <- list("Dx" = df_dx,
                     "Dx_trial" = df_dx_trial,
                     "DX_treat" = df_dx_treat)
  
  wb <- createWorkbook()
  
  arialStyle <- createStyle(fontName = "Arial", fontSize = 8)
  
  for (sheet_name in names(write.list)) {
    addWorksheet(wb, sheet_name)
    
    writeData(wb, sheet_name, write.list[[sheet_name]])
    n_rows <- nrow(write.list[[sheet_name]])
    n_cols <- ncol(write.list[[sheet_name]])
    
    addStyle(
      wb,
      sheet = sheet_name,
      style = arialStyle,
      rows = 1:(n_rows + 1),
      cols = 1:n_cols,
      gridExpand = TRUE
    )
  }
  saveWorkbook(wb, file = outputfile, overwrite = TRUE)

  ##########################################
  ## leave-one-trial-out cross-validation ##
  ##########################################
  
  # Function to calculate the Cauchy Combination Test (CCT) p-value
  cauchy_combination <- function(p_values) {
    # Convert p-values to Z-scores
    z_scores <- qnorm(1 - p_values)
    
    # Calculate the combined Z-score
    combined_z <- sum(z_scores) / sqrt(length(p_values))
    
    # Convert combined Z-score to p-value
    combined_p_value <- 1 - pnorm(combined_z)
    
    return(combined_p_value)
  }
  
  # Function to compute MinP-CCT-MinP (MCM)
  mcm_combination <- function(p_values) {
    # Ensure p-values are between 0 and 1
    if (any(p_values <= 0 | p_values >= 1)) {
      stop("p-values should be in the interval (0, 1).")
    }
    
    # Compute the MinP p-value (smallest p-value)
    min_p_value <- min(p_values)
    
    # Compute the Cauchy Combination Test (CCT) p-value
    cct_p_value <- cauchy_combination(p_values)
    
    # Final p-value is the minimum of MinP and CCT p-values
    final_p_value <- min(min_p_value, cct_p_value)
    
    return(final_p_value)
  }
  
  
  
  # Function to compute the Harmonic Mean P-value
  harmonic_mean_pvalue <- function(p_values) {
    # Ensure p-values are between 0 and 1
    if (any(p_values <= 0 | p_values >= 1)) {
      stop("p-values should be in the interval (0, 1).")
    }
    
    # Calculate the harmonic mean of the p-values
    k <- length(p_values)
    harmonic_mean <- k / sum(1 / p_values)
    
    # Convert the harmonic mean to a p-value
    # Note: The harmonic mean p-value does not have a standard distribution
    # Therefore, we'll use the chi-squared approximation
    chi_squared_statistic <- -2 * sum(log(p_values))
    df <- 2 * k
    combined_p_value <- pchisq(chi_squared_statistic, df, lower.tail = FALSE)
    
    return(combined_p_value)
  }
  

  results_list <- list()
  
  for (LOSO in c('trial', 'treatment')) {
    print(LOSO)
    
    if (LOSO == 'trial') {
      trials <- unique(dfwide$trial[dfwide$trial!='HC'])
    } else if (LOSO == 'treatment') {
      trials <- unique(dfwide$treatment[dfwide$trial!='HC'])
    }
    
   
    for (j in 1:length(dptvars)) {
      covinterest <- dptvars[j]
      print(paste("dependent variable = ", covinterest))
      
      # pre-allocation
      
      mse_mat <- matrix(
        0,
        nrow = length(graphmeasures),
        ncol = 1,
        dimnames = list(graphmeasures, c("MSE"))
      )
      rrmse_mat <- matrix(
        0,
        nrow = length(graphmeasures),
        ncol = 1,
        dimnames = list(graphmeasures, c("RRMSE"))
      )
      pvalues_mat <- matrix(
        0,
        nrow = length(graphmeasures),
        ncol = 10,
        dimnames = list(
          graphmeasures,
          c(
            "min_P",
            "max_P",
            "med_P",
            'mean_P',
            'SD_P',
            'Fisher-P',
            'Stoufer-P',
            'weighted-z-P',
            'MCM-P',
            'harmonic-P'
          )
        )
      )
      
      for (i in 1:length(graphmeasures)) {
        graph = graphmeasures[i]
        
        print(paste("global graph measure = ", graph))
        
        
        predictions <- c()
        actuals <- c()
        pvalues <- c()
        samplesize <- c()
        for (trial_id in trials) {
          if (LOSO == 'trial') {
            training_data <- dfwide %>%
              filter(trial != trial_id)
            test_data <- dfwide %>%
              filter(trial == trial_id)
            
            
          } else if (LOSO == 'treatment') {
            training_data <- dfwide %>%
              filter(treatment != trial_id)
            test_data <- dfwide %>%
              filter(treatment == trial_id)
            
          }
          
          model_graph <- lm(formula = as.formula(paste(
            graph, "~ ", covinterest, " + age + seks"
          )), data = training_data)
          
          rows <- grep(paste0("^", covinterest),
                       rownames(summary(model_graph)$coefficients),
                       value = TRUE)
          
          pvalues <- c(pvalues, as.numeric(summary(model_graph)$coefficients[rows, c(4), drop =
                                                                               FALSE]))
          samplesize <- c(samplesize, as.numeric(length(training_data)))
          
          preds <- predict(model_graph, newdata = test_data)
          predictions <- c(predictions, preds)
          actuals <- c(actuals, test_data[[graph]])
          
          
        }
        # save pvalue range
        # Calculate Fisher's combined test statistic
        X2 <- -2 * sum(log(pvalues))
        # Calculate degrees of freedom (twice the number of p-values)
        df <- 2 * length(pvalues)
        # Calculate the combined p-value
        Fisher_pvalue <- pchisq(X2, df, lower.tail = FALSE)
        #Stouffer
        z_scores <- qnorm(1 - pvalues)
        comb_z <- sum(z_scores) / sqrt(length(pvalues))
        Stoufer_pvalue <- 1 - pnorm(comb_z)
        # weighted Z-value
        weighted_z <- sum(samplesize * z_scores) / sqrt(sum(samplesize ^ 2))
        weighted_z_pvalue <- 1 - pnorm(weighted_z)
        
        
        # Compute the combined p-value using MinP-CCT-MinP (MCM)
        MCM_pvalue <- mcm_combination(pvalues)
        # harmonic mean P-value
        harmonic_pvalue <- harmonic_mean_pvalue(pvalues)
        
        
        
        pvalues_mat[i, ] <- c(
          min(pvalues),
          max(pvalues),
          median(pvalues),
          mean(pvalues),
          sd(pvalues),
          Fisher_pvalue,
          Stoufer_pvalue,
          weighted_z_pvalue,
          MCM_pvalue,
          harmonic_pvalue
        )
        
        mse_mat[i] <- mean((predictions - actuals) ^ 2)
        rrmse_mat[i] <- sqrt(mean((predictions - actuals) ^ 2)) / mean(actuals)
        
        
        r_squared <- 1 - (sum((actuals - predictions) ^ 2) / sum((actuals -
                                                                    mean(actuals)) ^ 2))
        
      }
      results_list[[paste0("df_", j, '_', LOSO)]] <- cbind(graphmeasures,sapply(as.data.frame(cbind(pvalues_mat, rrmse_mat)), round, 3))
      rm(pvalues_mat, rrmse_mat, mse_mat)
    }
    
  }
  
  df_LOSO_trial_Dx <- results_list$df_1_trial
  df_LOSO_trt_Dx <- results_list$df_1_treatment
  
  
  
  
  
  
  
  write.list2 <- list("Dx_trial" = df_LOSO_trial_Dx, "Dx_treatment" = df_LOSO_trt_Dx)
  
  
  wb <- createWorkbook()
  
  arialStyle <- createStyle(fontName = "Arial", fontSize = 8)
  
  for (sheet_name in names(write.list2)) {
    addWorksheet(wb, sheet_name)
    
    writeData(wb, sheet_name, write.list2[[sheet_name]])
    n_rows <- nrow(write.list2[[sheet_name]])
    n_cols <- ncol(write.list2[[sheet_name]])
    
    addStyle(
      wb,
      sheet = sheet_name,
      style = arialStyle,
      rows = 1:(n_rows + 1),
      cols = 1:n_cols,
      gridExpand = TRUE
    )
  }
  saveWorkbook(wb, file = outputfile2, overwrite = TRUE)
  
  
}
