## C. Vriend - Amsterdam UMC - July '24
## perform regression analysis on pre-to-post-treatment CORE data and associate with percentage improvement
## Leave-one-sample-out validation (trial or treatment) and calculation of harmonic P-value across folds.

# clear variables
rm(list=ls())

#library(tidyverse)
library(dplyr)
library(readr)
library(moderndive)
library(reshape2)
library(lme4)
library(readxl)
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
library(tidyr)



# Define the function to perform operations
clean_dataframe <- function(df) {
  
  # multiply by 100 and round to two decimals
  df[c('B','B.1','95CI-','95CI+','95CI-.1','95CI+.1','SE','SE.1')] <- df[c('B','B.1','95CI-','95CI+','95CI-.1','95CI+.1','SE','SE.1')]*100
  df[c('95CI-','95CI+','95CI-.1','95CI+.1','SE','SE.1')] <- sapply(df[c('95CI-','95CI+','95CI-.1','95CI+.1','SE','SE.1')],round,2)
  # unite columns for table  
  df <- df %>%
    unite("CI_crude", c("95CI-","95CI+"), sep= " | ", 
          remove = TRUE) %>% unite("CI_adj", c("95CI-.1","95CI+.1"), sep= " | ", 
                                   remove = TRUE)
  df$BSE_crude <- sprintf("%2.3f [%2.3f]", df$B, df$SE)
  df$BSE_adj <- sprintf("%2.3f [%2.3f]", df$B.1, df$SE.1)
  
  df<-select(df, -c(B.1, B, SE,SE.1))
  df$Pvalues_crude <- sprintf("%2.3f", df[['P-value']])
  df$Pvalues_adj <- sprintf("%2.3f", df[['P-value.1']])
  # remove cols
  df<-select(df, -c('P-value','P-value.1'))
  
  
  df<-df %>% relocate("Pvalues_crude", .after = "BSE_crude") %>% 
    relocate("CI_crude", .before = "Pvalues_crude") %>% 
    relocate("Pvalues_adj", .after = "BSE_adj") %>% 
    relocate("CI_adj", .before = "Pvalues_adj")  
  return(df)
}

date <- Sys.Date()

basedir='/data/anw/anw-work/NP/projects/data_chris/CORE'
statsdir<-file.path(basedir,'stats')
clinfile<-file.path(basedir,'clin','OCDPTSD_trials4analyses.xlsx')

for (acq in c('multi', 'dwi', 'func')) {
  
graphfile<-file.path(statsdir,paste0('global_input_longMM_acq-', acq,'_CORE.csv'))

outputfile<-file.path(basedir,'stats',paste0("CORE_MixedModel_long_global_acq-",acq,"_",date,".xlsx"))
outputfile2 <- file.path(basedir,
                         'stats',
                         paste0("CORE_LOSO_long_global_acq-", acq, "_", date, ".xlsx"))

# Import datasets
dflong <-read_csv(graphfile)
dfclin <-read_excel(clinfile)


# drop column T1MRI
dfclin<-select(dfclin, -c(T1MRI))

# convert characters to factors
dflong[sapply(dflong, is.character)] <- lapply(dflong[sapply(dflong, is.character)], 
                                       as.factor)
dflong<-merge(dflong,dfclin, by='Subj')


# Assuming your data frame is named df
dflong <- dflong %>%
  mutate(Dx = case_when(
    grepl("ARRIBA|TIPICCO", Subj) ~ "OCD",
    grepl("PROSPER", Subj) ~ "PTSD",
    TRUE ~ NA_character_  # Fills with NA if no condition is met
  ))


dflong <- dflong %>% 
  mutate(clinsev = case_when (
    session == "T0" ~ baseline_sev,
    session == "T1" ~ posttrial_sev
  ))

dflong <- dflong %>%
  group_by( Subj) %>%
  mutate(
    clinsev_mean = mean(clinsev, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(Dx) %>%
  # z-transform separately per group
  mutate(
    clinsev_centered = clinsev - clinsev_mean,
    zclinsev_centered = scale(clinsev_centered, center = TRUE, scale = TRUE)[,1]
  ) %>%
  ungroup()




if (acq =='dwi'){
#  dflong <- dflong %>% mutate_at(vars('GE'), ~(. * 100) )
  graphmeasures=c('GE','Qyeo','PC_glob_Schaef','Sigma','Q','PC_glob_ind')
  
  
} else if (acq=='func') {
 # dflong <- dflong %>% mutate_at(vars('GE'), ~(. * 10) )
  graphmeasures=c('GE','Qyeo','PC_glob_Schaef','Sigma','Q','PC_glob_ind')
  
  
} else if (acq=='multi'){
  graphmeasures=c('eigenvector','eccentricity','BC')
  
}

dptvars<-c('zclinsev_centered')

# List to store resulting dataframes
results_list <- list()
results_list_treat <- list()
results_list_trial <- list()


for (j in 1:length(dptvars)) {

covinterest<-dptvars[j]
print(paste("dependent variable = ", covinterest))


# pre-allocation
temps_graph_coeff <- matrix(0, nrow = length(graphmeasures), ncol = 10,
                      dimnames=list(graphmeasures,rep(c("B","SE","P-value","95CI-","95CI+"),2)))

temps_graph_coeff_trial <- matrix(0, nrow = length(graphmeasures), ncol = 10,
                            dimnames=list(graphmeasures,rep(c("B","SE","P-value","95CI-","95CI+"),2)))
temps_graph_coeff_treat <- matrix(0, nrow = length(graphmeasures), ncol = 10,
                            dimnames=list(graphmeasures,rep(c("B","SE","P-value","95CI-","95CI+"),2)))


for (i in 1:length(graphmeasures)) {
  graph=graphmeasures[i]
  
  print(paste("global graph measure = ", graph))
   
  df <- dflong %>%
    group_by(Subj) %>%
    mutate(
      graph_mean = mean(.data[[graph]], na.rm = TRUE),
        ) %>% 
    ungroup() %>%
    mutate(
      graph_centered = .data[[graph]] - graph_mean,
      zgraph_centered = scale(graph_centered, center = TRUE, scale = TRUE)[,1],

    )
  
  model_graph_crude<-lm(formula = as.formula(paste("zgraph_centered ~", covinterest, "+ as.factor(Subj)")),data=df)
  model_graph_adj<-lm(formula = as.formula(paste("zgraph_centered ~", covinterest, "+ age + seks + as.factor(Subj)")),data=df)
  
  # exact same outcome due to subj as fixed effect
  #model_graph_crude_trial<-lm(formula = as.formula(paste("zgraph_centered ~", covinterest, "+ as.factor(Subj) + trial")),data=df)
  #model_graph_crude_treat<-lm(formula = as.formula(paste("zgraph_centered ~", covinterest, "+ as.factor(Subj) + treatment")),data=df)
  
  #model_graph_adj_trial<-lm(formula = as.formula(paste("zgraph_centered ~", covinterest, "+ as.factor(Subj) + age + seks + trial")),data=df)
  #model_graph_adj_treat<-lm(formula = as.formula(paste("zgraph_centered ~", covinterest, "+ as.factor(Subj) + age + seks + treatment")),data=df)
  
  
  # standard models
    rows_crude <- grep(paste0("^",covinterest), rownames(summary(model_graph_crude)$coefficients), value = TRUE)
    rows_adj <- grep(paste0("^",covinterest), rownames(summary(model_graph_adj)$coefficients), value = TRUE)

    coeff_graph_crude<-cbind(summary(model_graph_crude)$coefficients[rows_crude,c(1,2,4),drop=FALSE],
                            confint(model_graph_crude)[rows_crude, ,drop=FALSE])
    coeff_graph_adj<-cbind(summary(model_graph_adj)$coefficients[rows_adj,c(1,2,4),drop=FALSE],
                           confint(model_graph_adj)[rows_adj, ,drop=FALSE])

  #  -trial models  
 #   rows_crude <- grep(paste0("^",covinterest), rownames(summary(model_graph_crude_trial)$coefficients), value = TRUE)
 #   rows_adj <- grep(paste0("^",covinterest), rownames(summary(model_graph_adj_trial)$coefficients), value = TRUE)
    
 #   coeff_graph_crude_trial<-cbind(summary(model_graph_crude_trial)$coefficients[rows_crude,c(1,2,4),drop=FALSE],
  #                           confint(model_graph_crude_trial)[rows_crude, ,drop=FALSE])
   
  #  coeff_graph_adj_trial<-cbind(summary(model_graph_adj_trial)$coefficients[rows_crude,c(1,2,4),drop=FALSE],
  #                           confint(model_graph_adj_trial)[rows_crude, ,drop=FALSE])
    
    # -treatment models  
 #   rows_crude <- grep(paste0("^",covinterest), rownames(summary(model_graph_crude_treat)$coefficients), value = TRUE)

  #  coeff_graph_crude_treat<-cbind(summary(model_graph_crude_treat)$coefficients[rows_crude,c(1,2,4),drop=FALSE],
  #                                 confint(model_graph_crude_treat)[rows_crude, ,drop=FALSE])
  #  coeff_graph_adj_treat<-cbind(summary(model_graph_adj_treat)$coefficients[rows_crude,c(1,2,4),drop=FALSE],
   #                          confint(model_graph_adj_treat)[rows_crude, ,drop=FALSE])
    
    
 # round to 3 decimals
 coeff_graph_crude<-sapply(coeff_graph_crude,round,3)
 coeff_graph_adj<-sapply(coeff_graph_adj,round,3)
 #coeff_graph_crude_trial<-sapply(coeff_graph_crude_trial,round,3)
 #coeff_graph_crude_treat<-sapply(coeff_graph_crude_treat,round,3)
 #coeff_graph_adj_trial<-sapply(coeff_graph_adj_trial,round,3)
 #coeff_graph_adj_treat<-sapply(coeff_graph_adj_treat,round,3)

 temps_graph_coeff[i,1:5]=coeff_graph_crude
 temps_graph_coeff[i,6:10]=coeff_graph_adj
 
 #temps_graph_coeff_trial[i,1:5]=coeff_graph_crude_trial
 #temps_graph_coeff_trial[i,6:10]=coeff_graph_adj_trial
 #temps_graph_coeff_treat[i,1:5]=coeff_graph_crude_treat
 #temps_graph_coeff_treat[i,6:10]=coeff_graph_adj_treat
 

 rm(coeff_graph_crude,coeff_graph_crude_trial,coeff_graph_crude_treat,
    coeff_graph_adj,coeff_graph_adj_treat,coeff_graph_adj_trial,df)

}

all_coeff<-as.data.frame(temps_graph_coeff)
#all_coeff_trial<-as.data.frame(temps_graph_coeff_trial)
#all_coeff_treat<-as.data.frame(temps_graph_coeff_treat)
names(all_coeff) <- make.unique(names(all_coeff))
#names(all_coeff_trial) <- make.unique(names(all_coeff_trial))
#names(all_coeff_treat) <- make.unique(names(all_coeff_treat))

df_globalall <- tibble::rownames_to_column(all_coeff, "globalmeasure")
#df_globalall_trial <- tibble::rownames_to_column(all_coeff_trial, "globalmeasure")
#df_globalall_treat <- tibble::rownames_to_column(all_coeff_treat, "globalmeasure")


results_list[[paste0("df_", j)]] <- clean_dataframe(df_globalall)
#esults_list_trial[[paste0("df_", j)]] <-clean_dataframe(df_globalall_trial)
#results_list_treat[[paste0("df_", j)]] <-clean_dataframe(df_globalall_treat)


}

# Rename and save the dataframes outside the for loop
df_percchange <- results_list$df_1
#df_percchange_trial <- results_list_trial$df_1
#df_percchange_treat <- results_list_treat$df_1


### save output ###
library(openxlsx)


write.list <- list("perchange" = df_percchange,
                   "perchange_trial" = df_percchange_trial,
                   "perchange_treat" = df_percchange_treat)

wb<-createWorkbook()

arialStyle <-createStyle(fontName = "Arial", fontSize=8)

for (sheet_name in names(write.list)) {
  addWorksheet(wb, sheet_name)
  
  writeData(wb,sheet_name,write.list[[sheet_name]])
  n_rows<-nrow(write.list[[sheet_name]])
  n_cols<-ncol(write.list[[sheet_name]])
  
  addStyle(wb,sheet = sheet_name, style=arialStyle,
           rows=1:(n_rows + 1 ),
           cols=1:n_cols,gridExpand=TRUE)
}
saveWorkbook(wb,file=outputfile,overwrite=TRUE)

write.xlsx(write.list, file = outputfile)


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

# treatment outcome variables
dptvars<-c('zclinsev_centered')

results_list <- list()

for (LOSO in c('trial', 'treatment')) {
  print(LOSO)

  if (LOSO == 'trial') {
    trials <- unique(dflong$trial)
  } else if (LOSO == 'treatment') {
    trials <- unique(dflong$treatment)
  }


  for (j in 1:length(dptvars)) {
    covinterest <- dptvars[j]
    print(paste("dependent variable = ", covinterest))

    # pre-allocation

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

        df <- dflong %>%
        group_by(Subj) %>%
        mutate(
          graph_mean = mean(.data[[graph]], na.rm = TRUE),
        ) %>%
        ungroup() %>%
        mutate(
          graph_centered = .data[[graph]] - graph_mean,
          zgraph_centered = scale(graph_centered, center = TRUE, scale = TRUE)[,1],

        )
            
        df <- df %>%
              group_by( Subj) %>%
              mutate(
                clinsev_mean = mean(clinsev, na.rm = TRUE)
              ) %>%
              ungroup() %>%
              group_by(Dx) %>%
              # z-transform separately per group
              mutate(
                clinsev_centered = clinsev - clinsev_mean,
                zclinsev_centered = scale(clinsev_centered, center = TRUE, scale = TRUE)[,1]
              ) %>%
              ungroup()
           
            
      pvalues <- c()
      samplesize <- c()
      for (trial_id in trials) {
        if (LOSO == 'trial') {
          training_data <- df %>%
            filter(trial != trial_id)
          test_data <- df %>%
            filter(trial == trial_id)


        } else if (LOSO == 'treatment') {
          training_data <- df %>%
            filter(treatment != trial_id)
          test_data <- df %>%
            filter(treatment == trial_id)

        }


        model_graph <- lm(formula = as.formula(paste(
          "zgraph_centered", "~ ", covinterest, " + as.factor(Subj)")), data = training_data)

        rows <- grep(paste0("^", covinterest),
                     rownames(summary(model_graph)$coefficients),
                     value = TRUE)

        pvalues <- c(pvalues, as.numeric(summary(model_graph)$coefficients[rows, c(4), drop =
                                                                             FALSE]))
        samplesize <- c(samplesize, as.numeric(length(training_data)))

   
      }
      # save pvalue range
      # Calculate Fisher's combined test statistic
      X2 <- -2 * sum(log(pvalues))
      # Calculate degrees of freedom (twice the number of p-values)
      dfx <- 2 * length(pvalues)
      # Calculate the combined p-value
      Fisher_pvalue <- pchisq(X2, dfx, lower.tail = FALSE)
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

 
    }
    results_list[[paste0("df_", j, '_', LOSO)]] <- cbind(graphmeasures,sapply(as.data.frame(cbind(pvalues_mat)), round, 3))

        rm(pvalues_mat)
  }

}

# save to dataframes
df_LOSO_trial_percchange <- results_list$df_1_trial
df_LOSO_trt_percchange <- results_list$df_1_treatment


write.list2 <- list(
  "perchange_trial" = df_LOSO_trial_percchange,
  "perchange_treatment" = df_LOSO_trt_percchange
)

# save and format output
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