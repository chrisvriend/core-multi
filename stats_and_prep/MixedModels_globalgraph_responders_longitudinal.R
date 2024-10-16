## C. Vriend - Amsterdam UMC - July '24
## perform multivariate mixed model analysis on CORE data


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
library(rmcorr)
library(ggplot2)


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
  df<-select(df, -c('P-value','P-value.1'))
  
  
  df<-df %>% relocate("Pvalues_crude", .after = "BSE_crude") %>% 
    relocate("CI_crude", .before = "Pvalues_crude") %>% 
    relocate("Pvalues_adj", .after = "BSE_adj") %>% 
    relocate("CI_adj", .before = "Pvalues_adj")  
  return(df)
}

# Define the function to perform operations
clean_dataframe_mod <- function(df) {
  
  # multiply by 100 and round to two decimals
  df[c('B','95CI-','95CI+','SE')] <- df[c('B','95CI-','95CI+','SE')]*100
  df[c('95CI-','95CI+','SE')] <- sapply(df[c('95CI-','95CI+','SE')],round,2)
  # unite columns for table  
  df <- df %>%
    unite("CI_crude", c("95CI-","95CI+"), sep= " | ", 
          remove = TRUE)
  df$BSE_crude <- sprintf("%2.3f [%2.3f]", df$B, df$SE)

  df<-select(df, -c(B, SE))
  
  # calculate FDR correction ##
  
  #df['Pfdr_crude']<-p.adjust(df[['P-value']], method = "BH")
  #df['Pfdr_adj']<-p.adjust(df[['P-value.1']], method = "BH")
  
  # round to 3 decs
  #df[c('Pfdr_crude','Pfdr_adj')] <- sapply(df[c('Pfdr_crude','Pfdr_adj')],round,3)
  
  #df$Pvalues_crude <- sprintf("%2.3f (%2.3f)", df$Pfdr_crude, df[['P-value']])
  #df$Pvalues_adj <- sprintf("%2.3f (%2.3f)", df$Pfdr_adj, df[['P-value.1']])
  df$Pvalues_crude <- sprintf("%2.3f", df[['P-value']])
  # remove cols
  df<-select(df, -c('P-value'))
  
  
  df<-df %>% relocate("Pvalues_crude", .after = "BSE_crude") %>% 
    relocate("CI_crude", .before = "Pvalues_crude")
  return(df)
}

date <- Sys.Date()

basedir='/data/anw/anw-work/NP/projects/data_chris/CORE'
statsdir<-file.path(basedir,'stats')
clinfile<-file.path(basedir,'clin','OCDPTSD_trials4analyses.xlsx')

for (acq in c('multi', 'dwi', 'func')) {
  

graphfile<-file.path(statsdir,paste0('global_input_longMM_acq-', acq,'_CORE.csv'))
outputfile<-file.path(basedir,'stats',paste0("CORE_MixedModel_long_global_responders_acq-",acq,"_",date,".xlsx"))
outputfile2<-file.path(basedir,'stats',paste0("CORE_LOSO_long_global_responders_acq-",acq,"_",date,".xlsx"))

# Import datasets
dflong <-read_csv(graphfile)
dfclin <-read_excel(clinfile)


# drop column T1MRI
dfclin<-select(dfclin, -c(T1MRI))


if (acq =='dwi'){
#  dflong <- dflong %>% mutate_at(vars('GE'), ~(. * 100) )
  graphmeasures=c('GE','Qyeo','PC_glob_Schaef','Sigma','Q','PC_glob_ind')
  dflong <-select(dflong, -c(GE_norm,GCC,GCC_norm,lambda,lambda_norm))
  
  
} else if (acq=='func') {
 # dflong <- dflong %>% mutate_at(vars('GE'), ~(. * 10) )
  graphmeasures=c('GE','Qyeo','PC_glob_Schaef','Sigma','Q','PC_glob_ind')
  dflong <-select(dflong, -c(GE_norm,GCC,GCC_norm,lambda,lambda_norm))
 
  
} else if (acq=='multi'){
  graphmeasures=c('eigenvector','eccentricity','BC')
  
}

dfwide <- dflong %>%
  pivot_wider(names_from = session, values_from=graphmeasures)


# additional exclusions? 
#dflong <-dflong %>% filter(subjID!='sub-2261' | subjID!='sub-2274' | subjID!='sub-3117' | subjID!='sub-4177')


dfwide<-merge(dfwide,dfclin, by='Subj')
dflong<-merge(dflong,dfclin, by='Subj')


## plot 

ggplot(dflong,aes(x=responder,y=Sigma,fill=session)) + 
  geom_bar(stat="identity",position="dodge")


ggplot(dflong,aes(x=responder,y=eigenvector,fill=session)) + 
  geom_bar(stat="identity",position="dodge")



# convert characters to factors
dfwide[sapply(dfwide, is.character)] <- lapply(dfwide[sapply(dfwide, is.character)], 
                                               as.factor)


# convert characters to factors
dflong[sapply(dflong, is.character)] <- lapply(dflong[sapply(dflong, is.character)], 
                                               as.factor)

dptvars<-c('responder')

# List to store resulting dataframes
results_list <- list()
results_list_treat <- list()
results_list_trial <- list()


for (j in 1:length(dptvars)) {

covinterest<-dptvars[j]
print(paste("dependent variable = ", covinterest))


# pre-allocation
temps_graph_paired <- matrix(0, nrow = length(graphmeasures), ncol = 4,
                             dimnames=list(graphmeasures,rep(c("mean_diff","t","df","P-value"),1)))

temps_graph_coeff <- matrix(0, nrow = length(graphmeasures), ncol = 10,
                      dimnames=list(graphmeasures,rep(c("B","SE","P-value","95CI-","95CI+"),2)))

temps_graph_coeff_trial <- matrix(0, nrow = length(graphmeasures), ncol = 5,
                            dimnames=list(graphmeasures,rep(c("B","SE","P-value","95CI-","95CI+"),1)))
temps_graph_coeff_treat <- matrix(0, nrow = length(graphmeasures), ncol = 5,
                            dimnames=list(graphmeasures,rep(c("B","SE","P-value","95CI-","95CI+"),1)))



for (i in 1:length(graphmeasures)) {
  graph=graphmeasures[i]
  
  df_baseline <- dflong[dflong$session == 'T0', c("Subj", graph)]
  names(df_baseline)[2] <- paste0(graph,'_T0')
  
  # Merge baseline measure back into the main dataframe
  dflong <- merge(dflong, df_baseline, by = "Subj")
  
  
  print(paste("global graph measure = ", graph))
  
  model_ttest<-t.test(dfwide[[paste0(graph,"_T0")]],dfwide[[paste0(graph,"_T1")]],paired=TRUE,alternative="two.sided")
  # model_graph_crude<-lm(formula = as.formula(paste0(graph,"_T1"," ~ responder + ",graph,"_T0")),data=dfwide)
  # model_graph_adj<-lm(formula = as.formula(paste0(graph,"_T1"," ~ responder + ",graph,"_T0"," + age + seks")),data=dfwide)
  # model_graph_trial<-lm(formula = as.formula(paste0(graph,"_T1"," ~ responder + ",graph,"_T0"," + age + seks + trial")),data=dfwide)
  # model_graph_treat<-lm(formula = as.formula(paste0(graph,"_T1"," ~ responder + ",graph,"_T0"," + age + seks + treatment")),data=dfwide)
  # 
  
  model_graph_crude<-lmer(formula = as.formula(paste0(graph," ~ session * ", covinterest," + ",graph ,"_T0"," + (1|Subj)")),data=dflong)
  model_graph_adj<-lmer(formula = as.formula(paste0(graph," ~ session * ", covinterest," + ",graph ,"_T0"," + age + seks + (1|Subj)")),data=dflong)
  model_graph_trial<-lmer(formula = as.formula(paste0(graph," ~ session * ", covinterest," + ",graph ,"_T0"," + age + seks + trial + (1|Subj)")),data=dflong)
  model_graph_treat<-lmer(formula = as.formula(paste0(graph," ~ session * ", covinterest," + ",graph ,"_T0"," + age + seks + treatment + (1|Subj)")),data=dflong)
  
  
  # standard models
    temps_graph_paired[i,1:4]<-cbind(model_ttest$estimate,model_ttest$statistic,model_ttest$parameter,model_ttest$p.value)
    rows_crude<-grep(paste0('sessionT1:',covinterest), rownames(summary(model_graph_crude)$coefficients), value = TRUE)

    coeff_graph_crude<-cbind(summary(model_graph_crude)$coefficients[rows_crude,c(1,2,5),drop=FALSE],
                            confint(model_graph_crude)[rows_crude, ,drop=FALSE])
    coeff_graph_adj<-cbind(summary(model_graph_adj)$coefficients[rows_crude,c(1,2,5),drop=FALSE],
                           confint(model_graph_adj)[rows_crude, ,drop=FALSE])
     coeff_graph_trial<-cbind(summary(model_graph_trial)$coefficients[rows_crude,c(1,2,5),drop=FALSE],
                             confint(model_graph_trial)[rows_crude, ,drop=FALSE])
     coeff_graph_treat<-cbind(summary(model_graph_treat)$coefficients[rows_crude,c(1,2,5),drop=FALSE],
                                   confint(model_graph_treat)[rows_crude, ,drop=FALSE])
   
  
 # round to 3 decimals
     
 coeff_graph_crude<-sapply(coeff_graph_crude,round,3)
 coeff_graph_adj<-sapply(coeff_graph_adj,round,3)
 coeff_graph_trial<-sapply(coeff_graph_trial,round,3)
 coeff_graph_treat<-sapply(coeff_graph_treat,round,3)
 
 temps_graph_coeff[i,1:5]=coeff_graph_crude
 temps_graph_coeff[i,6:10]=coeff_graph_adj
 temps_graph_coeff_trial[i,1:5]=coeff_graph_trial
 temps_graph_coeff_treat[i,1:5]=coeff_graph_treat
# temps_graph_paired<-sapply(temps_graph_paired,round,3)

 rm(coeff_graph_crude,coeff_graph_trial,coeff_graph_treat)

}
all_ttest<-as.data.frame(temps_graph_paired)

all_coeff<-as.data.frame(temps_graph_coeff)
all_coeff_trial<-as.data.frame(temps_graph_coeff_trial)
all_coeff_treat<-as.data.frame(temps_graph_coeff_treat)
names(all_ttest) <- make.unique(names(all_ttest))
names(all_coeff) <- make.unique(names(all_coeff))
names(all_coeff_trial) <- make.unique(names(all_coeff_trial))
names(all_coeff_treat) <- make.unique(names(all_coeff_treat))

df_paired <- tibble::rownames_to_column(all_ttest, "globalmeasure")

df_globalall <- tibble::rownames_to_column(all_coeff, "globalmeasure")
df_globalall_trial <- tibble::rownames_to_column(all_coeff_trial, "globalmeasure")
df_globalall_treat <- tibble::rownames_to_column(all_coeff_treat, "globalmeasure")


results_list[[paste0("df_", j)]] <- clean_dataframe(df_globalall)
results_list_trial[[paste0("df_", j)]] <-clean_dataframe_mod(df_globalall_trial)
results_list_treat[[paste0("df_", j)]] <-clean_dataframe_mod(df_globalall_treat)


}

# Rename and save the dataframes outside the for loop
df_responder <- results_list$df_1
df_responder_trial <- results_list_trial$df_1
df_responder_treat <- results_list_treat$df_1


### save output ###
library(openxlsx)


write.list <- list("Paired t-test" = df_paired,
                   "responder" = df_responder,
                   "responder_trial" = df_responder_trial,
                   "responder_treat" = df_responder_treat)

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

#write.xlsx(write.list, file = outputfile)



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
dptvars <- c('responder')

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
    # 
    # mse_mat <- matrix(
    #   0,
    #   nrow = length(graphmeasures),
    #   ncol = 1,
    #   dimnames = list(graphmeasures, c("MSE"))
    # )
    # rrmse_mat <- matrix(
    #   0,
    #   nrow = length(graphmeasures),
    #   ncol = 1,
    #   dimnames = list(graphmeasures, c("RRMSE"))
    # )
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
      
      
 #     predictions <- c()
#      actuals <- c()
      pvalues <- c()
      samplesize <- c()
      for (trial_id in trials) {
        if (LOSO == 'trial') {
          training_data <- dflong %>%
            filter(trial != trial_id)
       #   test_data <- dflong %>%
      #      filter(trial == trial_id)
          
          
        } else if (LOSO == 'treatment') {
          training_data <- dflong %>%
            filter(treatment != trial_id)
       #   test_data <- dflong %>%
      #      filter(treatment == trial_id)
          
        }
        
        model_graph<-lmer(formula = as.formula(paste0(graph," ~ session * ", covinterest," + ",graph ,"_T0"," + (1|Subj)")),data=training_data)
     
        
        rows <- grep(paste0('sessionT1:',covinterest),
                     rownames(summary(model_graph)$coefficients),
                     value = TRUE)
      
        # 5t column because it is lmer
        pvalues <- c(pvalues, as.numeric(summary(model_graph)$coefficients[rows, c(5), drop =
                                                                             FALSE]))  
        # pvalues <- c(pvalues, as.numeric(summary(model_graph)$coefficients[rows, c(4), drop =
        #                                                                      FALSE]))
        samplesize <- c(samplesize, as.numeric(length(training_data)))
        
       # preds <- predict(model_graph, newdata = test_data)
      #  predictions <- c(predictions, preds)
      #  actuals <- c(actuals, test_data[[graph]])
        
        
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
      
      # mse_mat[i] <- mean((predictions - actuals) ^ 2)
      # rrmse_mat[i] <- sqrt(mean((predictions - actuals) ^ 2)) / mean(actuals)
      # 
      # 
      # r_squared <- 1 - (sum((actuals - predictions) ^ 2) / sum((actuals -
      #                                                             mean(actuals)) ^ 2))
      
    }
    results_list[[paste0("df_", j, '_', LOSO)]] <- cbind(graphmeasures,sapply(as.data.frame(cbind(pvalues_mat)), round, 3))
    rm(pvalues_mat)
  }
  
}

# save to dataframes
df_LOSO_trial_responder <- results_list$df_1_trial
df_LOSO_trt_responder <- results_list$df_1_treatment


write.list2 <- list(
  "responder_trial" = df_LOSO_trial_responder,
  "responder_treatment" = df_LOSO_trt_responder
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



