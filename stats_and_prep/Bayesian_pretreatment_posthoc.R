

# clear variables
rm(list = ls())

#library(tidyverse)
library(dplyr)
library(readr)
library(moderndive)
library(reshape2)
library(brms)
library(readxl)
library(tidyr)
library(openxlsx)
library(ggplot2)


# input variables

date <- Sys.Date()
basedir = '/data/anw/anw-work/NP/projects/data_chris/CORE'
session = 'ses-T0'

for (acq in c('multi', 'dwi', 'func')) {
  if (acq == 'multi') {
    graphfile <- file.path(basedir,
                           acq,
                           session,
                           'output',
                           paste0('CORE_ML_', session, '_global.csv'))
  }  else  {
    graphfile <- file.path(basedir,
                           acq,
                           session,
                           'graph',
                           paste0('graph_global_acq-', acq, '_CORE.csv'))
  }
  
  clinfile <- file.path(basedir, 'clin', 'OCDPTSD_trials4analyses.xlsx')
 
  
  # Import datasets
  dfwide <- read_csv(graphfile)
  dfclin <- read_excel(clinfile)
  
  # convert characters to factors
  dfwide[sapply(dfwide, is.character)] <- lapply(dfwide[sapply(dfwide, is.character)], as.factor)
  
  
  # additional exclusions?
  #dfwide <-dfwide %>% filter(subjID!='sub-2261' | subjID!='sub-2274' | subjID!='sub-3117' | subjID!='sub-4177')
  
  dfclin <- dfclin %>%
    mutate(subjID = paste0(Subj, '_', session))
  
  # drop column T1MRI
  dfclin <- select(dfclin, -c(T1MRI))
  
  
  dfwide <- merge(dfwide, dfclin, by = 'subjID')
  
  # select only OCD trials
  # dfwide <- dfwide %>% filter(trial != 'PROSPER-B' &
  #                              trial != 'PROSPER-C')
  #
  #####################################################
  
  if (nrow(dfwide) < 120) {
    outputfile <- file.path(
      basedir,
      'stats',
      paste0("CORE_MixedModel_global_acq-", acq, "_", date, "_OCD.xlsx")
    )
    outputfile2 <- file.path(basedir,
                             'stats',
                             paste0("CORE_LOSO_global_acq-", acq, "_", date, "_OCD.xlsx"))
    dptvars <- c('perc_improv', 'responder', 'dYBOCS_dCAPS')
    if (acq == 'multi') {
      dfwide <- dfwide %>% mutate_at(vars('eigenvector'), ~ (. * 100))
      dfwide <- dfwide %>% mutate_at(vars('eccentricity'), ~ (. * 100))
    }
  } else {
    dptvars <- c('perc_improv', 'responder')
    
    
  }
  
  ######
  if (acq == 'multi') {
    graphmeasures <- c('eigenvector', 'eccentricity')
    
    
  } else  {
    graphmeasures <- c('GE', 'Qyeo', 'PC_glob_Schaef', 'Sigma')
    # only to check c('Q','PC_glob_ind')
    
    
  }
  
  
  for (j in 1:length(dptvars)) {
    covinterest <- dptvars[j]
    print(paste("dependent variable = ", covinterest))
    
    posteriors <- list()
    pplus <- c()

    
    for (i in 1:length(graphmeasures)) {
      graph = graphmeasures[i]
      
      print(paste("global graph measure = ", graph))
    #  formula <- paste(graph, "~ ", covinterest, " + age + seks + (1|treatment)")
      formula <- paste(graph, "~ ", covinterest, " + age + seks")
      
      # Gaussian distribution
      fit <- brm(
        formula = formula,
        data = dfwide,
        family = student(),
        iter = 4000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        prior = c(set_prior("normal(0,1)",class= "b"),set_prior("normal(0,1",class="Intercept"))
      )
      
    
      # evaluate fit
      sumfit<-capture.output(summary(fit))
      checkfit<-capture.output(pp_check(fit))
      fitcheck<-c("Model Summary:\n",sumfit,"\nPosterior Predictive Check:\n",checkfit)
      outputtxt<-file.path(basedir,"stats","bayes_global",paste0('BayesianSummary-acq-',acq,'_',graph,'_',covinterest,'.txt'))
      writeLines(fitcheck,outputtxt)

      posterior <- as_draws_df(fit)
      posterior_interest<-posterior[[paste0('b_',covinterest)]]
      posteriors[i]<-as.data.frame(posterior_interest)
      p_pos<-mean(posterior_interest> 0)
      pplus<-c(pplus,p_pos)
        
     

      # plot posterior distribution
      plot_color <- ifelse(p_pos <= 0.05, "blue", ifelse(p_pos <= 0.1, "cyan", ifelse(p_pos <=0.15,"skyblue",ifelse(
        p_pos <= 0.85,'grey41',ifelse(p_pos<=0.9,"yellow4",
        ifelse(p_pos <= 0.95, "yellow", "#C9182B"))))))

      ggplot(data.frame(posterior_interest),
             aes(x = posterior_interest)) + geom_density(fill = plot_color, alpha = 0.7) +
        geom_vline(
          xintercept = 0,
          linetype = "solid",
          alpha = 1,
          size = 1,
          color = "red"
        ) + theme_classic() +
        geom_text(aes(x=max(posterior_interest)*1.2,y=max(density(data.frame(posterior_interest)$posterior_interest)$y)*1.1,label=paste("P+ =",round(p_pos,2))),hjust=1.5,vjust=0,size=6,color="black") +
        theme(axis.title.x = element_text(size = 12),
              axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())
             


      ggsave(file.path(basedir,"stats","bayes_global",paste0('Ridgeplot_acq-',acq,'_',graph,'_',covinterest,'.pdf')),width=4,height=3)
    #  ggsave(file.path(basedir,"stats","bayes_global",paste0('Ridgeplot_acq-',acq,'_',graph,'_',covinterest,'.png')),dpi=400)
      
    }
    # 
    # posteriors_comb<-as.data.frame(do.call(cbind,posteriors))
    # names(posteriors_comb)<-graphmeasures
    # posteriors_long<-posteriors_comb %>% pivot_longer(cols = c(as.character(graphmeasures)), names_to = "parameter",values_to ="value")
    # 
    # ggplot(posteriors_long,aes(x=value)) + geom_density(alpha=0.5,adjust =1.5) + 
    #   geom_vline(
    #         xintercept = 0,
    #         linetype = "solid",
    #         alpha = 1,
    #         size = 1,
    #         color = "red"
    #       ) + theme_classic() + facet_wrap(~ parameter,ncol=-1,scales="free_y")

  }
  
  
}