#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:15:07 2024

@author: cvriend
"""


import os 
import pandas as pd
#import seaborn as sns
#import matplotlib.pyplot as plt
#import statsmodels.api as sm
import numpy as np


basedir='/data/anw/anw-work/NP/projects/data_chris/CORE'
clindir=os.path.join(basedir,'clin')
statsdir=os.path.join(basedir,'stats')


# List of acquisition types
acq_list = ['dwi', 'func']

# Loop through each acquisition type
for acq in acq_list:
    
    # load data
    df_T0=pd.read_csv(os.path.join(basedir,acq,'ses-T0','graph',('graph_global_acq-' + acq + '_CORE.csv')))
    df_T1=pd.read_csv(os.path.join(basedir,acq,'ses-T1','graph',('graph_global_acq-' + acq + '_CORE.csv')))
    
    df_T0[['Subj','session']]=df_T0.subjID.str.split(pat="_ses-",expand=True)
    df_T1[['Subj','session']]=df_T1.subjID.str.split(pat="_ses-",expand=True)
    
    df_T0.drop(['subjID'],axis=1,inplace=True)
    df_T1.drop(['subjID'],axis=1,inplace=True)
    df_T0.set_index(['Subj'],inplace=True)
    df_T1.set_index(['Subj'],inplace=True)
    
    
    df_T0 = df_T0.loc[df_T0.index.isin(df_T1.index)]
    
    df_long = pd.concat([df_T0, df_T1]).reset_index()
    # set T2 to T1
    df_long['session']=df_long['session'].replace({'T2': 'T1'}).astype('category')
    
    ## add clinical data
    # df_long.set_index(['Subj'],inplace=True)
    
    # # load clinical data
    # df_clin=pd.read_excel(os.path.join(clindir,'OCD_trials4analyses.xlsx'))
    
    # df_clin['trial']=df_clin['trial'].astype('category')
    # df_clin['treatment']=df_clin['treatment'].astype('category')
    # df_clin.set_index('Subj',inplace=True)
    
    # df_long=df_long.join(df_clin).reset_index()
    
    
    df_long.set_index(['Subj','session'],inplace=True)
    df_long=df_long.round(6)
    df_long=df_long.sort_values(['Subj','session'])
    
    df_long.to_csv(os.path.join(statsdir,('global_input_longMM_acq-' + acq +'_CORE.csv')))

###################
## multilayer ##
##################

acq='multi'

# load data
df_T0=pd.read_csv(os.path.join(basedir,acq,'ses-T0','output',('CORE_ML_ses-T0_global.csv')))
df_T1=pd.read_csv(os.path.join(basedir,acq,'ses-T1','output',('CORE_ML_ses-T1_global.csv')))

df_T0[['Subj','session']]=df_T0.subjID.str.split(pat="_ses-",expand=True)
df_T1[['Subj','session']]=df_T1.subjID.str.split(pat="_ses-",expand=True)

df_T0.drop(['subjID'],axis=1,inplace=True)
df_T1.drop(['subjID'],axis=1,inplace=True)
df_T0.set_index(['Subj'],inplace=True)
df_T1.set_index(['Subj'],inplace=True)


df_T0 = df_T0.loc[df_T0.index.isin(df_T1.index)]

df_long = pd.concat([df_T0, df_T1]).reset_index()
# set T2 to T1
df_long['session']=df_long['session'].replace({'T2': 'T1'}).astype('category')

## add clinical data
# df_long.set_index(['Subj'],inplace=True)

# # load clinical data
# df_clin=pd.read_excel(os.path.join(clindir,'OCD_trials4analyses.xlsx'))

# df_clin['trial']=df_clin['trial'].astype('category')
# df_clin['treatment']=df_clin['treatment'].astype('category')
# df_clin.set_index('Subj',inplace=True)

# df_long=df_long.join(df_clin).reset_index()


df_long.set_index(['Subj','session'],inplace=True)
df_long=df_long.round(6)
df_long=df_long.sort_values(['Subj','session'])

df_long.to_csv(os.path.join(statsdir,('global_input_longMM_acq-' + acq +'_CORE.csv')))


