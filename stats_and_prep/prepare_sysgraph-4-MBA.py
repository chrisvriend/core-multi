#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:15:07 2024

@author: cvriend
"""


import os 
import pandas as pd


import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np



basedir='/data/anw/anw-work/NP/projects/data_chris/CORE'
acq='dwi'
session='ses-T0'

graphdir=os.path.join(basedir,acq,session,'graph')
clindir=os.path.join(basedir,'clin')
statsdir=os.path.join(basedir,'stats','MBA')



# Define the order of ROIs
roi_order = ['DMN', 'FPN', 'DAN', 'VAN', 'SMN', 'LIM', 'VIS','SUBC']


# load clinical data
df_clin=pd.read_excel(os.path.join(clindir,'OCDPTSD_trials4analyses.xlsx'))
df_clin['Subj']=df_clin['Subj'] + '_' + session
df_clin['trial']=df_clin['trial'].astype('category')
df_clin['treatment']=df_clin['treatment'].astype('category')
df_clin=df_clin.set_index(['Subj'])
df_clin.drop(['T1MRI','Med'],axis=1,inplace=True)

df_MBA=pd.read_csv(os.path.join(graphdir,('MBA_input_acq-' + acq +'_CORE.txt')),delim_whitespace=True,index_col=['Subj'])

df_MBA2=df_MBA.join(df_clin).drop(['posttrial_sev'],axis=1)
df_MBA2=df_MBA2.rename({'seks':'sex'},axis=1)
df_MBA2['perc_improv']=df_MBA2['perc_improv']*100
df_MBA2['treatment']=df_MBA2['treatment'].replace({'DGT-EMDR':'DGT','SFT-IMRS':'SFT'}).astype('category')


df_MBA2['responder']=df_MBA2['responder'].replace({0:'znon-resp',1:'aresp'}).astype('category')

df_MBA2=df_MBA2.round(4)
df_MBA2.index.name='Subj'


df_MBA2=df_MBA2.reset_index()


# Create a dictionary to map ROIs to their indices
roi_index = {roi: index for index, roi in enumerate(roi_order)}

# Function to filter each group
def filter_upper_triangle(group):
    return group[group.apply(lambda row: roi_index[row['ROI1']] <= roi_index[row['ROI2']], axis=1)]

# Function to filter each group, removing lower triangle and diagonal
def filter_upper_triangle_exclude_diagonal(group):
    return group[group.apply(lambda row: roi_index[row['ROI1']] < roi_index[row['ROI2']], axis=1)]



# Group by 'Subj' and apply the filter function
#filtered_df = df_MBA2.groupby('Subj').apply(filter_upper_triangle).reset_index(drop=True).set_index(['Subj'])

# delete lower triangle and diagonal 
filtered_df = df_MBA2.groupby('Subj').apply(filter_upper_triangle_exclude_diagonal).reset_index(drop=True)
filtered_df.set_index('Subj',inplace=True)


filtered_df.to_csv(os.path.join(statsdir,('MBA_input_acq-' + acq +'_CORE_upd.txt')),sep=' ')

