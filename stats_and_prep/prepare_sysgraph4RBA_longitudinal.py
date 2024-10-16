#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:15:07 2024

Updated by: cvriend
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm

# Base directory for data
basedir = '/data/anw/anw-work/NP/projects/data_chris/CORE'
sessionT1 = 'ses-T1'

# Define the order of ROIs
roi_order = ['DMN', 'FPN', 'DAN', 'VAN', 'SMN', 'LIM', 'VIS', 'SUBC']

# Load clinical data
clindir = os.path.join(basedir, 'clin')
df_clin = pd.read_excel(os.path.join(clindir, 'OCDPTSD_trials4analyses.xlsx'))
df_clin['subjid'] = df_clin['Subj']
df_clin['trial'] = df_clin['trial'].astype('category')
df_clin['treatment'] = df_clin['treatment'].astype('category')
df_clin = df_clin.set_index(['subjid']).drop(['Subj'], axis=1)
df_clin.drop(['T1MRI','Med'], axis=1, inplace=True)

# Define acquisition types to iterate over
acq_types = ['multi', 'dwi', 'func']

# Loop over acquisition types
for acq in acq_types:
    # Define directories
    graphdir = os.path.join(basedir, acq, sessionT1, 'graph')
    graphdirT0 = os.path.join(basedir, acq, 'ses-T0', 'graph')
    statsdir = os.path.join(basedir, 'stats', 'RBA')
    
    # Process ses-T1
    if acq == 'multi':
        outdir = os.path.join(basedir, 'multi', sessionT1, 'output')
        outdirT0 = os.path.join(basedir, 'multi', 'ses-T0', 'output')
        df_temp = pd.read_csv(os.path.join(outdir, f'CORE_ML_{sessionT1}_system.csv'))
        df_tempT0 = pd.read_csv(os.path.join(outdirT0, 'CORE_ML_ses-T0_system.csv'))
    else:
        df_temp = pd.read_csv(os.path.join(graphdir, f'graph_PC_subnetworks_acq-{acq}_CORE.csv'))
        df_tempT0 = pd.read_csv(os.path.join(graphdirT0, f'graph_PC_subnetworks_acq-{acq}_CORE.csv'))
    
    # Extract subject ID and time from the 'subjid' or 'subjID' column
    id_column = 'subjID' 
    df_temp[['Subj', 'time']] = df_temp[id_column].str.split('_ses-', expand=True)
    df_temp['time'] = df_temp['time'].replace({'T2': 'T1'}).astype('category')
    df_temp.set_index(['Subj'], inplace=True)
    df_temp.drop([id_column, 'time'], inplace=True, axis=1)
    
    df_tempT0[['Subj', 'time']] = df_tempT0[id_column].str.split('_ses-', expand=True)
    df_tempT0.set_index(['Subj'], inplace=True)
    df_tempT0.drop([id_column, 'time'], inplace=True, axis=1)
    
    # Align and subtract values between T1 and T0
    _, df_tempT0_aligned = df_temp.align(df_tempT0, join='inner')
    
    df_long = pd.DataFrame(columns=roi_order)
    for roi in roi_order:
        df_long[roi] = df_temp[roi] - df_tempT0_aligned[roi]
    
    df_long.reset_index(inplace=True)
    df_PC = pd.melt(df_long, id_vars=['Subj'], var_name='network', value_name='Y').set_index(['Subj'])
    df_PC['network'] = df_PC['network'] + '_' + acq
    
    # Join with clinical info
    df_PC = df_PC.join(df_clin).drop(['posttrial_sev'], axis=1)
    df_PC = df_PC.rename({'seks': 'sex'}, axis=1)
    df_PC['perc_improv'] = df_PC['perc_improv'] * 100
    df_PC['responder'] = df_PC['responder'].replace({0: 'znon-resp', 1: 'aresp'}).astype('category')
    df_PC = df_PC.round(4)
    
    # Save to file
    df_PC.to_csv(os.path.join(statsdir, f'RBA_PC_acq-{acq}_CORE_long.txt'), sep=' ')
    
    # Temporary plot for 'LIM' network and improvement percentage
    # LIM_eccentricity = df_PC[df_PC['network'] == f'LIM_{acq}']['Y']
    # perc_improv = df_PC[df_PC['network'] == f'LIM_{acq}']['perc_improv']
    # dfplot = pd.concat([LIM_eccentricity, perc_improv], axis=1)
    # sns.regplot(x='perc_improv', y='Y', data=dfplot, scatter=True, color='black')
    # plt.show()

# Note: Ensure that both 'dwi' and 'func' data exist before running analysis

