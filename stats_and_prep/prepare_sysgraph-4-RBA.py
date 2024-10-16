#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:15:07 2024

@author: cvriend
"""

import os
import pandas as pd
import numpy as np

# Define base directories
basedir = '/data/anw/anw-work/NP/projects/data_chris/CORE'
session = 'ses-T0'

# Define clinical and stats directories
clindir = os.path.join(basedir, 'clin')
statsdir = os.path.join(basedir, 'stats', 'RBA')

# Load clinical data
df_clin = pd.read_excel(os.path.join(clindir, 'OCDPTSD_trials4analyses.xlsx'))
df_clin['Subj'] = df_clin['Subj'] + '_' + session
df_clin['trial'] = df_clin['trial'].astype('category')
df_clin['treatment'] = df_clin['treatment'].astype('category')
df_clin = df_clin.set_index(['Subj']).drop(['T1MRI'], axis=1)

# List of acquisition types
acq_list = ['dwi', 'multi', 'func']

# Loop through each acquisition type
for acq in acq_list:
    # Define specific directories
    graphdir = os.path.join(basedir, acq, session, 'graph')
    outdir = os.path.join(basedir, 'multi', session, 'output') if acq == 'multi' else None
    
    if acq == 'multi':
        # Load subnetwork means for multi
        df_temp = pd.read_csv(os.path.join(outdir, f'CORE_ML_{session}_system.csv'), index_col=['subjID'])
        df_temp.index.name = 'Subj'
    else:
        # Load graph data for dwi and func
        df_temp = pd.read_csv(os.path.join(graphdir, f'graph_PC_subnetworks_acq-{acq}_CORE.csv'))
        df_temp = pd.melt(df_temp, id_vars=['subjID'], var_name='network', value_name='Y').set_index(['subjID'])
        df_temp['network'] = df_temp['network'] + '_' + acq
        df_temp.index.name = 'Subj'
    
    # Join with clinical data
    df_data = df_temp.join(df_clin).drop(['posttrial_sev','Med'], axis=1)
    
    # Rename columns and apply necessary transformations
    df_data = df_data.rename({'seks': 'sex'}, axis=1)
    df_data['perc_improv'] = df_data['perc_improv'] * 100
    df_data['responder'] = df_data['responder'].replace({0: 'znon-resp', 1: 'aresp'}).astype('category')

    if acq == 'multi':
        # Specific changes for multi acquisition
        df_data['treatment'] = df_data['treatment'].replace({'DGT-EMDR': 'DGT', 'SFT-IMRS': 'SFT'}).astype('category')

    # Round the numerical data
    df_data = df_data.round(4)

    # Save the processed data to a text file
    df_data.to_csv(os.path.join(statsdir, f'RBA_input_acq-{acq}_CORE_upd.txt'), sep=' ')

# Note: Merging data from different acquisitions can be done separately if required.
# Example:
# df_PC_dwi = pd.read_csv(os.path.join(statsdir, 'RBA_input_acq-dwi_CORE_upd.txt'), sep=' ', index_col=['Subj'])
# df_PC_func = pd.read_csv(os.path.join(statsdir, 'RBA_input_acq-func_CORE_upd.txt'), sep=' ', index_col=['Subj'])
# df_PC_all = pd.concat([df_PC_dwi, df_PC_func])

