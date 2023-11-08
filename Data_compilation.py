#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 11:00:54 2023

@author: martaviagonzalez
"""
#%%import pandas as pd
import numpy as np
import glob
import os as os
import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy import stats
#%% paths definition
path_py_mac="/Users/martaviagonzalez/Documents/Documents - MacBook Pro de MVIA/GitHub/All_Treatment"
path_py_wdws =""
path_data_mac="/Users/martaviagonzalez/Documents/GitHub/EU_Overview/Data/All/"
path_data_wdws="C:/Users/maria/Documents/GitHub/EU_Overview/Data/All/"
path_folder_mac="/Users/martaviagonzalez/Documents/GitHub/EU_Overview/Data/"
path_folder_wdws = "C:/Users/maria/Documents/GitHub/EU_Overview/Data/"
#
mac_or_wdws = 'mac' #Introduce OS here
#
if mac_or_wdws=='mac':
    path_py = path_py_mac
    path_data = path_data_mac
    path_folder = path_folder_mac
else:
    path_py = path_py_wdws 
    path_data = path_data_wdws 
    path_folder = path_folder_wdws
#%% Plots definitions
bp = dict(linestyle='-', linewidth=0.6)
bp_grey =dict(linestyle='-', linewidth=0.6, color='gainsboro')
mp = dict(marker='o', linewidth=0.6,markeredgecolor='black', markerfacecolor='black')
mdp = dict(color='k',  linewidth=0.6)
wp = dict(linestyle='-', linewidth=0.6)
#%% Import treatment
os.chdir(path_py)
from Treatment import *
trt = Basics(5)
trt.Hello()
print(trt.x)
#%% Import Composition files
os.chdir(path_data)
all_files=glob.glob(path_data_mac+'*Composition.txt')
chem_comp =pd.DataFrame(data=[pd.read_csv(i, sep='\t') for i in all_files], columns=['Chemical_composition'])
li_site_names = [j[-28:-25] for j in all_files]
#%%Tweaks on the li_site_names
# li_site_names[10]="BO" 
li_site_names = ['HPB','MI', 'CAO-NIC', 'IPR', 'NOA', 'INO', 'CGR', 'PD', 'ATOLL', 'SMR','BO', 'CMN']
chem_comp.index=li_site_names
#%% We import metadatafiles
os.chdir(path_folder)
metadata = pd.read_csv("Sites_metadata.txt", sep='\t')
#%% We reorder the chem_comp df with the same order as in metadata.
chem_comp=chem_comp.reindex(metadata['Acronym']) 
chem_comp.drop(labels='SMR', inplace=True) #We remove hyttyala at the moment
#%% Average plots compound
os.chdir(path_folder_mac + "Preliminar Plots/Chemical Composition/")
comp='Chl'
chem_compound,chem_compound_t, chem_dt=pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
for i in range(0,len(chem_comp)):
    df1 = pd.DataFrame(chem_comp.iloc[i][0])
    df1.reset_index(inplace=True, drop=True)
    chem_compound=pd.concat([chem_compound, df1[comp]], axis=1)
    chem_dt = pd.concat([chem_dt, pd.to_datetime(df1['Time (Local)'], dayfirst=True)], axis=0)
    chem_compound_t=pd.concat([chem_compound_t, df1[comp]], axis=0,)
chem_compound.columns=chem_comp.index
chem_compound_t.columns=['All sites']
chem_dt.columns=['All times']
chem_compound_t['Hour'] = chem_dt['All times'].dt.hour
chem_compound_t['Month'] = chem_dt['All times'].dt.month
chem_compound_2=chem_compound.copy(deep = True)
# chem_compound_2[['NOA', 'BO', 'CAO-NIC', 'PD', 'MI']].fillna(0, inplace=True)
chem_compound_2[['NOA','BO', 'CAO-NIC', 'PD', 'MI']]=np.nan#.fillna({1:0}, inplace=True)
# df[1].fillna(0, inplace=True)
#Plotting boxplots
fig, axs=plt.subplots(figsize=(10,10), nrows=2, ncols=2)
chem_compound_t.boxplot(ax=axs[0,0],column='All sites', fontsize=12,boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
chem_compound.boxplot(ax=axs[0,1], boxprops=bp, fontsize=10, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False, rot=90)
chem_compound_2.boxplot(ax=axs[0,1], boxprops=bp_grey, fontsize=10, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False, rot=90)
chem_compound_t.boxplot(by='Hour', column='All sites', fontsize=10,  ax=axs[1,0], boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
chem_compound_t.boxplot(by='Month', column='All sites', fontsize=10,  ax=axs[1,1], boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
axs[1,0].set_title('')
axs[1,1].set_title('')
fig.suptitle(comp, fontsize=16)
axs[0,1].legend(['Urban site', 'Non-Urban site'], loc='upper right', labelcolor=['k', 'grey'])
plt.savefig('Boxplots_'+comp+'.png')
#%%


























