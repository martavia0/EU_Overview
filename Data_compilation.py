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
path_py_wdws ="C:/Users/maria/Documents/Marta Via/1. PhD/F. Scripts/Python Scripts"
path_data_mac="/Users/martaviagonzalez/Documents/GitHub/EU_Overview/Data/All/"
path_data_wdws="C:/Users/maria/Documents/GitHub/EU_Overview/Data/All/"
path_folder_mac="/Users/martaviagonzalez/Documents/GitHub/EU_Overview/Data/"
path_folder_wdws = "C:/Users/maria/Documents/GitHub/EU_Overview/Data/"
#
mac_or_wdws = 'wdws' #Introduce OS here
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
all_files=glob.glob(path_data+'*Composition.txt')
chem_comp=pd.DataFrame()
chem_comp['Chemical_composition']=[pd.read_csv(i, sep='\t') for i in all_files]
# chem_comp =pd.DataFrame(data=[pd.read_csv(i, sep='\t') for i in all_files], columns=['Chemical_composition'])
li_site_names = [j[-28:-25] for j in all_files]
#%%Tweaks on the li_site_names
# li_site_names[10]="BO" 
# li_site_names = ['HPB', 'MI', 'CAO-NIC', 'IPR', 'NOA', 'INO', 'CGR', 'PD', 'ATOLL', 'SMR','BO', 'CMN'] #for mac ordering
li_site_names = ['ATOLL','BO', 'CAO-NIC', 'CGR', 'CMN', 'HPB', 'INO', 'IPR', 'MI', 'NOA', 'PD', 'SMR'] #for wdws ordering
chem_comp.index=li_site_names
#%% We import metadatafiles
os.chdir(path_folder)
metadata = pd.read_csv("Sites_metadata.txt", sep='\t')
metadata=metadata.sort_values('Acronym')
#%% We reorder the chem_comp df with the same order as in metadata.
chem_comp=chem_comp.reindex(metadata['Acronym']) 
chem_comp.drop(labels='SMR', inplace=True) #We remove hyttyala at the moment
#%% Average plots compound - Boxplots
os.chdir(path_folder + "Preliminar Plots/Chemical Composition/")
comp='Org'
li_dfs=[]
chem_compound,chem_compound_t, chem_dt=pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
for i in range(0,len(chem_comp)):
    df1 = pd.DataFrame(chem_comp.iloc[i][0])
    df1.reset_index(inplace=True, drop=True)
    df1['datetime']=pd.to_datetime(df1['Time (UTC)'], dayfirst=True)
    chem_compound=pd.concat([chem_compound, df1[comp]], axis=1)
    chem_dt = pd.concat([chem_dt, pd.to_datetime(df1['Time (Local)'], dayfirst=True)], axis=0)
    chem_compound_t=pd.concat([chem_compound_t, df1[comp]], axis=0,)
    df1.index=df1.datetime
    li_dfs.append(df1)
chem_compound.columns=chem_comp.index
chem_compound_t.columns=['All sites']
chem_dt.columns=['All times']
chem_compound_t['Hour'] = chem_dt['All times'].dt.hour
chem_compound_t['Month'] = chem_dt['All times'].dt.month
chem_compound_2=chem_compound.copy(deep = True)
mask_urban=(metadata['Type']!='UB')
non_urban=metadata['Acronym'].loc[mask_urban].to_list()
chem_compound_2[non_urban]=np.nan
#%%Plotting boxplots
fig, axs=plt.subplots(figsize=(10,10), nrows=2, ncols=2)
chem_compound_t.boxplot(ax=axs[0,0],column='All sites', fontsize=12,boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
chem_compound.boxplot(ax=axs[0,1], boxprops=bp, fontsize=10, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False, rot=90)
# chem_compound_2.boxplot(ax=axs[0,1], boxprops=bp_grey, fontsize=10, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False, rot=90)
chem_compound_t.boxplot(by='Hour', column='All sites', fontsize=10,  ax=axs[1,0], boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
chem_compound_t.boxplot(by='Month', column='All sites', fontsize=10,  ax=axs[1,1], boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
axs[1,0].set_title('')
axs[1,1].set_title('')
fig.suptitle(comp, fontsize=16)
axs[0,1].legend(['Urban site', 'Non-Urban site'], loc='upper right', labelcolor=['k', 'grey'])
plt.savefig('Boxplots_'+comp+'.png')
li_dfs=[j.drop('datetime', inplace=True, axis=1) for j in li_dfs]

#%% Filtering process
nr_comp=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
DL=pd.Series([0.148,0.024,0.012,0.284,0.011], index=nr_comp)
max_conc=pd.Series([30,20,20,10,5], index=nr_comp)

for j in range(0, len(li_dfs)):
    df=li_dfs[i]
    df=df.reindex(nr_comp, axis=1) 
    DL_mask = (df['Org']<=-DL['Org']) & (df['SO4']<=-DL['SO4']) & (df['NO3']<=-DL['NO3']) & (df['NH4']<=-DL['NH4'])& (df['Chl']<=-DL['Chl'])
    M_mask = (df['Org']<=-max_conc['Org']) & (df['SO4']<=-max_conc['SO4']) & (df['NO3']<=-max_conc['NO3']) & (df['NH4']<=-max_conc['NH4'])& (df['Chl']<=-max_conc['Chl'])
    df_f = df.loc[DL_mask]
    df_f = df_f.loc[M_mask]
    df_f.plot(subplots=True) 
#%% Maybe now I should put them in teh daily pattern
min_date = pd.date_range(start='01/01/2010', end = '31/12/2023', freq='D').strftime('%d/%m/%Y') #The complete time series
li_days=[]
for i in range(0,len(chem_comp)):
    df1=pd.DataFrame(chem_comp.iloc[i][0])
    df1.reset_index(inplace=True, drop=True)
    df1['datetime']=pd.to_datetime(df1['Time (UTC)'], dayfirst=True) 
    df1['date']=df1['datetime'].dt.date
    df1d=df1.groupby(df1['date']).mean()
    df1d['datetime']=pd.to_datetime(df1d.index, dayfirst=True)
    df1d.columns=['Chl_'+li_site_names[i],'NH4_'+li_site_names[i], 'NO3_'+li_site_names[i],'Org_'+li_site_names[i], 'SO4_'+li_site_names[i], 'datetime_'+li_site_names[i] ]
    li_days.append(df1d) #List of the datetimes
    

#%% Merging DFS
df=pd.DataFrame(pd.to_datetime(min_date, dayfirst=True), columns=['date'])
dates = pd.merge(df,li_days[0], how='outer',left_on='date', right_index=True, sort=True)
for i in range (1,len(li_days)):
    dates = pd.merge(dates,li_days[i], how='outer',left_on='date', right_index=True, sort=True)
dates=dates.drop(dates.index[-1])
dates.index=dates['date']
dates.plot(legend=False)
#%% PLot
dates_plot=dates.loc[:, dates.columns.str.startswith('datetime')]
colors=['grey', ]
dates_plot=dates_plot.notnull().astype('int')
dates_plot=dates_plot.replace(0, np.nan)
for i in range(0,len(dates_plot.columns)):
    dates_plot[dates_plot.columns[i]]=dates_plot[dates_plot.columns[i]]*(i+1)
fig, axs = plt.subplots(figsize=(9,6))
for m, col, c in zip(li_marker, dates_plot.columns, li_color):
    dates_plot[col].plot(marker=m, lw=0,legend=False, ax=axs, color=c, grid=True)
axs.set_yticks(range(0,len(dates_plot.columns)+1))
axs.set_yticklabels(['']+li_site_names[:-1])
axs.set_xlabel('Time (years)', fontsize=14)
axs.set_ylabel('Sites', fontsize=14)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='green', label='Regional background'), 
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   Line2D([0], [0],  color='sienna', label='Mountain'), 
                   Line2D([0], [0],  color='mediumpurple', label='Coastal'),
                   Line2D([0], [0], color='sienna', label='Mountain'), 
                   Line2D([0], [0], marker='D', color='grey', label='AMS'),
                   Line2D([0], [0], marker='s', color='grey', label='Q-ACSM'),
                   Line2D([0], [0], marker='o', color='grey', label='ToF-ACSM')]
axs.legend(handles=legend_elements, loc = (1.02,0.5))#'upper right')
#Still to do: plot by type of instrument  or type of site. 
fig.savefig('Site_availability.png')

#%%
for m, col in zip('xosd', df):
    df[col].plot(marker=m)
plt.legend()

li_marker=[]
li_color=[]
for i in range(0,len(metadata[:-1])):
    if metadata['Type.1'].iloc[i]=='AMS':
        li_marker.append('D')
    if metadata['Type.1'].iloc[i]=='Q':
        li_marker.append('s')
    if metadata['Type.1'].iloc[i]=='ToF':
        li_marker.append('o')
    #
    if metadata['Type'].iloc[i]=='UB':
        li_color.append('royalblue')
    if metadata['Type'].iloc[i]=='RB':
        li_color.append('green')
    if metadata['Type'].iloc[i]=='C':
        li_color.append('mediumpurple')
    if metadata['Type'].iloc[i]=='SU':
        li_color.append('darkorange')
    if metadata['Type'].iloc[i]=='M':
        li_color.append('sienna')



