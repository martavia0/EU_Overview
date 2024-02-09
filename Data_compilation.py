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
path_data_mac="/Users/martaviagonzalez/Documents/Documents - MacBook Pro de MVIA/IDAEA-CSIC/Overview/All/"
path_data_wdws="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/Overview/All/"
path_individual_wdws="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/Overview/Individual_plots/"
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
li_site_names = [j[-29:-25] for j in all_files]
chem_comp=pd.DataFrame()
chem_comp['Chemical_composition']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_files]
li_sites_names = ['ATOLL', 'BAQS', 'BCN', 'BIR', 'BO', 'CAO-NIC', 'CGR', 'CMN', 'CRE', 
                  'CRP', 'DEM', 'DUB', 'FKL','GRA', 'HEL','HPB', 'HTM', 'INO', 'IPR',  
                  'KOS', 'KRK', 'LON-MR', 'LON-NK', 'LYO', 'MAD-CIE', 'MAG','MAQS', 'MAR-LCP', 
                  'MEL', 'MET', 'MH', 'MI', 'MSC', 'MSY', 'NOA', 'PAR-GEN', 'PAR-BPE', 'PAR-HAL',
                  'PD','POI', 'PRE', 'PRG-SUCH','PUY', 'REN', 'RUG',  'SIRTA', 'SPC',  'STR',
                  'TAL','TAR','VIR', 'VLN', 'ZEP', 'ZUR'] 
chem_comp.index = li_sites_names
#%% We import metadatafiles
os.chdir(path_folder)
metadata = pd.read_csv("Sites_metadata.txt", sep='\t')
metadata=metadata.sort_values('Acronym')
#chem_comp.index = metadata['Acronym']
#%% Colors
cl_nrpm1=['forestgreen', 'red', 'blue', 'gold', 'fuchsia']
#%%
""" NR-PM1 COMPOUNDS! """

#%% INDIVIDUAL PLOTS I
nr_dfs=[]
chem_comp['Chemical_composition']=[pd.read_csv(i, sep='\t', na_values='np.nan', keep_default_na=True) for i in all_files]
os.chdir(path_individual_wdws)
for i in range(0,len(chem_comp)):
    print(i,li_sites_names[i])
    dfi=chem_comp.iloc[i][0]
    dfi=dfi.replace('', np.nan)
    dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed')
    dfi.drop(dfi['Time (UTC)'], errors='ignore')
    dfi.drop(dfi['Time (Local)'], errors='ignore')
    dfi.index = dfi['datetime']
    del dfi['datetime']
    dfi=dfi[['Org', 'SO4', 'NO3', 'NH4', 'Chl']]
    dfi.columns=['OA', 'Sulphate', 'Nitrate', 'Ammonium', 'Chloride']
    fig, axs=plt.subplots(figsize=(10,10))
    dfi.plot(subplots=True, color=cl_nrpm1, ax=axs, title =li_sites_names[i])
    nr_dfs.append(dfi)
    plt.savefig(li_sites_names[i]+'_NRPM1.png')
#%% INDIVIDUAL PLOTS II
nr_pm1=['OA', 'Sulphate', 'Nitrate', 'Ammonium', 'Chloride']
for i in range(0,len(nr_dfs)):
    dfi = nr_dfs[i].copy(deep=True)
    dfi['datetime']=pd.to_datetime(dfi.index, dayfirst=True)
    dfi_mean=dfi.mean(numeric_only=True)
    if dfi_mean.lt(0.0).any() :
        dfi_mean = abs(dfi_mean)
    dfi['Month']= dfi['datetime'].dt.month
    dfi['Year']=dfi['datetime'].dt.year
    dfi['Hour']=dfi['datetime'].dt.hour
    fig, axs=plt.subplots(ncols = 4, figsize=(12,3))
    dfi_mean.plot.pie(ax=axs[0], legend=False, autopct='%1.0f%%', colors=cl_nrpm1, title='Avg. conc.', startangle=90, labels=None, counterclock=True)
    dfi.groupby('Year').mean().plot( y=nr_pm1, ax=axs[1], legend=False, color=cl_nrpm1, marker = 'o')
    dfi.groupby('Month').mean().plot(marker='o', y=nr_pm1,  ax=axs[2], legend=False, color=cl_nrpm1)
    dfi.groupby('Hour').mean().plot( y=nr_pm1, ax=axs[3], legend=False, color=cl_nrpm1)
    plt.suptitle(li_sites_names[i])
    plt.savefig(li_sites_names[i]+'_NRPM1_means.png')


#%% Average plots compound - Boxplot
os.chdir(path_folder + "Preliminar Plots/Chemical Composition/")
comp='Org'
li_dfs=[]
chem_compound,chem_compound_t, chem_dt=pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
for i in range(0,len(chem_comp)):
    print(i, li_sites_names[i])
    df1 = pd.DataFrame(chem_comp.iloc[i][0])
    df1.reset_index(inplace=True, drop=True)
    df1['datetime']=pd.to_datetime(df1['Time (UTC)'], dayfirst=True, format='mixed')
    chem_compound=pd.concat([chem_compound, df1[comp]], axis=1)
    chem_dt = pd.concat([chem_dt, pd.to_datetime(df1['Time (UTC)'], dayfirst=True, format='mixed')], axis=0)
    chem_compound_t=pd.concat([chem_compound_t, pd.to_numeric(df1[comp])], axis=0)
    df1.index=df1.datetime
    li_dfs.append(df1)
chem_compound.columns=chem_comp.index
chem_compound_t.columns=['All sites']#­metadata['Acronym']
chem_dt.columns=['All times']

#%%
chem_dt['All_times'] = pd.to_datetime(chem_dt['All times'], dayfirst=True, utc=True)
chem_compound_t['Hour'] = chem_dt['All_times'].dt.hour
chem_compound_t['Month'] = chem_dt['All_times'].dt.month
chem_compound_2=chem_compound.copy(deep = True)
mask_nonurban=(metadata['Type']!='UB')
non_urban=metadata['Acronym'].loc[mask_nonurban].to_list()
chem_compound_2[non_urban]=np.nan
#%%Plotting boxplots
os.chdir(path_folder + "Preliminar Plots/Chemical Composition/")
fig, axs=plt.subplots(figsize=(10,10),nrows=2, ncols=2, width_ratios=[1,2])
chem_compound_t.boxplot(column = 'All sites', ax=axs[0,0], fontsize=12,boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
chem_compound.boxplot(ax=axs[0,1], boxprops=bp, fontsize=10, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False, rot=90)
# chem_compound_2.boxplot(ax=axs[0,1], boxprops=bp_grey, fontsize=10, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False, rot=90)
chem_compound_t.boxplot(by='Hour', column='All sites', fontsize=10,  ax=axs[1,0], boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
chem_compound_t.boxplot(by='Month', column='All sites', fontsize=10,  ax=axs[1,1], boxprops=bp, medianprops=mdp,meanprops=mp, whiskerprops=wp, showfliers=False)
axs[1,0].set_title('')
axs[1,1].set_title('')
fig.suptitle(comp, fontsize=16)
axs[0,1].legend(['Urban site', 'Non-Urban site'], loc='upper right', labelcolor=['k', 'grey'])
plt.savefig('Boxplots_'+comp+'.png')
# li_dfs=[j.drop('datetime', inplace=True, axis=1) for j in li_dfs]


#%% Maybe now I should put them in the daily pattern
min_date = pd.date_range(start='01/01/2009', end = '31/12/2023', freq='D').strftime('%d/%m/%Y') #The complete time series
li_days, li_nr, li_dfs=[], [], []
for i in range(0,len(chem_comp)):
    print(i, li_sites_names[i], li_site_names[i])
    df1=pd.DataFrame(chem_comp.iloc[i][0])
    df1.reset_index(inplace=True, drop=True)
    df1['datetime']=pd.to_datetime(df1['Time (UTC)'], dayfirst=True, format='mixed') 
    df1['date']=df1['datetime'].dt.date
    df1['Org'].astype(float)
    df1=df1.drop(columns=['Time (Local)', 'Time (UTC)', 'MSA', 'Seasalt', 'datetime'], axis=1, errors='ignore')
    df1d=df1.groupby(by=df1['date']).mean(numeric_only=True)
    df1d['datetime']=pd.to_datetime(df1d.index, dayfirst=True)
    df1d.columns=['Chl_'+li_sites_names[i],'NH4_'+li_sites_names[i], 'NO3_'+li_sites_names[i],'Org_'+li_sites_names[i], 'SO4_'+li_sites_names[i], 'datetime_'+li_sites_names[i] ]
    li_dfs.append(df1d)
    nr=pd.DataFrame({li_sites_names[i]: df1d.drop('datetime', errors='ignore').sum(axis=1, numeric_only=True)})
    li_days.append(df1d) #List of the datetimes
    li_nr.append(nr)
#%% Merging DFS
df=pd.DataFrame(pd.to_datetime(min_date, dayfirst=True), columns=['date'])
li_days[0]['datet'] = pd.to_datetime(li_days[0].index)
dates = pd.merge(df,li_days[0], how='outer',left_on='date', right_on='datet', sort=True)
for i in range (1,len(li_days)):
    li_days[i]['datet']= pd.to_datetime(li_days[i].index)
    dates = pd.merge(dates,li_days[i], how='outer',left_on='date', right_on='datet', sort=True, suffixes=('_'+li_sites_names[i], '_y'))
dates=dates.drop(dates.index[-1])
dates.index=dates['date']
# dates.plot(legend=False)

#%% PLot availability
metadata.drop(36, axis=0)
dates_plot=dates.loc[:, dates.columns.str.startswith('datetime')]
colors=['grey', ]

li_marker=[]
li_color=[]
metadata = metadata.sort_values('Acronym')
metadata.index =range(0,len(metadata))
for i in range(0,len(metadata)):
    print(i, metadata['Acronym'][i])
    if metadata['Type.1'].iloc[i]=='AMS' or metadata['Type.1'].iloc[i]=='Q-AMS' or  metadata['Type.1'].iloc[i]=='c-ToF-AMS' :
        li_marker.append('D')
    elif metadata['Type.1'].iloc[i]=='Q':
        li_marker.append('s')
    elif metadata['Type.1'].iloc[i]=='ToF':
        li_marker.append('o')
    #
    if metadata['Type'].iloc[i]=='UB':
        li_color.append('royalblue')
    elif metadata['Type'].iloc[i]=='RB':
        li_color.append('green')
    elif metadata['Type'].iloc[i]=='C':
        li_color.append('mediumpurple')
    elif metadata['Type'].iloc[i]=='SU':
        li_color.append('darkorange')
    elif metadata['Type'].iloc[i]=='M':
        li_color.append('sienna')
    elif metadata['Type'].iloc[i]=='A':
        li_color.append('hotpink')
    elif metadata['Type'].iloc[i]=='TR':
        li_color.append('darkcyan')
       
dates_plot=dates_plot.notnull().astype('int')
dates_plot=dates_plot.replace(0, np.nan)

li_color.append('royalblue')
li_marker.append('s')

li_sites_names.reverse()
# 
for i in range(0,len(dates_plot.columns)):
    # print(dates_plot.columns[i])
    dates_plot[dates_plot.columns[i]]=dates_plot[dates_plot.columns[i]]*54-(i)
fig, axs = plt.subplots(figsize=(8,12))
for m, col, c in zip(li_marker, dates_plot.columns, li_color):
    dates_plot[col].plot(marker=m, lw=0, legend=False, ax=axs, color=c, grid=True)
axs.set_yticks(range(0,len(dates_plot.columns)+1))
axs.set_yticklabels(['']+li_sites_names, fontsize=12)
axs.set_xticklabels(range(2008,2025,2), fontsize=14, rotation=0)

axs.set_xlabel('Time (years)', fontsize=14)
axs.set_ylabel('Sites', fontsize=14)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='green', label='Regional background'), 
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   Line2D([0], [0],  color='mediumpurple', label='Coastal'),
                   Line2D([0], [0], color='sienna', label='Mountain'), 
                   Line2D([0], [0], color='hotpink', label='Arctic'), 
                   Line2D([0], [0], color='darkcyan', label='Traffic'), 
                   Line2D([0], [0], marker='D', color='grey', label='AMS'),
                   Line2D([0], [0], marker='s', color='grey', label='Q-ACSM'),
                   Line2D([0], [0], marker='o', color='grey', label='ToF-ACSM')]
axs.legend(handles=legend_elements, loc = (1.03,0.65))#'upper right')
plt.title('Data availability')
fig.savefig('Site_availability.png')

#%%
li_sites_names.reverse()
# metadata=metadata.drop(35)
metadata=metadata.reset_index()
metadata['index']=metadata.index
li_sites_types=metadata['Type'].to_list()

sites=pd.DataFrame({'Name': li_sites_names, 'Type': li_sites_types})
#%%Plots all data for eaqch typeand in black the mean of them
comp = 'Org'
col_list=[]
types =['UB', 'RB', 'SU', 'C', 'M', 'A', 'TR']

ub, rb, su, c, m, a, tr= pd.DataFrame(), pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(), pd.DataFrame()
compound_dates=dates.loc[:, dates.columns.str.startswith(comp)]
for i in range(0,len(compound_dates.columns)):
    col_list = [col for col in compound_dates.columns if col.endswith(li_sites_names[i])]
    if li_sites_types[i]=='UB':
        ub=pd.concat([ub, compound_dates[col_list]], axis=1)
    if li_sites_types[i]=='RB':
        rb=pd.concat([rb, compound_dates[col_list]], axis=1)    
    if li_sites_types[i]=='SU':
        su=pd.concat([su, compound_dates[col_list]], axis=1) 
    if li_sites_types[i]=='C':
        c=pd.concat([c, compound_dates[col_list]], axis=1) 
    if li_sites_types[i]=='M':
        m=pd.concat([m, compound_dates[col_list]], axis=1)
    if li_sites_types[i]=='A':
        a=pd.concat([a, compound_dates[col_list]], axis=1)
    if li_sites_types[i]=='TR':
        tr=pd.concat([tr, compound_dates[col_list]], axis=1)
fig, axs=plt.subplots(nrows=7, sharex=True, figsize=(10,10))
ub.plot(ax=axs[0], legend=False)
ub.median(axis=1).plot(ax=axs[0], legend=False, c='k', lw=1.5)
rb.plot(ax=axs[1], legend=False)
rb.median(axis=1).plot(ax=axs[1], legend=False, c='k', lw=1.5)
su.plot(ax=axs[2], legend=False)
su.median(axis=1).plot(ax=axs[2], legend=False, c='k', lw=1.5)
c.plot(ax=axs[3], legend=False)
c.median(axis=1).plot(ax=axs[3], legend=False, c='k', lw=1.5)
m.plot(ax=axs[4], legend=False)
m.median(axis=1).plot(ax=axs[4], legend=False, c='k', lw=1.5)
a.plot(ax=axs[5], legend=False)
a.median(axis=1).plot(ax=axs[5], legend=False, c='k', lw=1.5)
tr.plot(ax=axs[6], legend=False)
tr.median(axis=1).plot(ax=axs[6], legend=False, c='k', lw=1.5)
for i in range(0,len(types)):
    axs[i].set_ylabel(types[i])


#%%
fig, axs=plt.subplots( figsize=(5,5))
ub.median().plot(kind='kde', ax=axs, lw=3, color='royalblue')
rb.median().plot(kind='kde', ax=axs, lw=3, c ='green')
su.median().plot(kind='kde', ax=axs, lw=3, c='darkorange')
# c.mean().plot(kind='kde', ax=axs)
m.median().plot(kind='kde', ax=axs, lw=3, color = 'sienna' )
# a.median().plot(kind='kde', ax=axs)
tr.mean().plot(kind='kde', ax=axs, lw=3, color ='darkcyan' )
axs.legend(['UB', 'RB', 'SU', 'M', 'TR']) 
plt.title(comp, fontsize=15)
axs.set_xlabel('Concentration ($μg·m^{-3}$)', fontsize=14)
axs.set_ylabel('Frequency', fontsize=14)
# colors=['darkorange', 'royalblue', 'green', 'sienna', 'darkcyan', 'hotpink', 'mediumpurple'])

#%%
'''     SUPERATIONS!!    '''

limit_who_25_daily = 15
limit_who_15_yearly= 5

li_sup_daily, li_ndays = [],[]
for i in range(0,len(li_nr)):
    dfi=li_nr[i]
    mask=dfi.iloc[:,0]>=limit_who_25_daily
    sup_daily = dfi.loc[mask].count()[0]
    ndays=len(dfi)
    li_sup_daily.append(sup_daily)
    li_ndays.append(ndays)
sup_daily= pd.DataFrame(data={'Daily WHO superations':li_sup_daily, 'Nb days accounted':li_ndays})
sup_daily.index=li_sites_names
sup_daily['Percentage superations'] = 100*sup_daily['Daily WHO superations']/sup_daily['Nb days accounted']
sup_daily['Type']=li_sites_types
sup_daily['Type_int']=sup_daily['Type']
sup_daily['Type_int']=sup_daily['Type_int'].replace('UB',0).replace('RB', 1).replace('SU', 2).replace('C', 3).replace('M', 4).replace('A',5).replace('TR', 6)


colors=['royalblue','green', 'darkorange', 'mediumpurple', 'sienna', 'hotpink', 'darkcyan']
fig, axs = plt.subplots(figsize=(8,10), ncols=2, width_ratios=[3,2])
sup_daily=sup_daily.sort_values(by='Percentage superations')
sup_daily['Percentage superations'].plot(kind='barh', ax=axs[0], color=[colors[i] for i in sup_daily['Type_int']], zorder=3)
axs[0].set_ylabel('Percentage of WHO PM$_{2.5}$ daily thresholds superation', fontsize=13)
axs[0].set_xlim(0,100)
axs[0].grid(axis='x', zorder=0)
axs[0].set_title('NR-PM$_1$ concentration')

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='green', label='Regional background'), 
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   Line2D([0], [0],  color='mediumpurple', label='Coastal'),
                   Line2D([0], [0], color='sienna', label='Mountain'), 
                   Line2D([0], [0], color='hotpink', label='Arctic'), 
                   Line2D([0], [0], color='darkcyan', label='Traffic')]

axs[0].legend(handles=legend_elements, loc = (1.05,0.752))#'upper right')

sup_pie=sup_daily.groupby('Type').mean()
sup_pie=sup_pie.sort_values('Percentage superations', ascending=False)
sup_pie.plot.pie(y='Percentage superations', ax=axs[1], legend=False,autopct='%2.0f%%', labels=None,pctdistance=0.7,fontsize=12, 
                 startangle=90, counterclock=False, ylabel='', colors=['darkorange', 'royalblue', 'green', 'mediumpurple', 'darkcyan','sienna','hotpink', ])

#%%
from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="constrained",  figsize=(10,10))

gs = GridSpec(2, 2)
ax1 = fig.add_subplot(gs[:, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1:, -1])

sup_daily['Percentage superations'].plot(kind='barh', ax=ax1, color=[colors[i] for i in sup_daily['Type_int']], zorder=3)
ax1.set_ylabel('Percentage of WHO PM$_{2.5}$ daily thresholds superation', fontsize=13)
ax1.set_xlim(0,100)
ax1.grid(axis='x', zorder=0)
ax1.set_title('NR-PM$_1$ concentration')

legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='green', label='Regional background'), 
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   Line2D([0], [0],  color='mediumpurple', label='Coastal'),
                   Line2D([0], [0], color='sienna', label='Mountain'), 
                   Line2D([0], [0], color='hotpink', label='Arctic'), 
                   Line2D([0], [0], color='darkcyan', label='Traffic')]

ax1.legend(handles=legend_elements, loc = (0.4,0.02))#'upper right')

sup_pie=sup_daily.groupby('Type').mean()
sup_pie=sup_pie.sort_values('Percentage superations', ascending=False)
sup_pie.plot.pie(y='Percentage superations', ax=ax2, legend=False,autopct='%2.0f%%', labels=None,pctdistance=0.7,fontsize=15, 
                 startangle=90, counterclock=False, ylabel='', colors=['darkorange', 'royalblue', 'green', 'mediumpurple', 'darkcyan','sienna','hotpink', ])
sup_bp=sup_daily.sort_values(by = 'Type_int')

boxplot = sup_bp.boxplot(column=['Percentage superations'], by='Type', ax=ax3, showmeans=True,
                            boxprops=bp, medianprops=mdp, meanprops=mp, whiskerprops=wp) 
ax3.set_title('')
plt.suptitle('')
plt.show()
#%% By years and types of site.
li_year=[]
df_year=pd.DataFrame()
day_count=pd.DataFrame()
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
df_year.index=range(2010, 2024)
day_count.index=range(2010, 2024)
for i in range(0, len(li_nr)):
    a=li_nr[i].copy(deep=True)
    a['dt']=li_days[i].loc[:,li_days[i].columns.str.startswith('datetime')]
    a['date']=pd.to_datetime(a['dt'], dayfirst=True)
    a['Year']=a['date'].dt.year    
    a['Month'] = a['date'].dt.month
    a['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in a.date]
    mask=a.iloc[:,0]>=limit_who_25_daily
    b = a.loc[mask]
    df_year[li_sites_names[i]]=b.groupby('Year').count()['dt']
    day_count[li_sites_names[i]] = a.groupby('Year').count()['dt']
 
    
df_year=df_year.T
day_count=day_count.T
df_year['Type'], day_count['Type'] = li_sites_types, li_sites_types
dft_year=df_year.groupby('Type').sum().T
dayt_count=day_count.groupby('Type').sum().T
dfplot = 100*dft_year / dayt_count
dfplot.sort_values(by='Type', axis=1, ascending=False, inplace=True)
colors_types=['royalblue', 'darkcyan', 'orange', 'green', 'saddlebrown', 'purple', 'hotpink']

fig, axs =plt.subplots(figsize=(8,4))
dfplot.plot(legend=True,color=colors_types, ax=axs, marker='o', zorder=3)
axs.set_xlabel('Years', fontsize=12)
axs.set_ylabel('Percentage of days with \nsuperations (%)', fontsize=12)
axs.grid(axis='y', zorder=0)
axs.grid(axis='x', zorder=0)

#%% Per each site, proportion of each sesason per superation days
li_year=[]
df_year=pd.DataFrame()
day_count=pd.DataFrame()
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
# df_year.index=range(2010, 2024)
# day_count.index=range(2010, 2024)
for i in range(0, len(li_nr)):
    a=li_nr[i].copy(deep=True)
    a['dt']=li_days[i].loc[:,li_days[i].columns.str.startswith('datetime')]
    a['date']=pd.to_datetime(a['dt'], dayfirst=True)
    a['Year']=a['date'].dt.year    
    a['Month'] = a['date'].dt.month
    a['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in a.date]
    mask=a.iloc[:,0]>=limit_who_25_daily
    b = a.loc[mask]
    df_year[li_sites_names[i]]=b.groupby('Season').count()['dt']
    day_count[li_sites_names[i]] = a.groupby('Season').count()['dt']

fig, axs=plt.subplots(figsize=(4,9))
df_plt=100*df_year /day_count
df_plot=100*df_plt/df_plt.sum()
df_plot=df_plot.T
df_plot=df_plot[['DJF', 'MAM', 'JJA', 'SON' ]]
df_plot=df_plot.iloc[::-1]
df_plot.plot(kind='barh', stacked=True,ax=axs, legend=False, zorder=3,
             color=['royalblue', 'yellowgreen', 'gold', 'orange'])
axs.grid('y', zorder=2)
axs.set_xlabel('Percentage of superations per season (%)')
axs.set_ylabel('Site')

axs.legend(loc=(-0.15,-0.08), ncol=4)
#%% Average days per season with superation
fig, axs=plt.subplots(figsize=(4,4))
df_plot.mean().plot(kind='bar',  yerr=df_plot.std(), color='grey', ax=axs, zorder=2)
axs.grid(axis='y', zorder=0)
axs.set_ylabel('Percentage of days with \nsuperation per season (%)')

#%% Average days per season per type of site
li_sites_types.reverse()
df_plot['Types']=li_sites_types
dfp=df_plot.groupby('Types').mean()
dfp=dfp.iloc[::-1]

fig, axs=plt.subplots(figsize=(6,3))
dfp.plot(kind='bar', zorder=3, color = ['royalblue', 'yellowgreen', 'gold', 'orange'], ax=axs)
axs.grid(axis='y', zorder=1)
axs.legend(loc=(1.01, 0.5))
axs.set_ylabel('Percentage of days \nof superation (%)')
axs.set_xlabel('Types of site')
#%%Composition on superation days (NR-PM1)
count_sup, li_sup=[], []
df_sup, df_sup_count=pd.DataFrame(), pd.DataFrame()
for i in range(0, len(li_nr)):
    print(i, li_sites_names[i])
    a=li_dfs[i].copy(deep=True)
    a['NR']=li_nr[i].iloc[:,0]
    a['dt']=li_days[i].loc[:,li_days[i].columns.str.startswith('datetime')]
    a['date']=a['dt'].dt.date
    b=a.groupby(a['date']).mean()
    mask =b['NR']>=limit_who_25_daily
    c = b.loc[mask].mean(axis=0, numeric_only=True)
    d = b.loc[mask].count(axis=0, numeric_only=True)
    c.drop('NR', inplace=True)
    c.index=['Chl', 'NH4', 'NO3', 'OA', 'SO4']
    li_sup.append(c)
    count_sup.append(100*d['NR']/len(b))
df_sup = pd.DataFrame(li_sup)
df_sup_count=pd.Series(count_sup)
df_sup.index = li_sites_names
df_sup = df_sup[['OA', 'SO4', 'NO3', 'NH4', 'Chl']]

color_nr = ['green', 'red', 'blue', 'gold', 'fuchsia']
fig, axs = plt.subplots(figsize=(10,4))
df_sup.plot(kind='bar', stacked=True, ax=axs, color = color_nr, zorder=7, width=0.9)
axs2=axs.twinx()
df_sup_count.plot(ax=axs2, marker='D', lw=0, color='k',zorder=8, markersize=3)
axs.legend(loc=(0.15, -0.5), ncols=5)

axs.set_ylabel('Absolute Concentration \n $(μg·m^{-3})$', fontsize=12)
axs2.set_ylabel('Percentage of days \nwith superation (%)', fontsize=12)
axs.set_xlabel('Site', fontsize=12)
fig.suptitle('Days with superation')
axs2.set_ylim(-2,100)
#%% Ordered by perc with superations
color_nr = ['green', 'red', 'blue', 'gold', 'fuchsia']
fig, axs = plt.subplots(figsize=(10,4))
df_sup_count.index=li_sites_names
df_sup['count']=df_sup_count
df_sup=df_sup.sort_values(by='count')
df_sup.plot(y=['OA', 'SO4', 'NO3', 'NH4', 'Chl'], kind='bar', stacked=True, ax=axs, color = color_nr, zorder=7, width=0.9)
axs2=axs.twinx()
df_sup['count'].plot(ax=axs2, marker='D', lw=0, color='k',zorder=8, markersize=3)
axs.legend(loc=(0.15, -0.5), ncols=5)
axs.set_ylabel('Absolute Concentration \n $(μg·m^{-3})$', fontsize=12)
axs2.set_ylabel('Percentage of days \nwith superation (%)', fontsize=12)
axs.set_xlabel('Site', fontsize=12)
fig.suptitle('Days with superation')
axs2.set_ylim(-2,100)
#%% In relative terms
fig, axs = plt.subplots(figsize=(10,4))
nr=['OA', 'SO4', 'NO3', 'NH4', 'Chl']
df_sup_count.index=li_sites_names
df_sup['count']=df_sup_count
df_sup=df_sup.sort_values(by='count')
df_sup['sum']=df_sup[nr].sum(axis=1)
df_sup_rel = pd.DataFrame({'OA':100*df_sup['OA']/df_sup['sum'], 'SO4':100*df_sup['SO4']/df_sup['sum'],
                           'NO3':100*df_sup['NO3']/df_sup['sum'], 'NH4':100*df_sup['NH4']/df_sup['sum'], 
                           'Chl':100*df_sup['Chl']/df_sup['sum']})
df_sup_rel.plot(y=nr, kind='bar', stacked=True, ax=axs, color = color_nr, zorder=7, width=0.9)
axs2=axs.twinx()
df_sup['count'].plot(ax=axs2, marker='D', lw=0, color='k',zorder=8, markersize=3)
axs.legend(loc=(0.15, -0.5), ncols=5)
axs.set_ylabel('Absolute Concentration \n $(μg·m^{-3})$', fontsize=12)
axs2.set_ylabel('Percentage of days \nwith superation (%)', fontsize=12)
axs.set_xlabel('Site', fontsize=12)
fig.suptitle('Days with superation')

axs2.set_ylim(-2,100)
#%% MAPPPPP
li_nr_avg=[]
for i in range(0,len(nr_dfs)):
    df=nr_dfs[i]
    if li_sites_names[i] =='ZEP':
        df['nr'] = df['OA']+df['Sulphate']+df['Nitrate']+df['Chloride']
    elif li_sites_names[i] =='FKL' or li_sites_names[i]=='MH' or li_sites_names[i]=='BIR':
        df['nr'] = df['OA']+df['Sulphate']+df['Nitrate']+df['Ammonium']
    else:
        df['nr'] = df['OA']+df['Sulphate']+df['Nitrate']+df['Ammonium']+df['Chloride']
    li_nr_avg.append(df['nr'].mean())
    print(i, li_sites_names[i], df['nr'].mean())
metadata['NR']=li_nr_avg


plt.rcParams['lines.markersize'] ** 2
import geopandas as gpd
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
fig, axs = plt.subplots(figsize=(10,8))
countries.head()
countries.plot(color="lightgrey", ax=axs)
md=metadata.copy(deep=True)
countries.boundary.plot(ax=axs, color='k', zorder=0)
mask_ub=md['Type']=='UB'
mask_rb=md['Type']=='RB'
mask_su=md['Type']=='SU'
mask_c=md['Type']=='C'
mask_m=md['Type']=='M'
mask_tr=md['Type']=='TR'
mask_a=md['Type']=='A'
md['NR2']=md['NR']*15

md.loc[mask_ub].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='royalblue', ax=axs,linewidths=2, zorder=4, label = 'Urban background')
md.loc[mask_rb].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='green', ax=axs,linewidths=2, zorder=6, label='Regional Background')
md.loc[mask_su].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='darkorange', ax=axs,linewidths=2,zorder=5, label = 'Suburban'  )
md.loc[mask_c].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='purple', ax=axs,linewidths=2, zorder=7, label = 'Coastal' )
md.loc[mask_m].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='saddlebrown', ax=axs,linewidths=2, zorder=8, label = 'Mountain' )
md.loc[mask_tr].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='darkcyan', ax=axs,linewidths=2, zorder=4, label = 'Traffic')
md.loc[mask_a].plot.scatter(x="Lon", y="Lat", s='NR2', c='hotpink', ax=axs, zorder=9,linewidths=2,label = 'Arctic' )

# plt.legend(*sc.legend_elements("sizes"))
axs.set_xlim(-20, 40)
axs.set_ylim(33, 75)
axs.set_ylabel('')
axs.set_xlabel('')
plt.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), labelspacing=1)
#%%
plt.rcParams['lines.markersize'] ** 2
import geopandas as gpd
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
antarct=countries[countries.name == 'Antarctica']

fig, axs = plt.subplots(figsize=(5,5))
# countries.head()
antarct.plot(color="lightgrey", ax=axs)
antarct.boundary.plot(ax=axs, color='k', zorder=0)
md.loc[mask_a].plot.scatter(x="Lon", y="Lat", s='NR2', c='hotpink', ax=axs, zorder=9,linewidths=2,legend=False)
# axs.set_xlim(-80,-40)
# axs.set_ylim(-90,-60 )
axs.set_ylabel('')
axs.set_xlabel('')
#%%
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.geometry import Point

crs = ccrs.SouthPolarStereo()
crs_proj4 = crs.proj4_init
w = world.to_crs(crs_proj4)

fig, ax = plt.subplots(subplot_kw=dict(projection=crs))
w.plot(ax=ax, facecolor='sandybrown', edgecolor='black')

# Show the plot
plt.show()
#%%
# Import libraries
import os
import matplotlib.pyplot as plt
import geopandas as gpd
import earthpy as et

# Get the data & set working dir
data = et.data.get_data('spatial-vector-lidar')
os.chdir(os.path.join(et.io.HOME, 'earth-analytics', "data"))



#%%
'''
**************** SOURCE APPORTIONMENT *********************'''
#%%
os.chdir(path_data)
# hpb=pd.read_csv('HPB_PMF_Output_TS.txt', sep='\t')

# hpb['datetime']=pd.to_datetime(hpb['PMF Time (UTC)'], dayfirst=True)#, format="%d/m/%Y %H:%M")
#%% IMPORT ALL FILES
os.chdir(path_data)
all_files=glob.glob(path_data+'*Output_TS.txt')
li_site_names_oa = [j[-21:-18] for j in all_files]
li_site_names_oa_all = [j[-31:-12] for j in all_files]

oa_sa=pd.DataFrame()
oa_sa['OA_SA']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_files]
li_sites_names_oa = ['ATOLL', 'BCN', 'BIR', 'BO', 'CAO-NIC', 'CGR', 'CMN', 'CRP', 'DEM', 'DUB', 
                      'GRA', 'HPB', 'INO', 'IPR',  'KOS', 'KRK', 'LON-MR', 'LON-NK', 'MAG', 'MAR-LCP',
                      'MEL','MI', 'PD', 'PRE', 'RUG', 'SIRTA',  'SMR', 'SPC', 'TAR','VLN', 'ZUR'] 
oa_sa.index = li_sites_names_oa
li_oasa=[]
for i in range(0,len(oa_sa)):
    print(i, li_sites_names_oa[i])
    oai=oa_sa.loc[li_sites_names_oa[i]][0]
    oai['datetime'] = pd.to_datetime(oai['PMF Time (UTC)'], dayfirst=True)
    li_oasa.append(oai)
#%% Factors adjustments

factors_names = pd.DataFrame([dfi.columns[2:-1] for dfi in oa_sa['OA_SA']])
factors_names.index=li_sites_names_oa
colors_oasa=factors_names.replace({'HOA': 'grey', 'COA': 'mediumpurple', 'Amine-OA': 'skyblue',
                                   'BBOA': 'saddlebrown', 'LO-OOA': 'yellowgreen','MO-OOA':'darkgreen', 
                                   'OOA': 'green', 'Total OOA':'green', 'OOA_BB': 'olivedrab', 'OOA_BBaq':'olive',
                                   'HOA1': 'grey', 'HOA2': 'dimgrey', 'CSOA': 'rosybrown', 'Wood': 'sienna', 
                                   'Peat': 'sienna', 'Coal': 'sienna', 'POA': 'darkkhaki', 'CCOA': 'sandybrown',
                                   '58-OA': 'hotpink', 'ShInd-OA': 'purple', 'seasalt':'darkcyan',
                                   'BBOA1': 'saddlebrown', 'BBOA2': 'saddlebrown', 'LOA':'yellow'})
colors_oasa.iloc[29] = pd.Series(['yellow', 'darkkhaki', 'saddlebrown', 'grey', 'darkgreen', 'yellowgreen', 'None'])
#%% SA: INDIVIDUAL PLOTS I
from matplotlib.patches import Rectangle
os.chdir(path_individual_wdws)
oa=[]
for i in range(0,len(oa_sa)):
    print(i, li_sites_names_oa[i])
    oai=oa_sa.loc[li_sites_names_oa[i]][0]
    colours= [x for x in colors_oasa.loc[li_sites_names_oa[i]].tolist() if x is not None]
    oai_m=oai.iloc[:,2:-1].mean()
    fig, axs=plt.subplots(figsize=(5,5))
    oai_m.plot.pie(autopct='%1.0f%%',colors=colours, fontsize=12, startangle=90, counterclock=False)
    fig.suptitle(li_sites_names_oa[i])
    axs.add_patch(Rectangle((-0.26, -0.22),0.6, 0.4, facecolor='white', fill=True)) 
    oa_total = oai[[x for x in factors_names.loc[li_sites_names_oa[i]].tolist() if x is not None]].sum(axis=1).mean()
    axs.text(x=-0.24, y=-0.13, s='OA = '+str(oa_total.round(1))+'\n $μg·m^{-3}$', fontsize=12)
    plt.savefig(li_sites_names_oa[i]+'_OASA_means.png')
#%% Transform into daily
min_date = pd.date_range(start='01/01/2009', end = '31/12/2023', freq='D').strftime('%d/%m/%Y') #The complete time series
li_days, li_oad=[], []
for i in range(0,len(oa_sa)):
    print(i, li_sites_names_oa[i])
    df1=pd.DataFrame(li_oasa[i])
    # df1.reset_index(inplace=True, drop=True)
    df1['date']=df1['datetime'].dt.date
    # df1=df1.drop(columns=['Time (Local)', 'Time (UTC)', 'MSA', 'Seasalt', 'datetime'], axis=1, errors='ignore')
    df1d=df1.groupby(by=df1['date']).mean(numeric_only=True)
    df1d['datetime']=pd.to_datetime(df1d.index, dayfirst=True)
    li_days.append(df1d) #List of the datetimes
#%% Importing Crit poll meas
os.chdir(path_data)
all_gm_files=glob.glob(path_data+'*meteo.txt')
li_site_names = [j[-28:-22] for j in all_gm_files]
print(li_site_names)
gm=pd.DataFrame()
gm['GM']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_gm_files]
li_names_gm = ['ATOLL','BAQS', 'BCN', 'DEM','HEL', 'LON-MR','LON-NK', 
            'MAD-CIE', 'MRS-LCP', 'NOA', 'PRG-SUCH', 'SIRTA', 'ZUR'] 
gm.index = li_names_gm
#%% Joining nrpm1, oasa, gases, meteo
li_names = li_sites_names
min_date = pd.date_range(start='01/01/2009', end = '31/12/2023', freq='H').strftime('%d/%m/%Y') #The complete time series
li_all=[]
for i in range(0,len(chem_comp)):
    print(li_sites_names[i])
    dfi=chem_comp.iloc[i][0].copy(deep=True)
    dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed', errors='coerce')
    dfi['datehour']=dfi['datetime'].dt.round('H')
    dfi=dfi.groupby('datehour').mean(numeric_only=True)
    if (li_sites_names[i] in oa_sa.index)==True:
        print('IN OASA!')
        oai=oa_sa.loc[li_names[i]][0].copy(deep=True)
        oai['datetime']=pd.to_datetime(oai['PMF Time (UTC)'] ,dayfirst=True, errors='raise')
        oai['datehour']=oai['datetime'].dt.round('H')
        oai=oai.groupby('datehour').mean(numeric_only=True)
        dfi=pd.merge(left=dfi, right=oai, how='outer',left_on='datehour', right_on='datehour')
    if (li_names[i] in gm.index) == True:
        print('IN GM!')
        gmi=pd.DataFrame()
        gmi=gm.loc[li_names[i]][0].copy(deep=True)
        gmi['datetime']=pd.to_datetime(gmi['TIME UTC, end'], dayfirst=True,errors='raise')
        gmi['datehour']=gmi['datetime'].dt.round('H')
        gmi=gmi.groupby('datehour').mean(numeric_only=True)
        dfi=pd.merge(left=dfi, right=gmi, how='outer',left_on='datehour', right_on='datehour')
    
    # dfi.to_csv(li_names[i]+'_ALL.txt')
    li_all.append(dfi)
#%% Composition plot
li_means=[]
for i in range(0,len(li_all)):
    print(li_names[i])
    dfi = li_all[i]
    dfi['PM1_sum']=dfi['Org']+dfi['SO4']+dfi['NO3']+dfi['NH4']+dfi['Chl']
    if li_names[i] in gm.index:
        print(li_names[i] + 'has BC')
        dfi['PM1_sum'] = dfi['PM1_sum']+dfi['BC(ng/m3)']/1000.0
        dfi['BC']=dfi['BC(ng/m3)']/1000.0
    li_means.append(dfi.mean(numeric_only=True))
means=pd.DataFrame(li_means)

means_plot=pd.DataFrame()
means_plot.index=['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC', 'HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA']
for i in range(0,len(li_all)):
    if li_names[i] in gm.index and li_names[i] not in oa_sa.index:
        toadd=means[['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC']].iloc[i]
        toadd.index=['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
    elif li_names[i] in oa_sa.index and li_names[i] not in gm.index:
        toadd=means[['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA1', 'HOA2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl']].iloc[i]
        toadd.index=['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA1', 'HOA2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
    elif li_names[i] in oa_sa.index and li_names[i] in gm.index:
        toadd=means[['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA1', 'HOA2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl', 'BC']].iloc[i]
        toadd.index=['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA1', 'HOA2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl', 'BC']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
    else:
        toadd=means[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].iloc[i]
        toadd.index=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)        

means_plot=means_plot.T
means_plot = means_plot[[ 'SO4', 'NO3', 'NH4', 'Chl', 'BC', 'Org', 'HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','OOA','ShInd-OA', 'CSOA', 'CCOA', 
             'Coal', 'Peat', 'Wood', 'HOA1', 'HOA2',  'OOA_BBaq', 'OOA_BB', 'Amine-OA',]]

colors_all = ['red', 'blue', 'gold', 'fuchsia', 'black','green', 'grey', 'mediumorchid', 'saddlebrown', 'lightgreen', 'darkgreen','green',
              'purple','gainsboro', 'hotpink', 'hotpink', 'tan', 'sandybrown', 'grey', 'grey', 'forestgreen', 'mediumseagreen', 'skyblue' ]

fig, axs=plt.subplots(figsize=(10,5))
means_plot.plot(kind='bar', stacked=True, color=colors_all, ax=axs, width=0.85)
axs.set_xlabel('Sites', fontsize=13)
axs.set_ylabel('Concentration ($μg·m^{-3}$)', fontsize=13)
plt.legend(loc=(1.02,-0.24))
#%% Same but in abs terms
means_plot_mean = 100.0*means_plot.T/means_plot.sum(axis=1)

fig, axs=plt.subplots(figsize=(10,5))
means_plot_mean.T.plot(kind='bar', stacked=True, color=colors_all, ax=axs, width=0.85)
axs.set_xlabel('Sites', fontsize=13)
axs.set_ylabel('Absolute Concentration ($μg·m^{-3}$)', fontsize=13)
axs.set_ylim(0,100)
plt.legend(loc=(1.02,-0.24))

#%% Boxplots per site
factors=['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA']
dfs=[]
fig, axs = plt.subplots(nrows=5, figsize=(8,12), constrained_layout = True)
for j in range(0,len(factors)):
    factor=factors[j]
    li_factor=[]
    li_factor_name=[]
    for i in range(0,len(li_all)):
        dfi = li_all[i]
        dfi=dfi.reset_index()
        if (factor in dfi.columns) and (dfi[factor].mean != np.nan):
            li_factor.append(dfi[factor])
            li_factor_name.append(li_names[i])
    df_factor=pd.DataFrame(li_factor)
    df_factor = df_factor.T
    df_factor.columns=li_factor_name
    
    df_factor.boxplot(showfliers=False, ax=axs[j],showmeans=True, boxprops=bp, whiskerprops=wp, medianprops=mdp, meanprops=mp)
    axs[j].set_xticks(range(0,len(li_factor_name)))
    
    axs[j].set_xticklabels(labels=['']+df_factor.columns[:-1].tolist(),  rotation=90)
    axs[j].set_ylabel(factor, fontsize=12)
#%% Composition per type of site
means_plot_c = means_plot.copy(deep=True)
means_plot_c['Type']=metadata['Type']
means_plot_c['OOA']=means_plot_c['OOA'] + means_plot_c['OOA_BBaq']+means_plot_c['OOA_BB'] 
means_plot_c['HOA']=means_plot_c['HOA1'] +means_plot_c['HOA2'] + means_plot_c['HOA']
means_plot_c['CCOA']=means_plot_c['CCOA'] +means_plot_c['Coal']

means_plot_c.drop(['OOA_BB', 'OOA_BBaq', 'Coal', 'HOA1', 'HOA2'], inplace=True, axis=1, errors='ignore')
comp_bysite=means_plot_c.groupby('Type').mean(numeric_only=True)

colors_all = ['red', 'blue', 'gold', 'fuchsia', 'black','limegreen', 'grey', 'mediumorchid', 'saddlebrown', 'lightgreen', 'darkgreen','green',
              'purple','gainsboro', 'hotpink', 'tan', 'sandybrown', 'skyblue' ]

comp_bysite=comp_bysite.iloc[::-1]
fig, axs=plt.subplots(figsize=(5,5))
comp_bysite.plot(kind='bar', stacked=True, color=colors_all, ax=axs)
axs.set_ylabel('Absolute Concentration ($μg·m^{-3}$', fontsize=12)
axs.set_xlabel('Type of site', fontsize=12)

plt.legend(loc=(1.05, -0.0))
#%% Time series of OA factors


#%% INDIVIDUAL PLOTS I: Time series
oa_dfs=[]
os.chdir(path_individual_wdws)
for i in range(0,len(chem_comp)):
    if li_names[i] in oa_sa.index:
        print(li_names[i])            
        dfi=li_all[i].copy(deep=True)
        dfi.drop('OA', inplace=True, axis=1)
        dfi['datetime']=pd.to_datetime(dfi.index, dayfirst=True)
        dfi.index = dfi['datetime']
        del dfi['datetime']
        factors=[col for col in dfi.columns if col.endswith('OA')]
        dfi=dfi[factors]
        dfi.columns=factors
        fig, axs=plt.subplots(figsize=(10,10))
        dfi.plot(subplots=True, ax=axs, title =li_names[i])
        oa_dfs.append(dfi)
        plt.savefig(li_sites_names[i]+'_factors_TS.png')
    #%% INDIVIDUAL PLOTS II: Cycles
for i in range(0,len(li_all)):
    if li_names[i] in oa_sa.index:
        dfi = pd.DataFrame(li_all[i]).copy(deep=True)
        print(li_names[i])
        factors=[col for col in dfi.columns if col.endswith('OA')]
        factors = factors[:-1]
        dfi.drop('OA', inplace=True, errors='ignore', axis=1)
        dfi['datetime']=pd.to_datetime(dfi.index, dayfirst=True)
        dfi['Month']= dfi['datetime'].dt.month
        dfi['Year']=dfi['datetime'].dt.year
        dfi['Hour']=dfi['datetime'].dt.hour
        dfi_mean = dfi[factors].mean(numeric_only=True)
        fig, axs=plt.subplots(ncols = 4, figsize=(12,3))
        dfi_mean.plot.pie(ax=axs[0], legend=False, autopct='%1.0f%%',  title='Avg. conc.', startangle=90, labels=None, counterclock=True)
        dfi.groupby('Year').mean(numeric_only=True).plot( y=factors, ax=axs[1],   marker = 'o')
        dfi.groupby('Month').mean(numeric_only=True).plot(marker='o', y=factors,  legend=False, ax=axs[2], )
        dfi.groupby('Hour').mean(numeric_only=True).plot( y=factors, ax=axs[3], legend=False)
        plt.suptitle(li_sites_names[i])
        plt.savefig(li_sites_names[i]+'Factor_means.png')
#%% INDIVIDUAL PLOTS III: Intercomp OA
interc = pd.DataFrame()
os.chdir(path_individual_wdws)
for i in range(0,len(li_all)): 
    fig, axs=plt.subplots(figsize=(5,5))
    dfi =li_all[i]
    cols=[col for col in dfi.columns if col.endswith('OA')]
    dfi['OA']=dfi[cols].sum(numeric_only=True, axis=1)
    pl=dfi.plot(x='Org', y='OA', kind='scatter',  ax=axs, c=dfi.index, title=li_names[i])
    axs.set_xlim(0,dfi['Org'].quantile(0.9))
    axs.set_ylim(0,dfi['OA'].quantile(0.9))
    axs.text(x=0,y=0, s='y = '+str(trt.slope(dfi['OA'], dfi['Org'])[0])+' · x + '+ 
             str(trt.slope(dfi['OA'], dfi['Org'])[1])+'\nR$^2$ = '+str(trt.R2(dfi['OA'], dfi['Org'])))
    plt.savefig(li_names[i]+'_OA_closure.png')

