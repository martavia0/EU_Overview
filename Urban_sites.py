# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 12:00:12 2024

@author: Marta Via
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
os.chdir(path_py)
from Treatment import *
trt = Basics(5)
#%% paths definition
path_py_wdws ="C:/Users/maria/Documents/Marta Via/1. PhD/F. Scripts/Python Scripts"
path_data_wdws="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/Overview/Selected/"
path_individual_wdws="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/Overview/Selected_plots/"
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
#%% Colors
nr_colors=['green', 'red', 'blue', 'gold', 'fuchsia']

#%% Import Composition files
os.chdir(path_data)
all_files=glob.glob(path_data+'*Composition.txt')
li_site_names = [j[-29:-25] for j in all_files]
chem_comp=pd.DataFrame()
chem_comp['Chemical_composition']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_files]
li_names = ['ATOLL', 'BAQS', 'BCN', 'BO', 'CAO-NIC', 'DEM', 'DUB', 'HEL','INO','KRK', 'LON-MR', 'LON-NK', 'MAD-CIE',
             'MAQS', 'MAR-LCP',  'NOA', 'PRG-SUCH', 'SIRTA', 'TAR', 'ZUR'] 
chem_comp.index = li_names
#%% We import metadatafiles
os.chdir(path_folder)
metadata = pd.read_csv("Sites_metadata_selected.txt", sep='\t')
metadata=metadata.sort_values('Acronym')

#%% Dayly index now
min_date = pd.date_range(start='01/01/2009', end = '31/12/2023', freq='D').strftime('%d/%m/%Y') #The complete time series
li_days, li_nr, li_dfs=[], [], []
for i in range(0,len(chem_comp)):
    print(i, li_names[i])
    df1=pd.DataFrame(chem_comp.iloc[i][0])
    df1.reset_index(inplace=True, drop=True)
    df1['datetime']=pd.to_datetime(df1['Time (UTC)'], dayfirst=True, format='mixed') 
    df1['date']=df1['datetime'].dt.date
    df1['Org'].astype(float)
    df1=df1.drop(columns=['Time (Local)', 'Time (UTC)', 'MSA', 'Seasalt', 'datetime'], axis=1, errors='ignore')
    df1d=df1.groupby(by=df1['date']).mean(numeric_only=True)
    df1d['datetime']=pd.to_datetime(df1d.index, dayfirst=True)
    df1d.columns=['Chl_'+li_names[i],'NH4_'+li_names[i], 'NO3_'+li_names[i],'Org_'+li_names[i], 'SO4_'+li_names[i], 'datetime_'+li_names[i] ]
    li_dfs.append(df1d)
    nr=pd.DataFrame({li_names[i]: df1d.drop('datetime', errors='ignore').sum(axis=1, numeric_only=True)})
    li_days.append(df1d) #List of the datetimes
    li_nr.append(nr)
#%% Merging DFS
df=pd.DataFrame(pd.to_datetime(min_date, dayfirst=True), columns=['date'])
li_days[0]['datet'] = pd.to_datetime(li_days[0].index)
dates = pd.merge(df,li_days[0], how='outer',left_on='date', right_on='datet', sort=True)
for i in range (1,len(li_days)):
    li_days[i]['datet']= pd.to_datetime(li_days[i].index)
    dates = pd.merge(dates,li_days[i], how='outer',left_on='date', right_on='datet', sort=True, suffixes=('_'+li_names[i], '_y'))
dates=dates.drop(dates.index[-1])
dates.index=dates['date']
# dates.drop('date', inplace=True)
#%% PLot availability

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
    elif metadata['Type'].iloc[i]=='SU':
        li_color.append('darkorange')
    elif metadata['Type'].iloc[i]=='M':
        li_color.append('sienna')
    elif metadata['Type'].iloc[i]=='TR':
        li_color.append('grey')
dates_plot=dates_plot.notnull().astype('int')
dates_plot=dates_plot.replace(0, np.nan)
li_color.append('royalblue')

# li_names.reverse()
dates_plot.columns=li_names

for i in range(0,len(dates_plot.columns)):
    # print(dates_plot.columns[i])
    dates_plot[dates_plot.columns[i]]=dates_plot[dates_plot.columns[i]]*20-(i)
fig, axs = plt.subplots(figsize=(10,7))
for m, col, c in zip(li_marker, dates_plot.columns, li_color):
    dates_plot[col].plot(marker=m, lw=0, legend=False, ax=axs, color=c, grid=True)
axs.set_yticks(range(0,len(dates_plot.columns)+1))
axs.set_yticklabels(['']+li_names, fontsize=12)
# axs.set_xticklabels(range(2011,2025,2), fontsize=14, rotation=0)

axs.set_xlabel('Time (years)', fontsize=14)
axs.set_ylabel('Sites', fontsize=14)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   Line2D([0], [0], color='grey', label='Traffic'), 
                   Line2D([0], [0], marker='D', color='grey', label='AMS'),
                   Line2D([0], [0], marker='s', color='grey', label='Q-ACSM'),
                   Line2D([0], [0], marker='o', color='grey', label='ToF-ACSM')]
axs.legend(handles=legend_elements, loc = (1.03,0.45))#'upper right')
plt.title('Data availability')
fig.savefig('Site_availability.png')
#%% Site data organisation
# li_names.reverse()
li_dfi=[]
for i in range(0,len(chem_comp)):
    print(li_names[i])
    dfi = pd.DataFrame(chem_comp.iloc[i][0])
    dfi.drop('Time (Local)', inplace=True, errors='ignore')
    dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed')
    li_dfi.append(dfi)    
#%% Mean site composition
means= pd.DataFrame([i.mean(numeric_only=True) for i in li_dfi])
means=means[['Org', 'SO4', 'NO3', 'NH4', 'Chl']]
means.index=li_names
means['nr']=means['Org']+means['SO4']+means['NO3']+means['NH4']+ means['Chl']
means.sort_values(by='nr', inplace=True)
means.drop('nr',axis=1, inplace=True)
fig, axs=plt.subplots(figsize=(6,3))
means.plot(kind='bar', stacked=True, ax=axs, color=nr_colors, zorder=3)
axs.grid(axis='y', zorder=1)
axs.set_ylabel('Mean composition ($μg·m^{-3}$)')
axs.set_xlabel('Sites')
plt.savefig('ChemicalComposition_sites.png')
#%% Importing SA
os.chdir(path_data)
all_sa_files=glob.glob(path_data+'*Output_TS.txt')
li_site_names = [j[-23:-18] for j in all_sa_files]
oasa=pd.DataFrame()
oasa['SA']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True, infer_datetime_format=True) for i in all_sa_files]
li_names_sa = ['ATOLL','BCN', 'BO', 'CAO-NIC', 'DEM', 'DUB','INO','KRK', 
            'LON-MR', 'LON-NK', 'MAR-LCP', 'SIRTA', 'TAR', 'ZUR'] 
oasa.index = li_names_sa
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
min_date = pd.date_range(start='01/01/2009', end = '31/12/2023', freq='H').strftime('%d/%m/%Y') #The complete time series
li_all=[]
for i in range(0,len(chem_comp)):
    print(li_names[i])
    dfi=chem_comp.iloc[i][0].copy(deep=True)
    dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed', errors='coerce')
    dfi['datehour']=dfi['datetime'].dt.round('H')
    dfi=dfi.groupby('datehour').mean(numeric_only=True)
    if (li_names[i] in oasa.index)==True:
        print('IN OASA!')
        oai=oasa.loc[li_names[i]][0].copy(deep=True)
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
    dfi.to_csv(li_names[i]+'_ALL.txt')
    li_all.append(dfi)
#%%
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
#%% Composition plot
means_plot=pd.DataFrame()
means_plot.index=['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC', 'HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA']
for i in range(0,len(li_all)):
    if li_names[i] in gm.index and li_names[i] not in oasa.index:
        toadd=means[['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC']].iloc[i]
        toadd.index=['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
    elif li_names[i] in oasa.index and li_names[i] not in gm.index:
        toadd=means[['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA-1', 'HOA-2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl']].iloc[i]
        toadd.index=['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA-1', 'HOA-2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
    elif li_names[i] in oasa.index and li_names[i] in gm.index:
        toadd=means[['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA-1', 'HOA-2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl', 'BC']].iloc[i]
        toadd.index=['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'Coal', 'Peat', 'Wood', 'HOA-1', 'HOA-2', 'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
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
             'Coal', 'Peat', 'Wood', 'HOA-1', 'HOA-2',  'OOA_BBaq', 'OOA_BB', 'Amine-OA',]]

colors_all = ['red', 'blue', 'gold', 'fuchsia', 'black','green', 'grey', 'mediumorchid', 'saddlebrown', 'lightgreen', 'darkgreen','green',
              'purple','gainsboro', 'hotpink', 'hotpink', 'tan', 'sandybrown', 'grey', 'grey', 'forestgreen', 'mediumseagreen', 'skyblue' ]

fig, axs=plt.subplots(figsize=(7,5))
means_plot.plot(kind='bar', stacked=True, color=colors_all, ax=axs, width=0.85)
axs.set_xlabel('Sites', fontsize=13)
axs.set_ylabel('Concentration ($μg·m^{-3}$)', fontsize=13)
plt.legend(loc=(1.02,-0.24))
#%% PM1, OA intercomp!
for i in range(0,len(li_names)):
    dfi = li_all[i]
    print(li_names[i], dfi.columns)
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
    if li_names =='ZUR':
        col_list = ['Org', 'SO4', 'NO3', 'NH4', 'Chl']
    if li_names[i] in gm.index and li_names != 'ZUR':
        dfi['BC_ug']=dfi['BC(ng/m3)']/1000.0
        col_list = ['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC_ug']
    dfi['PM1_ACSM'] = dfi[col_list].sum(axis=1)
    dfi.plot(kind='scatter', x='PM1(μg/m3)', y='PM1_ACSM', ax=axs[0,0], color=dfi.index)
    if li_names[i] != 'LON-NK':
        dfi.plot(kind='scatter', x='PM2.5(μg/m3)', y='PM1_ACSM', ax=axs[1,0], color=dfi.index)
    if li_names[i] in oasa.index:
        col_list = [col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]
        dfi['OA_app'] = dfi[col_list].sum(axis=1)
        dfi.plot(kind='scatter', x='Org', y='OA_app', ax=axs[1,1], color=dfi.index)
    plt.suptitle(li_names[i])



