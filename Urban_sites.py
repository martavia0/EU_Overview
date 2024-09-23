# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 12:00:12 2024

@author: Marta Via
"""

# %% paths definition
path_py_wdws =r"C:/Users/maria/Documents/Marta Via/1. PhD/F. Scripts/Python Scripts/"
path_py_mac ="/Users/martaviagonzalez/Documents/GitHub/EU_Overview/"

path_data_wdws=r"C:/Users/marta/Documents/IDAEA-CSIC/Overview/Selected"
path_data_mac="/Users/martaviagonzalez/Documents/Documents - MacBook Pro de MVIA/Work/IDAEA-CSIC/Overview/Selected/"

path_individual_wdws=r"C:/Users/marta/Documents/IDAEA-CSIC/Overview/Selected_plots/"
path_individual_mac="/Users/martaviagonzalez/Documents/Documents - MacBook Pro de MVIA/Work/IDAEA-CSIC/Overview/Selected_plots/"

path_folder_wdws = r"C:\Users\marta\Documents\IDAEA-CSIC\Overview"
path_folder_mac = "/Users/martaviagonzalez/Documents/GitHub/EU_Overview/"
#
mac_or_wdws = 'wdws' #Introduce OS here
#
if mac_or_wdws=='mac':
    path_py = path_py_mac
    path_data = path_data_mac
    path_folder = path_folder_mac
    path_individual=path_individual_mac
else:
    path_py = path_py_wdws 
    path_data = path_data_wdws 
    path_folder = path_folder_wdws
    path_individual=path_individual_wdws
#%%import pandas as pd
import numpy as np
import glob
import os as os
import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy import stats
# os.chdir(path_py)
# from Treatment import *
# trt = Basics(5)
#%% Colors
nr_colors=['green', 'red', 'blue', 'gold', 'fuchsia']
#%% Boxprops
bp = dict(linestyle='-', linewidth=0.6)
bp_grey =dict(linestyle='-', linewidth=0.6, color='gainsboro')
mp = dict(marker='o', linewidth=0.1,markeredgecolor='black', markerfacecolor='black', markersize=3)
mdp = dict(color='k',  linewidth=0.6)
wp = dict(linestyle='-', linewidth=0.6)
#%% Import Composition files
os.chdir(path_data)
all_files=glob.glob('*Composition.txt')
li_site_names = [j[-29:-25] for j in all_files]
chem_comp=pd.DataFrame()
chem_comp['Chemical_composition']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_files]
# li_names = ['PAR-BPE', 'HPB','MAR-LCP', 'MI', 'CAO-NIC', 'LYO', 'SIRTA', 'HEL', 'LON-MR', 
#             'BAQS', 'BCN', 'TAR','GRA','PRG-SUCH', 'DEM','MAQS', 'NOA', 'PAR-HAL', 'INO', 
#             'MAD-CIE', 'DUB', 'POI', 'KRK', 'ZUR', 'PD', 'ATOLL', 'STR','TAL','REN', 'CRE', 
#             'LON-NK', 'MET', 'VLN', 'PAR-GEN', 'BO']
# chem_comp.index = li_names
li_names_good = ['ATOLL', 'BAQS', 'BCN', 'BO', 'CAO-NIC', 'CRE', 'DEM', 'DUB', 'GRA','HEL',
                 'HPB', 'INO','KRK', 'LON-MR', 'LON-NK', 'LYO', 'MAD-CIE', 'MAQS', 'MAR-LCP', 
                 'MET', 'MI', 'NOA', 'PAR-GEN', 'PAR-BPE','PAR-HAL', 'PD', 'POI', 'PRG-SUCH', 
                 'REN', 'SIRTA','STR', 'TAL', 'TAR','VLN', 'ZUR'] 
# chem_comp=chem_comp.reindex(li_names_good)
chem_comp.index=li_names_good
li_names=li_names_good
#%% We import metadatafiles
os.chdir(r"C:\Users\marta\Documents\IDAEA-CSIC\Overview\Selected\Data")
metadata = pd.read_csv("Sites_metadata_selected.txt", sep='\t')
metadata=metadata.sort_values('Acronym')
li_sites_types = metadata["Type"]
#%%
for i in range(0,len(chem_comp)):
    print(len(chem_comp.iloc[i][0]))
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
    # print(i, metadata['Acronym'][i], metadata['Type'][i])
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
    elif metadata['Type'].iloc[i]=='RB':
        li_color.append('green')
    print(i, metadata['Acronym'][i], metadata['Type'][i], li_color[i])
dates_plot=dates_plot.notnull().astype('int')
dates_plot=dates_plot.replace(0, np.nan)
li_color.append('royalblue')

li_names.reverse()
dates_plot.columns=li_names

for i in range(0,len(dates_plot.columns)):
    # print(dates_plot.columns[i])
    dates_plot[dates_plot.columns[i]]=dates_plot[dates_plot.columns[i]]*35-(i)
# fig, axs = plt.subplots(figsize=(6,9))
fig, axs = plt.subplots(figsize=(42, 63))

for m, col, c in zip(li_marker, dates_plot.columns, li_color):
    dates_plot[col].plot(marker=m, lw=0, legend=False, ax=axs, color=c, grid=True, markersize=55)
axs.set_yticks(range(0,len(dates_plot.columns)+1))
axs.set_yticklabels(['']+li_names, fontsize=70)
# axs.set_xticklabels( fontsize=70, rotation=0)
axs.tick_params(axis='x', labelsize=60)

axs.set_xlabel('Time (years)', fontsize=70)
axs.set_ylabel('Sites', fontsize=70)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   Line2D([0], [0], color='grey', label='Traffic'), 
                   Line2D([0], [0], color='green', label='Remote background'), 
                   Line2D([0], [0], marker='D', color='grey', label='AMS'),
                   Line2D([0], [0], marker='s', color='grey', label='Q-ACSM'),
                   Line2D([0], [0], marker='o', color='grey', label='ToF-ACSM')]
axs.legend(handles=legend_elements, loc = (1.03,0.45), fontsize=70)#'upper right')
plt.title('Data availability\n', fontsize=90)
fig.savefig('Site_availability.png')
#%% Site data organisation
li_names.reverse()
li_dfi=[]
for i in range(0,len(chem_comp)):
    print(li_names[i])
    dfi = pd.DataFrame(chem_comp.iloc[i][0])
    dfi.drop('Time (Local)', inplace=True, errors='ignore')
    dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed')
    li_dfi.append(dfi)    

#%% Theil-Senn in a loop Absolute and relative
from scipy import stats

comps=['NR-PM1', 'Org', 'SO4', 'NO3', 'NH4', 'Chl']
long_ts = ['ATOLL', 'BCN', 'CRE', 'DUB', 'HEL', 'HPB', 'INO', 'LYO', 'MAD-CIE', 'MAR-LCP', 'MET', 'NOA', 'POI', 'SIRTA', 'TAL', 'ZUR']
colors=['orange', 'royalblue', 'royalblue', 'royalblue', 'grey', 'green', 'orange', 
                 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'orange', 'royalblue', 'royalblue']

fig, axs=plt.subplots(figsize=(10,8), nrows=6, ncols=2, sharex=True)
for j in range(0,len(comps)):
    mk, mkr=[],[]
    print(comps[j])
    for i in range(0,len(li_dfi)):
        dfi=li_dfi[i]
        dfi['NR-PM1']=dfi[comps[1:]].sum(axis=1)
        dfid=dfi.groupby(dfi['datetime']).mean(numeric_only=True)
        res = stats.theilslopes(y=dfid[comps[j]], x=range(len(dfid)), alpha=0.95)
        if j!=0:
            dfidr=100.0*dfid[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].divide(dfid[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].sum(axis=1), axis=0)
            resr = stats.theilslopes(y=dfidr[comps[j]], x=range(len(dfid)), alpha=0.95)
            mkr.append(resr)
        mk.append(res)
    mk_df=pd.DataFrame(mk)
    mkr_df=pd.DataFrame(mkr)
    mk_df.to_csv('Mann_Kendall_results_'+comps[j]+'.txt', sep='\t')
#%%
#Selection of only the years of highest availability  
    if j!=0:
        mkr_df.to_csv('Mann_Kendall_results_Relative_'+comps[j]+'.txt', sep='\t')
        mkr_df['Type']=metadata['Type']
        mkr_df['Acr']=li_names
        df_mkr=mkr_df[mkr_df['Acr'].isin(long_ts)]
        df_mkr=df_mkr.sort_values( 'Acr', ascending=True)
        df_mkr.index=range(0,len(df_mkr))

    mk_df['Type']=metadata['Type']
    mk_df['Acr']=li_names
    df_mk=mk_df[mk_df['Acr'].isin(long_ts)]

    df_mk=df_mk.sort_values( 'Acr', ascending=True)
    df_mk.index=range(0,len(df_mk))
comp_labels = ['NR-PM$_1$', 'OA', 'SO$_4^{2-}$', 'NO$_3^{-}$', 'NH$_4^{+}$', 'Cl$^-$']
# Theil-Senn Plot
    df_mk['x']=df_mk.index  
    df_mk.plot.scatter(y='slope',x='x', ax=axs[j,0], marker='s',color=colors, s=50, zorder=7)#,
    axs[j,0].errorbar(y=df_mk['slope'],x=df_mk['x'],yerr=[df_mk['low_slope'], df_mk['high_slope']],fmt='o', color='grey', zorder=1)
    axs[j,0].set_ylabel(comp_labels[j])
    axs[j,0].set_xticks(range(0,len(df_mk)))
    axs[j,0].set_xticklabels(df_mk['Acr'], rotation=90)
    axs[j,0].grid()

    if j!=0:
        df_mkr['x']=df_mkr.index
        df_mkr.plot.scatter(y='slope',x='x', ax=axs[j,1], marker='s',color=colors, s=50, zorder=7)#,
        axs[j,1].errorbar(y=df_mkr['slope'],x=df_mkr['x'],yerr=[df_mkr['low_slope'], df_mkr['high_slope']],fmt='o', color='grey', zorder=1)
        axs[j,1].set_ylabel('')
        axs[j,1].set_xticks(range(0,len(df_mkr)))
        axs[j,1].set_xticklabels(df_mkr['Acr'], rotation=90)
        axs[j,1].grid()

axs[5,0].set_xlabel('\nSites', fontsize=14)
axs[5,1].set_xlabel('\nSites', fontsize=14)
fig.text(-0.03, 0.5, 'Theil-Senn slopes', va='center', rotation='vertical', fontsize=14)
fig.text(0.3, 0.93, 'Absolute', ha='center',fontsize=14)
fig.text(0.7, 0.93, 'Relative', ha='center',fontsize=14)

legend_elements = [Line2D([0], [0], color='royalblue', label='UB', ),
                   Line2D([0], [0], color='darkorange', label='SU'), 
                   Line2D([0], [0], color='grey', label='TR'), 
                   Line2D([0], [0], color='green', label='RB')]
plt.legend(handles=legend_elements, loc = (0.82,6.02))#'upper right')
fig.delaxes(axs[0,1])

fig.tight_layout()

#%% Mean site composition
means= pd.DataFrame([i.mean(numeric_only=True) for i in li_dfi])
means=means[['Org', 'SO4', 'NO3', 'NH4', 'Chl']]
means.index=li_names
means['nr']=means['Org']+means['SO4']+means['NO3']+means['NH4']+ means['Chl']
li_color2= li_color[:-1]
means['color']=li_color2
means.sort_values(by='nr', inplace=True)
means.drop('nr',axis=1, inplace=True)
means_rel = 100.0*means[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].divide(means[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].sum(axis=1), axis=0)
df_nr=pd.DataFrame()

for i in range(0,len(li_dfs)):
    df=pd.DataFrame(li_dfs[i].copy(deep=True))
    df.columns=['Chl', 'NH4', 'NO3', 'Org', 'SO4', 'datetime', 'datet']
    df['nr'] = df[['Chl', 'NH4', 'NO3', 'Org', 'SO4']].sum(numeric_only=True, axis=1)
    df_nr=pd.concat([df_nr, df['nr']], axis=1, ignore_index=True)
df_nr.columns=li_names
df_nr = df_nr[means.index.tolist()]

fig, axs=plt.subplots(figsize=(12,8), nrows=2, sharex=True, layout="constrained" )
positions=range(0,35)
bp_colors =[dict(linestyle='-', linewidth=0.6, color=li_color[i]) for i in range(0,len(li_color))]
bp_dict = df_nr.boxplot(showfliers=False, showmeans=True, ax=axs[0],positions=positions,return_type='dict',
                        medianprops=mdp, meanprops=mp, whiskerprops=wp, zorder=4,  patch_artist = True)
for i, (box, color) in enumerate(zip(bp_dict['boxes'], means['color'])):
    box.set_facecolor(color)
legend_elements = [Line2D([0], [0], color='royalblue', label='UB', ),
                   Line2D([0], [0], color='darkorange', label='SU'), 
                   Line2D([0], [0], color='grey', label='TR'), 
                   Line2D([0], [0], color='green', label='RB')]
axs[0].legend(handles=legend_elements, loc = (1.02,0.65), fontsize=11)#'upper right')

axs[0].set_xticklabels(df_nr.columns, fontsize=10, rotation=90)
axs[0].set_ylabel('NR-PM$_1$ concentration \n ($μg·m^{-3}$)', fontsize=18)
axs[0].set_xlabel('Sites', fontsize=12)
axs[0].set_ylim(0,90)
axs[0].set_xlim(0,36)
axs[0].text(x=0.1, y=17, s='WHO PM$_{2.5}$ limit', color='red')
x1, y1 = [-10, 40], [15, 15]
axs[0].plot(x1, y1, zorder=1, color='red', lw=1)
means_rel.plot(kind='bar', stacked=True, ax=axs[1], color=nr_colors, zorder=3,legend=False)
axs[1].grid(axis='y', zorder=1)
# axs[1].set_ylabel('NR-PM$_1$ mean composition \n($μg·m^{-3}$)', fontsize=12)
axs[1].set_ylabel('NR-PM$_1$ relative \ncomposition (%)', fontsize=18)
axs[1].set_xlabel('\n\nSites')
handles_comp = [Line2D([0], [0], color='green', label='$OA$', ),
                Line2D([0], [0], color='red', label='$SO_4^{2-}$'), 
                Line2D([0], [0], color='blue', label='$NO_3^{-}$'), 
                Line2D([0], [0], color='goldenrod', label='$NH_4^{+}$'),
                Line2D([0], [0], color='fuchsia', label='$Cl^{-}$')]
axs[1].legend(handles=handles_comp, loc=(1.02, 0.45), fontsize=11)
plt.savefig('NRlevels_Relcomp.png')

#%% Importing SA
os.chdir(path_data)
all_sa_files=glob.glob('*Output_TS.txt')
li_site_names = [j[-23:-18] for j in all_sa_files]
oasa=pd.DataFrame()
oasa['SA']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True, infer_datetime_format=True) for i in all_sa_files]
# li_names_sa_wrong=['TAR','MI', 'MAR-LCP', 'ZUR', 'VLN', 'CAO-NIC','PD', 'GRA', 'DEM','BO', 
#                     'DUB', 'KRK', 'HPB', 'LON-NK', 'BCN', 'ATOLL', 'INO', 'LON-MR', 'SIRTA'] 
# li_names_sa_wrong =['ATOLL', 'BCN', 'BO', 'CAO', 'DEM', 'DUB', 'GRA', 'HPB', 'INO', 'KRK', 
                    # 'LONMR', 'LONNK', 'ARLCP', 'MI', 'PD', 'SIRTA', 'TAR', 'VLN', 'ZUR']       
oasa.index = li_site_names

li_names_sa = ['ATOLL','BCN', 'BO', 'CAO-NIC', 'DEM', 'DUB','GRA', 'HPB', 'INO','KRK', 
               'LON-MR', 'LON-NK', 'MAR-LCP', 'MI', 'PD', 'SIRTA', 'TAR', 'VLN', 'ZUR'] 
# oasa=oasa.reset_index(li_names_gm)
oasa.index = li_names_sa
#%% Importing Crit poll meas
os.chdir(path_data)
all_gm_files=glob.glob('*meteo.txt')
li_site_names = [j[-28:-22] for j in all_gm_files]
gm=pd.DataFrame()
gm['GM']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_gm_files]
li_names_gm_wrong =   ['LON-NK', 'SIRTA', 'BAQS', 'MAR-LCP', 'BCN', 'DEM', 'PRG-SUCH', 
                       'NOA', 'ATOLL', 'LON-MR', 'HEL', 'MAD-CIE']
gm.index = li_names_gm_wrong
li_names_gm= ['ATOLL','BAQS', 'BCN', 'DEM','HEL', 'LON-MR','LON-NK', 
               'MAD-CIE', 'MRS-LCP', 'NOA', 'PRG-SUCH', 'SIRTA'] 
gm=gm.reindex(li_names_gm)
#%% Importing BC
os.chdir(path_data)
all_bc_files=glob.glob('*BC.txt')
li_site_names = [j[-12:-7] for j in all_bc_files]
bc=pd.DataFrame()
li_names_bc_wrong = ['MAR-LCP', 'HEL', 'GRA', 'LON-MR', 'SIRTA', 'BAQS', 'BCN', 'LON-NK', 
               'NOA', 'MAD-CIE', 'INO', 'MI', 'ZUR', 'DEM', 'ATOLL']
bc['BC']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_bc_files]
bc.index = li_names_bc_wrong
li_names_bc = ['ATOLL','BAQS', 'BCN', 'DEM', 'GRA', 'HEL', 'INO', 'ISP', 'LON-MR','LON-NK', 
               'MAD-CIE', 'MI', 'MRS-LCP', 'NOA', 'SIRTA', 'ZUR']
bc=bc.reindex(li_names_bc) 
#%% Joining nrpm1, oasa, gases, meteo
min_date = pd.date_range(start='01/01/2009', end = '31/12/2023', freq='H').strftime('%d/%m/%Y') #The complete time series
li_all=[]
for i in range(0,len(chem_comp)):
    print(li_names[i])
    dfi=chem_comp.iloc[i][0].copy(deep=True)
    dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed', errors='coerce')
    dfi['datehour']=dfi['datetime'].dt.round('H')
    dfi_d=dfi.groupby('datehour').mean(numeric_only=True)
    dfi_d=dfi.reset_index()
    if (li_names[i] in oasa.index)==True:
        print('IN OASA!')
        oai=oasa.loc[li_names[i]][0].copy(deep=True)
        oai['datetime']=pd.to_datetime(oai['PMF Time (UTC)'] ,dayfirst=True, errors='raise')
        oai['datehour']=oai['datetime'].dt.round('H')
        oai=oai.groupby('datehour').mean(numeric_only=True)
        oai = oai.reset_index()
        dfi_d=pd.merge(left=dfi_d, right=oai, how='outer',left_on='datehour', right_on='datehour')
        print(dfi_d.mean(numeric_only=True))
    if (li_names[i] in gm.index) == True:
        print('IN GM!')
        gmi=pd.DataFrame()
        gmi=gm.loc[li_names[i]][0].copy(deep=True)
        gmi['datetime']=pd.to_datetime(gmi['TIME UTC, end'], dayfirst=True,errors='raise')
        gmi['datehour']=gmi['datetime'].dt.round('H')
        gmi=gmi.groupby('datehour').mean(numeric_only=True)
        gmi = gmi.reset_index()
        dfi_d=pd.merge(left=dfi_d, right=gmi, how='outer',left_on='datehour', right_on='datehour')
    if (li_names[i] in bc.index) == True:
        print('IN BC!')
        bci=pd.DataFrame()
        bci=bc.loc[li_names[i]][0].copy(deep=True)
        bci['datetime2']=pd.to_datetime(bci['datetime'], dayfirst=True, errors='raise')
        bci['datehour']=bci['datetime2'].dt.round('H')
        bci=bci.groupby('datehour').mean(numeric_only=True)
        bci = bci.reset_index()
        dfi_d=pd.merge(left=dfi_d, right=bci, how='outer',left_on='datehour', right_on='datehour')
    dfi_d.to_csv(li_names[i]+'_ALL.txt')
    li_all.append(dfi_d)
#%% Rearranging BC?
# li_means=[]
# for i in range(0,len(li_all)):
#     print(li_names[i], '\n')
#     dfi = li_all[i]
#     print(dfi.columns)
#     dfi['PM1_sum']=dfi['Org']+dfi['SO4']+dfi['NO3']+dfi['NH4']+dfi['Chl']
#     if li_names[i] in gm.index:
#         print(li_names[i] + 'has BC')
#         if li_names[i] == 'DEM':
#             dfi['PM1_sum'] = dfi['PM1_sum']+dfi['BC']/1000.0
#             dfi['BC']=dfi['BC']/1000.0
#         elif li_names[i] == 'SIRTA':
#             dfi['PM1_sum'] = dfi['PM1_sum']+dfi['BC(ng/m3)_x']/1000.0
#             dfi['BC']=dfi['BC(ng/m3)_x']/1000.0
#         else: 
#             dfi['PM1_sum'] = dfi['PM1_sum']+dfi['BC(ng/m3)']/1000.0
#             dfi['BC']=dfi['BC(ng/m3)']/1000.0
#     li_means.append(dfi.mean(numeric_only=True))
# means=pd.DataFrame(li_means)
#%% Composition plot
means_plot=pd.DataFrame()
plt.rcParams.update({'font.size': 18})

means_plot.index=['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC', 'HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA', 'Amine-OA','CCOA', 'CSOA', 'ShInd-OA', 'OOA_BB', 'OOA_BBaq']
for i in range(0,len(li_all)):
    dfi=li_all[i]
    print(li_names[i])
    # GM and BC
    if (li_names[i] not in oasa.index) and (li_names[i] in gm.index):
        if 'BC' in dfi.columns:
            toadd=dfi[['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC']].mean()
            toadd.columns=['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC']
        else:
            toadd=dfi[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].mean()
            toadd.columns=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
        print(1, toadd.columns)
    #Only GM
    elif (li_names[i] in oasa.index) and (li_names[i] in gm.index):
        if 'BC' in dfi.columns:
            toadd=dfi[[ 'SO4', 'NO3', 'NH4', 'Chl', 'BC']+[col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]].mean()
            toadd.columns=[ 'SO4', 'NO3', 'NH4', 'Chl', 'BC']+[col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]
        else: 
            toadd=dfi[[ 'SO4', 'NO3', 'NH4', 'Chl']+[col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]].mean()
            toadd.columns=[ 'SO4', 'NO3', 'NH4', 'Chl']+[col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot, toadd, how='outer', left_index=True, right_index=True)
        print(2, toadd.columns)
    # Only OA
    elif (li_names[i] in oasa.index) and (li_names[i] not in gm.index): 
        if 'BC' in dfi.columns:
            col_oa_list = ['SO4', 'NO3', 'NH4', 'Chl', 'BC']+[col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]
        else: 
            col_oa_list = ['SO4', 'NO3', 'NH4', 'Chl']+[col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]
        toadd=dfi[col_oa_list].mean()
        toadd.columns=col_oa_list
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot, toadd, how='outer', left_index=True, right_index=True) 
    else:
        if 'BC' in dfi.columns:
            col_oa_list=['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC']
        else:
            col_oa_list=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
        toadd=dfi[col_oa_list].mean()
        toadd.columns=col_oa_list
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
        print(4, toadd.columns)        

means_plot=means_plot.T

means_plot = means_plot[['Org','SO4','NO3', 'NH4', 'Chl', 'BC', 'HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','OOA', 'Amine-OA', 'CCOA', 'CSOA', 'ShInd-OA', 'OOA_BB', 'OOA_BBaq']]

colors_all = ['green', 'red', 'blue', 'gold', 'fuchsia', 'black','grey', 'mediumorchid', 'saddlebrown', 'lightgreen', 'darkgreen','green',
              'skyblue', 'hotpink','rosybrown', 'purple','sandybrown', 'tan']

fig, axs=plt.subplots(figsize=(16,8))
means_plot.plot(kind='bar', stacked=True, color=colors_all, ax=axs, width=0.75)
axs.set_xlabel('Sites', fontsize=22)
axs.set_xticks(range(0,len(means_plot.index)))
axs.set_xticklabels(means_plot.index, fontsize=18)

axs.set_ylabel('Concentration ($μg·m^{-3}$)', fontsize=22)
plt.legend(loc=(1.02,-0.1), fontsize=18)
plt.savefig("Compostition_plot.pdf")
#%% PM1, OA intercomp!
os.chdir(path_individual)
for i in range(0,len(li_names)):
    dfi = li_all[i]
    col_list=[]
    print(li_names[i])
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
    if ("BC(ng/m3)" in dfi.columns) and ("BC" not in dfi.columns):
        dfi['BC_ug']=dfi['BC(ng/m3)']/1000.0
        col_list = ['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC_ug']  
    elif ("BC(ng/m3)" not in dfi.columns) and ("BC" in dfi.columns):
        dfi['BC_ug']=dfi['BC']/1000.0
        col_list = ['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC_ug']  
    elif ("BC(ng/m3)" in dfi.columns) and ("BC" in dfi.columns):
        dfi['BC_ug']=dfi['BC(ng/m3)']/1000.0
        col_list = ['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BC_ug']  
        dfi.plot(kind='scatter', x='BC_ug', y='BC', ax=axs[1,0], color=dfi.index)
    elif ("BC(ng/m3)" not in dfi.columns) and ("BC" not in dfi.columns):
        col_list = ['Org', 'SO4', 'NO3', 'NH4', 'Chl']
    dfi['PM1_ACSM'] = dfi[col_list].sum(axis=1)
    if "PM1(μg/m3)" in dfi.columns:
        dfi.plot(kind='scatter', x='PM1(μg/m3)', y='PM1_ACSM', ax=axs[0,0], color=dfi.index)
    # if li_names[i] == 'LON-NK':
        # dfi.plot(kind='scatter', x='PM2.5(μg/m3)', y='PM1_ACSM', ax=axs[0,1], color=dfi.index)
    if li_names[i] in oasa.index:
        col_oa_list = [col for col in dfi.columns if (col.endswith('OA') or col.endswith('BB') or col =='Peat' or col =='Wood') ]
        dfi['OA_app'] = dfi[col_oa_list].sum(axis=1)
        dfi.plot(kind='scatter', x='Org', y='OA_app', ax=axs[1,1], color=dfi.index)
    print(dfi['PM1_ACSM'].mean())
    plt.suptitle(li_names[i])
    os.chdir(path_individual)
    fig.savefig(li_names[i]+'_mass_closure.png')


#%% Color arrangement and pies
factors_names=[]#pd.DataFrame()
for i in range(0,len(li_names_sa)):
    oasai=pd.DataFrame(oasa.iloc[i].tolist()[0])
    factors_names.append(oasai.columns[2:])
factors_names=pd.DataFrame(factors_names)
factors_names.index=li_names_sa
colors_oasa=factors_names.replace({'HOA': 'grey', 'COA': 'mediumpurple', 'Amine-OA': 'skyblue',
                                   'BBOA': 'saddlebrown', 'LO-OOA': 'yellowgreen','MO-OOA':'darkgreen', 
                                   'OOA': 'green', 'Total OOA':'green', 'OOA_BB': 'olivedrab', 'OOA_BBaq':'olive','LOA':'yellow',
                                   'HOA1': 'grey', 'HOA2': 'dimgrey', 'CSOA': 'rosybrown', 'WoodOA': 'sienna', 
                                   'PeatOA': 'sienna', 'CoalOA': 'rosybrown', 'POA': 'darkkhaki', 'CCOA': 'sandybrown',
                                   '58-OA': 'hotpink', 'ShInd-OA': 'purple', 'seasaltOA':'darkcyan','BBOA1': 'saddlebrown', 'BBOA2': 'saddlebrown'})
colors_oasa.loc['VLN'] = pd.Series(['slategrey', 'darkkhaki', 'saddlebrown', 'grey', 'darkgreen', 'yellowgreen', 'None'])
oasa_limeans=[]
fig, axs=plt.subplots(nrows=4, ncols=5, figsize=(10,8))
rows=range(0,4)
cols=range(0,5)
matrix_idx = [(i, j) for i in rows  for j in cols]
for k in range(0,len(oasa)):
    factors=oasa.copy(deep=True)
    oasa_i=pd.DataFrame(oasa.iloc[k][0]).mean(numeric_only=True)
    colors_i=colors_oasa.loc[li_names_sa[k]].iloc[0:len(oasa_i)].tolist()                
    colors_i = [i for i in colors_i if i is not None]
    print(li_names_sa[k], colors_i)
    oasa_i.plot.pie(ax=axs[matrix_idx[k]], title=li_names_sa[k], fontsize=8, colors=colors_i, 
                    autopct='%1.0f%%', startangle=90,counterclock=False, ylabel='')
    oasa_limeans.append(oasa_i)
fig.delaxes(axs[3,4])
#%% Redoing the pies plot in bars.
fig, axs=plt.subplots(figsize=(8,3))
oasa_df=pd.DataFrame(oasa_limeans)
oasa_df_norm = 100.0*oasa_df.divide(oasa_df.sum(axis=1), axis=0)
oasa_df_norm.index=li_names_sa
oasa_df_norm=oasa_df_norm[['LO-OOA', 'MO-OOA', 'OOA', 'OOA_BB', 'OOA_BBaq',
                           'HOA', 'HOA1', 'HOA2', 'BBOA', 'PeatOA', 'WoodOA', 'CoalOA',
                           'COA', 'Amine-OA', 'ShInd-OA', 'LOA ', 'POA ', 'CSOA']]
colors_OA= ['yellowgreen','darkgreen','green', 'olivedrab', 'olive', 
            'grey', 'dimgrey','silver', 'saddlebrown', 'darkkhaki', 'goldenrod', 'rosybrown', 
            'purple', 'skyblue', 'darkcyan', 'hotpink', 'yellow','steelblue']
oasa_df_norm.plot(kind='bar', stacked=True, ax=axs, color=colors_OA, zorder=7, width=0.83)
axs.legend(loc=(0.01, -0.91), ncol=5)
axs.grid(zorder=3)
axs.set_xlabel('Sites', fontsize=13)
axs.set_ylabel('Relative concentrations \n(%)', fontsize=13)
axs.set_title('OA sources', fontsize=14)
#%% Site diel plots.
oasadiel, oasadiel_std =[],[]
factors_names=[]
colors_oa.index= li_names_sa
# oasa.iloc[3].drop(axis=0,ArithmeticErrorinplace=True)
for i in range(0,len(li_names_sa)):
    print(i, li_names_sa[i])
    oasai=pd.DataFrame(oasa.iloc[i].tolist()[0])
    oasai['datetime']=pd.to_datetime(oasai['PMF Time (UTC)'], dayfirst=True, format='mixed', errors='coerce')
    oasai['Hour']=oasai['datetime'].dt.hour
    oasai.drop(labels='PMF Time (Local)', axis=1, inplace=True)
    oasai_diel=oasai.groupby(oasai['Hour']).mean(numeric_only=True)
    oasai_diel_std=oasai.groupby(oasai['Hour']).std(numeric_only=True)
    oasadiel.append(oasai_diel)
    oasadiel_std.append(oasai_diel)
    factors_names.append(oasai_diel.columns.tolist())
        
factors_names=pd.DataFrame(factors_names)
# factors_names3.drop(labels=['datetime', 'Hour'], axis=1, inplace=True, errors='ignore')
colors_oa=factors_names.replace({'HOA': 'grey', 'COA': 'mediumpurple', 'Amine-OA': 'skyblue',
                                   'BBOA': 'saddlebrown', 'LO-OOA': 'yellowgreen','MO-OOA':'darkgreen', 
                                   'OOA': 'green', 'OOA_BB': 'olivedrab', 'OOA_BBaq':'olive','LOA':'hotpink',
                                   'HOA1': 'dimgrey', 'HOA2': 'silver', 'CSOA': 'steelblue', 'WoodOA': 'goldenrod', 
                                   'PeatOA': 'darkkhaki', 'CoalOA': 'rosybrown', 'POA ': 'yellow', 'LOA ': 'hotpink', 
                                   '58-OA': 'hotpink', 'ShInd-OA': 'purple', 'seasaltOA':'darkcyan','BBOA1': 'saddlebrown', 'BBOA2': 'saddlebrown'})
# colors_oa.index=li_names_sa
colors_oa.index= li_names_sa
rows=range(0,4)
cols=range(0,5)
matrix_idx = [(i, j) for i in rows  for j in cols]
fig, axs=plt.subplots(nrows=4, ncols=5, figsize=(10,8), sharex=True, tight_layout=True)

for k in range(0,len(oasadiel)):
    diel_k=oasadiel[k].copy(deep=True)
    diel_k_norm = 100.0*diel_k.divide(diel_k.sum(axis=1), axis=0)
    print(li_names_sa[k], matrix_idx[k])
    colors_i=colors_oa.loc[li_names_sa[k]].iloc[0:len(oasa_i)].tolist()                
    colors_i = [i for i in colors_i if i is not None]
    diel_k_norm.plot(ax=axs[matrix_idx[k]], legend=False, color=colors_i, title=li_names_sa[k])
    axs[matrix_idx[k]].set_xlabel('')
    axs[matrix_idx[k]].set_xticks(range(0,25,3))
    axs[matrix_idx[k]].set_xticklabels(range(0,25,3))
    axs[matrix_idx[k]].grid()
axs[2,4].set_xticks(range(0,25,3))#
axs[2,4].set_xticklabels(range(0,25,3))#
fig.text(x=0.05, y=0.45, s='Normalised Concentration (%)', fontsize=14, va='center', rotation = 'vertical')
fig.text(x=0.55, y=0.05, s='Local time (h)', fontsize=14, ha='center')

fig.delaxes(axs[3,4])

legend_elements = [Line2D([0], [0], color='grey', label='HOA'), Line2D([0], [0], color='saddlebrown', label='BBOA' ),
                   Line2D([0], [0], color='yellowgreen', label='LO-OOA'), Line2D([0], [0], color='darkgreen', label='MO-OOA'),
                   Line2D([0], [0], color='green', label='OOA'), Line2D([0], [0], color='mediumorchid', label='COA'), 
                   Line2D([0], [0], color='rosybrown', label='Coal OA'),Line2D([0], [0], color='darkkhaki', label='Peat OA'), 
                   Line2D([0], [0], color='goldenrod', label='Wood OA'), Line2D([0], [0], color='skyblue', label='Amine-OA'),
                   Line2D([0], [0], color='purple', label='ShIndOA'), Line2D([0], [0], color='yellow', label='POA'), 
                   Line2D([0], [0], color='hotpink', label='LOA'), Line2D([0], [0], color='steelblue', label='CSOA')]
plt.legend(loc=((-3.7,-1)), handles=legend_elements,ncol=5)

#%% Site monthly plots.
oasamonthly, oasamonthly_std =[],[]
factors_names3=[]
for i in range(0,len(li_names_sa)):
    print(li_names_sa[i])
    oasai=pd.DataFrame(oasa.iloc[i].tolist()[0])
    oasai['datetime']=pd.to_datetime(oasai['PMF Time (UTC)'], dayfirst=True, format='mixed')
    oasai['Month']=oasai['datetime'].dt.month
    oasai.drop(labels='PMF Time (Local)', axis=1, inplace=True)
    oasai.drop(labels='PMF Time (UTC)', axis=1, inplace=True)    
    oasai_monthly=oasai.groupby(oasai['Month']).mean(numeric_only=True)
    print(oasai_monthly)
    oasai_monthly_std=oasai.groupby(oasai['Month']).std()
    oasamonthly.append(oasai_monthly)
    oasamonthly_std.append(oasai_monthly)
    factors_names3.append(oasai_monthly.columns[:-2].tolist())
factors_names=pd.DataFrame(factors_names3)
#%%
# factors_names3.drop(labels=['datetime', 'Hour'], axis=1, inplace=True, errors='ignore')
colors_oa=factors_names.replace({'HOA': 'grey', 'COA': 'mediumpurple', 'Amine-OA': 'skyblue',
                                   'BBOA': 'saddlebrown', 'LO-OOA': 'yellowgreen','MO-OOA':'darkgreen', 
                                   'OOA': 'green', 'OOA_BB': 'olivedrab', 'OOA_BBaq':'olive','LOA':'hotpink',
                                   'HOA1': 'dimgrey', 'HOA2': 'silver', 'CSOA': 'steelblue', 'WoodOA': 'goldenrod', 
                                   'PeatOA': 'darkkhaki', 'CoalOA': 'rosybrown', 'POA ': 'yellow', 'LOA ': 'hotpink', 
                                   '58-OA': 'hotpink', 'ShInd-OA': 'purple', 'seasaltOA':'darkcyan','BBOA1': 'saddlebrown', 'BBOA2': 'saddlebrown'})
colors_oa.index=li_names_sa
rows=range(0,4)
cols=range(0,5)
matrix_idx = [(i, j) for i in rows  for j in cols]
fig, axs=plt.subplots(nrows=4, ncols=5, figsize=(10,8), sharex=True, tight_layout=True)
# dfidr=100.0*dfid[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].divide(dfid[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].sum(axis=1), axis=0)

for k in range(0,len(oasadiel)):
    monthly_k=oasamonthly[k].copy(deep=True)
    monthly_k_norm = 100.0*monthly_k.divide(monthly_k.sum(axis=1), axis=0)
    # monthly_k_norm.drop(['Month'], axis=1, inplace=True)
    print(li_names_sa[k], matrix_idx[k])
    colors_i=colors_oa.loc[li_names_sa[k]].iloc[0:len(oasa_i)].tolist()                
    colors_i = [i for i in colors_i if i is not None]
    monthly_k_norm.plot(ax=axs[matrix_idx[k]], legend=False, color=colors_i, title=li_names_sa[k])
    axs[matrix_idx[k]].set_xlabel('')
    axs[matrix_idx[k]].set_xticks(range(0,13,1))
    axs[matrix_idx[k]].set_xticklabels(['', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], minor=True)
    axs[matrix_idx[k]].grid( alpha=0.7, axis='x')
axs[2,4].set_xticks(range(0,13,1))#
axs[2,4].set_xticklabels(['', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])#
fig.text(x=-0.05, y=0.5, s='Normalised Concentration (%)', fontsize=14, va='center', rotation = 'vertical')
fig.text(x=0.45, y=-0.05, s='Local time (h)', fontsize=14, ha='center')

fig.delaxes(axs[3,4])

legend_elements = [Line2D([0], [0], color='grey', label='HOA'), Line2D([0], [0], color='saddlebrown', label='BBOA' ),
                   Line2D([0], [0], color='yellowgreen', label='LO-OOA'), Line2D([0], [0], color='darkgreen', label='MO-OOA')]
plt.legend(loc=((-3.7,-1)), handles=legend_elements,ncol=5)
#%% 
m_all = pd.DataFrame()
for i in range(0,len(oasamonthly)):
    monthly_k=oasamonthly[i].copy(deep=True)
    monthly_k_norm = 100.0*monthly_k.divide(monthly_k.sum(axis=1), axis=0)
    m_all = pd.concat([m_all, monthly_k_norm], axis=1)
m_all.drop('CSOA', axis=1, inplace=True)
m_all.drop('OOA', axis=1, inplace=True)

m_all=m_all.T
m_all['factor'] = m_all.index
m_all_mean = m_all.groupby(by=m_all['factor']).mean().T
m_all_p25 = m_all.groupby(by=m_all.index).quantile(0.25, numeric_only=True).T
m_all_p75 = m_all.groupby(by=m_all['factor']).quantile(0.75, numeric_only=True).T
factors= m_all.columns
factors_all = ['COA', 'HOA', 'BBOA', 'LO-OOA', 'MO-OOA']
colors = ['saddlebrown', 'mediumpurple', 'grey', 'yellowgreen', 'darkgreen', 'green']
m_all_p25.index = range(0,12) 
m_all_p25 = m_all_p25[factors_all]
m_all_p75.index = range(0,12) 
m_all_p75 = m_all_p75[factors_all]
m_all_mean.index = range(0,12) 
m_all_mean = m_all_mean[factors_all]



fig, axs = plt.subplots(figsize=(6,4))
for j in range(0, len(m_all_mean.columns)):
    axs.plot(m_all_mean.iloc[:,j], color=colors[j], zorder=3)
    axs.fill_between(m_all_p25.index, m_all_p25.iloc[:,j], m_all_p75.iloc[:,j], color=colors[j], alpha=0.5, zorder=6)
axs.set_xticks(range(0,12))
axs.grid(axis='x', zorder=9)
axs.set_xticklabels( ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], fontsize=14)
axs.set_ylabel('Normalised Concentration (%)', fontsize=12)
axs.set_xlabel('Month', fontsize=12)

legend_elements = [Line2D([0], [0], color='grey', label='HOA'), 
                   Line2D([0], [0], color='saddlebrown', label='BBOA' ),
                   Line2D([0], [0], color='mediumpurple', label='COA'),
                   Line2D([0], [0], color='yellowgreen', label='LO-OOA'),
                   Line2D([0], [0], color='darkgreen', label='MO-OOA')]
plt.legend(loc=(-0.1,-0.3), handles=legend_elements, ncol = 5)
#%%
m_all_mean.plot.bar(stacked=True, color= colors)
#%%
m_all_HOA= m_all[m_all['factor']=='HOA']

#%% Mean compounds by type
means_bytype = means_plot[['SO4', 'NO3', 'NH4', 'Chl', 'BC', 'Org', 'HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA',
                           'OOA']]#,'ShInd-OA', 'CSOA', 'CCOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA']]
means_bytype['OOA']=means_bytype['MO-OOA']+means_bytype['LO-OOA']
metadata2=metadata.copy(deep=True)
metadata2.index=metadata['Acronym']
means_bytype['Type']=metadata2['Type']
means_bytype2= means_bytype.groupby(means_bytype['Type']).mean()
means_bytype2_desv= means_bytype.groupby(means_bytype['Type']).std()
means_bytype2=means_bytype2.reindex(['UB', 'SU', 'TR', 'RB'])
means_bytype2_desv=means_bytype2_desv.reindex(['UB', 'SU', 'TR', 'RB'])

fig, axs=plt.subplots(figsize=(8,4))
means_bytype2.T.plot(kind='bar', ax=axs, width=0.8, yerr=means_bytype2_desv.T, 
                     color=['royalblue', 'orange', 'grey', 'green'])
axs.set_ylabel('Mean Concentratizon ($μg·m^{-3}$)', fontsize=13)
axs.set_xlabel('Main NR-PM$_1$ compounds and sources', fontsize=13)
axs.legend(loc='upper left')
#%% MAPPPPP
li_nr_avg=[]
for i in range(0,len(li_dfs)):
    df=li_dfs[i]
    df['nr'] = df.sum(numeric_only=True, axis=1)
    li_nr_avg.append(df['nr'].mean())
    print(i, li_names[i], df['nr'].mean())
metadata['NR']=li_nr_avg

plt.rcParams['lines.markersize'] **2
plt.rcParams.update({'font.size': 20})

import geopandas as gpd
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
fig, axs = plt.subplots(figsize=(22,16))
countries.head()
countries.plot(color="lightgrey", ax=axs)
md=metadata.copy(deep=True)
countries.boundary.plot(ax=axs, color='k', zorder=0)
mask_ub=md['Type']=='UB'
mask_rb=md['Type']=='RB'
mask_su=md['Type']=='SU'
mask_tr=md['Type']=='TR'
md['NR2']=md['NR']*30

md.loc[mask_ub].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='royalblue', ax=axs,linewidths=2, zorder=4, label = 'Urban background')
md.loc[mask_su].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='darkorange', ax=axs,linewidths=2,zorder=5, label = 'Suburban'  )
md.loc[mask_rb].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='green', ax=axs,linewidths=2, zorder=6, label='Regional Background')
md.loc[mask_tr].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='grey', ax=axs,linewidths=2, zorder=4, label = 'Traffic')

# plt.legend(*sc.legend_elements("sizes"))
axs.set_xlim(-18, 40)
axs.set_ylim(33, 68)
axs.set_ylabel('')
axs.set_xlabel('')
plt.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), labelspacing=2, fontsize=18)
#%% MAPPPPP SAAAAA
li_nr_avg=[]
li_name=[]
for i in range(0,len(li_dfs)):
    if li_names[i] in li_names_sa:
        df=li_dfs[i]
        df['nr'] = df.sum(numeric_only=True, axis=1)
        li_nr_avg.append(df['nr'].mean())
        li_name.append(li_names[i])

plt.rcParams['lines.markersize'] **2
plt.rcParams.update({'font.size': 20})

import geopandas as gpd
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
fig, axs = plt.subplots(figsize=(16,12))
countries.head()
countries.plot(color="lightgrey", ax=axs)
md=metadata.copy(deep=True)

filtered_df = md[md['Acronym'].isin(li_name)]
md = filtered_df

countries.boundary.plot(ax=axs, color='k', zorder=0)
mask_ub=md['Type']=='UB'
mask_rb=md['Type']=='RB'
mask_su=md['Type']=='SU'
mask_tr=md['Type']=='TR'
md['NR2']=md['NR']*35

md.loc[mask_ub].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='royalblue', ax=axs,linewidths=2, zorder=4, label = 'Urban background')
md.loc[mask_su].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='darkorange', ax=axs,linewidths=2,zorder=5, label = 'Suburban'  )
md.loc[mask_rb].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='green', ax=axs,linewidths=2, zorder=6, label='Regional Background')
md.loc[mask_tr].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='grey', ax=axs,linewidths=2, zorder=4, label = 'Traffic')

# plt.legend(*sc.legend_elements("sizes"))
axs.set_xlim(-18, 40)
axs.set_ylim(33, 68)
axs.set_ylabel('')
axs.set_xlabel('')
plt.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), labelspacing=2, fontsize=18)

%% Diels and monthlys
fig, axs=plt.subplots(nrows=5, ncols=2, figsize=(10,10))
factors=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
for j in range(0,len(factors)):
    factor=factors[j]
    print(factor)
    df_h = pd.DataFrame()
    df_m = pd.DataFrame()
    for i in range(0,len(li_dfs)):
        if li_names[i] in li_names[i]:
            print(li_names[i])
            dfi=pd.DataFrame(chem_comp.iloc[i][0])
            dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed')
            dfi['Month']=dfi['datetime'].dt.month
            dfi['Hour']=dfi['datetime'].dt.hour
            df_h=pd.concat([df_h, dfi[factor].groupby(by=dfi['Hour']).mean(numeric_only=True)], axis=1)
            df_m=pd.concat([df_m, dfi[factor].groupby(by=dfi['Month']).mean(numeric_only=True)], axis=1)
        # df_h.columns = li_names
    # df_m.columns = li_names

    df_m_n=df_m/df_m.sum(axis=0)
    df_h_n=df_h/df_h.sum(axis=0)
    df_m_n.plot(ax=axs[j,0],  color='grey', legend=False)
    df_m_n.mean(axis=1).plot(ax=axs[j,0],  color=nr_colors[j], legend=False)
    df_h_n.plot(ax=axs[j,1],  color='grey', legend=False)
    df_h_n.mean(axis=1).plot(ax=axs[j,1],  color=nr_colors[j], legend=False)
    axs[j,0].set_ylim(-0.05, 0.5)
    axs[j, 0].set_ylabel(factors[j], fontsize=12)
axs[4,0].set_xlabel('Monthly cycle', fontsize=12)
axs[4,1].set_xlabel('Diel cycle', fontsize=12)

# fig.suptitle('\nGroup 4: High Cl$^-$' , fontsize=14)#$NO_3 ^{-}$ ~ $SO_4^{2-}$'
#%% Diels and monthlys per type
fig, axs=plt.subplots(nrows=5, ncols=2, figsize=(10,10))
factors=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
factors_name=['OA', 'SO$_4^{3-}$', 'NO$_3^-$', 'NH$_4^+}', 'Cl$^-$']

for j in range(0,len(factors)):
    factor=factors[j]
    print(factor)
    df_h= pd.DataFrame()
    df_m=pd.DataFrame()
    for i in range(0,len(li_dfs)):
        print(li_names[i])
        dfi=pd.DataFrame(chem_comp.iloc[i][0])
        dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format= 'mixed')
        dfi['Month']=dfi['datetime'].dt.month
        dfi['Hour']=dfi['datetime'].dt.hour
        df_h=pd.concat([df_h, dfi[factor].groupby(by=dfi['Hour']).mean(numeric_only=True)], axis=1)
        df_m=pd.concat([df_m, dfi[factor].groupby(by=dfi['Month']).mean(numeric_only=True)], axis=1)

    df_m.columns=li_names
    df_h.columns=li_names
    df_m_n=df_m/df_m.sum(axis=0)
    df_h_n=df_h/df_h.sum(axis=0)


    mask_ub = [li_names[i] for i in range(0,len(li_names)) if metadata['Type'].iloc[i]=='UB' ]
    mask_su = [li_names[i] for i in range(0,len(li_names)) if metadata['Type'].iloc[i]=='SU' ]
    mask_tr = [li_names[i] for i in range(0,len(li_names)) if metadata['Type'].iloc[i]=='TR' ]
    mask_rb = [li_names[i] for i in range(0,len(li_names)) if metadata['Type'].iloc[i]=='RB' ]

    #URBAN
    df_m_n[mask_ub].mean(axis=1).plot(ax=axs[j,0],  color='royalblue', legend=False)
    std_m_urban = df_m_n[mask_ub].std(axis=1)
    axs[j,0].fill_between(df_m_n.index,df_m_n[mask_ub].mean(axis=1)+std_m_urban,
                          df_m_n[mask_ub].mean(axis=1)-std_m_urban, color='royalblue', alpha=0.3)
    df_h_n[mask_ub].mean(axis=1).plot(ax=axs[j,1],  color='royalblue', legend=False)
    std_h_urban = df_h_n[mask_ub].std(axis=1)
    axs[j,1].fill_between(df_h_n.index,df_h_n[mask_ub].mean(axis=1)+std_h_urban,
                          df_h_n[mask_ub].mean(axis=1)-std_h_urban, color='royalblue', alpha=0.3)
    #SUBURBAN
    df_m_n[mask_su].mean(axis=1).plot(ax=axs[j,0],  color='orange', legend=False)
    std_m_su = df_m_n[mask_su].std(axis=1)
    axs[j,0].fill_between(df_m_n.index,df_m_n[mask_su].mean(axis=1)+std_m_su,
                          df_m_n[mask_su].mean(axis=1)-std_m_su, color='orange', alpha=0.3)
    df_h_n[mask_su].mean(axis=1).plot(ax=axs[j,1],  color='orange', legend=False)
    std_h_su = df_h_n[mask_su].std(axis=1)
    axs[j,1].fill_between(df_h_n.index,df_h_n[mask_su].mean(axis=1)+std_h_su,
                          df_h_n[mask_su].mean(axis=1)-std_h_su, color='orange', alpha=0.3)
    #TRAFFIC
    df_m_n[mask_tr].mean(axis=1).plot(ax=axs[j,0],  color='grey', legend=False)
    std_m_tr = df_m_n[mask_su].std(axis=1)
    axs[j,0].fill_between(df_m_n.index,df_m_n[mask_tr].mean(axis=1)+std_m_tr,
                          df_m_n[mask_tr].mean(axis=1)-std_m_tr, color='grey', alpha=0.3)
    df_h_n[mask_tr].mean(axis=1).plot(ax=axs[j,1],  color='grey', legend=False)
    std_h_tr = df_h_n[mask_tr].std(axis=1)
    axs[j,1].fill_between(df_h_n.index,df_h_n[mask_tr].mean(axis=1)+std_h_tr,
                          df_h_n[mask_tr].mean(axis=1)-std_h_tr, color='grey', alpha=0.3)
    #REGIONAL
    df_m_n[mask_rb].mean(axis=1).plot(ax=axs[j,0],  color='green', legend=False)
    std_m_rb = df_m_n[mask_rb].std(axis=1)
    axs[j,0].fill_between(df_m_n.index,df_m_n[mask_rb].mean(axis=1)+std_m_rb,
                      df_m_n[mask_rb].mean(axis=1)-std_m_rb, color='green', alpha=0.3)
    df_h_n[mask_rb].mean(axis=1).plot(ax=axs[j,1],  color='green', legend=False)
    std_h_rb = df_h_n[mask_rb].std(axis=1)
    axs[j,1].fill_between(df_h_n.index,df_h_n[mask_rb].mean(axis=1)+std_h_rb,
                          df_h_n[mask_rb].mean(axis=1)-std_h_rb, color='grey', alpha=0.3)
    
    axs[j,0].set_ylim(-0.05, 0.5)
    axs[j,1].set_ylim(-0.005, 0.10)

    axs[j, 0].set_ylabel(factors_name[j])
axs[4,0].set_xlabel('Monthly cycle')
axs[4,1].set_xlabel('Diel cycle')
fig.legend(loc=(1.1, 0.5))
# fig.suptitle('\nGroup 4: High Cl$^-$' , fontsize=14)#$NO_3 ^{-}$ ~ $SO_4^{2-}$'
#%% Diels and monthlys
fig, axs=plt.subplots(nrows=5, ncols=2, figsize=(10,10))
factors=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
for j in range(0,len(factors)):
    factor=factors[j]
    print(factor)
    df_h = pd.DataFrame()
    df_m = pd.DataFrame()
    for i in range(0,len(li_dfs)):
        if li_names[i] in li_names[i]:
            print(li_names[i])
            dfi=pd.DataFrame(chem_comp.iloc[i][0])
            dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed' )
            dfi['Month']=dfi['datetime'].dt.month
            dfi['Hour']=dfi['datetime'].dt.hour
            df_h=pd.concat([df_h, dfi[factor].groupby(by=dfi['Hour']).mean(numeric_only=True)], axis=1)
            df_m=pd.concat([df_m, dfi[factor].groupby(by=dfi['Month']).mean(numeric_only=True)], axis=1)
        # df_h.columns = li_names
    # df_m.columns = li_names

    df_m_n=df_m/df_m.sum(axis=0)
    df_h_n=df_h/df_h.sum(axis=0)
    df_m_n.plot(ax=axs[j,0],  color='grey', legend=False)
    df_m_n.mean(axis=1).plot(ax=axs[j,0],  color=nr_colors[j], legend=False)
    df_h_n.plot(ax=axs[j,1],  color='grey', legend=False)
    df_h_n.mean(axis=1).plot(ax=axs[j,1],  color=nr_colors[j], legend=False)
    axs[j,0].set_ylim(-0.05, 0.5)
    axs[j,1].set_ylim(-0.05, 0.2)
    axs[j, 0].set_ylabel(factors[j])
axs[4,0].set_xlabel('Monthly cycle')
axs[4,1].set_xlabel('Diel cycle')

# fig.suptitle('\nGroup 4: High Cl$^-$' , fontsize=14)#$NO_3 ^{-}$ ~ $SO_4^{2-}$'
#%% Individual diel and monthly plots I
factors=['Org', 'SO4', 'NO3', 'NH4', 'Chl']
diel_comp, monthly_comp, diel_comp_std, monthly_comp_std = [], [], [], []
for j in range(0,len(factors)):
    factor=factors[j]
    print(factor)
    df_h = pd.DataFrame()
    df_m = pd.DataFrame()
    for i in range(0,len(li_dfs)):
        if li_names[i] in li_names[i]:
            print(li_names[i])
            dfi=pd.DataFrame(chem_comp.iloc[i][0])
            dfi['datetime']=pd.to_datetime(dfi['Time (UTC)'], dayfirst=True, format='mixed' )
            dfi['Month']=dfi['datetime'].dt.month
            dfi['Hour']=dfi['datetime'].dt.hour
            df_h=pd.concat([df_h, dfi[factor].groupby(by=dfi['Hour']).mean(numeric_only=True)], axis=1)
            df_h_std=pd.concat([df_h_std, dfi[factor].groupby(by=dfi['Hour']).std(numeric_only=True)], axis=1)
            df_m=pd.concat([df_m, dfi[factor].groupby(by=dfi['Month']).mean(numeric_only=True)], axis=1)
            df_m_std=pd.concat([df_m_std, dfi[factor].groupby(by=dfi['Month']).std(numeric_only=True)], axis=1)

    diel_comp.append(df_h)
    diel_comp_std.append(df_h_std)
    monthly_comp.append(df_m)
    monthly_comp_std.append(df_m_std)

#%%Individual diel and monthly plots II
fig, axs=plt.subplots(nrows=5, ncols=7, figsize=(20,10), sharex=True, constrained_layout=True)
axes = [(j,i) for j in range(0,5) for i in range(0,7)]
for j in range(0,len(chem_comp)):
    for i in range(0,len(factors)):
        # print(color_nr[i], i, factors[i])
        if i == 0:
            ax2= axs[axes[j]].twinx()
            ax2.plot(diel_comp[i].iloc[:, j], color=color_nr[i])
            ax2.fill_between(range(0,24), diel_comp[i].iloc[:, j]-diel_comp_std[i].iloc[:, j], 
                                      diel_comp[i].iloc[:, j]+diel_comp_std[i].iloc[:, j],
                                      color=color_nr[i], alpha=0.3)
            ax2.spines['right'].set_color('green')
            ax2.yaxis.label.set_color('green')
            ax2.tick_params(axis='y', colors='green')
            if j %7 ==6:
                ax2.set_ylabel('OA\n($μg·m^{-3}$)', color='green')

        else: 
            axs[axes[j]].plot(diel_comp[i].iloc[:, j], color=color_nr[i])
            axs[axes[j]].fill_between(range(0,24), diel_comp[i].iloc[:, j]-diel_comp_std[i].iloc[:, j], 
                                      diel_comp[i].iloc[:, j]+diel_comp_std[i].iloc[:, j],
                                      color=color_nr[i], alpha=0.5)
        axs[axes[j]].set_xticks(range(0,24))
        axs[axes[j]].set_xticklabels(['0','','','','','','6', '','','','','','12', '','','','','','18', '','','','','',])
        axs[axes[j]].xaxis.set_major_locator(plt.MultipleLocator(6))
        axs[axes[j]].grid(which='major', axis='x', linestyle='--')
        axs[axes[j]].set_title(li_names_good[j])
        if j%7 ==0:
            axs[axes[j]].set_ylabel('Composition\n($μg·m^{-3}$)', color='k')
# plt.subplots_adjust(hspace=0.05)
#%%Individual diel and monthly plots II
fig, axs=plt.subplots(nrows=5, ncols=7, figsize=(20,10), sharex=True, constrained_layout=True)
axes = [(j,i) for j in range(0,5) for i in range(0,7)]
for j in range(0,len(chem_comp)):
    for i in range(0,len(factors)):
        # print(color_nr[i], i, factors[i])
        if i == 0:
            ax2= axs[axes[j]].twinx()
            ax2.plot(monthly_comp[i].iloc[:, j], color=color_nr[i])
            ax2.fill_between(range(1,13), monthly_comp[i].iloc[:, j]-monthly_comp_std[i].iloc[:, j], 
                                      monthly_comp[i].iloc[:, j]+monthly_comp_std[i].iloc[:, j],
                                      color=color_nr[i], alpha=0.3)
            ax2.spines['right'].set_color('green')
            ax2.yaxis.label.set_color('green')
            ax2.tick_params(axis='y', colors='green')
            if j %7 ==6:
                ax2.set_ylabel('OA\n($μg·m^{-3}$)', color='green')

        else: 
            axs[axes[j]].plot(monthly_comp[i].iloc[:, j], color=color_nr[i])
            axs[axes[j]].fill_between(range(1,13), monthly_comp[i].iloc[:, j]-monthly_comp_std[i].iloc[:, j], 
                                      monthly_comp[i].iloc[:, j]+monthly_comp_std[i].iloc[:, j],
                                      color=color_nr[i], alpha=0.5)
        axs[axes[j]].set_xticks(range(0,12))
        axs[axes[j]].set_xticklabels(['0','','','3','','','6', '','','9','','',])
        axs[axes[j]].xaxis.set_major_locator(plt.MultipleLocator(3))
        axs[axes[j]].grid(which='major', axis='x', linestyle='--')
        axs[axes[j]].set_title(li_names_good[j])
        if j%7 ==0:
            axs[axes[j]].set_ylabel('Composition\n($μg·m^{-3}$)', color='k')
# plt.subplots_adjust(hspace=0.05)
#%%
no3_like = ['MAQS',' STR', 'REN', 'DUB', 'MAD-CIE', 'DEM', 'PAR-HAL', 'POI', 'MET', 'SIRTA', 'LON-NK', 'LYO', 'TAL', 'PAR-GEN', 'CRE', 'LON-MR', 'VLN', 'ATOLL', 'BO', 'PAR-BPE', 'IPR', 'GRA', 'MI']
so4_like = ['TAR', 'MAR-LCP', 'BCN', 'CAO-NIC', 'NOA', 'PD']
no3so4_like = ['HEL', 'DUB', 'PRG-SUCH',' ZUR', 'INO']
chl_like = ['KRK']
#%% Superations arrangement

'''     SUPERATIONS!!    '''

limit_who_25_daily = 15
limit_who_15_yearly= 5

li_sup_daily, li_ndays, li_name, li_type = [],[], [],[]
for i in range(0,len(li_nr)):
    print(li_nr[i].columns[0])
    dfi=li_nr[i]
    mask=dfi.iloc[:,0]>=limit_who_25_daily
    sup_daily = dfi.loc[mask].count()[0]
    ndays=len(dfi)
    li_sup_daily.append(sup_daily)
    li_ndays.append(ndays)
    li_name.append(li_nr[i].columns[0])
    li_type.append(metadata['Type'].iloc[i])
sup_daily= pd.DataFrame(data={'Daily WHO superations':li_sup_daily, 'Nb days accounted':li_ndays})
sup_daily.index=li_name
sup_daily['Percentage superations'] = 100*sup_daily['Daily WHO superations']/sup_daily['Nb days accounted']

sup_daily['Type']=li_type
sup_daily['Type_int']=li_type
sup_daily['Type_int']=sup_daily['Type_int'].replace('UB',0).replace('RB', 1).replace('SU', 2).replace('C', 3).replace('M', 4).replace('A',5).replace('TR', 6)

colors=['royalblue','green', 'darkorange', 'mediumpurple', 'sienna', 'hotpink', 'darkcyan']
fig, axs = plt.subplots(figsize=(10,18), ncols=2, width_ratios=[3,2])
sup_daily=sup_daily.sort_values(by='Percentage superations')
sup_daily['Percentage superations'].plot(kind='barh', ax=axs[0], color=[colors[i] for i in sup_daily['Type_int']], zorder=3)
axs[0].set_ylabel('Percentage of WHO PM$_{2.5}$ daily thresholds superation', fontsize=20)
axs[0].set_xlim(0,100)
axs[0].grid(axis='x', zorder=0)
axs[0].set_title('NR-PM$_1$ concentration')

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   # Line2D([0], [0],  color='mediumpurple', label='Coastal'),
                   # Line2D([0], [0], color='sienna', label='Mountain'), 
                   # Line2D([0], [0], color='hotpink', label='Arctic'), 
                   Line2D([0], [0], color='darkcyan', label='Traffic'),
                   Line2D([0], [0], color='green', label='Regional background')]

axs[0].legend(handles=legend_elements, loc = (1.05,0.85))#'upper right')

sup_pie=sup_daily.groupby('Type').mean()
sup_pie=sup_pie.sort_values('Percentage superations', ascending=False)
sup_pie.plot.pie(y='Percentage superations', ax=axs[1], legend=False, autopct='%2.0f%%', labels=None, pctdistance=0.7,fontsize=14, 
                 startangle=90, counterclock=False, ylabel='', colors=['grey', 'orange', 'royalblue', 'green' ])

#%%
from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="constrained",  figsize=(15,15))

gs = GridSpec(2, 2)
ax1 = fig.add_subplot(gs[:, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1:, -1])

sup_daily['Percentage superations'].plot(kind='barh', ax=ax1, color=[colors[i] for i in sup_daily['Type_int']], zorder=3)
ax1.set_ylabel('Percentage of WHO PM$_{2.5}$ daily thresholds exceedances', fontsize=18)
ax1.set_xlim(0,100)
ax1.grid(axis='x', zorder=0)
ax1.set_title('NR-PM$_1$ concentration', fontsize=18)

legend_elements = [Line2D([0], [0], color='royalblue', label='Urban background', ),
                   Line2D([0], [0], color='darkorange', label='Suburban'), 
                   Line2D([0], [0], color='darkcyan', label='Traffic'),
                   Line2D([0], [0], color='green', label='Regional background')]

ax1.legend(handles=legend_elements, loc = (0.25,0.02))#'upper right')

sup_pie=sup_daily.groupby('Type').mean()
sup_pie=sup_pie.sort_values('Percentage superations', ascending=False)
sup_pie.plot.pie(y='Percentage superations', ax=ax2, legend=False,autopct='%2.0f%%', labels=None,pctdistance=0.7,fontsize=20, 
                 startangle=90, counterclock=False, ylabel='', colors=['grey', 'darkorange', 'royalblue', 'green'  ])
sup_bp=sup_daily.sort_values(by = 'Type_int')
positions = [3,1,2,0]
boxplot = sup_bp.boxplot(column=['Percentage superations'], by='Type', ax=ax3, showmeans=True,
                            boxprops=bp, medianprops=mdp, meanprops=mp, whiskerprops=wp, positions = positions) 
ax3.set_title('Percentage of WHO PM$_{2.5}$ daily \nthresholds exceedances', fontsize=18)
ax3.set_xlabel("Type of site", fontsize=18)
plt.suptitle('')
fig.text(x=0.02, y=0.88, s="(a)")
fig.text(x=0.52, y=0.88, s="(b)")
fig.text(x=0.52, y=0.52, s="(c)")

plt.show()
#%% By years and types of site.
li_year=[]
li_sites_names = li_names
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
ratio=df_year.sum(axis=1).divide(day_count.sum(axis=1))
# ratio.plot()

df_year['Type'], day_count['Type'] = li_type, li_type
dft_year=df_year.groupby('Type').sum().T
dayt_count=day_count.groupby('Type').sum().T
dfplot = 100*dft_year / dayt_count
dfplot.sort_values(by='Type', axis=1, ascending=False, inplace=True)
colors_types=['royalblue', 'grey', 'orange', 'green', ]

fig, axs =plt.subplots(figsize=(12,8), nrows=2, tight_layout=True)

dfplot.plot(legend=True,color=colors_types, ax=axs[0], marker='o', zorder=3, fontsize=11)
axs[0].set_xlabel('Years', fontsize=12)
axs[0].set_ylabel('Percentage of days with \nexcedances (%)', fontsize=15)
axs[0].grid(axis='y', zorder=0)
axs[0].grid(axis='x', zorder=0)
axs[0].legend(loc=(1.02,0.4), ncol=1, fontsize=11, title="Type")
axs[0].text(x=0.68,y=53, s='(a)', fontsize=11)
# Per each site, proportion of each sesason per superation days
li_year=[]
df_seas=pd.DataFrame()
df_season=[]
df_ndays=[]
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
    print(i, b.groupby('Season').count()['dt'])
    df_season.append(b.groupby('Season').count()['dt'])
    df_ndays.append(a.groupby('Season').count()['dt'])

df_seas=pd.DataFrame(df_season)
df_seas.index= li_name
day_count=pd.DataFrame(df_ndays)
day_count.index= li_name

df_plt=100*df_seas /day_count
df_plot=100*df_plt.T/df_plt.sum(axis=1)
# df_plot=df_plt

# df_plot=df_plot[['DJF', 'MAM', 'JJA', 'SON' ]]
# df_plot=df_plot.iloc[::-1]
# fig, axs=plt.subplots(figsize=(9,4))
df_plot.T.plot(kind='bar', stacked=True, ax=axs[1], legend=False, zorder=3,color=['royalblue', 'yellowgreen', 'gold', 'orange'])
axs[1].grid('y', zorder=2)
axs[1].set_ylabel('Seasonal distribution of days\n with exceedances (%)', fontsize=15)
axs[1].set_xlabel('Site', fontsize=12)

axs[1].legend(loc=(1.02,0.35), ncol=1, fontsize=11, title="Season")
axs[1].text(x=0.1,y=110, s='(b)', fontsize=11)

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
    c.drop('nr', inplace=True)

    c.index=['Chl', 'NH4', 'NO3', 'OA', 'SO4']
    li_sup.append(c)
    count_sup.append(100*d['NR']/len(b))
df_sup = pd.DataFrame(li_sup)
df_sup_count=pd.Series(count_sup)
df_sup.index = li_sites_names
df_sup = df_sup[['OA', 'SO4', 'NO3', 'NH4', 'Chl']]

color_nr = ['green', 'red', 'blue', 'gold', 'fuchsia']
fig, axs = plt.subplots(figsize=(14,6))
df_sup.plot(kind='bar', stacked=True, ax=axs, color = color_nr, zorder=7, width=0.9)
axs2=axs.twinx()
df_sup_count.plot(ax=axs2, marker='D', lw=0, color='k',zorder=8, markersize=5)

axs.set_ylabel('Absolute Concentration \n $(μg·m^{-3})$', fontsize=16)
axs2.set_ylabel('Percentage of days \nwith exceedance (%)', fontsize=16)
axs.set_xlabel('Site', fontsize=16)
handles_comp = [Line2D([0], [0], color='green', label='$OA$', ),
                Line2D([0], [0], color='red', label='$SO_4^{2-}$'), 
                Line2D([0], [0], color='blue', label='$NO_3^{-}$'), 
                Line2D([0], [0], color='goldenrod', label='$NH_4^{+}$'),
                Line2D([0], [0], color='fuchsia', label='$Cl^{-}$')]
axs.legend(handles=handles_comp, loc=(0.15, -0.5), ncols=5)
fig.suptitle('Days with exceedance')
axs2.set_ylim(-2,100)
#%% Ordered by perc with superations
color_nr = ['green', 'red', 'blue', 'gold', 'fuchsia']
fig, axs = plt.subplots(figsize=(20,8))
df_sup_count.index=li_sites_names
df_sup['count']=df_sup_count
df_sup=df_sup.sort_values(by='count')
df_sup.plot(y=['OA', 'SO4', 'NO3', 'NH4', 'Chl'], kind='bar', stacked=True, ax=axs, color = color_nr, zorder=7, width=0.9)
axs2=axs.twinx()
df_sup['count'].plot(ax=axs2, marker='D', lw=0, color='k',zorder=8, markersize=8)
axs.legend(loc=(0.25, -0.45), ncols=5, fontsize=20)
axs.set_ylabel('Absolute Concentration \n $(μg·m^{-3})$', fontsize=24)
axs2.set_ylabel('Percentage of days \nwith superation (%)', fontsize=24)
axs.set_xlabel('Site', fontsize=24)
axs.set_xticks(range(0,len(df_sup.index)))
axs.set_xticklabels(df_sup.index, fontsize=22)
axs.set_yticklabels(range(0,50,5), fontsize=22)
axs2.set_yticks(range(0,100,10))
axs2.set_yticklabels(range(0,100,10), fontsize=22)

fig.suptitle('Days with superation', fontsize=28)
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

