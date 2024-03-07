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
#%% Boxprops
bp = dict(linestyle='-', linewidth=0.6)
bp_grey =dict(linestyle='-', linewidth=0.6, color='gainsboro')
mp = dict(marker='o', linewidth=0.1,markeredgecolor='black', markerfacecolor='black', markersize=3)
mdp = dict(color='k',  linewidth=0.6)
wp = dict(linestyle='-', linewidth=0.6)
#%% Import Composition files
os.chdir(path_data)
all_files=glob.glob(path_data+'*Composition.txt')
li_site_names = [j[-29:-25] for j in all_files]
chem_comp=pd.DataFrame()
chem_comp['Chemical_composition']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_files]
li_names = ['ATOLL', 'BAQS', 'BCN', 'BO', 'CAO-NIC', 'CRE', 'DEM', 'DUB', 'GRA','HEL','HPB', 'INO','KRK', 
            'LON-MR', 'LON-NK', 'LYO', 'MAD-CIE', 'MAQS', 'MAR-LCP', 'MET', 'MI', 'NOA', 'PAR-GEN', 'PAR-BPE','PAR-HAL',
            'PD', 'POI', 'PRG-SUCH', 'REN', 'SIRTA','STR', 'TAL', 'TAR','VLN', 'ZUR'] 
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
fig, axs = plt.subplots(figsize=(6,9))
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
                   Line2D([0], [0], color='green', label='Remote background'), 
                   Line2D([0], [0], marker='D', color='grey', label='AMS'),
                   Line2D([0], [0], marker='s', color='grey', label='Q-ACSM'),
                   Line2D([0], [0], marker='o', color='grey', label='ToF-ACSM')]
axs.legend(handles=legend_elements, loc = (1.03,0.45))#'upper right')
plt.title('Data availability')
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

fig, axs=plt.subplots(figsize=(8,8), nrows=2, sharex=True, layout="constrained" )
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
axs[0].legend(handles=legend_elements, loc = (1.03,0.45))#'upper right')

axs[0].set_xticklabels(df_nr.columns, fontsize=10, rotation=90)
axs[0].set_ylabel('NR-PM$_1$ concentration \n ($μg·m^{-3}$)', fontsize=12)
axs[0].set_xlabel('Sites', fontsize=12)
axs[0].set_ylim(0,90)
axs[0].set_xlim(0,36)
axs[0].text(x=0.1, y=16, s='WHO PM$_{2.5}$ limit', color='red')
x1, y1 = [-10, 40], [15, 15]
axs[0].plot(x1, y1, zorder=1, color='red', lw=1)
means_rel.plot(kind='bar', stacked=True, ax=axs[1], color=nr_colors, zorder=3,legend=False)
axs[1].grid(axis='y', zorder=1)
# axs[1].set_ylabel('NR-PM$_1$ mean composition \n($μg·m^{-3}$)', fontsize=12)
axs[1].set_ylabel('NR-PM$_1$ relative \ncomposition (%)', fontsize=12)
axs[1].set_xlabel('\n\nSites')
axs[1].legend( loc=(1.03, 0.15))
plt.savefig('NRlevels_Relcomp.png')

#%% Importing SA
os.chdir(path_data)
all_sa_files=glob.glob(path_data+'*Output_TS.txt')
li_site_names = [j[-23:-18] for j in all_sa_files]
oasa=pd.DataFrame()
oasa['SA']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True, infer_datetime_format=True) for i in all_sa_files]
li_names_sa = ['ATOLL','BCN', 'BO', 'CAO-NIC', 'DEM', 'DUB','GRA', 'HPB', 'INO','KRK', 
               'LON-MR', 'LON-NK', 'MAR-LCP', 'MI', 'PD', 'SIRTA', 'TAR', 'VLN', 'ZUR'] 
oasa.index = li_names_sa
#%% Importing Crit poll meas
os.chdir(path_data)
all_gm_files=glob.glob(path_data+'*meteo.txt')
li_site_names = [j[-28:-22] for j in all_gm_files]
gm=pd.DataFrame()
gm['GM']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_gm_files]
li_names_gm = ['ATOLL','BAQS', 'BCN', 'DEM','HEL', 'LON-MR','LON-NK', 
               'MAD-CIE', 'MRS-LCP', 'NOA', 'PRG-SUCH', 'SIRTA'] 
gm.index = li_names_gm
#%% Importing BC
os.chdir(path_data)
all_bc_files=glob.glob(path_data+'*BC.txt')
li_site_names = [j[-12:-7] for j in all_bc_files]
print(li_site_names)
bc=pd.DataFrame()
bc['BC']=[pd.read_csv(i, sep='\t', na_values='null', keep_default_na=True) for i in all_bc_files]
li_names_bc = ['ATOLL','BAQS', 'BCN', 'DEM', 'GRA', 'HEL', 'INO', 'ISP', 'LON-MR','LON-NK', 
               'MAD-CIE', 'MI', 'MRS-LCP', 'NOA', 'SIRTA', 'ZUR'] 
bc.index = li_names_bc
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
    if (li_names[i] in bc.index) == True:
        print('IN BC!')
        bci=pd.DataFrame()
        bci=bc.loc[li_names[i]][0].copy(deep=True)
        bci['datetime2']=pd.to_datetime(bci['datetime'], dayfirst=True, errors='raise')
        bci['datehour']=bci['datetime2'].dt.round('H')
        bci=bci.groupby('datehour').mean(numeric_only=True)
        dfi=pd.merge(left=dfi, right=bci, how='outer',left_on='datehour', right_on='datehour')
    dfi.to_csv(li_names[i]+'_ALL.txt')
    li_all.append(dfi)
#%% Rearranging BC?
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
                     'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl']].iloc[i]
        toadd.index=['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl']
        toadd.name=li_names[i]
        means_plot = pd.merge(means_plot,toadd, how='outer',left_index=True, right_index=True)
    elif li_names[i] in oasa.index and li_names[i] in gm.index:
        toadd=means[['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
                     'SO4', 'NO3', 'NH4', 'Chl', 'BC']].iloc[i]
        toadd.index=['HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA','ShInd-OA', 'CSOA', 'CCOA', 
                     'OOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA',
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
                          'OOA_BBaq', 'OOA_BB', 'Amine-OA',]]

colors_all = ['red', 'blue', 'gold', 'fuchsia', 'black','green', 'grey', 'mediumorchid', 'saddlebrown', 'lightgreen', 'darkgreen','green',
              'purple','gainsboro', 'hotpink', 'tan', 'sandybrown', 'skyblue' ]

fig, axs=plt.subplots(figsize=(7,4))
means_plot.plot(kind='bar', stacked=True, color=colors_all, ax=axs, width=0.85)
axs.set_xlabel('Sites', fontsize=13)
axs.set_ylabel('Concentration ($μg·m^{-3}$)', fontsize=13)
plt.legend(loc=(1.02,-0.24))
#%%
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
                                   'PeatOA': 'sienna', 'CoalOA': 'sienna', 'POA': 'darkkhaki', 'CCOA': 'sandybrown',
                                   '58-OA': 'hotpink', 'ShInd-OA': 'purple', 'seasaltOA':'darkcyan','BBOA1': 'saddlebrown', 'BBOA2': 'saddlebrown' })
colors_oasa.loc['VLN'] = pd.Series(['slategrey', 'darkkhaki', 'saddlebrown', 'grey', 'darkgreen', 'yellowgreen', 'None'])

fig, axs=plt.subplots(nrows=4, ncols=5, figsize=(8,6))
rows=range(0,4)
cols=range(0,5)
matrix_idx = [(i, j) for i in rows  for j in cols]
for k in range(0,len(oasa)):
    factors=oasa.copy(deep=True)
    oasa_i=pd.DataFrame(oasa.iloc[k][0]).mean(numeric_only=True)
    colors_i=colors_oasa.loc[li_names_sa[k]].iloc[0:len(oasa_i)].tolist()                
    colors_i = [i for i in colors_i if i is not None]
    print(li_names_sa[k], colors_i)
    oasa_i.plot.pie(ax=axs[matrix_idx[k]], title=li_names_sa[k], fontsize=7, colors=colors_i, 
                    autopct='%1.0f%%', startangle=90,counterclock=False)
fig.delaxes(axs[3,4])
#%%
means_bytype = means_plot[['SO4', 'NO3', 'NH4', 'Chl', 'BC', 'Org', 'HOA', 'COA', 'BBOA', 'LO-OOA', 'MO-OOA',
                           'OOA','ShInd-OA', 'CSOA', 'CCOA', 'OOA_BBaq', 'OOA_BB', 'Amine-OA']]
metadata2=metadata.copy(deep=True)
metadata2.index=metadata['Acronym']
means_bytype['Type']=metadata2['Type']
means_bytype2= means_bytype.groupby(means_plot['Type']).mean()
means_bytype2.plot()
#%% MAPPPPP
li_nr_avg=[]
for i in range(0,len(li_dfs)):
    df=li_dfs[i]
    df['nr'] = df.sum(numeric_only=True, axis=1)
    li_nr_avg.append(df['nr'].mean())
    print(i, li_names[i], df['nr'].mean())
metadata['NR']=li_nr_avg


plt.rcParams['lines.markersize'] ** 2
import geopandas as gpd
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
fig, axs = plt.subplots(figsize=(11,8))
countries.head()
countries.plot(color="lightgrey", ax=axs)
md=metadata.copy(deep=True)
countries.boundary.plot(ax=axs, color='k', zorder=0)
mask_ub=md['Type']=='UB'
mask_rb=md['Type']=='RB'
mask_su=md['Type']=='SU'
mask_tr=md['Type']=='TR'
md['NR2']=md['NR']*5

md.loc[mask_ub].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='royalblue', ax=axs,linewidths=2, zorder=4, label = 'Urban background')
md.loc[mask_su].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='darkorange', ax=axs,linewidths=2,zorder=5, label = 'Suburban'  )
md.loc[mask_rb].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='green', ax=axs,linewidths=2, zorder=6, label='Regional Background')
md.loc[mask_tr].plot.scatter(x="Lon", y="Lat", s='NR2', alpha=0.5, c='grey', ax=axs,linewidths=2, zorder=4, label = 'Traffic')

# plt.legend(*sc.legend_elements("sizes"))
axs.set_xlim(-18, 40)
axs.set_ylim(33, 68)
axs.set_ylabel('')
axs.set_xlabel('')
plt.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), labelspacing=1)

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
    axs[j,0].set_ylim(-0.05, 0.3)
    axs[j, 0].set_ylabel(factors[j])
axs[4,0].set_xlabel('Monthly cycle')
axs[4,1].set_xlabel('Diel cycle')

# fig.suptitle('\nGroup 4: High Cl$^-$' , fontsize=14)#$NO_3 ^{-}$ ~ $SO_4^{2-}$'
#%%

#%%
no3_like = ['MAQS',' STR', 'REN', 'DUB', 'MAD-CIE', 'DEM', 'PAR-HAL', 'POI', 'MET', 'SIRTA', 'LON-NK', 'LYO', 'TAL', 'PAR-GEN', 'CRE', 'LON-MR', 'VLN', 'ATOLL', 'BO', 'PAR-BPE', 'IPR', 'GRA', 'MI']
so4_like = ['TAR', 'MAR-LCP', 'BCN', 'CAO-NIC', 'NOA', 'PD']
no3so4_like = ['HEL', 'DUB', 'PRG-SUCH',' ZUR', 'INO']
chl_like = ['KRK']
#%%Heatmap diel and monthly


