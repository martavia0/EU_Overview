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
path_py_mac="/Users/martaviagonzalez/Documents/Documents - MacBook Pro de MVIA/GitHub/All_Treatment"
os.chdir(path_py_mac)
from Treatment import *
#%% Import treatment
trt = Basics(5)
trt.Hello()
print(trt.x)
#%% Import Composition files
path_data_mac="/Users/martaviagonzalez/Documents/GitHub/EU_Overview/Data/All/"
os.chdir(path_data_mac)
all_files=glob.glob(path_data_mac+'*Composition.txt')
li_files = [pd.read_csv(i, sep='\t') for i in all_files]
li_site_names = []
#%%
