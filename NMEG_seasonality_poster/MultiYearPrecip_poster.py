import sys
sys.path.append( '/home/greg/current/NMEG_utils/py_modules/' )

import load_nmeg as ld
import transform_nmeg as tr
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb as pdb
import datetime as dt

# Cold season = Nov - Feb
# Spring = March - June
# Summer = July - October

PRISM_datapath = ('/home/greg/sftp/eddyflux/Ancillary_met_data/PRISM_daily/')

PRISM_filelist = os.listdir(PRISM_datapath)

startYear = 2007
endYear = 2013

PRISM_idx = pd.date_range(str(startYear) + '-01-01',
        str(endYear + 1) + '-01-01', freq = '1D')
PRISM_df = pd.DataFrame(index = PRISM_idx)

data = pd.DataFrame()
for i in range(startYear, endYear + 1):
    fName = 'PRISM_DailyPrecip_{0}.csv'.format(i)
    fData = ld.loadPRISMfile(PRISM_datapath + fName)
    data = data.append(fData)
    
PRISM_df = data.reindex(PRISM_idx)

# Add water year and season columns
PRISM_df = tr.add_WY_cols(PRISM_df)

# Rename columns and export
PRISM_df = PRISM_df.rename(columns={'MCon': 'US-Vcm', 'PPine': 'US-Vcp',
    'JSav' : 'US-Wjs', 'PJ' : 'US-Mpj', 'PJ_girdle' : 'US-Mpg', 
    'GLand' : 'US-Seg', 'New_GLand' : 'US-Sen', 'SLand' : 'US-Ses'})

PRISM_df.to_csv(datapath + 'PRISM_daily.csv', na_rep='NaN')


# Put dictionary keys in a list
keyList = ['MCon', 'PPine', 'JSav', 'PJC', 'PJG',
    'GLand', 'NewGLand', 'SLand']
# List AF site names in same order
siteNames = ['US-Vcm', 'US-Vcp', 'US-Wjs', 'US-Mpj', 'US-Mpg',
    'US-Seg', 'US-Sen', 'US-Ses']


# Sum precip data
grouped = PRISM_df.groupby(['year_w', 'season'])
wyrSeasonalSum = grouped.aggregate(np.sum)
wyrSeasonalSum.doy_w = None

ipdb.set_trace()

wyrSeasonalSum.plot()
wyrSeasonalSum.plot(kind='bar', stacked=True)

grouped = PRISM_df.groupby(['year_w'])
wyrAnnualSum = grouped.aggregate(np.sum)



PrecipDaily.to_csv('dailyprecip.csv', na_rep='NaN')
wyrPrecipSum.to_csv('wyearPrecipSums.csv', na_rep='Nan') 
