import sys
sys.path.append( '/home/greg/current/NMEG_utils/py_modules/' )

import load_nmeg as ld
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb as pdb
import datetime as dt

datapath = ('/home/greg/sftp/eddyflux/Ameriflux_files/provisional/')
fileList = os.listdir(datapath)

startYear = 2007
endYear = 2013

# Create empty dataframe with date range for 2007 to 2013
newidx = pd.date_range(str(startYear) + '-01-01',
        str(endYear + 1) + '-01-01', freq = '30T')
df = pd.DataFrame(index = newidx)

# Put dictionary keys in a list
keyList = ['MCon', 'PPine', 'JSav', 'PJC', 'PJG',
    'GLand', 'NewGLand', 'SLand']
# List AF site names in same order
siteNames = ['US-Vcm', 'US-Vcp', 'US-Wjs', 'US-Mpj', 'US-Mpg',
    'US-Seg', 'US-Sen', 'US-Ses']
# Create the empty dictionary to fill
dfDict = {}
PrecipDF = df.copy()

for i, site in enumerate(siteNames):
    # Get the files in the folder and an indexed dataframe
    siteFileList = [s for s in fileList if site in s]
    siteFileList = [s for s in siteFileList if 'with_gaps' in s]
    siteDf = pd.DataFrame()
    # Loop through each year and fill the dataframe
    for j in range(startYear, endYear + 1):
        fName = '{0}_{1}_with_gaps.txt'.format(site, j)
        
        if fName in siteFileList:
            yearDf = ld.load_aflx_file(datapath + fName, j)
            #siteDf = pd.merge(siteDf, yearDf,
            #        left_index=True, right_index=True,
            #how='left')
            siteDf = siteDf.append(yearDf)

    # Put the new multiyear precip value in a dataframe
    idxyrs = siteDf.index.year > 2005;
    siteDf = siteDf.iloc[idxyrs, :]
    siteDf = siteDf.reindex(newidx)
    PrecipDF[site] = siteDf.P


PrecipDaily = PrecipDF.resample('1D', how='sum')
wy = PrecipDaily.index + dt.timedelta(days=61)
PrecipDaily['year_w'] = wy.year
PrecipDaily['doy_w'] = wy.dayofyear

grouped = PrecipDaily.groupby('year_w')
wyrPrecipSum = grouped.aggregate(np.sum)

PrecipDaily.to_csv('dailyprecip_with_gaps.csv', na_rep='NaN')
wyrPrecipSum.to_csv('wyearPrecipSums.csv', na_rep='Nan') 
