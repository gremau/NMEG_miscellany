import load_nmeg as ld
import transform_nmeg as tr
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ipdb as ipdb

datapath = ('/home/greg/sftp/eddyflux/Ameriflux_files/provisional/'
            '2007-2014_Reichstein/')
fileList = os.listdir(datapath)

startyear = 2007
endyear = 2014

# List AF site names in same order
siteNames = [ 'US-Vcm', 'US-Vcp', 'US-Wjs', 'US-Mpj', 'US-Sen', 'US-Ses' ]

for i, site in enumerate( siteNames ):
    # Get the multi-year ameriflux dataframe
    site_df = ld.get_multiyr_aflx( site, datapath, gapfilled=True, 
                                   startyear=startyear, endyear=endyear )

    # Create a daily dataframe these are pretty much the defaults
    site_df_resamp = tr.resample_aflx( site_df, 
            freq='1D', c_fluxes=[ 'GPP', 'RE', 'FC' ], 
            le_flux=[ 'LE' ], avg_cols=[ 'TA', 'RH', 'Rg', 'RNET' ], 
            precip_col='PRECIP' , tair_col='TA' )

    # Export
    site_df_resamp.to_csv( '../processed_data/' + site 
            + '_biederman_synth.csv',
            na_rep = '-9999')


# Function for plotting fixed and original data
#def plotfixed(df, var1):
#    # Plot filled series over original data
#    plt.figure()
#    plt.plot(df.index, df[var1], '-r', df.index, df[var1 + '_fixed'], '-b')
#    #plt.plot(df.index, df[var2], 'om')
#    plt.legend(['Orig ' + var1, 'Fixed ' + var1, 'Lin'])
#    plt.title(var1 + ' data filling')
