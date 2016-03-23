import sys
sys.path.append( '/home/greg/current/NMEG_utils/py_modules/' )
import load_nmeg as ld
import transform_nmeg as tr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb as pdb

remake_files = False

datapath = ('/home/greg/current/NMEG_utils/processed_data/')

startYear = 2007
endYear = 2013

def make_files() :
    
    # Path to ameriflux files
    afpath = ('/home/greg/sftp/eddyflux/Ameriflux_files_GM/'
            '2007-2013_Reichstein/')

    # Load measurement specific dataframes that will contain a
    # multi-year column for each site
    FC_df = ld.get_multiyr_var_allsites('FC', afpath, startYear, endYear)
    GPP_df = ld.get_multiyr_var_allsites('GPP', afpath, startYear, endYear)
    RE_df = ld.get_multiyr_var_allsites('RE', afpath, startYear, endYear)

    # Integrate C fluxes
    FC_int_df = ld.integrate_C_flux(FC_df)
    GPP_int_df = ld.integrate_C_flux(GPP_df)
    RE_int_df = ld.integrate_C_flux(RE_df)

    # Resample to daily sums
    FC_int_daily = FC_int_df.resample('1D', how='sum')
    GPP_int_daily = GPP_int_df.resample('1D', how='sum')
    RE_int_daily = RE_int_df.resample('1D', how='sum')

    # Add water year and season columns
    FC_int_daily = tr.add_WY_cols(FC_int_daily)
    GPP_int_daily = tr.add_WY_cols(GPP_int_daily)
    RE_int_daily = tr.add_WY_cols(RE_int_daily)

    # output to csv
    FC_int_daily.to_csv(datapath + 'FC_int_daily.csv', na_rep='NaN')
    GPP_int_daily.to_csv(datapath + 'GPP_int_daily.csv', na_rep='NaN')
    RE_int_daily.to_csv(datapath + 'RE_int_daily.csv', na_rep='NaN')

    print('new files made')

# Rebuild files in needed (make sure AF directory is available)
if remake_files :
    make_files()

FCdf = pd.read_csv(datapath + 'FC_int_daily.csv',
        index_col=0, parse_dates=True)
GPPdf = pd.read_csv(datapath + 'GPP_int_daily.csv',
        index_col=0, parse_dates=True)
REdf = pd.read_csv(datapath + 'RE_int_daily.csv',
        index_col=0, parse_dates=True)

# Make a plot showing bimodal NEE pattern at a site
subset = FCdf[FCdf.year_w==2008]
plt.plot(subset.index, subset['US-Seg'])
plt.show()

# Now we need to sum and resample some columns



