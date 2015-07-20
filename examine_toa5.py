'''
A script for resampling 5 minute JSav TOA5 files to 30 minute files. and
have 5 minute frequency timestamps. 

This script depends on the NMEG_utils repository found on Greg Maurer's
github. You will need to ensure that this repository is accessable and the
correct path is added by python (edit line 19 of the script as appropriate).

Then:

1. Navigate to this directory in a terminal.
2. Run the 'resample_toa5'script (type 'python resample_toa5.py')
3. Copy the resampled files in ./resampled_toa5 to the parent directory
    ( /JSav/toa5 )
'''

import sys, os
# append the path to the NMEG_utils repository (greg's github)
sys.path.append(r'C:\Users\greg\Documents\GitHub\NMEG_utils\py_modules')

import load_nmeg as ld
import pdb as pdb

data_path = r'C:\Research_Flux_Towers\Flux_Tower_Data_by_Site' \
        r'\SLand\Raw_data_from_cards\archived_ts_data2\\'

file_list = os.listdir(data_path)

toa5_files = [ f for f in file_list if 'TOA5_' in f ]

new = ld.load_toa5_file(data_path + toa5_files[0])

print(new.columns)

#new.Ux.plot()
#new.Uy.plot()
#new.Uz.plot()
#new.Ts.plot()
#new.diag_csat.plot()


