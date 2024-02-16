import os
import pandas as pd
import glob
import sys
import numpy as np

__author__		= "Andreas Lawaetz"
__copyright__	= "Copyright 2022"
__version__		= "0.0.1"
__credits__		= ["Andreas Lawaetz"]
__maintainer__	= "Andreas Lawaetz"
__email__		= "acl58@bath.ac.uk"
__status__		= "Production"

##################################################################################
#
#	supplyCapTerm.py
#
#
#	Copyright (c) Andreas Lawaetz 2022
#
#	Permission is hereby granted, free of charge, to any person obtaining a copy
#	of this software and associated documentation files (the "Software"), to deal
#	in the Software without restriction, including without limitation the rights
#	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#	copies of the Software, and to permit persons to whom the Software is
#	furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#	THE SOFTWARE.
#
##################################################################################

### This script takes the peak files from peak_folder_0
### It then sorts rows where a peak has been called but where the peak has been flagged (z < z_ii)
### If the peak is found in Cap-seq or Term-seq datasets, then the row is added to the equivalent file
### in peak_folder_3


input_file_0 = sys.argv[1]
z_ii = int(sys.argv[2])
capterm = sys.argv[3]


df_in = pd.read_csv(input_file_0, sep = '\t', header = 0)

if df_in['strand'][0] == 0:
    strand_id = 'fwd'
else:
    strand_id = 'rev'

capterm_data = np.load(capterm, allow_pickle='TRUE').item()

df = df_in[(df_in['peak'] != '.') & (df_in['shift_cutoff'] < z_ii) & (df_in['shift_cutoff'] < z_ii)]

new_df = pd.DataFrame()
for index, row in df.iterrows():
    if row['peak'] == 'START':
        peak_id = 'starts'
    if row['peak'] == 'STOP':
        peak_id = 'stops'
    
    if len([i for i in capterm_data[f'{peak_id}_{strand_id}'] if abs(i - row['location']) <= 5]) != 0:
        new_df = pd.concat([new_df, df[df.index == index]])


new_df['shift_cutoff'] = [z_ii] * len(new_df)
new_df['shift_cutoff_small'] = [z_ii] * len(new_df)



df_3 = pd.concat([df_in, new_df])
df_3 = df_3.drop_duplicates(['location', 'strand', 'end'], keep = 'last')
if strand_id == 'fwd':
    df_3 = df_3.sort_values(['location', 'end'], axis = 0, ascending = True)
else:
    df_3 = df_3.sort_values(['location', 'end'], axis = 0, ascending = False)

df_3.to_csv(input_file_0, sep = '\t', index = False)




