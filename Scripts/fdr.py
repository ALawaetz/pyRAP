import os
import pandas as pd
import glob
import statistics as stat
import random
# from Scripts import folders
import sys
from multiprocessing import Process, Pool

__author__		= "Andreas Lawaetz"
__copyright__	= "Copyright 2022"
__version__		= "0.0.1"
__credits__		= ["Andreas Lawaetz"]
__maintainer__	= "Andreas Lawaetz"
__email__		= "acl58@bath.ac.uk"
__status__		= "Production"

##################################################################################
#
#	options_unique.py
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


##################################################################################
### Peaks are re-evaluated using this false discovery rate module
### Peaks are discarded if their read value is less than what is expected by chance
### Read value expected by chance is calculated by taking the sum of all read values surrounding a peak (25 bp upstream and 
### 25 bp downstream) and distributing those reads at random across a hypothetical region of 50 bp. 
### We then count the number of reads (peak height) "mapping" to each position in this 50 bp region.
### A threshold value (z_fdr) is calculated as 
### z_fdr = mean(peak height) + fdr * standard_deviation(peak height)
### fdr is defined by the user using the --fdr flag
### For help type
### python3 pipeline.py --help

### Peaks are discarded if their read value is less than z_fdr


input_folder = sys.argv[1]
z_score = float(sys.argv[2])
output_folder = sys.argv[3]
processors = int(sys.argv[4])

def fdr_local(file):
    wig = pd.read_csv(file, sep = '\t', header = 0)
    peaks_df = wig[wig['peak'] != '.'].copy()
    indexes = peaks_df.index.tolist()
    values = peaks_df['value'].tolist()
    peak_value_index = []

    for v, i in zip(values, indexes):
        peak_value_index.append((v, i))

    z_heights = []
    discard = pd.DataFrame()

    for pvi in peak_value_index:
        start = wig['location'][pvi[1]] - 50
        end = wig['location'][pvi[1]] + 50
        reads = wig['value'][(wig['location'] >= start) & (wig['location'] <= end) & (wig['value'] != 0)]
        sum_reads = int(sum(reads))

        ### randomly distribute reads across gene
        gene_length = abs(end - start)
        genomic_positions = [i for i in range(gene_length)]
        random_reads_allocation = []
        for i in range(sum_reads):
            random_reads_allocation.append(random.choice(genomic_positions))

        ### bin heights; what is the distribution of heights
        heights = []
        for i in range(0, gene_length):
            heights.append(random_reads_allocation.count(i))

        try:
            z_fdr = stat.mean(heights) + z_score * stat.stdev(heights)
        except:
            print('review code in fdr module')

        if pvi[0] < z_fdr:
            discard = pd.concat([discard, wig[wig.index == pvi[1]]])
            wig = wig[wig.index != pvi[1]]
            z_heights.append(z_fdr)
    discard['z_fdr_height'] = z_heights
    wig.to_csv(os.path.join(output_folder, os.path.basename(file)), sep = '\t', index = False)
    discard.to_csv(os.path.join(output_folder, os.path.basename(file) + '_discard'), sep = '\t', index = False)


if __name__ == '__main__':
    with Pool(processors) as p:
        p.map(fdr_local, sorted(glob.glob(input_folder + '*csv')))

