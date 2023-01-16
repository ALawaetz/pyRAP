import os
import pandas as pd
import glob
import statistics as stat
import random
from Scripts import folders
import sys

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
### Create false discovery rate z-score by randomly distributing reads mapping to regions
### Regions are between peaks (peak[n]:peak[n + 1], peak[n + 1]:peak[n + 2],..., peak[i -1]:peak[i]), i == total number of peaks.
### For each region it is calculated the height necessary to reach the defined z-score.
### Low confidence peaks (z-score less than z_Start) are discarded if their value is less than the FDR determined height.

def fdr_local(file, z_score, output_folder):
    wig = pd.read_csv(file, sep = '\t', header = 0)
    peaks_df = wig[wig['peak'] != '.'].copy()
    indexes = peaks_df.index.tolist()
    peaks = peaks_df['peak'].tolist()
    values = peaks_df['value'].tolist()
    peak_value_index = []

    for p, v, i in zip(peaks, values, indexes):
        peak_value_index.append((p, v, i))

    z_heights = []
    discard = pd.DataFrame()
    for pvi in peak_value_index:
        start = wig['location'][pvi[2]] - 50
        end = wig['location'][pvi[2]] + 50
        reads = wig['value'][(wig['location'] >= start) & (wig['location'] <= end) & (wig['value'] != 0)]
        if stat.median(reads) < 1:
            sum_reads = int(sum(reads) / min(reads))
        else:
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
        if stat.median(reads) < 1:
            themin = min(reads)
            heights = [i * themin for i in heights]
        else:
            pass
        # for i in list(set(random_reads_allocation)):
        #     heights.append(random_reads_allocation.count(i))
        try:
            z_fdr = stat.mean(heights) + z_score * stat.stdev(heights)
        except:
            print('review code in fdr module')
        if pvi[1] < z_fdr:
            discard = pd.concat([discard, wig[wig.index == pvi[2]]])
            wig = wig[wig.index != pvi[2]]
            z_heights.append(z_fdr)
    discard['z_fdr_height'] = z_heights
    wig.to_csv(os.path.join(output_folder, os.path.basename(file)), sep = '\t', index = False)
    discard.to_csv(os.path.join(output_folder, os.path.basename(file) + '_discard'), sep = '\t', index = False)



def fdr(file, z_score, input_folder, output_folder_1, output_folder_2, z_start, FDR_tuples_folder):
    print('Calculating FDR')
    df = pd.read_csv(file, sep = '\t', header = 0)
    peak_dict = {}
    peak_dict['startfwd'] = sorted(list(set(df['start'][df['strand'] == '+'].tolist())), reverse = True)
    peak_dict['stopfwd'] = sorted(list(set(df['end'][df['strand'] == '+'].tolist())))
    peak_dict['startrev'] = sorted(list(set(df['end'][df['strand'] == '-'].tolist())))
    peak_dict['stoprev'] = sorted(list(set(df['start'][df['strand'] == '-'].tolist())), reverse = True)

    pooled_peaks = {}
    for key, value in peak_dict.items():
        new = []
        for item in value:
            if len(new) == 0:
                new.append(int(item))
            else:
                if abs(new[-1] - item) <= 10:
                    new.remove(new[-1])
                    new.append(int(item))
                else:
                    new.append(int(item))
        pooled_peaks[key] = sorted(new)

    double_pooled = {}
    double_pooled['fwd'] = sorted(list(set(pooled_peaks['startfwd'] + pooled_peaks['stopfwd'])))
    double_pooled['rev'] = sorted(list(set(pooled_peaks['startrev'] + pooled_peaks['stoprev'])))

    wiglist = []
    for file in sorted(glob.glob(input_folder + '*.csv')):
        wiglist.append(os.path.basename(file))


    for z in range(0, len(glob.glob(input_folder + '*.csv')), 2):
        wig = pd.read_csv(input_folder + wiglist[z], sep = '\t', header = 0)
        wig = pd.concat([wig, pd.read_csv(input_folder + wiglist[z + 1], sep = '\t', header = 0)])

        ### determine height to create z-score above 10 by allocating reads across genes at random
        z_tuples = []
        # n = 0
        for key, value in double_pooled.items():
            if key == 'fwd':
                strand = 0
            else:
                strand = 1
            for i in range(0, len(value) - 1):
                start = value[i]
                end = value[i + 1]
                reads = int(sum(wig['value'][(wig['location'] >= start) & (wig['location'] <= end) & (wig['strand'] == strand)]))
                ### randomly distribute reads across gene
                gene_length = abs(end - start)
                genomic_positions = [i for i in range(gene_length)]
                random_reads_allocation = []
                for i in range(reads):
                    random_reads_allocation.append(random.choice(genomic_positions))

                ### bin heights; what is the distribution of heights
                heights = []
                for i in range(0, gene_length):
                    heights.append(random_reads_allocation.count(i))
                # for i in list(set(random_reads_allocation)):
                #     heights.append(random_reads_allocation.count(i))
                try:
                    z_tuples.append((start, end, strand, stat.mean(heights) + z_score * stat.stdev(heights)))
                except:
                    z_tuples.append((start, end, strand, 0))
                # n += 1
                # if n % 500 == 0:
                #     print(n)
                else:
                    pass

        df = pd.DataFrame()
        start = []
        end = []
        strand = []
        height = []
        for item in z_tuples:
            start.append(item[0])
            end.append(item[1])
            strand.append(item[2])
            height.append(item[3])

        df['start'] = start
        df['end'] = end
        df['strand'] = strand
        df['height'] = height

        wig.to_csv(os.path.join(output_folder_1, '{}_combined_wigs.csv'.format(wiglist[z])), sep = '\t', index = False)
        df.to_csv(os.path.join(FDR_tuples_folder, '{}_z_tuples.csv'.format(wiglist[z])), sep = '\t', index = False)

    ## takes about an hour to reach this point pr. condition


    ## remove peaks with low FDR z-score
    print('Removing peaks with low FDR')
    for file in sorted(glob.glob(output_folder_1 + '*_combined_wigs.csv')):
        df = pd.read_csv(file, sep = '\t', header = 0)
        df = df[(df['peak'] != '.')] #& (df['shift_cutoff'] < z_start)]
        tuples = pd.read_csv(FDR_tuples_folder + os.path.basename(file)[:-18] + '_z_tuples.csv', sep = '\t', header = 0)
        if len(tuples) == 0:
            print('''No annotations made.
            Consider expanding search area or lowering cutoff values.
            For help type "python3 pipeline.py --help"''')
            sys.exit()
        multiple_starts = pd.read_csv(folders.peak_folder_0 + 'multiple_starts', sep = '\t', header = 0)['multiple_starts'].tolist()
        multiple_stops = pd.read_csv(folders.peak_folder_0 + 'multiple_stops', sep = '\t', header = 0)['multiple_stops'].tolist()
        discard = pd.DataFrame()
        indexlist = []
        for index, row in df.iterrows():
            n = 0
            if row['end'] == 5 and row['location'] in multiple_starts:
                pass
            elif row['end'] == 3 and row['location'] in multiple_stops:
                pass
            else:
                heights = tuples['height'][(row['location'] >= tuples['start']) & (row['location'] <= tuples['end']) & (row['strand'] == tuples['strand'])].tolist()
                if row['location'] < min(tuples['start'][tuples['strand'] == row['strand']].tolist()):
                    heights = tuples['height'][tuples['start'] == min(tuples['start'][tuples['strand'] == row['strand']].tolist())].tolist()
                elif row['location'] > max(tuples['start'][tuples['strand'] == row['strand']].tolist()):
                    heights = tuples['height'][tuples['start'] == max(tuples['start'][tuples['strand'] == row['strand']].tolist())].tolist()
                append_to_indexlist = 'Yes'
                for height in heights:
                    if row['value'] > height:
                        append_to_indexlist = 'No'
                        break
                    else:
                        pass
                if append_to_indexlist == 'Yes':
                    discard = pd.concat([discard, df[df.index == index]])
                    indexlist.append(index)
                else:
                    pass

                # for i, r in tuples.iterrows():
                #     if row['location'] >= r['start'] and row['location'] <= r['end'] and row['strand'] == r['strand']:
                #         if row['value'] >= r['height']:
                #             if file == '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main/Pipeline/Peak_folder_2/WT_LB_exp_DeLoughery_fwd.csv_combined_wigs.csv':
                #                 print('???')
                #                 print(row['location'])
                #                 print(r['height'])
                #                 print('???')
                #             break
                #         else:
                #             n += 1
                #             if n == 2:
                #                 # unhash this line if you wanna save discarded peaks in seperate dataframe for later analysis
                #                 discard = pd.concat([discard, df[df.index == index]])
                #                 indexlist.append(index)
                #                 break
                #             else:
                #                 pass
                #     else:
                #         pass

        df = pd.read_csv(output_folder_1 + os.path.basename(file), sep = '\t', header = 0)
        for item in indexlist:
            df = df[df.index != item]
        fwd = df[df['strand'] == 0]
        rev = df[df['strand'] == 1]
        fwd = fwd.sort_values(['location'], axis = 'index').reset_index(drop = True)
        rev = rev.sort_values(['location'], axis = 'index', ascending = False).reset_index(drop = True)
        fwd.to_csv(os.path.join(output_folder_2, os.path.basename(file)[:-25] + 'fwd.csv'), sep = '\t', index = False)
        rev.to_csv(os.path.join(output_folder_2, os.path.basename(file)[:-25] + 'rev.csv'), sep = '\t', index = False)
        os.remove(file)
        # os.remove(output_folder_1 + os.path.basename(file)[:-18] + '_z_tuples.csv')
        discard.to_csv(os.path.join(input_folder, '{}_discard.txt'.format(os.path.basename(file))), sep = '\t', index = False)
