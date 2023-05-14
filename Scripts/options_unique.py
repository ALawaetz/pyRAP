import os
import pandas as pd
import glob

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


### Peaks are sorted based on whether they appear in a single sample or in multiple samples.
### Unique peaks (peaks appearing in one sample only) are discarded if their z-score is below 
### the threshold set by the z_ii flag
### Greater stringency can be applied using the --unique flag so that peaks appearing in only N number of samples 
### are defined as unique.
### For help type
### python3 pipeline.py --help

### Uniqueness is defined as a peak that only appears in N number of samples (--unique flag, default is 1) +/- 5 bp.



class unique:

    def __init__(self, file, peak, strand, unique):
        self.df = pd.read_csv(file, sep = '\t', header = 0)
        self.df = self.df.sort_values(['strand', peak])
        self.total_list = self.df[peak][self.df['strand'] == strand].tolist()
        self.attributes = self.df['attributes'][self.df['strand'] == strand].tolist()
        self.unique_list = []
        self.multiple_list = []
        self.unique = unique
        del self.df


    ### Function to determine if a peak is unique (more than 5 bp from next)
    def unique_and_multiple(self):
        for i in range(0, len(self.total_list)):
            attributes_list = []
            if self.total_list[i] in self.unique_list or self.total_list[i] in self.multiple_list:
                pass ### this is for unique start that are connected to multiple stops so it will not be added to the list twice.
            elif i == len(self.total_list) - 1 and abs(self.total_list[i] - self.total_list[i - 1]) > 5:
                self.unique_list.append(self.total_list[i])
            elif i == len(self.total_list) - 1 and abs(self.total_list[i] - self.total_list[i - 1]) <= 5:
                self.multiple_list.append(self.total_list[i])
            elif i == 0 and abs(self.total_list[i] - self.total_list[i + 1]) > 5:
                self.unique_list.append(self.total_list[i])
            elif i != 0 and i != len(self.total_list) - 1 and abs(self.total_list[i] - self.total_list[i + 1]) > 5 and abs(self.total_list[i] - self.total_list[i - 1]) > 5:
                self.unique_list.append(self.total_list[i])
            else:
                attributes_list.append(self.attributes[i])
                n = 0
                while i + n + 1 < len(self.total_list) and abs(self.total_list[i + n] - self.total_list[i + n + 1]) <= 5:
                    attributes_list.append(self.attributes[i + n + 1])
                    n += 1
                n = 0
                while i - n - 1 >= 0 and abs(self.total_list[i - n] - self.total_list[i - n - 1]) <= 5:
                    attributes_list.append(self.attributes[i - n - 1])
                    n += 1
                else:
                    if len(list(set(attributes_list))) <= self.unique:
                        self.unique_list.append(self.total_list[i])
                    else:
                        self.multiple_list.append(self.total_list[i])

# ### test if there is any bidirectional starts or stops at the excact same position
def bidirectional(list1, list2):
    global answer
    answer = 'no'
    for item in list1:
        if item in list2:
            answer = 'yes'
            break
        else:
            answer = 'no'


        ### If 'yes' is printed we have to be aware of bidirectional start and stop sites
        ### when reevalutating peaks and when making venn diagrams further down the pipeline

class saver:

    def __init__(self, list1, list2, destination, name):
        self.list1 = list1
        self.list2 = list2
        self.combined = self.list1 + self.list2
        self.destination = destination
        self.name = name
        self.path = destination + '/' + name

    def to_disk(self):
        peak = pd.DataFrame()
        peak[self.name] = self.combined
        peak.to_csv(os.path.join(self.destination, self.name), sep = '\t', index = False)




class reevaluate:

    def __init__(self, unique_starts, unique_stops, multiple_starts, multiple_stops, peak_folder, z_ii, output_folder):
        self.unique_starts = pd.read_csv(unique_starts, sep = '\t', header = 0)[os.path.basename(unique_starts)].tolist()
        self.unique_stops = pd.read_csv(unique_stops, sep = '\t', header = 0)[os.path.basename(unique_stops)].tolist()
        self.multiple_starts = pd.read_csv(multiple_starts, sep = '\t', header = 0)[os.path.basename(multiple_starts)].tolist()
        self.multiple_stops = pd.read_csv(multiple_stops, sep = '\t', header = 0)[os.path.basename(multiple_stops)].tolist()
        self.peak_folder = peak_folder
        self.z_ii = z_ii
        self.output_folder = output_folder


    def func(self):
        for file in glob.glob(self.peak_folder + '/*csv'):
            new_peaks = []
            df = pd.read_csv(file, sep = '\t', header = 0)
            peak = df['peak'].tolist()
            location = df['location'].tolist()
            cutoff = df['shift_cutoff'].tolist()
            cutoff_small = df['shift_cutoff_small'].tolist()
            del df
            for i in range(0, len(peak)):
                if peak[i] == '.':
                    new_peaks.append('.')
                elif cutoff[i] >= self.z_ii or cutoff_small[i] >= self.z_ii:
                    new_peaks.append(peak[i])
                elif peak[i] == 'START' and location[i] in self.multiple_starts or peak[i] == 'STOP' and location[i] in self.multiple_stops:
                    new_peaks.append(peak[i])
                else:
                    new_peaks.append('.')
            df = pd.read_csv(file, sep = '\t', header = 0)
            df['peak'] = new_peaks
            df.to_csv(os.path.join(self.output_folder, os.path.basename(file)), sep = '\t', index = False)
            del peak, location, cutoff, new_peaks, df
