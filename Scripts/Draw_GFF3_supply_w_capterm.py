import pandas as pd
import numpy as np
import os
import glob
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
#	Draw_GFF3_supply_w_capterm.py
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


### This script is used to verify low-confidence peaks that might otherwise be discarded with options_unique.py and fdr.py modules. 
### Coordinates from alternative sources (e.g. BioCyc) are feeded into the pipeline with the flags --compare_starts and compare_stops.
### Using the flag --search_distance you tell pyRAP the disrepancy allowed to still pass as a match (e.g. coordinates found
### within 5 bp (default) of each other is considered a match). A peak verified this way will be used to make annotations.
### A second source of coordinates can be used simultaneously with the flags --operon_starts, --operon_stops and --search_distance_0 (default is 30 bp)



def density_full(csv, start, end):
    csv_list = csv['value'][(csv['location'] >= start) & (csv['location'] <= end)].tolist()
    try:
        the_dens = len([c for c in csv_list if c != 0])/len(csv_list)
    except ZeroDivisionError:
        the_dens = 0
    return the_dens

def density(csv, start, end):
    csv_datapoints = csv['value'][(csv['location'] > start) & (csv['location'] < end)].tolist()
    csv_datalist = []
    for i in range(0, len(csv_datapoints), 100):
        csv_datalist.append(csv_datapoints[i:i + 100])
    thegroup = []
    for csv_list in csv_datalist:
        thegroup.append(len([c for c in csv_list if c != 0])/len(csv_list))
    try:
        themin = min(thegroup)
    except ValueError:
        themin = 0
    return themin

def density_immediate(csv, peak, end, strand):
    if end == 5 and strand == 0 or end == 3 and strand == 1:
        csv_datapoints = csv['value'][(csv['location'] >= peak) & (csv['location'] < peak + 20)].tolist()
    if end == 5 and strand == 1 or end == 3 and strand == 0:
        csv_datapoints = csv['value'][(csv['location'] > peak - 20) & (csv['location'] <= peak)].tolist()
    dens = len([c for c in csv_datapoints if c != 0])/len(csv_datapoints)
    return dens

def overlap_gff(gff, start, end, strand):
    conditions = (start <= gff['start']) & (end > gff['start']) | (start >= gff['start']) & (start < gff['end']) | (start == gff['start']) & (end == gff['end'])
    overlap = gff[conditions]
    if strand == 0:
        strand = '+'
    else:
        strand = '-'
    overlap = overlap[overlap['strand'] == strand]
    return len(overlap)

def cap_term_inbetween(cap_term, start, end):
    if len([c for c in cap_term if c > start + 5 and c < end - 5]) == 0:
        inbetween = 'No'
    else:
        inbetween = 'Yes'
    return inbetween

def genes_inbetween(old_annotation, start, end):
    genes = old_annotation[(old_annotation['start'] > start) & (old_annotation['end'] < end)]
    CDS = genes[genes['feature'] == 'CDS']
    thelist = [len(genes), len(CDS)]
    return thelist

def corresponding_peak(corresponding_list, peak, end, strand):
    thelist = sorted(list(set(corresponding_list + [peak])))
    if end == 5 and strand == 0 or end == 3 and strand == 1:
        try:
            the_corresponding_peak = thelist[thelist.index(peak) + 1]
        except IndexError:
            the_corresponding_peak = None
    else:
        if thelist.index(peak) != 0:
            the_corresponding_peak = thelist[thelist.index(peak) - 1]
        else:
            the_corresponding_peak = None
    return the_corresponding_peak

def nearest_old_annotation(old_annotation, peak, end, strand):
    is_zero = 'No'
    if end == 5 and strand == '+':
        nearest_gene = old_annotation[(old_annotation['strand'] == '+') & (old_annotation['end'] >= peak) & (old_annotation['source'] != 'prediction')]
        nearest_gene = nearest_gene.sort_values(['start'])
        nearest_gene = nearest_gene.reset_index(drop = True)
        if len(nearest_gene) == 0:
            is_zero = 'Yes'
        else:
            go_to = nearest_gene['end'][0]
            short_to = nearest_gene['start'][0]
            if peak > nearest_gene['start'][0]:
                dens_to = nearest_gene['end'][0]
            else:
                dens_to = nearest_gene['start'][0]
    elif end == 5 and strand == '-':
        nearest_gene = old_annotation[(old_annotation['strand'] == '-') & (old_annotation['start'] <= peak) & (old_annotation['source'] != 'prediction')]
        nearest_gene = nearest_gene.sort_values(['end'], ascending = False)
        nearest_gene = nearest_gene.reset_index(drop = True)
        if len(nearest_gene) == 0:
            is_zero = 'Yes'
        else:
            go_to = nearest_gene['start'][0]
            short_to = nearest_gene['end'][0]
            if peak < nearest_gene['end'][0]:
                dens_to = nearest_gene['start'][0]
            else:
                dens_to = nearest_gene['end'][0]
    elif end == 3 and strand == '+':
        nearest_gene = old_annotation[(old_annotation['strand'] == '+') & (old_annotation['start'] <= peak) & (old_annotation['source'] != 'prediction')]
        nearest_gene = nearest_gene.sort_values(['end'], ascending = False)
        nearest_gene = nearest_gene.reset_index(drop = True)
        if len(nearest_gene) == 0:
            is_zero = 'Yes'
        else:
            go_to = nearest_gene['start'][0]
            short_to = nearest_gene['end'][0]
            if peak < nearest_gene['end'][0]:
                dens_to = nearest_gene['start'][0]
            else:
                dens_to = nearest_gene['end'][0]
    elif end == 3 and strand == '-':
        nearest_gene = old_annotation[(old_annotation['strand'] == '-') & (old_annotation['end'] >= peak) & (old_annotation['source'] != 'prediction')]
        nearest_gene = nearest_gene.sort_values(['start'])
        nearest_gene = nearest_gene.reset_index(drop = True)
        if len(nearest_gene) == 0:
            is_zero = 'Yes'
        else:
            go_to = nearest_gene['end'][0]
            short_to = nearest_gene['start'][0]
            if peak > nearest_gene['start'][0]:
                dens_to = nearest_gene['end'][0]
            else:
                dens_to = nearest_gene['start'][0]
    else:
        pass
    if is_zero == 'Yes':
        return None
    else:
        return [dens_to, go_to, short_to]

def highest_peak(csv, peak, end, strand, search_distance):
    # csv = csv[(csv['location'] > peak - 50) & (csv['location'] < peak + 50) & (csv['strand'] == strand) & (csv['end'] == end)].copy()
    if end == 5 and strand == 0 or end == 3 and strand == 1:
        csv = csv[(csv['location'] > peak - search_distance) & (csv['location'] < peak + 100) & (csv['strand'] == strand) & (csv['end'] == end)].copy()
    else:
        csv = csv[(csv['location'] < peak + search_distance) & (csv['location'] > peak - 100) & (csv['strand'] == strand) & (csv['end'] == end)].copy()
    try:
        themax = max(csv['value'].tolist())
    except:
        themax = 0
    if csv['value'].tolist().count(themax) / len([i for i in csv['value'].tolist() if i != 0]) > 0.34:
        go = 'No'
    else:
        go = 'Yes'
    try:
        corresponding_max = max(csv['value'][csv['end'] != end].tolist())
    except:
        corresponding_max = 0
    return [themax, go, corresponding_max]

def same_operon(csv, peak1, end1, peak2, end2, strand):
    operon1 = csv['operon'][(csv['location'] == peak1) & (csv['end'] == end1) & (csv['strand'] == strand)].iloc[0]
    operon2 = csv['operon'][(csv['location'] == peak2) & (csv['end'] == end2) & (csv['strand'] == strand)].iloc[0]
    if operon1 == operon2:
        is_same = 'Yes'
    else:
        is_same = 'No'
    return is_same

def overlap_half(gff, genome_name, condition):
    half_df = gff[gff['phase'] != '.'].copy()
    new_annotations = pd.DataFrame()
    for index, row in half_df.iterrows():
        conditions = (row['start'] <= half_df['start']) & (row['end'] > half_df['start']) | (row['start'] >= half_df['start']) & (row['start'] < half_df['end']) | (row['start'] == half_df['start']) & (row['end'] == half_df['end'])
        overlap = half_df[conditions].copy()
        overlap = overlap[(overlap['strand'] == row['strand']) & (overlap['phase'] != row['phase'])]
        if len(overlap) > 0:
            gff = gff[gff.index != index]
            if row['phase'] == 'start_site_resolved_only' and row['strand'] == '+':
                for end in overlap['end'].tolist():
                    if end == row['end']:
                        pass
                    else:
                        mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [row['start']], 'end': [end], 'score': [1000000], 'strand': [row['strand']], 'phase': ['.'], 'attributes': [condition]}
                        new_row = pd.DataFrame(mydata)
                        new_annotations = pd.concat([new_annotations, new_row])
            elif row['phase'] == 'start_site_resolved_only' and row['strand'] == '-':
                for start in overlap['start'].tolist():
                    if start == row['start']:
                        pass
                    else:
                        mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [row['end']], 'score': [1000000], 'strand': [row['strand']], 'phase': ['.'], 'attributes': [condition]}
                        new_row = pd.DataFrame(mydata)
                        new_annotations = pd.concat([new_annotations, new_row])
            else:
                pass
    gff = pd.concat([gff, new_annotations])
    gff = gff.drop_duplicates().reset_index(drop = True)

    return gff



class missing_peaks:
    def __init__(self, csv, gff, cap_term, genome_name, condition, old_annotation, search_distance, gff_folder, number):
        self.csv_file = csv
        self.csv = pd.read_csv(csv, sep = '\t', header = 0)
        self.gff = pd.read_csv(gff, sep = '\t', header = 0)
        self.cap_term = np.load(cap_term, allow_pickle='TRUE').item()
        self.genome_name = genome_name
        self.condition = condition
        self.old_annotation = pd.read_csv(old_annotation, sep = '\t', header = 0)
        self.search_distance = search_distance
        self.missing_peaks = pd.DataFrame({'location': [], 'value': [], 'strand': [], 'end': [], 'peak': [], 'shift_cutoff': [], 'shift_cutoff_small': []})
        self.gff_folder = gff_folder

        peaks = list(set(self.gff['start'].tolist() + self.gff['end'].tolist()))

        ### find peaks using cappable-seq and term-seq
        csv_peak = self.csv[self.csv['peak'] != '.'].copy()
        for index, row in csv_peak.iterrows():
            thehighest_peak = highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)
            if row['end'] == 5 and row['strand'] == 0 and row['value'] == thehighest_peak[0] and thehighest_peak[1] == 'Yes':
                if row['location'] not in peaks:# and len([c for c in self.cap_term['starts_fwd'] if abs(c - row['location']) <= search_distance]) != 0:
                    self.missing_peaks = pd.concat([self.missing_peaks, csv_peak[csv_peak.index == index]])
            elif row['end'] == 5 and row['strand'] == 1 and row['value'] == thehighest_peak[0] and thehighest_peak[1] == 'Yes':
                if row['location'] not in peaks:# and len([c for c in self.cap_term['starts_rev'] if abs(c - row['location']) <= search_distance]) != 0:
                    self.missing_peaks = pd.concat([self.missing_peaks, csv_peak[csv_peak.index == index]])
            else:
                pass
            ### If no strand information for term-seq
            if 'stops_all' in list(self.cap_term.keys()):
                if row['end'] == 3 and row['value'] == thehighest_peak[0] and thehighest_peak[1] == 'Yes':
                    if row['location'] not in peaks:# and len([c for c in self.cap_term['stops_all'] if abs(c - row['location']) <= 1]) != 0:
                        self.missing_peaks = pd.concat([self.missing_peaks, csv_peak[csv_peak.index == index]])
                else:
                    pass
            else:
                if row['end'] == 3 and row['strand'] == 0 and row['value'] == thehighest_peak[0] and thehighest_peak[1] == 'Yes':
                    if row['location'] not in peaks:# and len([c for c in self.cap_term['stops_fwd'] if abs(c - row['location']) <= search_distance]) != 0:
                        self.missing_peaks = pd.concat([self.missing_peaks, csv_peak[csv_peak.index == index]])
                elif row['end'] == 3 and row['strand'] == 1 and row['value'] == thehighest_peak[0] and thehighest_peak[1] == 'Yes':
                    if row['location'] not in peaks:# and len([c for c in self.cap_term['stops_rev'] if abs(c - row['location']) <= search_distance]) != 0:
                        self.missing_peaks = pd.concat([self.missing_peaks, csv_peak[csv_peak.index == index]])
                else:
                    pass

        if number == 1:
            self.missing_peaks.to_csv(os.path.join(folders.missing_peaks, os.path.basename(self.csv_file)), sep = '\t', index = False)
        try:
            other_missing = pd.read_csv(folders.missing_peaks + os.path.basename(csv), sep = '\t', header = 0)
        except FileNotFoundError:
            other_missing = pd.DataFrame()
        self.missing_peaks = pd.concat([self.missing_peaks, other_missing])
        self.missing_peaks = self.missing_peaks.drop_duplicates()
        self.missing_peaks = self.missing_peaks.sort_values(['location'])
        starts_fwd = sorted(list(set(self.gff['start'][self.gff['strand'] == '+'].tolist())))
        starts_rev = sorted(list(set(self.gff['end'][self.gff['strand'] == '-'].tolist())))
        stops_fwd = sorted(list(set(self.gff['end'][self.gff['strand'] == '+'].tolist())))
        stops_rev = sorted(list(set(self.gff['start'][self.gff['strand'] == '-'].tolist())))
        self.missing_starts_fwd_list = sorted(list(set(self.missing_peaks['location'][(self.missing_peaks['strand'] == 0) & (self.missing_peaks['end'] == 5)].tolist())))
        self.missing_starts_rev_list = sorted(list(set(self.missing_peaks['location'][(self.missing_peaks['strand'] == 1) & (self.missing_peaks['end'] == 5)].tolist())))
        self.missing_stops_fwd_list = sorted(list(set(self.missing_peaks['location'][(self.missing_peaks['strand'] == 0) & (self.missing_peaks['end'] == 3)].tolist())))
        self.missing_stops_rev_list = sorted(list(set(self.missing_peaks['location'][(self.missing_peaks['strand'] == 1) & (self.missing_peaks['end'] == 3)].tolist())))
        self.starts_fwd = starts_fwd + self.missing_starts_fwd_list
        self.starts_rev = starts_rev + self.missing_starts_rev_list
        self.stops_fwd = stops_fwd + self.missing_stops_fwd_list
        self.stops_rev = stops_rev + self.missing_stops_rev_list


    def supply_lone_peaks(self):
        missing_annotations = pd.DataFrame({'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []})
        for index, row in self.missing_peaks.iterrows():
            mydata = {}
            ### by adding a start site to the list of stop site and order the list
            ### we can make an annotation to the nearest corresponding peak.
            ### for forward start sites
            if row['end'] == 5 and row['strand'] == 0 and corresponding_peak(self.stops_fwd, row['location'], row['end'], row['strand']) != None:
                start = row['location']
                end = corresponding_peak(self.stops_fwd, row['location'], row['end'], row['strand'])
                if density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 and same_operon(self.csv, row['location'], row['end'], end, 3, row['strand']) == 'Yes':
                    if density(self.csv, start, end) >= 0.1 and cap_term_inbetween(self.cap_term['stops_fwd'], start, end) == 'No' or genes_inbetween(self.old_annotation, start, end)[1] == 1 and density(self.csv, start, end) >= 0.01 or genes_inbetween(self.old_annotation, start, end)[0] <= 1 and density(self.csv, start, end) >= 0.1:
                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [self.condition]}
            ### for reverse start sites
            elif row['end'] == 5 and row['strand'] == 1 and corresponding_peak(self.stops_rev, row['location'], row['end'], row['strand']) != None:
                end = row['location']
                start = corresponding_peak(self.stops_rev, row['location'], row['end'], row['strand'])
                if density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 and same_operon(self.csv, row['location'], row['end'], start, 3, row['strand']) == 'Yes':
                    if density(self.csv, start, end) >= 0.1 and cap_term_inbetween(self.cap_term['stops_rev'], start, end) == 'No' or genes_inbetween(self.old_annotation, start, end)[1] == 1 and density(self.csv, start, end) >= 0.01 or genes_inbetween(self.old_annotation, start, end)[0] <= 1 and density(self.csv, start, end) >= 0.1:
                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [self.condition]}
            ### for fwd stop sites
            elif row['end'] == 3 and row['strand'] == 0 and corresponding_peak(self.starts_fwd, row['location'], row['end'], row['strand']) != None:
                end = row['location']
                start = corresponding_peak(self.starts_fwd, row['location'], row['end'], row['strand'])
                if density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 and same_operon(self.csv, row['location'], row['end'], start, 5, row['strand']) == 'Yes':
                    if density(self.csv, start, end) >= 0.1 and cap_term_inbetween(self.cap_term['starts_fwd'], start, end) == 'No' or genes_inbetween(self.old_annotation, start, end)[1] == 1 and density(self.csv, start, end) >= 0.01 or genes_inbetween(self.old_annotation, start, end)[0] <= 1 and density(self.csv, start, end) >= 0.1:
                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [self.condition]}
            ### for rev stop sites
            elif row['end'] == 3 and row['strand'] == 1 and corresponding_peak(self.starts_rev, row['location'], row['end'], row['strand']) != None:
                start = row['location']
                end = corresponding_peak(self.starts_rev, row['location'], row['end'], row['strand'])
                if density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 and same_operon(self.csv, row['location'], row['end'], end, 5, row['strand']) == 'Yes':
                    if density(self.csv, start, end) >= 0.1 and cap_term_inbetween(self.cap_term['starts_rev'], start, end) == 'No' or genes_inbetween(self.old_annotation, start, end)[1] == 1 and density(self.csv, start, end) >= 0.01 or genes_inbetween(self.old_annotation, start, end)[0] <= 1 and density(self.csv, start, end) >= 0.1:
                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [self.condition]}
            else:
                pass
            try:
                newrow = pd.DataFrame(data = mydata)
                missing_annotations = pd.concat([missing_annotations, newrow])
            except:
                pass

        df = pd.concat([self.gff, missing_annotations])
        df = df.drop_duplicates()
        df = df.sort_values(['start', 'end'])
        df.to_csv(os.path.join(self.gff_folder, self.condition + '.gff3'), sep = '\t', index = False)


    def supply_remaining_peaks(self):
        missing_annotations = pd.DataFrame({'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []})
        for index, row in self.missing_peaks.iterrows():
            # mydata = {}
            ### by adding a start site to the list of stop site and order the list
            ### we can make an annotation to the nearest corresponding peak.
            ### for forward start sites
            if row['end'] == 5 and row['strand'] == 0 and corresponding_peak(self.stops_fwd, row['location'], row['end'], row['strand']) != None:
                start = row['location']
                end = corresponding_peak(self.stops_fwd, row['location'], row['end'], row['strand'])
                if overlap_gff(self.gff, start, end, row['strand']) == 0 and cap_term_inbetween(self.cap_term['stops_fwd'], start, end) == 'No' and same_operon(self.csv, row['location'], row['end'], end, 3, row['strand']) == 'Yes' and density_full(self.csv, start, end) > 0.05 and density(self.csv, start, end) >= 0.1:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [self.condition]}
            ### for reverse start sites
            elif row['end'] == 5 and row['strand'] == 1 and corresponding_peak(self.stops_rev, row['location'], row['end'], row['strand']) != None:
                end = row['location']
                start = corresponding_peak(self.stops_rev, row['location'], row['end'], row['strand'])
                if overlap_gff(self.gff, start, end, row['strand']) == 0 and cap_term_inbetween(self.cap_term['stops_rev'], start, end) == 'No' and same_operon(self.csv, row['location'], row['end'], start, 3, row['strand']) == 'Yes' and density_full(self.csv, start, end) > 0.05 and density(self.csv, start, end) >= 0.1:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [self.condition]}
            ### for fwd stop sites
            elif row['end'] == 3 and row['strand'] == 0 and corresponding_peak(self.starts_fwd, row['location'], row['end'], row['strand']) != None:
                end = row['location']
                start = corresponding_peak(self.starts_fwd, row['location'], row['end'], row['strand'])
                if overlap_gff(self.gff, start, end, row['strand']) == 0 and cap_term_inbetween(self.cap_term['starts_fwd'], start, end) == 'No' and same_operon(self.csv, row['location'], row['end'], start, 5, row['strand']) == 'Yes' and density_full(self.csv, start, end) > 0.05 and density(self.csv, start, end) >= 0.1:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [self.condition]}
            ### for rev stop sites
            elif row['end'] == 3 and row['strand'] == 1 and corresponding_peak(self.starts_rev, row['location'], row['end'], row['strand']) != None:
                start = row['location']
                end = corresponding_peak(self.starts_rev, row['location'], row['end'], row['strand'])
                if overlap_gff(self.gff, start, end, row['strand']) == 0 and cap_term_inbetween(self.cap_term['starts_rev'], start, end) == 'No' and same_operon(self.csv, row['location'], row['end'], end, 5, row['strand']) == 'Yes' and density_full(self.csv, start, end) > 0.05 and density(self.csv, start, end) >= 0.1:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [start], 'end': [end], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [self.condition]}
            else:
                pass
            try:
                newrow = pd.DataFrame(data = mydata)
                missing_annotations = pd.concat([missing_annotations, newrow])
            except:
                pass

        df = pd.concat([self.gff, missing_annotations])
        df = df.drop_duplicates()
        df = df.sort_values(['start', 'end'])
        df.to_csv(os.path.join(self.gff_folder, self.condition + '.gff3'), sep = '\t', index = False)

    def supply_half_peaks(self):
        missing_annotations = pd.DataFrame({'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []})
        for index, row in self.missing_peaks.iterrows():
            mydata = {}
            ### we adding a start site to the list of stop site and order the list
            ### we can make an annotation to the nearest corresponding peak.
            ### for forward start sites
            if row['end'] == 5 and row['strand'] == 0:
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 5, '+')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if density(self.csv, row['location'], go_to[0]) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and len([c for c in self.cap_term['starts_fwd'] if abs(c - row['location']) <= self.search_distance]) != 0:# or overlap_gff(self.gff, row['location'], go_to[1], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [row['location']], 'end': [go_to[1]], 'score': [1000000], 'strand': ['+'], 'phase': ['start_site_resolved_only'], 'attributes': [self.condition]}
            ### for reverse start sites
            elif row['end'] == 5 and row['strand'] == 1:
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 5, '-')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if density(self.csv, go_to[0], row['location']) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and len([c for c in self.cap_term['starts_rev'] if abs(c - row['location']) <= self.search_distance]) != 0:# or overlap_gff(self.gff, go_to[1], row['location'], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [go_to[1]], 'end': [row['location']], 'score': [1000000], 'strand': ['-'], 'phase': ['start_site_resolved_only'], 'attributes': [self.condition]}
            ### for fwd stop sites
            elif row['end'] == 3 and row['strand'] == 0:
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 3, '+')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if density(self.csv, go_to[0], row['location']) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and len([c for c in self.cap_term['stops_fwd'] if abs(c - row['location']) <= self.search_distance]) != 0:# or overlap_gff(self.gff, go_to[1], row['location'], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [go_to[1]], 'end': [row['location']], 'score': [1000000], 'strand': ['+'], 'phase': ['stop_site_resolved_only'], 'attributes': [self.condition]}
            ### for rev stop sites
            elif row['end'] == 3 and row['strand'] == 1:
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 3, '-')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if density(self.csv, row['location'], go_to[0]) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and len([c for c in self.cap_term['stops_rev'] if abs(c - row['location']) <= self.search_distance]) != 0:# or overlap_gff(self.gff, row['location'], go_to[1], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [row['location']], 'end': [go_to[1]], 'score': [1000000], 'strand': ['-'], 'phase': ['stop_site_resolved_only'], 'attributes': [self.condition]}
            else:
                pass
            try:
                newrow = pd.DataFrame(data = mydata)
                missing_annotations = pd.concat([missing_annotations, newrow])
            except:
                pass


        df = pd.concat([self.gff, missing_annotations], ignore_index = True)
        df = df.drop_duplicates()
        try:
            if self.gff['strand'][0] == '+':
                start_duplicated_list = df.duplicated(subset = ['start', 'strand'], keep = False).tolist()
                stop_duplicated_list = df.duplicated(subset = ['end', 'strand'], keep = False).tolist()
            else:
                start_duplicated_list = df.duplicated(subset = ['end', 'strand'], keep = False).tolist()
                stop_duplicated_list = df.duplicated(subset = ['start', 'strand'], keep = False).tolist()
            df['start_duplicated'] = start_duplicated_list
            df['stop_duplicated'] = stop_duplicated_list
        except:
            print('''No annotations made for {}
            Consider analysing a larger region or lowering cutoff values
            for help type "python3 pipeline.py --help"'''.format(os.path.basename(self.csv_file)))

        index_list = []
        try:
            for index, row in df.iterrows():
                if row['phase'] == 'start_site_resolved_only' and row['start_duplicated'] == True:
                    index_list.append(index)
                elif row['phase'] == 'stop_site_resolved_only' and row['stop_duplicated'] == True:
                    index_list.append(index)
        except:
            pass

        for item in index_list:
            df = df[df.index != item]

        try:
            df = df.drop(['start_duplicated', 'stop_duplicated'], axis = 1)
            df = df.sort_values(['start', 'end'])
        except:
            pass
        df.to_csv(os.path.join(self.gff_folder, self.condition + '.gff3'), sep = '\t', index = False)

    def supply_half_peaks_operons(self):
        missing_annotations = pd.DataFrame({'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []})
        for index, row in self.missing_peaks.iterrows():
            mydata = {}
            ### we adding a start site to the list of stop site and order the list
            ### we can make an annotation to the nearest corresponding peak.
            ### for forward start sites
            if row['end'] == 5 and row['strand'] == 0:
                if len([c for c in self.cap_term['starts_fwd'] if abs(c - row['location']) <= self.search_distance]) != 0:
                    operon_yes = 'Yes'
                else:
                    operon_yes = 'No'
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 5, '+')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if operon_yes == 'Yes' and density(self.csv, row['location'], go_to[0]) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and operon_yes == 'Yes':# or overlap_gff(self.gff, row['location'], go_to[1], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [row['location']], 'end': [go_to[1]], 'score': [1000000], 'strand': ['+'], 'phase': ['start_site_resolved_only'], 'attributes': [self.condition]}
            ### for reverse start sites
            elif row['end'] == 5 and row['strand'] == 1:
                if len([c for c in self.cap_term['starts_rev'] if abs(c - row['location']) <= self.search_distance]) != 0:
                    operon_yes = 'Yes'
                else:
                    operon_yes = 'No'
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 5, '-')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if operon_yes == 'Yes' and density(self.csv, go_to[0], row['location']) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and operon_yes == 'Yes':# or overlap_gff(self.gff, go_to[1], row['location'], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [go_to[1]], 'end': [row['location']], 'score': [1000000], 'strand': ['-'], 'phase': ['start_site_resolved_only'], 'attributes': [self.condition]}
            ### for fwd stop sites
            elif row['end'] == 3 and row['strand'] == 0:
                if len([c for c in self.cap_term['stops_fwd'] if abs(c - row['location']) <= self.search_distance]) != 0:
                    operon_yes = 'Yes'
                else:
                    operon_yes = 'No'
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 3, '+')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if operon_yes == 'Yes' and density(self.csv, go_to[0], row['location']) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and operon_yes == 'Yes':# or overlap_gff(self.gff, go_to[1], row['location'], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [go_to[1]], 'end': [row['location']], 'score': [1000000], 'strand': ['+'], 'phase': ['stop_site_resolved_only'], 'attributes': [self.condition]}
            ### for rev stop sites
            elif row['end'] == 3 and row['strand'] == 1:
                if len([c for c in self.cap_term['stops_rev'] if abs(c - row['location']) <= self.search_distance]) != 0:
                    operon_yes = 'Yes'
                else:
                    operon_yes = 'No'
                go_to = nearest_old_annotation(self.old_annotation, row['location'], 3, '-')
                if go_to == None or highest_peak(self.csv, row['location'], row['end'], row['strand'], self.search_distance)[2] >= row['value']:
                    continue
                if operon_yes == 'Yes' and density(self.csv, row['location'], go_to[0]) > 0.1 and density_immediate(self.csv, row['location'], row['end'], row['strand']) >= 0.1 or abs(row['location'] - go_to[2]) < 50 and operon_yes == 'Yes':# or overlap_gff(self.gff, row['location'], go_to[1], row['strand']) == 0:
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [row['location']], 'end': [go_to[1]], 'score': [1000000], 'strand': ['-'], 'phase': ['stop_site_resolved_only'], 'attributes': [self.condition]}
            else:
                pass
            try:
                newrow = pd.DataFrame(data = mydata)
                missing_annotations = pd.concat([missing_annotations, newrow])
            except:
                pass

        df = pd.concat([self.gff, missing_annotations], ignore_index = True)
        df = df.drop_duplicates()
        try:
            if self.gff['strand'][0] == '+':
                start_duplicated_list = df.duplicated(subset = ['start', 'strand'], keep = False).tolist()
                stop_duplicated_list = df.duplicated(subset = ['end', 'strand'], keep = False).tolist()
            else:
                start_duplicated_list = df.duplicated(subset = ['end', 'strand'], keep = False).tolist()
                stop_duplicated_list = df.duplicated(subset = ['start', 'strand'], keep = False).tolist()
            df['start_duplicated'] = start_duplicated_list
            df['stop_duplicated'] = stop_duplicated_list
        except:
            print('''No annotations made for {}
            Consider analysing a larger region or lowering cutoff values
            for help type "python3 pipeline.py --help"'''.format(os.path.basename(self.csv_file)))

        index_list = []
        try:
            for index, row in df.iterrows():
                if row['phase'] == 'start_site_resolved_only' and row['start_duplicated'] == True:
                    index_list.append(index)
                elif row['phase'] == 'stop_site_resolved_only' and row['stop_duplicated'] == True:
                    index_list.append(index)
        except:
            pass

        for item in index_list:
            df = df[df.index != item]

        try:
            df = df.drop(['start_duplicated', 'stop_duplicated'], axis = 1)
            df = df.sort_values(['start', 'end'])
        except:
            pass

        df = overlap_half(df, self.genome_name, self.condition)
        df.to_csv(os.path.join(self.gff_folder, self.condition + '.gff3'), sep = '\t', index = False)
