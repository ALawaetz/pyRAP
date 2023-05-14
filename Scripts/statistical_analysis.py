import pandas as pd
import os
import glob
import sys
import statistics as stat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys



# # ################################################################################
# # ################################################################################
# ######################### Statistical analysis #################################
# ### All figures are put into the analysis folder or appropriate folders therein.
# ### Figures can also be found in the word document 'Resolving transcript boundaries reveils RNase dependent sRNA processing'
# ### which is found in the base_folder
def lighten_color(color, amount=0.5): ### reference for function: https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
# #

# # ### Figure Read density 2
def read_density(input_folder, list_of_conditions, analysis_folder):
    tracklist = []
    for file in sorted(glob.glob(input_folder + '*.csv')):
        tracklist.append(file)
    merge_dict = {}
    n = 0
    for i in range(0, len(tracklist), 2):
        values = pd.read_csv(tracklist[i], sep = ',', usecols = ['value'])['value'].tolist()
        values += pd.read_csv(tracklist[i + 1], sep = ',', usecols = ['value'])['value'].tolist()
        read_dens = len([i for i in values if i != 0]) / len(values)
        merge_dict[list_of_conditions[n]] = read_dens
        n += 1
    data = []
    for key, value in merge_dict.items():
        data.append(value)

    plt.rcdefaults()
    fig, ax = plt.subplots()
    condition = list_of_conditions
    y_pos = np.arange(len(condition))
    ax.grid(axis = 'x', linestyle = '--', linewidth = 0.5)
    ax.set_axisbelow(True)
    ax.barh(y_pos, data, color = lighten_color('orange', 0.5))
    ax.set_yticks(y_pos, labels=condition)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_title('Read density')
    plt.savefig(os.path.join(analysis_folder, 'Read_density.png'), bbox_inches='tight')
    plt.close()



def collect_all_peaks(list_of_files):
    for file in list_of_files:
        df = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

        starts_fwd = df['start'][(df['strand'] == '+') & (df['phase'] == '.')].tolist()
        starts_rev = df['end'][(df['strand'] == '-') & (df['phase'] == '.')].tolist()

        stops_fwd = df['end'][(df['strand'] == '+') & (df['phase'] == '.')].tolist()
        stops_rev = df['start'][(df['strand'] == '-') & (df['phase'] == '.')].tolist()

        half_resolved = df[df['phase'] != '.'].copy()
        if len(half_resolved) != 0:
            for index, row in half_resolved.iterrows():
                if row['strand'] == '+' and row['attributes'].find('start_site_resolved_only') != -1:
                    starts_fwd += [row['start']]
                elif row['strand'] == '-' and row['attributes'].find('start_site_resolved_only') != -1:
                    starts_rev += [row['end']]
                elif row['strand'] == '+' and row['attributes'].find('stop_site_resolved_only') != -1:
                    stops_fwd += [row['end']]
                elif row['strand'] == '-' and row['attributes'].find('stop_site_resolved_only') != -1:
                    stops_rev += [row['start']]
                else:
                    print('review code in count peaks')

        list_of_lists = [starts_fwd, starts_rev, stops_fwd, stops_rev]
        namelist = ['starts_fwd', 'starts_rev', 'stops_fwd', 'stops_rev']
        adict = {}
        for item, name in zip(list_of_lists, namelist):
            new_list = []
            sorted_unique_list = sorted(list(set(item)))
            for s in sorted_unique_list:
                if len([i for i in sorted_unique_list if i >= s and i <= s + 5]) > 1:
                    pass
                else:
                    new_list.append(s)
            adict[name] = new_list


    return [sorted(list(set(adict['starts_fwd']))), sorted(list(set(adict['starts_rev']))), sorted(list(set(adict['stops_fwd']))), sorted(list(set(adict['stops_rev'])))]


def sort_peaks_by_cutoff(in_file, z_min, not_file, z_max, peaks, end):
    in_csv = pd.read_csv(in_file, sep = '\t', header = 0)
    not_csv = pd.read_csv(not_file, sep = '\t', header = 0)
    sorted_peaks = []
    for p in peaks:
        if len(in_csv[(in_csv['location'] == p) & (in_csv['shift_cutoff'] >= z_min) & (in_csv['end'] == end) | (in_csv['location'] == p) & (in_csv['shift_cutoff_small'] >= z_min) & (in_csv['end'] == end)]) == 1:
            if len(not_csv[(not_csv['location'] >= p - 5) & (not_csv['location'] <= p + 5) & (not_csv['shift_cutoff'] > z_max) & (not_csv['end'] == end) | (not_csv['location'] >= p - 5) & (not_csv['location'] <= p + 5) & (not_csv['shift_cutoff_small'] > z_max) & (not_csv['end'] == end)]) == 0:
                sorted_peaks.append(p)
    return sorted_peaks


####################### Number of start and stop sites #########################
def plot_n_peaks(list_of_files, analysis_folder):
    dict_starts = {}
    dict_stops = {}
    condition = []
    for file in list_of_files:
        all_peaks = collect_all_peaks([file])

        all_starts = all_peaks[0]
        all_starts += all_peaks[1]
        dict_starts[os.path.basename(file)] = len(all_starts)

        all_stops = all_peaks[2]
        all_stops += all_peaks[3]
        dict_stops[os.path.basename(file)] = len(all_stops)

        condition.append(os.path.basename(file)[:-5])

    n_starts = []
    for key, item in dict_starts.items():
        n_starts.append(item)

    n_stops = []
    for key, item in dict_stops.items():
        n_stops.append(item)

    plt.rcdefaults()
    fig, ax = plt.subplots()
    y_pos = np.arange(len(condition))
    height_start = n_starts
    height_stop = n_stops
    ax.grid(axis = 'x', linestyle = '--', linewidth = 0.5)
    ax.set_axisbelow(True)
    ax.barh(y_pos, height_start, color=lighten_color('orange', 0.5))
    ax.barh(y_pos, height_stop, color = 'none', edgecolor = 'blue')
    ax.set_yticks(y_pos, labels=condition)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_title('# start and stop sites')
    ax.legend(['start', 'stop'], loc="upper right")
    plt.savefig(os.path.join(analysis_folder, 'Number_starts_and_stops'), bbox_inches='tight')
    plt.close()
    for c, sta, sto in zip(condition, height_start, height_stop):
        print('{}: {} starts and {} stops'.format(c, sta, sto))


def agreement(compare_list, rend_list, name_list, agree_distance):

    agree = {}

    for compare, rend, name in zip(compare_list, rend_list, name_list):
        agree_list = []
        for i in compare:
            if len([r for r in rend if r >= i - agree_distance and r <= i + agree_distance]) != 0:
                agree_list.append(i)
            else:
                pass
        agree[name] = agree_list

    return agree

def disagreement(compare_list, rend_list, name_list, agree_distance):

    disagree = {}

    for compare, rend, name in zip(compare_list, rend_list, name_list):
        disagree_list = []
        for i in compare:
            if len([r for r in rend if r >= i - agree_distance and r <= i + agree_distance]) != 0:
                pass
            else:
                disagree_list.append(i)
        disagree[name] = disagree_list

    return disagree


############################### Coverage #######################################
### To estimate how much of the transcriptome is covered by Rend-seq
### we asses the fraction of operon start and stop sites that are covered.
### Cap-seq and term-seq is only in one condition and hence don't cover entire transcriptome.
### Operon start and stop sites dont include all start and stop sites (e.g. processed transcripts and alternative start sites)
### But provides an estimate of total coverage.
def plot_coverage(operon_starts_file, operon_stops_file, rend_file, agree_distance, analysis_folder):
    operon_starts = pd.read_csv(operon_starts_file, sep = ',', header = 0)
    operon_stops = pd.read_csv(operon_stops_file, sep = ',', header = 0)

    operon_starts_fwd = operon_starts['location'][operon_starts['strand'] == '+'].tolist()
    operon_starts_rev = operon_starts['location'][operon_starts['strand'] == '-'].tolist()
    operon_stops_fwd = operon_stops['location'][operon_stops['strand'] == '+'].tolist()
    operon_stops_rev = operon_stops['location'][operon_stops['strand'] == '-'].tolist()

    rend_peaks = collect_all_peaks([rend_file])
    rend_starts_fwd = rend_peaks[0]
    rend_starts_rev = rend_peaks[1]
    rend_stops_fwd = rend_peaks[2]
    rend_stops_rev = rend_peaks[3]

    operon_list = [operon_starts_fwd, operon_starts_rev, operon_stops_fwd, operon_stops_rev]
    rend_list = [rend_starts_fwd, rend_starts_rev, rend_stops_fwd, rend_stops_rev]
    name_list = ['operon_starts_fwd', 'operon_starts_rev', 'operon_stops_fwd', 'operon_stops_rev']

    agree = agreement(operon_list, rend_list, name_list, agree_distance)

    percentage_agreed_starts = (len(agree['operon_starts_fwd']) + len(agree['operon_starts_rev'])) / (len(operon_starts_fwd) + len(operon_starts_rev)) * 100
    percentage_agreed_stops = (len(agree['operon_stops_fwd']) + len(agree['operon_stops_rev'])) / (len(operon_stops_fwd) + len(operon_stops_rev)) * 100

    plt.rcdefaults()
    fig, ax = plt.subplots()
    condition = ['Starts', 'Stops']
    y_pos = np.arange(len(condition))
    height = [percentage_agreed_starts, percentage_agreed_stops]
    ax.set_yticks(y_pos, labels=condition)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_title('Percentage of operons covered by Rend-seq')
    ax.grid(axis = 'x', linestyle = '--', linewidth = 0.5)
    ax.set_axisbelow(True)
    ax.barh(y_pos, height, color=lighten_color('orange', 0.5), height= 0.5)
    ax.set_xlim([0, 100])
    ax.set(xlabel='%')
    plt.savefig(os.path.join(analysis_folder, 'Coverage'), bbox_inches='tight')
    plt.close()
    print('{} % of operon starts is covered by Rend-seq'.format(int(percentage_agreed_starts)))
    print('{} % of operon stops is covered by Rend-seq'.format(int(percentage_agreed_stops)))




#################### Agreement with Cap/Term seq ###############################
def plot_agreement_w_CapTerm(cap_file, term_file, rend_file, analysis_folder):
    cap_starts = pd.read_csv(cap_file, sep = ',', header = 0)
    term_stops = pd.read_csv(term_file, sep = ',', header = 0)

    cap_starts_fwd = cap_starts['location'][cap_starts['strand'] == '+'].tolist()
    cap_starts_rev = cap_starts['location'][cap_starts['strand'] == '-'].tolist()
    term_stops_fwd = term_stops['location'][term_stops['strand'] == '+'].tolist()
    term_stops_rev = term_stops['location'][term_stops['strand'] == '-'].tolist()

    rend_peaks = collect_all_peaks([rend_file])
    rend_starts_fwd = rend_peaks[0]
    rend_starts_rev = rend_peaks[1]
    rend_stops_fwd = rend_peaks[2]
    rend_stops_rev = rend_peaks[3]

    cap_term_list = [cap_starts_fwd, cap_starts_rev, term_stops_fwd, term_stops_rev]
    rend_list = [rend_starts_fwd, rend_starts_rev, rend_stops_fwd, rend_stops_rev]
    name_list = ['starts_fwd', 'starts_rev', 'stops_fwd', 'stops_rev']

    CapTerm_in_Rend_0 = agreement(cap_term_list, rend_list, name_list, 0)
    CapTerm_in_Rend_5 = agreement(cap_term_list, rend_list, name_list, 5)
    CapTerm_in_Rend_10 = agreement(cap_term_list, rend_list, name_list, 10)

    Rend_in_CapTerm_0 = agreement(rend_list, cap_term_list, name_list, 0)
    Rend_in_CapTerm_5 = agreement(rend_list, cap_term_list, name_list, 5)
    Rend_in_CapTerm_10 = agreement(rend_list, cap_term_list, name_list, 10)

    percentage_Cap_agreed_starts_0 = (len(CapTerm_in_Rend_0['starts_fwd']) + len(CapTerm_in_Rend_0['starts_rev'])) / (len(cap_starts_fwd) + len(cap_starts_rev)) * 100
    percentage_Cap_agreed_starts_5 = (len(CapTerm_in_Rend_5['starts_fwd']) + len(CapTerm_in_Rend_5['starts_rev'])) / (len(cap_starts_fwd) + len(cap_starts_rev)) * 100
    percentage_Cap_agreed_starts_10 = (len(CapTerm_in_Rend_10['starts_fwd']) + len(CapTerm_in_Rend_10['starts_rev'])) / (len(cap_starts_fwd) + len(cap_starts_rev)) * 100

    percentage_Cap_agreed_stops_0 = (len(CapTerm_in_Rend_0['stops_fwd']) + len(CapTerm_in_Rend_0['stops_rev'])) / (len(term_stops_fwd) + len(term_stops_rev)) * 100
    percentage_Cap_agreed_stops_5 = (len(CapTerm_in_Rend_5['stops_fwd']) + len(CapTerm_in_Rend_5['stops_rev'])) / (len(term_stops_fwd) + len(term_stops_rev)) * 100
    percentage_Cap_agreed_stops_10 = (len(CapTerm_in_Rend_10['stops_fwd']) + len(CapTerm_in_Rend_10['stops_rev'])) / (len(term_stops_fwd) + len(term_stops_rev)) * 100


    percentage_Rend_agreed_starts_0 = (len(Rend_in_CapTerm_0['starts_fwd']) + len(Rend_in_CapTerm_0['starts_rev'])) / (len(rend_starts_fwd) + len(rend_starts_rev)) * 100
    percentage_Rend_agreed_starts_5 = (len(Rend_in_CapTerm_5['starts_fwd']) + len(Rend_in_CapTerm_5['starts_rev'])) / (len(rend_starts_fwd) + len(rend_starts_rev)) * 100
    percentage_Rend_agreed_starts_10 = (len(Rend_in_CapTerm_10['starts_fwd']) + len(Rend_in_CapTerm_10['starts_rev'])) / (len(rend_starts_fwd) + len(rend_starts_rev)) * 100

    percentage_Rend_agreed_stops_0 = (len(Rend_in_CapTerm_0['stops_fwd']) + len(Rend_in_CapTerm_0['stops_rev'])) / (len(rend_stops_fwd) + len(rend_stops_rev)) * 100
    percentage_Rend_agreed_stops_5 = (len(Rend_in_CapTerm_5['stops_fwd']) + len(Rend_in_CapTerm_5['stops_rev'])) / (len(rend_stops_fwd) + len(rend_stops_rev)) * 100
    percentage_Rend_agreed_stops_10 = (len(Rend_in_CapTerm_10['stops_fwd']) + len(Rend_in_CapTerm_10['stops_rev'])) / (len(rend_stops_fwd) + len(rend_stops_rev)) * 100


    labels = [r'$\Delta$ 0bp', r'$\Delta$ 5bp', r'$\Delta$ 10bp']
    Cap_in_Rend = [percentage_Cap_agreed_starts_0, percentage_Cap_agreed_starts_5, percentage_Cap_agreed_starts_10]
    Rend_in_Cap = [percentage_Rend_agreed_starts_0, percentage_Rend_agreed_starts_5, percentage_Rend_agreed_starts_10]

    x = np.arange(len(labels))
    width = 0.35

    plt.rcdefaults()
    fig, ax = plt.subplots()
    ax.bar(x - width/2, Cap_in_Rend, width, label='Cap_in_Rend', color=lighten_color('orange', 0.5))
    ax.bar(x + width/2, Rend_in_Cap, width, label='Rend_in_Cap', color='blue')
    ax.set_ylabel('%')
    ax.set_title('Agreement w. Cap-seq')
    ax.set_xticks(x, labels)
    ax.legend()
    ax.grid(axis = 'y', linestyle = '--', linewidth = 0.5)
    ax.set_axisbelow(True)
    ax.set_ylim([0, 100])
    plt.savefig(os.path.join(analysis_folder, 'Agreement_CapSeq'), bbox_inches='tight')
    plt.close()

    Cap_in_Rend = [percentage_Cap_agreed_stops_0, percentage_Cap_agreed_stops_5, percentage_Cap_agreed_stops_10]
    Rend_in_Cap = [percentage_Rend_agreed_stops_0, percentage_Rend_agreed_stops_5, percentage_Rend_agreed_stops_10]
    fig, ax = plt.subplots()
    ax.bar(x - width/2, Cap_in_Rend, width, label='Term_in_Rend', color=lighten_color('orange', 0.5))
    ax.bar(x + width/2, Rend_in_Cap, width, label='Rend_in_Term', color='blue')
    ax.set_ylabel('%')
    ax.set_title('Agreement w. Term-seq')
    ax.set_xticks(x, labels)
    ax.legend()
    ax.grid(axis = 'y', linestyle = '--', linewidth = 0.5)
    ax.set_axisbelow(True)
    ax.set_ylim([0, 100])
    plt.savefig(os.path.join(analysis_folder, 'Agreement_TermSeq'), bbox_inches='tight')
    plt.close()

    # print('{} % of Cap-seq starts is covered by Rend-seq with 0bp discrepancy'.format(int(percentage_Cap_agreed_starts_0)))
    print('{} % of Cap-seq starts is covered by Rend-seq with 5bp discrepancy'.format(int(percentage_Cap_agreed_starts_5)))
    # print('{} % of Cap-seq starts is covered by Rend-seq with 10bp discrepancy'.format(int(percentage_Cap_agreed_starts_10)))
    # print('{} % of Term-seq stops is covered by Rend-seq with 0bp discrepancy'.format(int(percentage_Cap_agreed_stops_0)))
    print('{} % of Term-seq stops is covered by Rend-seq with 5bp discrepancy'.format(int(percentage_Cap_agreed_stops_5)))
    # print('{} % of Term-seq stops is covered by Rend-seq with 10bp discrepancy'.format(int(percentage_Cap_agreed_stops_10)))

    # print('{} % of Rend-seq starts is found in Cap-seq with 0bp discrepancy'.format(int(percentage_Rend_agreed_starts_0)))
    print('{} % of Rend-seq starts is found in Cap-seq with 5bp discrepancy'.format(int(percentage_Rend_agreed_starts_5)))
    # print('{} % of Rend-seq starts is found in Cap-seq with 10bp discrepancy'.format(int(percentage_Rend_agreed_starts_10)))
    # print('{} % of Rend-seq stops is found in Term-seq with 0bp discrepancy'.format(int(percentage_Rend_agreed_stops_0)))
    print('{} % of Rend-seq stops is found in Term-seq with 5bp discrepancy'.format(int(percentage_Rend_agreed_stops_5)))
    # print('{} % of Rend-seq stops is found in Term-seq with 10bp discrepancy'.format(int(percentage_Rend_agreed_stops_10)))



def distance(stat_file, distance_folder, feature_file):
    df = pd.read_csv(stat_file, sep = '\t', header = 0)
    df = df[df['phase'] == '.']
    fwd = df[df['strand'] == '+'].copy()
    fwd = fwd.drop_duplicates(subset=['feature', 'start', 'strand'])
    rev = df[df['strand'] == '-'].copy()
    rev = rev.drop_duplicates(subset=['feature', 'end', 'strand'])
    df = pd.concat([fwd, rev])
    feature_df = pd.read_csv(feature_file, sep = '\t', header = 0)
    feature = feature_df['feature'].tolist()
    dist_id = feature_df['dist_id'].tolist()
    sRNA_type = feature_df['sRNA_type'].tolist()
    titles = feature_df['title'].tolist()

    for feat, id, type, title in zip(feature, dist_id, sRNA_type, titles):
        if type != 'na':
            try:
                data = df[id][(df['feature'] == feat) & (df['sRNA_type'] == type)]
                plot = plt.hist(data, bins=range(-100, 100, 1))
                plt.xlabel('nt')
                plt.title(title)
                plt.savefig(os.path.join(distance_folder, title))
                plt.close()
            except:
                print('The combination of {}, and {} is not found in annotations'.format(feat, type))
        else:
            try:
                data = df[id][(df['feature'] == feat)]
                plot = plt.hist(data, bins=range(-100, 100, 1))
                plt.xlabel('nt')
                plt.title(title)
                plt.savefig(os.path.join(distance_folder, title))
                plt.close()
            except:
                print('The combination of {}, and {} is not found in annotations'.format(feat, type))


### do as above
### provide a csv file denoting what conditions should be compared to each other
def differentially_expressed_peaks(gff_folder, condition_compare_file, csv_folder, z_min, z_max, agree_distance, density_cutoff, DEG_Overlap_folder): #, differentially_expressed_genes_All_folder, differentially_expressed_genes_HighConfidence_folder, DEG_HIGH_Dens_folder, DEG_Overlap_lessStringent_folder):
    compare_df = pd.read_csv(condition_compare_file, sep = '\t', header = 0)
    peaks_in = compare_df['peaks_in'].tolist()
    peaks_not_in = compare_df['peaks_not_in'].tolist()

    for in_file, not_file in zip(peaks_in, peaks_not_in):

        in_peaks = collect_all_peaks([gff_folder + '/' + in_file + '.gff3'])
        not_peaks = collect_all_peaks([gff_folder + '/' + not_file + '.gff3'])
        df = pd.read_csv(gff_folder + '/' + in_file + '.gff3', sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

        phase = []
        z = 1

        namelist = ['starts_fwd', 'starts_rev', 'stops_fwd', 'stops_rev']
        disagreem = disagreement(in_peaks, not_peaks, namelist, agree_distance)
        df_All_DEP = pd.DataFrame()
        for key, item in disagreem.items():
            if z == 1:
                sorted_in_peaks = item
                start_fwd_df = pd.DataFrame()
                for start_fwd in sorted_in_peaks:
                    start_fwd_df = pd.concat([start_fwd_df, df[(df['start'] == start_fwd) & (df['strand'] == '+')]])
                start_fwd_df = start_fwd_df.drop_duplicates(subset = ['start', 'strand'])
                df_All_DEP = pd.concat([df_All_DEP, start_fwd_df])
                phase += ['unique_start'] * len(start_fwd_df)
            elif z == 2:
                sorted_in_peaks = item
                start_rev_df = pd.DataFrame()
                for start_rev in sorted_in_peaks:
                    start_rev_df = pd.concat([start_rev_df, df[(df['end'] == start_rev) & (df['strand'] == '-')]])
                start_rev_df = start_rev_df.drop_duplicates(subset = ['start', 'strand'])
                df_All_DEP = pd.concat([df_All_DEP, start_rev_df])
                phase += ['unique_start'] * len(start_rev_df)
            elif z == 3:
                sorted_in_peaks = item
                stop_fwd_df = pd.DataFrame()
                for stop_fwd in sorted_in_peaks:
                    stop_fwd_df = pd.concat([stop_fwd_df, df[(df['end'] == stop_fwd) & (df['strand'] == '+')]])
                stop_fwd_df = stop_fwd_df.drop_duplicates(subset = ['start', 'strand'])
                df_All_DEP = pd.concat([df_All_DEP, stop_fwd_df])
                phase += ['unique_stop'] * len(stop_fwd_df)
            elif z == 4:
                sorted_in_peaks = item
                stop_rev_df = pd.DataFrame()
                for stop_rev in sorted_in_peaks:
                    stop_rev_df = pd.concat([stop_rev_df, df[(df['start'] == stop_rev) & (df['strand'] == '-')]])
                stop_rev_df = stop_rev_df.drop_duplicates(subset = ['start', 'strand'])
                df_All_DEP = pd.concat([df_All_DEP, stop_rev_df])
                phase += ['unique_stop'] * len(stop_rev_df)
            else:
                pass
            z += 1
        df_All_DEP['phase'] = phase
        # try:
        #     df_All_DEP = df_All_DEP.sort_values(['feature', 'start', 'end'])
        #     df_All_DEP = df_All_DEP.drop_duplicates(['start', 'end', 'strand'])
        #     df_All_DEP.to_csv(os.path.join(differentially_expressed_genes_All_folder, 'annotations_in_{}_not_in_{}'.format(in_file, not_file)), sep = '\t', index = False, header = None)
        # except:
        #     print('''No peaks were found in {} that were not also found in {}
        #     Consider lowering cutoff values with --z_min_compare and --z_max_compare
        #     Type python3 pipeline.py --help for more information'''.format(in_file, not_file))

        new_df = pd.DataFrame()
        sorted_in_peaks = []
        phase = []
        z = 1
        for peaks in in_peaks:
            if z == 1:
                sorted_in_peaks = sort_peaks_by_cutoff(csv_folder + '/' + in_file + '_fwd.csv', z_min, csv_folder + '/' + not_file + '_fwd.csv', z_max, peaks, 5)
                for start_fwd in sorted_in_peaks:
                    new_df = pd.concat([new_df, df[(df['start'] == start_fwd) & (df['strand'] == '+')]])
                    phase += ['unique_start'] * len(df[(df['start'] == start_fwd) & (df['strand'] == '+')])
            elif z == 2:
                sorted_in_peaks = sort_peaks_by_cutoff(csv_folder + '/' + in_file + '_rev.csv', z_min, csv_folder + '/' + not_file + '_rev.csv', z_max, peaks, 5)
                for start_rev in sorted_in_peaks:
                    new_df = pd.concat([new_df, df[(df['end'] == start_rev) & (df['strand'] == '-')]])
                    phase += ['unique_start'] * len(df[(df['end'] == start_rev) & (df['strand'] == '-')])
            elif z == 3:
                sorted_in_peaks = sort_peaks_by_cutoff(csv_folder + '/' + in_file + '_fwd.csv', z_min, csv_folder + '/' + not_file + '_fwd.csv', z_max, peaks, 3)
                for stop_fwd in sorted_in_peaks:
                    new_df = pd.concat([new_df, df[(df['end'] == stop_fwd) & (df['strand'] == '+')]])
                    phase += ['unique_stop'] * len(df[(df['end'] == stop_fwd) & (df['strand'] == '+')])
            elif z == 4:
                sorted_in_peaks = sort_peaks_by_cutoff(csv_folder + '/' + in_file + '_rev.csv', z_min, csv_folder + '/' + not_file + '_rev.csv', z_max, peaks, 3)
                for stop_rev in sorted_in_peaks:
                    new_df = pd.concat([new_df, df[(df['start'] == stop_rev) & (df['strand'] == '-')]])
                    phase += ['unique_stop'] * len(df[(df['start'] == stop_rev) & (df['strand'] == '-')])
            else:
                pass
            z += 1

        new_df['phase'] = phase

        # try:
        #     new_df = new_df.sort_values(['feature', 'start', 'end'])
        #     new_df = new_df.drop_duplicates(['start', 'end', 'strand'])
        #     new_df.to_csv(os.path.join(differentially_expressed_genes_HighConfidence_folder, 'annotations_in_{}_not_in_{}'.format(in_file, not_file)), sep = '\t', index = False, header = None)
        # except:
        #     print('''No high confidence peaks were found in {} that were not also found in {}
        #     Consider lowering cutoff values with --z_min_compare and --z_max_compare
        #     Type python3 pipeline.py --help for more information'''.format(in_file, not_file))

        if len(new_df) != 0:
            csv_fwd = pd.read_csv(csv_folder + '/' + in_file + '_fwd.csv', sep = '\t', header = 0)
            csv_rev = pd.read_csv(csv_folder + '/' + in_file + '_rev.csv', sep = '\t', header = 0)
            csv = pd.concat([csv_fwd, csv_rev])
            high_density = pd.DataFrame()
            for index, row in new_df.iterrows():
                if row['strand'] == '+':
                    thestrand = 0
                else:
                    thestrand = 1
                values = csv['value'][(csv['location'] >= row['start']) & (csv['location'] <= row['end']) & (csv['strand'] == thestrand)].tolist()
                density = len([v for v in values if v != 0]) / len(values)
                if density < density_cutoff:
                    pass
                else:
                    high_density = pd.concat([high_density, new_df[new_df.index == index]])

            # if len(high_density) != 0:
            #     high_density = high_density.drop_duplicates(['start', 'end', 'strand'])
            #     high_density.to_csv(os.path.join(DEG_HIGH_Dens_folder, 'annotations_in_{}_not_in_{}'.format(in_file, not_file)), sep = '\t', index = False, header = None)
            # else:
            #     print('''No HIGH density annotations were found in {} that were not also found in {}
            #     Consider lowering cutoff values with --density_compare
            #     Type python3 pipeline.py --help for more information'''.format(in_file, not_file))


            not_gff = pd.read_csv(gff_folder + '/' + not_file + '.gff3', sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            with_overlap = pd.DataFrame()
            for index, row in new_df.iterrows():
                if row['strand'] == '+' and row['phase'] == 'unique_start':
                    overlap = not_gff[(not_gff['start'] < row['start']) & (row['start'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
                elif row['strand'] == '-' and row['phase'] == 'unique_start':
                    overlap = not_gff[(not_gff['start'] < row['end']) & (row['end'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
                elif row['strand'] == '+' and row['phase'] == 'unique_stop':
                    overlap = not_gff[(not_gff['start'] < row['end']) & (row['end'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
                elif row['strand'] == '-' and row['phase'] == 'unique_stop':
                    overlap = not_gff[(not_gff['start'] < row['start']) & (row['start'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
                else:
                    print('review code DEG overlap')

                if len(overlap) > 0:
                    with_overlap = pd.concat([with_overlap, new_df[new_df.index == index]])

            if len(with_overlap) != 0:
                with_overlap = with_overlap.drop_duplicates(['start', 'end', 'strand'])
                with_overlap.to_csv(os.path.join(DEG_Overlap_folder, 'annotations_in_{}_not_in_{}'.format(in_file, not_file)), sep = '\t', index = False, header = None)
            else:
                print('''No overlap annotations were found in {} that were not also found in {}
                Consider lowering cutoff values with --density_compare
                Type python3 pipeline.py --help for more information'''.format(in_file, not_file))

            # with_overlap = pd.DataFrame()
            # for index, row in df_All_DEP.iterrows():
            #     if row['strand'] == '+' and row['phase'] == 'unique_start':
            #         overlap = not_gff[(not_gff['start'] < row['start']) & (row['start'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
            #     elif row['strand'] == '-' and row['phase'] == 'unique_start':
            #         overlap = not_gff[(not_gff['start'] < row['end']) & (row['end'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
            #     elif row['strand'] == '+' and row['phase'] == 'unique_stop':
            #         overlap = not_gff[(not_gff['start'] < row['end']) & (row['end'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
            #     elif row['strand'] == '-' and row['phase'] == 'unique_stop':
            #         overlap = not_gff[(not_gff['start'] < row['start']) & (row['start'] < not_gff['end']) & (not_gff['strand'] == row['strand']) & (not_gff['phase'] == '.')].copy()
            #     else:
            #         print('review code DEG overlap')

            #     if len(overlap) > 0:
            #         with_overlap = pd.concat([with_overlap, df_All_DEP[df_All_DEP.index == index]])

            # if len(with_overlap) != 0:
            #     with_overlap = with_overlap.drop_duplicates(['start', 'end', 'strand'])
            #     with_overlap.to_csv(os.path.join(DEG_Overlap_lessStringent_folder, 'annotations_in_{}_not_in_{}'.format(in_file, not_file)), sep = '\t', index = False, header = None)
            # else:
            #     print('''No overlap annotations were found in {} that were not also found in {}
            #     Consider lowering cutoff values with --density_compare
            #     Type python3 pipeline.py --help for more information'''.format(in_file, not_file))
        else:
            print('''No overlap annotations were found in {} that were not also found in {}
                Consider lowering cutoff values with --density_compare
                Type python3 pipeline.py --help for more information'''.format(in_file, not_file))
