import pandas as pd
import os
from multiprocessing import Process, Pool
import folders
import sys
import glob

input_folder = sys.argv[1]
output_folder = sys.argv[2]
operon_file = sys.argv[3]
genome_name = sys.argv[4]
primary_operons = sys.argv[5]
processors = int(sys.argv[6])


### Annotations are made by connecting start and stop sites that belong to the same operon.
### Operons are made in the script "Make_RendSeq_operon_annotation.py" using the flags
### --operon_cutoff and --guide
### Operon coordinates can also be given using the --operon flag
### A hybrid option is also possible where pyRAP uses two operon files; the first being either of the two
### options described above and the second operon file is given using the --alternative_operons flag
### For help type
### python3 pipeline.py --help

### This script attempts to assign every peak to an operon.
### pyRAP first tries to assign peaks to operons made using the --operon_cutoff flag or 
### supplied using the --operon flag.
### Unasigned peaks will subsequently be attempted to assign to operons supplied with 
### the --alternative_operons flag (if used).
### Peaks are assigned to operons with which they overlap.
### If a start site doesn't overlap with an operon and no stop site is found between the peak and the 
### downstream operon, then the peak will be assigned to the downstream operon.
### If a stop site doesn't overlap with an operon and no start site is found between the peak and the 
### upstream operon, then the peak will be assigned to the upstream operon.

### Finally, annotations are made between start and stop sites assigned to the same operon.



def YourCode(file):
    df = pd.read_csv(file, sep = '\t', header = 0, usecols = ['location', 'strand', 'peak'])
    df = df.sort_values(['location'], ascending = True, axis = 0)
    thestrand = df['strand'][0]
    if thestrand == 0:
        thestrand = '+'
    else:
        thestrand = '-'
    locations = df['location'][df['peak'] != '.'].tolist()
    peaks = df['peak'][df['peak'] != '.'].tolist()
    del df
    
    if primary_operons == 'None':
        operonfile = glob.glob(folders.operon_folder + '*_FwdandRev.gff3')
        df_operon = pd.read_csv(operonfile[0], sep = '\t', usecols = ['start', 'end', 'strand', 'attributes'])
        df_operon = df_operon.sort_values(['start'], ascending = True, axis = 0)
    else:
        df_operon = pd.read_csv(primary_operons, sep = '\t', names = ['seqID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        df_operon = df_operon.sort_values(['start'], ascending = True, axis = 0)
    

    df_operon = df_operon[df_operon['strand'] == thestrand]

    operon_attributes = []

    for location, peak in zip(locations, peaks):
        overlap = df_operon[(df_operon['start'] <= location) & (df_operon['end'] >= location)].copy()
        overlap = overlap.reset_index(drop = True)
        if len(overlap) == 0 or len(overlap) > 2:
            if peak == 'START':
                operon_attributes.append('Orphan_start')
            else:
                operon_attributes.append('Orphan_stop')
        elif len(overlap) == 1:
            operon_attributes.append(overlap['attributes'][0])
        else:
            if thestrand == '+':
                if peak == 'START':
                    operon_attributes.append(overlap['attributes'][1])
                else:
                    operon_attributes.append(overlap['attributes'][0])
            else:
                if peak == 'START':
                    operon_attributes.append(overlap['attributes'][0])
                else:
                    operon_attributes.append(overlap['attributes'][1])



    revised_operons = []
    if thestrand == '-':
        operon_attributes = operon_attributes[::-1]
    else:
        pass

    
    for i in range(len(operon_attributes)):
        if operon_attributes[i] == 'Orphan_start':
            n = 0
            try:
                while operon_attributes[i + n] == 'Orphan_start':
                    n += 1
                    if operon_attributes[i + n] == 'Orphan_stop':
                        revised_operons.append('Cuckoo')
                    elif operon_attributes[i + n] == 'Orphan_start':
                        pass
                    else:
                        revised_operons.append(operon_attributes[i + n])
            except:
                revised_operons.append('Cuckoo')
        elif operon_attributes[i] == 'Orphan_stop':
            n = 0
            try:
                while operon_attributes[i - n] == 'Orphan_stop':
                    n += 1
                    if operon_attributes[i - n] == 'Orphan_start':
                        revised_operons.append('Cuckoo')
                    elif operon_attributes[i - n] == 'Orphan_stop':
                        pass
                    else:
                        revised_operons.append(operon_attributes[i - n])
            except:
                revised_operons.append('Cuckoo')
        else:
            revised_operons.append(operon_attributes[i])

    
    if thestrand == '-':
        revised_operons = revised_operons[::-1]
    

    df_rend = pd.read_csv(file, sep = '\t', header = 0)
    df_rend = df_rend.sort_values(['location'], ascending = True, axis = 0)
    df_rend = df_rend[df_rend['peak'] != '.']
    if len(revised_operons) != 0:
        df_rend['operon'] = revised_operons
    else:
        df_rend['operon'] = ['.'] * len(df_rend)

    ### read the file once more, but now include all rows
    df_rend_full = pd.read_csv(file, sep = '\t', header = 0)

    ### concat the two dataframes and then drop the duplicates and keep the last added
    df_rend = pd.concat([df_rend_full, df_rend])
    df_rend = df_rend.fillna('.')
    df_rend = df_rend.drop_duplicates(subset = ['location', 'end'], keep = 'last')

    if thestrand == '+':
        df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = True)
    else:
        df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = False)

    df_rend.to_csv(file, sep = '\t', index = False)

    del df_rend_full
    del df_rend          


    ### Draw GFF3 annotation
    name = file
    myheader = {'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []}
    GFF3 = pd.DataFrame(data = myheader)

    ### Load dataframe in df_full and omit all rows not containing peaks and save as df
    ### we will iterate through each row in df
    ### and to asses the read density between peaks we cross reference to df_full
    df_full = pd.read_csv(file, sep = '\t', usecols = ['location', 'value', 'strand', 'peak', 'operon'])

    thestrand = df_full['strand'][0]
    values = df_full['value'].tolist()
    full_location = df_full['location'].tolist()

    df = df_full[df_full['peak'] != '.']

    df_location = df['location'].tolist()
    df_peak = df['peak'].tolist()
    df_operon = df['operon'].tolist()
    del df, df_full

    for i in range(0, len(df_peak)):
        if df_peak[i] == 'START':
            read_dens = 1
            moving_dens = 1
            n = 0
            while i < len(df_peak) - 1 and read_dens >= 0.2 and moving_dens >= 0.1 and df_operon[i] == df_operon[i + n + 1]:

                ### also make sure that the first 25 bp downstream a start and upstream a stop has density above 0.2
                index_start = full_location.index(df_location[i + n])
                index_end = full_location.index(df_location[i + n + 1])
                condition_values = values[index_start: index_end + 1]
                if abs(df_location[i + n] - df_location[i+ n + 1]) > 100:
                    for c in range(0, len(condition_values), 100):
                        for_condition_values = condition_values[c: c + 100]
                        non_zero_values = [z for z in for_condition_values if z != 0]
                        moving_dens = len(non_zero_values)/len(for_condition_values)
                        if moving_dens < 0.1:
                            break
                        else:
                            pass
                else:
                    non_zero_values = [z for z in condition_values if z != 0]
                    read_dens = len(non_zero_values)/len(condition_values)

                if read_dens >= 0.2 and moving_dens >= 0.1 and df_peak[i + n + 1] == 'STOP':
                    if thestrand == 0:
                        mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i]], 'end': [df_location[i + n + 1]], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [name]}
                    else:
                        mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i + n + 1]], 'end': [df_location[i]], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [name]}
                    newrow = pd.DataFrame(data = mydata)
                    GFF3 = pd.concat([GFF3, newrow])
                else:
                    pass
                n += 1
                if i + n == len(df_peak) - 1:
                    break
                else:
                    pass
        else:
            pass

    ### now some peaks will still not have been connected in transcript because they been assigned to different operons
    ### Those such peaks will be connected to the nearest corresponding start or stop if the read density is above 0.2
    ### and only if one of the peaks is assigned to an operon with the ending 'a_half' (peaks that falled outside both
    ### rend-seq and bsg operons). This means that peaks from different operons will still not be connected.

    GFF3_remain = pd.DataFrame(data = myheader)

    ### Load dataframe in df_full and omit all rows not containing peaks and save as df
    ### we will iterate through each row in df
    ### and to asses the read density between peaks we cross reference to df_full

    startlist = GFF3['start'][GFF3['strand'] == '+'].tolist()
    startlist = startlist + GFF3['end'][GFF3['strand'] == '-'].tolist()
    endlist = GFF3['end'][GFF3['strand'] == '+'].tolist()
    endlist = endlist + GFF3['start'][GFF3['strand'] == '-'].tolist()
    for i in range(0, len(df_peak)):
        if df_peak[i] == 'START' and df_location[i] not in startlist:
            read_dens = 1
            moving_dens = 1
            n = 0
            ### make a empty list to append operons. When connecting starts to stops
            ### We can move from operon X to operon X__and_a_half but not to operon X+1
            move_operon = []
            move_operon.append(df_operon[i])
            while i < len(df_peak) - 1 and read_dens >= 0.2 and moving_dens >= 0.1:

                index_start = full_location.index(df_location[i + n])
                index_end = full_location.index(df_location[i + n + 1])
                condition_values = values[index_start: index_end + 1]
                if abs(df_location[i + n] - df_location[i+ n + 1]) > 200:
                    for c in range(0, len(condition_values), 200):
                        for_condition_values = condition_values[c: c + 200]
                        non_zero_values = [z for z in for_condition_values if z != 0]
                        moving_dens = len(non_zero_values)/len(for_condition_values)
                        if moving_dens < 0.1:
                            break
                        else:
                            pass
                else:
                    non_zero_values = [z for z in condition_values if z != 0]
                    read_dens = len(non_zero_values)/len(condition_values)

                if read_dens >= 0.2 and moving_dens >= 0.1 and df_peak[i + n + 1] == 'STOP':
                    if '_a_half' in df_operon[i] or '_a_half' in df_operon[i + n + 1]:
                        if thestrand == 0:
                            mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i]], 'end': [df_location[i + n + 1]], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [name]}
                        else:
                            mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i + n + 1]], 'end': [df_location[i]], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [name]}
                        newrow = pd.DataFrame(data = mydata)
                        GFF3_remain = pd.concat([GFF3_remain, newrow])
                        break
                    else:
                        move_operon.append(df_operon[i + n + 1])
                        if len(list(set(move_operon))) >= 2:
                            break
                        else:
                            pass
                else:
                    pass
                n += 1
                if i + n == len(df_peak) - 1:
                    break
                else:
                    pass
        elif df_peak[i] == 'STOP' and df_location[i] not in endlist:
            read_dens = 1
            moving_dens = 0.2
            n = 0
            move_operon = []
            move_operon.append(df_operon[i])
            while i > 0 and read_dens >= 0.2 and moving_dens >= 0.1:

                index_start = full_location.index(df_location[i - n - 1])
                index_end = full_location.index(df_location[i - n])
                condition_values = values[index_start: index_end + 1]
                if abs(df_location[i - n] - df_location[i - n - 1]) > 200:
                    for c in range(0, len(condition_values), 200):
                        for_condition_values = condition_values[c: c + 200]
                        non_zero_values = [z for z in for_condition_values if z != 0]
                        moving_dens = len(non_zero_values)/len(for_condition_values)
                        if moving_dens < 0.1:
                            break
                        else:
                            pass
                else:
                    non_zero_values = [z for z in condition_values if z != 0]
                    read_dens = len(non_zero_values)/len(condition_values)


                if read_dens >= 0.2 and moving_dens >= 0.1 and df_peak[i - n - 1] == 'START':
                    if '_a_half' in df_operon[i] or '_a_half' in df_operon[i - n - 1]:
                        if thestrand == 0:
                            mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i - n - 1]], 'end': [df_location[i]], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [name]}
                        else:
                            mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i]], 'end': [df_location[i - n - 1]], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [name]}
                        newrow = pd.DataFrame(data = mydata)
                        GFF3_remain = pd.concat([GFF3_remain, newrow])
                        break
                    else:
                        move_operon.append(df_operon[i - n - 1])
                        if len(list(set(move_operon))) >= 2:
                            break
                        else:
                            pass
                else:
                    pass
                n += 1
                if i - n == 0:
                    break
                else:
                    pass
        else:
            pass

    GFF3 = pd.concat([GFF3, GFF3_remain])
    GFF3 = GFF3.drop_duplicates()
    GFF3 = GFF3.sort_values(['start'], axis = 0)
    GFF3 = GFF3.reset_index(drop = True)
    ### omit all annotations less than 15 nucleotides
    droplist = []
    startlist = GFF3['start'].tolist()
    endlist = GFF3['end'].tolist()
    i = 0
    for s, e in zip(startlist, endlist):
        if e - s < 15:
            droplist.append(i)
        else:
            pass
        i += 1
    GFF3 = GFF3.drop(GFF3.index[droplist])
    GFF3.to_csv(os.path.join(output_folder + os.path.basename(file) + '.gff3'), sep = '\t', index = False)

    del GFF3, GFF3_remain


if __name__ == '__main__':
    with Pool(processors) as p:
        p.map(YourCode, sorted(glob.glob(input_folder + '/*csv')))
