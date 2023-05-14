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



def YourCode(file):
    df = pd.read_csv(file, sep = '\t', header = 0, usecols = ['location', 'strand', 'peak']) 
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
        overlap = df_operon[(df_operon['start'] <= location) & (df_operon['end'] >= location)]
        if len(overlap) == 0:
            if peak == 'START':
                operon_attributes.append('Orphan_start')
            else:
                operon_attributes.append('Orphan_stop')
        elif len(overlap) == 1:
            operon_attributes.append(overlap['attributes'][0])
        elif len(overlap) == 2:
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
        else:
            print('Check your operon file. There are one or more instances where more than two operons overlap on the same strand')
            sys.exit()

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
    

    df_rend = pd.read_csv(file, sep = '\t', header = 0)
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

    if thestrand == 0:
        df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = True)
    else:
        df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = False)

    df_rend.to_csv(file, sep = '\t', index = False)

    del df_rend_full
                    


        




    
    





if __name__ == '__main__':
    with Pool(processors) as p:
        p.map(YourCode, sorted(glob.glob(input_folder + '*csv')))