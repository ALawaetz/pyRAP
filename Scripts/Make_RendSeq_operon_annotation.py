import pandas as pd
import os
from multiprocessing import Process, Pool
import folders
import sys

base_name = sys.argv[1]
genome_name = sys.argv[2]
operon_cutoff = float(sys.argv[3])


### pyRAP creates annotations by connection start sites and stop sites within operons.
### Peaks therefore need to be assigned to operons before annotations can be made.
### Operon data can be provided to the pipeline by the user by using the --operon flag
### For help type:
### python3 pipeline.py --help

### If no operons are provided, pyRAP creates its own operons by using the Rend-seq file set using the --guide flag
### Operons in this context are defined as regions where the read density between peaks is above the 
### threshold set by the --operon_cutoff flag.
### pyRAP calculates read density between neighboring peaks.
### Read density is calculated as the fraction of non-zero values to the distance 
### between the two peaks multiplied by two (as each genomic coordinate has two read values; a 5’ read and a 3’ read).
### Two peaks belong to seperate operons when the read density between the two peaks is below 
### the threshold set by the --operon_cutoff flag

def YourCode(file):
    if base_name in file:
        name = file
        myheader = {'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []}
        GFF3_operon = pd.DataFrame(data = myheader)

        ### Load dataframe in df_full and omit all rows not containing peaks and save as df
        ### Iterate through each row in df
        ### To asses the read density between peaks we cross reference to df_full
        df_full = pd.read_csv(folders.peak_folder_0 + file, sep = '\t', header = 0)
        df = df_full[df_full['peak'] != '.']
        df = df.reset_index(drop = True)

        try:
            thestrand = max(df['strand'])
        except:
            print('No operons could be made. Try analysing a larger region of genome. Type ctrl+c to exit.')
            sys.exit()

        ### Operons will be named by numbers starting with 1
        operon_name = 1
        for index, row in df.iterrows():
            if row['peak'] == 'START' and index < len(df.index) - 1:
                read_dens = 0
                if index > 0:
                    ### if the read density between this start and the upstream peak is less than the cutoff (--operon_cutoff flag)
                    ### then we have an operon start
                    start = row['location']
                    upstream = df['location'][index - 1]
                    if thestrand == 0:
                        conditions = (df_full['location'] >= upstream) & (df_full['location'] <= start)
                    else:
                        conditions = (df_full['location'] <= upstream) & (df_full['location'] >= start)

                    read_dens = len(df_full[conditions][df_full[conditions]['value'] != 0])/len(df_full[conditions])

                if read_dens < operon_cutoff or index == 0:
                    ### go forth until the density between peaks is less than the cutoff
                    start = row['location']
                    downstream = df['location'][index + 1]
                    if thestrand == 0:
                        conditions = (df_full['location'] >= start) & (df_full['location'] <= downstream)
                    else:
                        conditions = (df_full['location'] <= start) & (df_full['location'] >= downstream)

                    read_dens = len(df_full[conditions][df_full[conditions]['value'] != 0])/len(df_full[conditions])
                    ### if read density between the start peak and the next peak is above the cutoff then see if read density
                    ### between the next two peaks is above the cutoff and so on until read density is below the cutoff
                    if read_dens >= operon_cutoff:
                        n = 0
                        while read_dens >= operon_cutoff and index + n + 2 < len(df.index):
                            n += 1
                            start = df['location'][index + n]
                            stop = df['location'][index + n + 1]
                            if thestrand == 0:
                                conditions = (df_full['location'] >= start) & (df_full['location'] <= stop)
                            else:
                                conditions = (df_full['location'] <= start) & (df_full['location'] >= stop)

                            read_dens = len(df_full[conditions][df_full[conditions]['value'] != 0])/len(df_full[conditions])
                        ### now we are moving out of an operon. If we are at a stop, then draw the operon to here
                        if df['peak'][index + n] == 'STOP':
                            if thestrand == 0:
                                mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [row['location']], 'end': [df['location'][index + n]], 'score': [1000000], 'strand': [row['strand']], 'phase': ['.'], 'attributes': ['ID=Rend-seq_operon-fwd-{}'.format(operon_name)]}
                            else:
                                mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df['location'][index + n]], 'end': [row['location']], 'score': [1000000], 'strand': [row['strand']], 'phase': ['.'], 'attributes': ['ID=Rend-seq_operon-rev-{}'.format(operon_name)]}
                            newrow = pd.DataFrame(data = mydata)
                            GFF3_operon = pd.concat([GFF3_operon, newrow], ignore_index = True)
                            operon_name += 1
                        ### If you are not at a stop, then go back untill you reach a stop and draw the operon to there.
                        else:
                            while df['peak'][index + n] != 'STOP':
                                n -= 1
                                if index + n < 0:
                                    break
                                else:
                                    pass
                            if n > 0:
                                if thestrand == 0:
                                    mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [row['location']], 'end': [df['location'][index + n]], 'score': [1000000], 'strand': [row['strand']], 'phase': ['.'], 'attributes': ['ID=Rend-seq_operon-fwd-{}'.format(operon_name)]}
                                else:
                                    mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df['location'][index + n]], 'end': [row['location']], 'score': [1000000], 'strand': [row['strand']], 'phase': ['.'], 'attributes': ['ID=Rend-seq_operon-rev-{}'.format(operon_name)]}
                                newrow = pd.DataFrame(data = mydata)
                                GFF3_operon = pd.concat([GFF3_operon, newrow], ignore_index = True)
                                operon_name += 1
                            else:
                                pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
        GFF3_operon.to_csv(os.path.join(folders.operon_folder + file + '_GFF3_operon_RendSeq.gff3'), sep = '\t', index = False)
    else:
        pass

if __name__ == '__main__':
    with Pool(2) as p:
        p.map(YourCode, os.listdir(folders.peak_folder_0))
