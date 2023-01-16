import pandas as pd
import os
from multiprocessing import Process, Pool
import folders
import sys

base_name = sys.argv[1]
genome_name = sys.argv[2]
operon_cutoff = float(sys.argv[3])



## Before we can analyse our wig files we need an operon annotation file in order
## to know what starts sites connect to what stop sites. We will use WT LB condition
## from Lalanne et al. 2018 because it has the highest read count and the broadest
## expression profile.

## An operon starts when the read density between a start and the upstream peak is less than 0.8
## An operon ends when the read density between a stop and the downstream peak is less than 0.8
## Read density is the fraction of non-zero values in a window to the total number of positions in that window.
def YourCode(file):
# for file in os.listdir():
    if base_name in file:
        name = file
        myheader = {'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []}
        GFF3_operon = pd.DataFrame(data = myheader)

        ### Load dataframe in df_full and omit all rows not containing peaks and save as df
        ### we will iterate through each row in df
        ### and to asses the read density between peaks we cross reference to df_full
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
                    ### if the read density between this start and the upstream peak is less than 0.8 then we have an operon start
                    start = row['location']
                    upstream = df['location'][index - 1]
                    if thestrand == 0:
                        conditions = (df_full['location'] >= upstream) & (df_full['location'] <= start)
                    else:
                        conditions = (df_full['location'] <= upstream) & (df_full['location'] >= start)

                    read_dens = len(df_full[conditions][df_full[conditions]['value'] != 0])/len(df_full[conditions])

                if read_dens < operon_cutoff or index == 0:
                    ### go forth until the density between peaks is < 0.5
                    start = row['location']
                    downstream = df['location'][index + 1]
                    if thestrand == 0:
                        conditions = (df_full['location'] >= start) & (df_full['location'] <= downstream)
                    else:
                        conditions = (df_full['location'] <= start) & (df_full['location'] >= downstream)

                    read_dens = len(df_full[conditions][df_full[conditions]['value'] != 0])/len(df_full[conditions])
                    ### if read density between the start peak and the next peak is above 0.8 then see if read density
                    ### between the next two peaks is above 0.8 and so on until read density is below 0.8
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
