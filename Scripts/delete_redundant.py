import pandas as pd
import os
from multiprocessing import Process, Pool
import folders
import sys
import glob

inputfolder = sys.argv[1]
outputfolder = sys.argv[2]
processors = int(sys.argv[3])


def YourCode(file):
    df = pd.read_csv(file, sep = '\t', header = 0)
    halfs = df[df['phase'] != '.'].copy()
    for indexx, row in halfs.iterrows():
        if row['phase'] == 'start_site_resolved_only':
            if row['strand'] == '+':
                if row['start'] in df['start'][(df['strand'] == '+') & (df['phase'] == '.')].tolist():
                    df = df[df.index != indexx]
                else:
                    pass
            else:
                if row['end'] in df['end'][(df['strand'] == '-') & (df['phase'] == '.')].tolist():
                    df = df[df.index != indexx]
                else:
                    pass
        else:
            if row['strand'] == '+':
                if row['end'] in df['end'][(df['strand'] == '+') & (df['phase'] == '.')].tolist():
                    df = df[df.index != indexx]
                else:
                    pass
            else:
                if row['start'] in df['start'][(df['strand'] == '-') & (df['phase'] == '.')].tolist():
                    df = df[df.index != indexx]
                else:
                    pass
    
    df.to_csv(os.path.join(outputfolder, os.path.basename(file)), sep = '\t', index = False)


if __name__ == '__main__':
    with Pool(processors) as p:
        p.map(YourCode, sorted(glob.glob(inputfolder + '/*.gff3')))