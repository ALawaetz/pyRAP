import os
import pandas as pd
import sys
import glob

thefolder = sys.argv[1]


for file in glob.glob(thefolder + '/*.gff3'):
    df = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    df['attributes'] = [os.path.basename(file)[:-5]] * len(df)
    df.to_csv(file, sep = '\t', index = False)
