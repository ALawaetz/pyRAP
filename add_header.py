import os
import pandas as pd
import glob

thefolder = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main/Pipeline/Annotation_files/'

for file in glob.glob(thefolder + '*.gff3'):
    df = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    df.to_csv(os.path.join(thefolder + os.path.basename(file)), sep = '\t', index = False)
