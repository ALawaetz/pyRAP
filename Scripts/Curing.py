import os, glob
import pandas as pd
import sys
import folders

wig_folder = sys.argv[1]
genome_size = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])


# ### Wig files are two dimensional datasets containing an x-coordinate, denoting the position in the genome and a y-coordinate which is the number of reads
# ### mapping to this position. The Wig files are produced from aligned files (bowtie standard output) by doing a 'pile up' to count the total number of reads
# ### mapping to each position in the genome (pipeline: https://github.com/jblalanne/Rend_seq_core_scripts). Positions with no reads are omitted in original Wig files.
#
# ############################# Curing ###########################################
# ### add missing zeros to genomic locations with no expression (necessary for downstream analysis)
# ## and add columns names to location and value

# print('adding missing zeros')

# os.chdir(wig_folder)
for file in glob.glob(wig_folder + '/*.wig'):
    with open(file, 'r') as f:
        text = f.read().splitlines(True)
        f.close()
        n = 0
        for item in text:
            alist = item.split('\t')
            try:
                float(alist[0])
                break
            except ValueError:
                n += 1

    df = pd.read_csv(file, sep = '\t', names = ('location', 'value'), skiprows = n)
    df = df[(df['location'] >= start) & (df['location'] <= end)]
    df = df.set_index('location').reindex(range(start, end + 1)).fillna(0).reset_index()
    df.to_csv(os.path.join(folders.cured_folder, os.path.basename(file) + '.csv'), index = False)

os.chdir(folders.cured_folder)

## Wig files contain in their file name information about their strand identity (forward or reverse) and their end identity (5' or 3').
## Add strand identifier (forward or reverse) as column in files
## 0 = +
## 1 = -

### print('adding strand information')
for file in glob.glob(folders.cured_folder + '*_[3-5]f*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['strand'] = [0]*len(df.index)
    df.to_csv(file, index = False)
    del df
for file in glob.glob(folders.cured_folder + '*_[3-5]r*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['strand'] = [1]*len(df.index)
    df.to_csv(file, index = False)
    del df
for file in glob.glob(folders.cured_folder + '*_[3-5]_f*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['strand'] = [0]*len(df.index)
    df.to_csv(file, index = False)
    del df
for file in glob.glob(folders.cured_folder + '*_[3-5]_r*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['strand'] = [1]*len(df.index)
    df.to_csv(file, index = False)
    del df

## Add end (5' or 3') identifier as column in files

### print('''adding 'end information''')
for file in glob.glob(folders.cured_folder + '*_3[f,r]*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['end'] = [3]*len(df.index)
    df.to_csv(file, index = False)
    del df
for file in glob.glob(folders.cured_folder + '*_5[f,r]*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['end'] = [5]*len(df.index)
    df.to_csv(file, index = False)
    del df
for file in glob.glob(folders.cured_folder + '*_3_[f,r]*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['end'] = [3]*len(df.index)
    df.to_csv(file, index = False)
    del df
for file in glob.glob(folders.cured_folder + '*_5_[f,r]*'):
    df = pd.read_csv(file, sep = ',', header = 0)
    df['end'] = [5]*len(df.index)
    df.to_csv(file, index = False)
    del df

### Wig files come in groups of four. Each condition has a 5' forward file, a 3' forward file, a 5' reverse file, and a 3' reverse file.

### sort folder alphabetically
### a folder sorted alphabetically will list each condition grouped together and in the order 3f, 3r, 5f, 5r
tracklist = sorted(os.listdir())
#
# ### merge 5' and 3' files from each condition into one file according to strand
# ### and name the merged file after the first file of the merger
### print('merge files from each condition according to strand')
i = 0
while i < len(tracklist):
    df = pd.read_csv(folders.cured_folder + tracklist[i], sep = ',', header = 0)
    df = pd.concat([df, pd.read_csv(folders.cured_folder + tracklist[i+2], sep = ',', header = 0)])
    df.to_csv('3and5wigs_' + tracklist[i], index = False)
    df = pd.read_csv(folders.cured_folder + tracklist[i+1], sep = ',', header = 0)
    df = pd.concat([df, pd.read_csv(folders.cured_folder + tracklist[i+3], sep = ',', header = 0)])
    df.to_csv('3and5wigs_' + tracklist[i+1], index = False)
    i += 4
    del df

### delete the non-merged files
for file in os.listdir():
    if file.find('3and5wigs') != 0:
        os.remove(file)
    else:
        pass

# ### sort rows by location in the direction of transcription, meaning that
# ### on the forward strand numbers go up and on the reverse strand numbers go down.
### print('sorting rows by location')
for file in os.listdir():
    df = pd.read_csv(file, sep = ',', header = 0)
    thestrand = df['strand'][0]
    if thestrand == 0:
        df = df.sort_values(['location', 'end'], axis = 0, ascending = True)
        df.to_csv(file, index = False)
    else:
        df = df.sort_values(['location', 'end'], axis = 0, ascending = False)
        df.to_csv(file, index = False)
    del df

os.chdir(folders.base_folder)
os.system('touch temp.txt')
