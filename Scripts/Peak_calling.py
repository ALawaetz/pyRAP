import pandas as pd
import os
from multiprocessing import Process, Pool
import statistics as stat
import folders
import sys

z_start = float(sys.argv[1])
z_stop = float(sys.argv[2])
z_small = float(sys.argv[3])


### In brief, the script goes through the genome and finds any values above 5 standard deviations of the mean in the 25 or
### 50 bp window downstream start sites or upstream stop sites.
### compute mean value for each position (i) in 25 and 50 bp windows (= 50 and 100 rows, since each position has a 5' and a 3' value)
### moving across the genome in 1 bp steps. i is center except at the ends of genome. 5' values we compare to downstream positions
### and 3' values we compare to upstream positions.
### For 3' values the window excludes the 10 immediate upstream values. This is because stop sites are often relatively broad meaning that
### the end of a transcript is not confined to one excat position but rather spread out on 5 bp approximately. That makes stop sites more difficult to detect
### than start sites that most often can be localised to one excat position. By ommiting in the calculation the 10 immediately upstream values of a potential
### stop site we can keep a high cutoff (5 or 6 standard deviations above the mean) which yields few false positives.
### For all positions cutoff values are determined within two window sizes (50 and 100) as some peaks are detected using big windows and other peaks
### only in small windows, e.g. if two peaks are positioned within 50 bp of each other.
### Start sites generally produce higher peaks than stop sites because as mentioned above start sites are typically confined to one positon as opposed to stop sites.
### Therefore we use higher cutoff values (6 and 7 standard deviations above the mean for small and big windows, respectively) than for stop
### sites (5 and 6 standard deviations above the mean for small and big windows, respectively).
### Lastly, a peak is only called, if it is the highest 5' or 3' value within +/- 5 bp. As rows in dataset are ordered by location and end value (5' or 3') this is determined
### as the maximum value of every second value.

def YourCode(file):
    if '.csv' in file and 'GFF3' not in file:
        df_rend = pd.read_csv(folders.cured_folder + file, sep = ',', header = 0)
        values = df_rend['value'].tolist()
        endlist = df_rend['end'].tolist()
        locationlist = df_rend['location'].tolist()
        del df_rend

        window_size = 100
        window_small = 50
        changing_mean = []
        cutoff = []
        cutoff_small = []
        for i in range(0, len(values)):
            if values[i] == 0:
                changing_mean.append('.')
                cutoff.append(0)
                cutoff_small.append(0)
                continue
            ### At the end of genome, i can't be centered in window
            elif i + window_size >= len(values):
                thegroup = values[i - window_size + 1: i + 1]
                themax = max(values[i - 20 + 2:i + 2][::2])
                theSD = z_start * stat.stdev(thegroup)
                thegroup_small = values[i - window_small + 1: i + 1]
                theSD_small = z_start * stat.stdev(thegroup_small)
            ### At the beginning of genome, i can't be centered in window
            elif i - window_size < 0:
                thegroup = values[i: i + window_size]
                themax = max(values[i:i + 20][::2])
                theSD = z_start * stat.stdev(thegroup)
                thegroup_small = values[i: i + window_small]
                theSD_small = z_start * stat.stdev(thegroup_small)
            ### upstream windows for 5' values
            elif endlist[i] == 5:
                thegroup = values[i: i + window_size]
                themax = max(values[i - 10: i + 10][::2])
                theSD = z_start * stat.stdev(thegroup)
                thegroup_small = values[i: i + window_small]
                theSD_small = z_small * stat.stdev(thegroup_small)
            ### downstream windows for 3' windows
            else:
                thegroup = values[i - window_size: i - 10] + [values[i]]
                themax = max(values[i - 10: i + 10][::2])
                theSD = z_stop * stat.stdev(thegroup)
                thegroup_small = values[i - window_small: i - 10] + [values[i]]
                theSD_small = z_small * stat.stdev(thegroup_small)
            if values[i] > stat.mean(thegroup) + theSD or values[i] > stat.mean(thegroup_small) + theSD_small:
                try:
                    noise_values = [v for v in values[i - 100 : i + 100] if v != 0]
                    noise_count = noise_values.count(values[i]) / len(noise_values)
                except:
                    noise_count = 0
                if values[i] == themax and noise_count <= 0.34:
                    if endlist[i] == 5:
                        changing_mean.append('START')
                    else:
                        changing_mean.append('STOP')
                    cutoff.append(round(abs(values[i] - stat.mean(thegroup))/stat.stdev(thegroup), 1))
                    try:
                        cutoff_small.append(round(abs(values[i] - stat.mean(thegroup_small))/stat.stdev(thegroup_small), 1))
                    except:
                        cutoff_small.append(round(abs(values[i] - stat.mean(thegroup_small))/0.001, 1))
                else:
                    changing_mean.append('.')
                    cutoff.append(round(abs(values[i] - stat.mean(thegroup))/stat.stdev(thegroup), 1))
                    cutoff_small.append(round(abs(values[i] - stat.mean(thegroup_small))/stat.stdev(thegroup_small), 1))
            else:
                changing_mean.append('.')
                try:
                    cutoff.append(round(abs(values[i] - stat.mean(thegroup))/stat.stdev(thegroup), 1))
                except:
                    cutoff.append(0)
                try:
                    cutoff_small.append(round(abs(values[i] - stat.mean(thegroup_small))/stat.stdev(thegroup_small), 1))
                except:
                    cutoff_small.append(0)


        df_rend = pd.read_csv(folders.cured_folder + file, sep = ',', header = 0)
        df_rend['peak'] = changing_mean
        df_rend['shift_cutoff'] = cutoff
        df_rend['shift_cutoff_small'] = cutoff_small

        df_rend.to_csv(os.path.join(folders.peak_folder_0, file), sep = '\t', index = False)

        del changing_mean, cutoff, cutoff_small, values, endlist
    else:
        pass
#
if __name__ == '__main__':
    with Pool(8) as p:
        p.map(YourCode, sorted(os.listdir(folders.cured_folder)))
