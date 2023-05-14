import os
import glob
import time

time_start = time.time()

script = os.getcwd() + '/pipeline.py'

### S. aureus
### Edit variable names to match the paths on your computer
wigs = '/Users/andreas/Bacillus/Bioinformatics/S_aureus_wigs' # -f {}
chrom = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/S_aureus_chromsize.txt' # -c {}
old_annotation = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/aureoWiki_w_predicted_miscRNA.gff3' # --old_annotation {}
output_folder = '/Users/andreas/Bacillus/Bioinformatics/pyRAP_output_7May2023/Thursday_14_May/S_aureus/' # --output_folder {}
guide = 'S_aureus_WT'
names = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/names_of_conditions_S_aureus.txt' # -n {}
operons = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/aureoOperons.gff3' # -o {}
operon_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/start_sites_s_aureus.csv'
operon_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/stop_sites_s_aureus.csv'
feature_distance = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_distance_S_aureus.txt'
compare_peaks_across_samples = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_peaks_across_samples_S_aureus.txt'


### If you only wanna analyse a small region of the genome specify location here
genome_start = 2385394 # --genome_start {}
genome_end = 2387825 # --genome_end {}


os.system('''python3 {} \
-f {} \
-c {} \
-g {} \
-n {} \
--operons {} \
--old_annotation {} \
--output_folder {} \
--statistics \
--operon_starts {} \
--operon_stops {} \
--unique 1 \
--feature_distance {} \
--compare_peaks_across_samples {} \
--source_exclude prediction \
--fdr 2 \
--z_start 7 \
--z_stop 5 \
--z_ii 6 \
--search_distance_0 50 \
--output All
'''.format(script, wigs, chrom, guide, names, operons, old_annotation, output_folder, operon_starts, operon_stops, feature_distance, compare_peaks_across_samples))





time_stop = time.time()

print('Time to execute: {} hours'.format((time_stop - time_start) / 60 / 60))