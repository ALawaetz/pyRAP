import os
import time

time_start = time.time()

homedir = os.getcwd()
script = homedir + '/pipeline.py'

### Example
### This script works without edit to analyse WT and pnpA conditions from the study by Lalanne et al.
### Simply navigate to the github directory and type in terminal 
### python3 Example_Terminal_input.py

### To analyse other wig files, modify paths accordingly.

wigs = f'{homedir}/Wig_files' # -f {}
chrom = f'{homedir}/Chrom_size_files/B_subtilis_chromsize.txt' # -c {}
guide = 'WT_LB_25s'
names = f'{homedir}/Txt_files_examples/-n/wig_names_Bacillus_Lalanne.txt' # -n {}
operons = f'{homedir}/Alternative_operons/BSG_operons_extended.csv' # --operons {}
old_annotation = f'{homedir}/Alternative_annotation_files/LevsGFF_w_predicted_miscRNA.gff3' # --old_annotation {}
output_folder = f'{homedir}/pyRAP_output/' # --output_folder {}
compare_starts = f'{homedir}/Alternative_coordinates/CappableSeq_Bacillus.csv'
compare_stops = f'{homedir}/Alternative_coordinates/comprehensive_termination_atlas_Bacillus.csv'
operon_starts = f'{homedir}/Alternative_coordinates/operon_start_sites_B_subtilis.csv' # --operon_starts {}
operon_stops = f'{homedir}/Alternative_coordinates/operon_stop_sites_B_subtilis.csv' # --operon_stops {}
feature_distance = f'{homedir}/Txt_files_examples/--feature_distance/Bacillus_feature_distance.txt'

### If you only want to analyse a small region of the genome specify location here
### and remember to add the flags in the terminal input below.
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
--compare_starts {} \
--compare_stops {} \
--operon_starts {} \
--operon_stops {} \
--unique 1 \
--fdr 2 \
--compare_condition {} \
--feature_distance {} \
--source_exclude prediction \
--genome_start 2385394 \
--genome_end 2387825 \
--output All
'''.format(script, wigs, chrom, guide, names, operons, old_annotation, output_folder, compare_starts, compare_stops, operon_starts, operon_stops, guide, feature_distance))


time_stop = time.time()

print('Time to execute: {} hours'.format((time_stop - time_start) / 60 / 60))
