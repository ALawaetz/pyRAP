import os
import glob

script = os.getcwd() + '/pipeline.py'

### Staphylococcus aureus
wigs = '/Users/andreas/Bacillus/Bioinformatics/S_aureus_wigs' # -f {}
chrom = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/S_aureus_chromsize.txt' # -c {}
guide = 'S_aureus_WT'
names = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/names_of_conditions_S_aureus.txt' # -n {}
operons = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/aureoOperons.gff3' # -o {}
old_annotation = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/aureoWiki_w_predicted_miscRNA.gff3' # --old_annotation {}
old_annotation_stat = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/aureoWiki.gff3'
output_folder = '/users/andreas/bacillus/Bioinformatics/pyRAP_test_27Nov22/S_aureus_7Jan/' # --output_folder {}
compare_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/start_sites_s_aureus.csv'
compare_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/stop_sites_s_aureus.csv'
feature_distance = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_distance_S_aureus.txt'
compare_peaks_across_samples = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_peaks_across_samples_S_aureus.txt'

### If you only wanna analyse a small region of the genome specify location here. Add flags to os.system input below.
genome_start = 0 # --genome_start {}
genome_end = 5000 # --genome_end {}

os.system('''python3 {} \
-f {} \
-c {} \
-g {} \
-n {} \
-o {} \
--old_annotation {} \
--old_annotation_stat {} \
--output_folder {} \
--statistics \
--operon_starts {} \
--operon_stops {} \
--fdr 8 \
--unique 1 \
--search_distance_0 50 \
--z_start 7 \
--z_stop 5 \
--z_ii 6 \
--output All \
--compare_condition {} \
--feature_distance {} \
--compare_peaks_across_samples {} \
'''.format(script, wigs, chrom, guide, names, operons, old_annotation, old_annotation_stat, output_folder, compare_starts, compare_stops, guide, feature_distance, compare_peaks_across_samples))
