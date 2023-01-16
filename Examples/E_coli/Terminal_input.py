import os
import glob

script = os.getcwd() + '/pipeline.py'

### E. coli
wigs = '/Users/andreas/Bacillus/Bioinformatics/E_coli_ALL_wig_files' # -f {}
chrom = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/E_coli_chromsize.txt' # -c {}
old_annotation = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/E_coli_NC_000913_2_w_predicted_miscRNA.gff3' # --old_annotation {}
old_annotation_stat = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/E_coli_NC_000913_2_duplicates_removed.gff3' # --old_annotation {}
output_folder = '/users/andreas/bacillus/Bioinformatics/pyRAP_test_27Nov22/E_coli_7Jan/' # --output_folder {}
guide = 'E_coli_WT'
names = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/E_coli_wig_names.txt' # -n {}
operons = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/E_coli_operons_NC_000913_2_regulonDB.gff3' # -o {}
compare_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/start_sites_e_coli_NC_000913_2.csv'
compare_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/Term_seq_stop_sites_e_coli_M63_LB_OD_0_4_POOL_NC_000913_2.csv'
operon_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_start_sites_e_coli.csv' # --operon_starts {}
operon_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_stop_sites_e_coli.csv' # --operon_stops {}
feature_distance = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_distance_E_coli.txt'
compare_peaks_across_samples = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_peaks_across_samples_E_coli.txt'

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
--compare_starts {} \
--compare_stops {} \
--operon_starts {} \
--operon_stops {} \
--unique 1 \
--output All \
--compare_condition {} \
--feature_distance {} \
--compare_peaks_across_samples {} \
'''.format(script, wigs, chrom, guide, names, operons, old_annotation, old_annotation_stat, output_folder, compare_starts, compare_stops, operon_starts, operon_stops, guide, feature_distance, compare_peaks_across_samples))
