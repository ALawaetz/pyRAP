import os
import glob

script = os.getcwd() + '/pipeline.py'


### Bacillus subtilis Lalanne study
### Edit variable names to match the paths on your computer
wigs = '/Users/andreas/Bacillus/Bioinformatics/Bacillus_wig_files_ALL_Lalanne' # -f {}
chrom = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/B_subtilis_chromsize.txt' # -c {}
guide = 'WT_LB_25s'
names = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/Wig_names_Lalanne.txt' # -n {}
operons = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/BSG_operons_no_header.txt' # -o {}
old_annotation = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/LevsGFF_w_predicted_miscRNA.gff3'
old_annotation_stat = '/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Options_parser_folder/Nicolas_et_al_annotation.txt' # --old_annotation {}
output_folder = '/users/andreas/bacillus/Bioinformatics/pyRAP_test_27Nov22/B_subtilis_7Jan_Lalanne/' # --output_folder {}
compare_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/CappableSeq_Bacillus.csv'
compare_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/comprehensive_termination_atlas_Bacillus.csv'
operon_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_start_sites_B_subtilis.csv' # --operon_starts {}
operon_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_stop_sites_B_subtilis.csv' # --operon_stops {}
feature_distance = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/Bacillus_feature_distance.txt'
compare_peaks_across_samples = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_peaks_across_samples_bacillus_Lalanne.txt'

### If you only wanna analyse a small region of the genome specify location here
genome_start = 1030000 # --genome_start {}
genome_end = 1042706 # --genome_end {}


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
