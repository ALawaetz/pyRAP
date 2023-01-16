import os
import glob
import time

timestart = time.time()

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

try:
    os.system('rm {}{}'.format('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up', '/end_pipe.txt'))
except:
    pass

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


if len(glob.glob('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up/' + '*end_pipe.txt')) == 1:
    pass
else:
    os.wait()

os.system('rm {}{}'.format('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up', '/end_pipe.txt'))

## Bacillus subtilis Deloughery study
## Edit variable names to match the paths on your computer
wigs = '/Users/andreas/Bacillus/Bioinformatics/Bacillus_wig_files_ALL_Deloughery' # -f {}
chrom = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/B_subtilis_chromsize.txt' # -c {}
guide = 'WT_Bsub_LB_exp'
names = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/wig_names_Deloughery.txt' # -n {}
operons = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/BSG_operons_no_header.txt' # -o {}
old_annotation = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/LevsGFF_w_predicted_miscRNA.gff3'
old_annotation_stat = '/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Options_parser_folder/Nicolas_et_al_annotation.txt' # --old_annotation {}
output_folder = '/users/andreas/bacillus/Bioinformatics/pyRAP_test_27Nov22/B_subtilis_7Jan_Deloughery/' # --output_folder {}
compare_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/CappableSeq_Bacillus.csv'
compare_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/comprehensive_termination_atlas_Bacillus.csv'
operon_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_start_sites_B_subtilis.csv' # --operon_starts {}
operon_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_stop_sites_B_subtilis.csv' # --operon_stops {}
feature_distance = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/Bacillus_feature_distance.txt'
compare_peaks_across_samples = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/Compare_peaks_across_samples_bacillus_Deloughery.txt'

### If you only wanna analyse a small region of the genome specify location here
genome_start = 1030600 # --genome_start {}
genome_end = 1038900 # --genome_end {}


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
--z_stop 5 \
--z_ii 6 \
--compare_condition {} \
--feature_distance {} \
--compare_peaks_across_samples {} \
'''.format(script, wigs, chrom, guide, names, operons, old_annotation, old_annotation_stat, output_folder, compare_starts, compare_stops, operon_starts, operon_stops, guide, feature_distance, compare_peaks_across_samples))

if len(glob.glob('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up/' + '*end_pipe.txt')) == 1:
    pass
else:
    os.wait()

os.system('rm {}{}'.format('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up', '/end_pipe.txt'))

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

### If you only wanna analyse a small region of the genome specify location here
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

if len(glob.glob('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up/' + '*end_pipe.txt')) == 1:
    pass
else:
    os.wait()

os.system('rm {}{}'.format('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up', '/end_pipe.txt'))

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

### If you only wanna analyse a small region of the genome specify location here
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


if len(glob.glob('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up/' + '*end_pipe.txt')) == 1:
    pass
else:
    os.wait()

os.system('rm {}{}'.format('/Users/andreas/Bacillus/Bioinformatics/AllWigFiles/To_Lauren/Tidy_up', '/end_pipe.txt'))

timestop = time.time()

print('time: {}'.format(timestop - timestart))
