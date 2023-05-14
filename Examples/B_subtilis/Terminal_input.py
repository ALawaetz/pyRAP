import os
import glob
import time

time_start = time.time()

script = os.getcwd() + '/pipeline.py'




### Bacillus subtilis WT samples from Lalanne and Deloughery study
### Edit variable names to match the paths on your computer
wigs = '/Users/andreas/Bacillus/Bioinformatics/B_subtilis_wigs_All_WT' # -f {}
chrom = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/B_subtilis_chromsize.txt' # -c {}
guide = 'WT_LB_25s'
names = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/wig_names_Bacillus_WT_Lalanne_and_DeLoughery.txt' # -n {}
operons = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/BSG_operons_extended.csv' # --operons {}
old_annotation = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/LevsGFF_w_predicted_miscRNA.gff3' # --old_annotation {}
output_folder = '/Users/andreas/Bacillus/Bioinformatics/pyRAP_output_7May2023/Thursday_14_May/B_subtilis_WT_Lalanne_and_DeLoughery/' # --output_folder {}
compare_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/CappableSeq_Bacillus.csv'
compare_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/comprehensive_termination_atlas_Bacillus.csv'
operon_starts = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_start_sites_B_subtilis.csv' # --operon_starts {}
operon_stops = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/operon_stop_sites_B_subtilis.csv' # --operon_stops {}
feature_distance = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/Bacillus_feature_distance.txt'

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
--compare_starts {} \
--compare_stops {} \
--operon_starts {} \
--operon_stops {} \
--unique 1 \
--fdr 2 \
--compare_condition {} \
--feature_distance {} \
--source_exclude prediction \
--output All
'''.format(script, wigs, chrom, guide, names, operons, old_annotation, output_folder, compare_starts, compare_stops, operon_starts, operon_stops, guide, feature_distance))


if len([i for i in os.listdir(output_folder) if i != '.DS_Store']) == 2:
    pass
else:
    os.wait()



### E coli
### Edit variable names to match the paths on your computer
wigs = '/Users/andreas/Bacillus/Bioinformatics/E_coli_ALL_wig_files' # -f {}
chrom = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-20230507/chrom_size_files/E_coli_chromsize.txt' # -c {}
old_annotation = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-20230507/Alternative_annotation_files/EcoCyc_w_predicted_miscRNA_k12v2.gff3' # --old_annotation {}
output_folder = '/Users/andreas/Bacillus/Bioinformatics/pyRAP_output_7May2023/Thursday_14_May/E_coli/' # --output_folder {}
guide = 'E_coli_WT'
names = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/E_coli_wig_names.txt' # -n {}
operons = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-20230507/Alternative_operons/E_coli_operons_NC_000913_2_regulonDB.gff3' # -o {}
compare_starts = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-20230507/Alternative_coordinates/start_sites_e_coli_NC_000913_2.csv'
compare_stops = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-20230507/Alternative_coordinates/Term_seq_stop_sites_e_coli_M63_LB_OD_0_4_POOL_NC_000913_2.csv'
operon_starts = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-20230507/Alternative_coordinates/operon_start_sites_e_coli.csv' # --operon_starts {}
operon_stops = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-20230507/Alternative_coordinates/operon_stop_sites_e_coli.csv' # --operon_stops {}
feature_distance = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_distance_E_coli_EcoCyc.txt'
compare_peaks_across_samples = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/compare_peaks_across_samples_E_coli.txt'


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
--compare_starts {} \
--compare_stops {} \
--operon_starts {} \
--operon_stops {} \
--unique 1 \
--fdr 2 \
--compare_condition {} \
--feature_distance {} \
--compare_peaks_across_samples {} \
--source_exclude prediction \
--output All
'''.format(script, wigs, chrom, guide, names, operons, old_annotation, output_folder, compare_starts, compare_stops, operon_starts, operon_stops, guide, feature_distance, compare_peaks_across_samples))


if len([i for i in os.listdir(output_folder) if i != '.DS_Store']) == 2:
    pass
else:
    os.wait()


